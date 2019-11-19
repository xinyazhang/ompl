/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2017, Xinya Zhang, University of Texas at Austin
*  Copyright (c) 2008, Willow Garage, Inc.
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Willow Garage nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Author: Xinya Zhang
 *         Modified from RRT codebase by Ioan Sucan */

#include "ompl/geometric/planners/rrt/ReRRT.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/tools/config/SelfConfig.h"
#include <limits>
#include <ostream>
#include <unordered_map>

ompl::geometric::ReRRT::ReRRT(const base::SpaceInformationPtr &si) : base::Planner(si, "Re-RRT")
{
    specs_.approximateSolutions = true;
    specs_.directed = true;

    goalBias_ = 0.05;
    maxDistance_ = 0.0;
    lastGoalMotion_ = nullptr;
    knearest_ = 1;

    Planner::declareParam<double>("range", this, &ReRRT::setRange, &ReRRT::getRange, "0.:1.:10000.");
    Planner::declareParam<double>("goal_bias", this, &ReRRT::setGoalBias, &ReRRT::getGoalBias, "0.:.05:1.");
}

ompl::geometric::ReRRT::~ReRRT()
{
    freeMemory();
}

void ompl::geometric::ReRRT::clear()
{
    Planner::clear();
    sampler_.reset();
    freeMemory();
    if (nn_)
        nn_->clear();
    lastGoalMotion_ = nullptr;
}

void ompl::geometric::ReRRT::setup()
{
    Planner::setup();
    tools::SelfConfig sc(si_, getName());
    sc.configurePlannerRange(maxDistance_);

    if (!nn_)
        nn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Motion*>(this));
    nn_->setDistanceFunction(std::bind(&ReRRT::distanceFunction, this, std::placeholders::_1, std::placeholders::_2));
}

void ompl::geometric::ReRRT::setSampleSet(const Eigen::Ref<const Eigen::MatrixXd> Q)
{
    predefined_samples_ = Q;
}

void ompl::geometric::ReRRT::setSampleSetEdges(const Eigen::Ref<const Eigen::MatrixXi> QB,
                                               const Eigen::Ref<const Eigen::MatrixXi> QE,
                                               const Eigen::Ref<const Eigen::MatrixXi> QEB)
{
    predefined_sample_tree_bases_ = QB;
    predefined_sample_edges_ = QE;
    predefined_sample_edge_bases_ = QEB;
}

void ompl::geometric::ReRRT::setSampleSetFlags(const Eigen::Ref<const Eigen::Matrix<uint32_t, -1, 1>> QF)
{
    pds_flags_ = QF;
}

void ompl::geometric::ReRRT::getSampleSetConnectivity(Eigen::SparseMatrix<int>& C)
{
    C.resize(1, predefined_samples_.rows());
    C.setZero();
    C.setFromTriplets(connectivity_tup_.begin(), connectivity_tup_.end());
}

void ompl::geometric::ReRRT::getCompactGraph(Eigen::Matrix<int64_t, -1, 1>& nouveau_vertex_id,
	                                     Eigen::MatrixXd& nouveau_vertices,
	                                     Eigen::Matrix<int64_t, -1, 2>& edges) const
{
    if (predefined_samples_.rows() <= 0) {
	OMPL_ERROR("%s: getCompactGraph requires setSampleSet", getName().c_str());
	return ;
    }
    std::unordered_set<Motion*> motions; // to store motions we are going to keep
    {
	std::vector<Motion*> all_motions;
	nn_->list(all_motions);
	for (auto m : all_motions) {
	    // Skip motion that's not in PDS
	    if (m->is_nouveau)
		continue;
	    // Add motions from PDS to root into motions
	    do {
		motions.emplace(m);
		m = m->parent;
	    } while (m);
	}
    }
    std::vector<Motion*> nouveau_motions;
    for (const auto& m : motions) {
	if (m->motion_type != Motion::MOTION_OF_SAMPLE)
	    continue;
	if (!m->is_nouveau)
	    continue;
	nouveau_motions.emplace_back(m);
    }
    std::sort(nouveau_motions.begin(),
	      nouveau_motions.end(),
	      [](const Motion* lhs, const Motion* rhs) {
	         return lhs->motion_index < rhs->motion_index;
	      });
    nouveau_vertex_id.resize(nouveau_motions.size());
    nouveau_vertices.resize(nouveau_motions.size(), predefined_samples_.cols());
    size_t rowi = 0;
    std::vector<double> reals;
    for (const auto& current : nouveau_motions) {
	// Vid
	nouveau_vertex_id(rowi) = current->motion_index;
	// V
	si_->getStateSpace()->copyToReals(reals, current->state);
	nouveau_vertices.row(rowi) = Eigen::Map<Eigen::VectorXd>(reals.data(), reals.size(), 1);
	rowi++;
    }
    // std::cerr << "nouveau motions filled\n";
    rowi = 0;
    // Number of edges:
    // 1. Tree structure, E=V-1
    // 2. We do not track the root (initial state) so E=V'
    size_t nedge = 0;
    for (const auto& current : motions) {
	Motion* parent = current->parent;
	if (parent)
	    nedge++;
    }
    edges.setZero(nedge, 2);
    // std::cerr << "edges resized to " << motions.size() << ", 2" << std::endl;
    for (const auto& current : motions) {
	Motion* parent = current->parent;
	if (!parent)
	    continue ;
	int64_t cindex = current->motion_index;
	if (current->motion_type == Motion::MOTION_OF_START)
	    cindex = -(cindex + 1);
	int64_t pindex = parent->motion_index;
	if (parent->motion_type == Motion::MOTION_OF_START)
	    pindex = -(pindex + 1);
	// E
	if (cindex < -2 || pindex < -2)
	    throw std::runtime_error("invalid index: cindex "+std::to_string(cindex)+" pindex " +std::to_string(pindex));
	edges.row(rowi) << cindex, pindex;
	// std::cerr << "edges(" << rowi << ")" << std::endl;
	rowi++;
    }
}

void ompl::geometric::ReRRT::freeMemory()
{
    if (nn_)
    {
        std::vector<Motion*> motions;
        nn_->list(motions);
        for (unsigned int i = 0 ; i < motions.size() ; ++i)
        {
            if (motions[i]->state)
                si_->freeState(motions[i]->state);
            delete motions[i];
        }
    }
}

namespace {

std::ostream&
    print_state(std::ostream& fout,
	const ompl::base::State* state,
	const ompl::base::SpaceInformationPtr& si)
{
        std::vector<double> reals;
        si->getStateSpace()->copyToReals(reals, state);
        std::copy(reals.begin(), reals.end(), std::ostream_iterator<double>(fout, " "));
        return fout;
}

}

ompl::base::PlannerStatus ompl::geometric::ReRRT::solve(const base::PlannerTerminationCondition &ptc)
{
    checkValidity();
    base::Goal                 *goal   = pdef_->getGoal().get();
    base::GoalSampleableRegion *goal_s = dynamic_cast<base::GoalSampleableRegion*>(goal);

    while (const base::State *st = pis_.nextStart())
    {
        Motion *motion = new Motion(si_);
	motion->motion_type = Motion::MOTION_OF_START;
	motion->motion_index = nn_->size();
        si_->copyState(motion->state, st);
	std::cout << getName() << ": Initialize by adding state ";
	print_state(std::cout, st, si_) << std::endl;
        nn_->add(motion);
    }
    // OMPL_INFORM("%s: Actual class of Goal: %s", getName().c_str(), typeid(*goal).name());

    if (nn_->size() == 0)
    {
        OMPL_ERROR("%s: There are no valid initial states!", getName().c_str());
        return base::PlannerStatus::INVALID_START;
    }

    if (!sampler_)
        sampler_ = si_->allocStateSampler();

    motion_of_start_size_ = nn_->size();
    OMPL_INFORM("%s: Starting planning with %u states already in datastructure", getName().c_str(), nn_->size());
    OMPL_INFORM("%s: Max Distance %f", getName().c_str(), maxDistance_);

    Motion *solution  = nullptr;
    Motion *approxsol = nullptr;
    double  approxdif = std::numeric_limits<double>::infinity();
    Motion *rmotion   = new Motion(si_);
    base::State *rstate = rmotion->state;
    base::State *xstate = si_->allocState();
    base::State *restate = si_->allocState();
    size_t nsample_created = 0;
    size_t nsample_injected = 0;
    bool sat = false;
    // Predefined sample set support
    ssize_t iteration = 0;
    connectivity_tup_.clear();
    ssize_t pds_size = predefined_samples_.rows();
    if (pds_size > 0 && knearest_ > 1) {
	throw std::runtime_error("ReRRT: setSampleSet is incompatitable with setKNearest");
    }
    bool enable_predefined_samples = pds_size > 0;
    bool enable_pds_edges = enable_predefined_samples && predefined_sample_edges_.rows() > 0;
    OMPL_INFORM("%s: PDS %s", getName().c_str(), enable_predefined_samples ? "Enabled" : "Disabled");
    OMPL_INFORM("%s: PDS Edge %s", getName().c_str(), enable_pds_edges ? "Enabled" : "Disabled");
    double retraction_ratio = 1.0;
    bool use_retracted_sample = getOptionReal("retract", retraction_ratio);
    OMPL_INFORM("%s: Retraction enabled %s ration ratio %f", getName().c_str(), use_retracted_sample ? "True" : "False", retraction_ratio);
    long bloom_limit = -1;
    bool has_bloom_limit = getOptionInt("bloom_limit", bloom_limit);
    OMPL_INFORM("%s: bloom_limit enabled %s value %ld", getName().c_str(), has_bloom_limit ? "True" : "False", bloom_limit);
    long ec_limit = -1;
    bool has_ec_limit = getOptionInt("ec_limit", ec_limit);
    OMPL_INFORM("%s: edge connection limit enabled %s value %ld",
		getName().c_str(),
		has_ec_limit ? "True" : "False",
		ec_limit);
    long pca_after = -1;

    while (ptc() == false && !sat)
    {
	bool injecting = (nsample_created >= sample_injection_ && nsample_injected < samples_to_inject_.size());
	if (enable_predefined_samples) {
	    // skip already added PDS
	    while (enable_pds_edges &&
		   iteration < pds_size &&
		   (pds_flags_(iteration) & _PDS_FLAG_ALREADY_IN_TREE)) {
                // OMPL_INFORM("\t%s: Skip PDS %ld since it's already in tree",
                //            getName().c_str(),
                //            iteration);
		iteration++;
            }
	    // Early termination
	    if (iteration >= predefined_samples_.rows())
		break;
	    // Override sampler
	    si_->getStateSpace()->copyFromEigen3(rstate, predefined_samples_.row(iteration));
#if 0
	    std::cout << "PDS [" << iteration << "]\t";
	    print_state(std::cout, rstate, si_);
	    std::cout << std::endl;
#endif
	} else if (injecting) {
	    si_->getStateSpace()->copyFromReals(rstate, samples_to_inject_[nsample_injected]);
#if 0
	    std::cout << getName() << ": Injecting state ";
	    print_state(std::cout, rstate, si_) << std::endl;
#endif
	    nsample_injected++;
	} else {
	    /* sample random state (with goal biasing) */
	    if (goal_s && rng_.uniform01() < goalBias_ && goal_s->canSample())
		goal_s->sampleGoal(rstate);
	    else
		sampler_->sampleUniform(rstate);
	}

	std::vector<Motion*> nmotions;
	if (knearest_ == 1) {
	    nmotions.emplace_back(nn_->nearest(rmotion));
	} else {
	    nn_->nearestK(rmotion, knearest_, nmotions);
	}
#if 0
	OMPL_INFORM("%s: Find %d neighbors", getName().c_str(), nmotions.size());
#endif
        /* find closest state in the tree */
	for (auto nmotion: nmotions) {
	    base::State *dstate = rstate;

	    /* find state to add */
	    if (!enable_predefined_samples && !injecting) {
		double d = si_->distance(nmotion->state, rstate);
		if (d > maxDistance_) {
		    si_->getStateSpace()->interpolate(nmotion->state, rstate, maxDistance_ / d, xstate);
		    dstate = xstate;
		}
	    }

	    bool to_create_motion = true;
	    std::pair<base::State*, double> lastValid(restate, 0.0);
	    bool directly_connected = si_->checkMotion(nmotion->state, dstate, lastValid);
	    if (!directly_connected) {
		if (lastValid.first != nullptr && lastValid.second < 1.0 - 1e-4 && lastValid.second > 1e-4) {
		    if (use_retracted_sample) {
			si_->getStateSpace()->interpolate(nmotion->state, lastValid.first, retraction_ratio, dstate);
		    } else {
			si_->copyState(dstate, lastValid.first);
		    }
		} else
		    to_create_motion = false;
		if (injecting) {
#if 0
		    OMPL_INFORM("%s: Retraction performed, state: %p, tau: %f", getName().c_str(), lastValid.first, lastValid.second);
		    print_state(std::cout, dstate, si_);
		    std::cout << std::endl;
#endif
		}
	    }
#if 0
	    if (injecting) {
		print_state(std::cout, dstate, si_);
		std::cout << std::endl;
	    }
#endif
	    if (to_create_motion) {
		/* create a motion */
		Motion *motion = new Motion(si_);
		si_->copyState(motion->state, dstate);
		motion->parent = nmotion;
		if (enable_predefined_samples) {
		    motion->is_nouveau = !directly_connected;
		    motion->motion_type = Motion::MOTION_OF_SAMPLE;
		    motion->motion_index = iteration;
		    // OMPL_INFORM("CT[%ld] := ..., directly connected = %d", iteration, (int)directly_connected);
		} else {
		    motion->is_nouveau = true;
		    motion->motion_type = Motion::MOTION_OF_SAMPLE;
		    motion->motion_index = nn_->size() - motion_of_start_size_;
		}
#if 0
		std::cout.precision(17);
		std::cout << "[" << iteration << "]" << (directly_connected ? "Direct" : "Indirect") << " RDT: Edge Creation " << std::endl
		          << "\t<" << nmotion->motion_index << "> ";
		print_state(std::cout, nmotion->state, si_);
		std::cout << "\n\t<" << motion->motion_index << "> ";
		print_state(std::cout, dstate, si_);
		std::cout << std::endl;

		if (injecting) {
		    std::cout << "Re-RRT: Creating Edge From " << std::endl << "\t";
		    print_state(std::cout, nmotion->state, si_) << std::endl << "\tTo\n\t";
		    print_state(std::cout, dstate, si_) << std::endl;
		}
#endif
		nn_->add(motion);
		double dist = 0.0;
		sat = goal->isSatisfied(motion->state, &dist);
		if (sat)
		{
		    approxdif = dist;
		    solution = motion;
		    break;
		}
		if (dist < approxdif)
		{
		    approxdif = dist;
		    approxsol = motion;
		}
		nsample_created++;
		if (directly_connected) {
                    if (enable_pds_edges)
                        addWholePdsTree(iteration, motion);
                    if (enable_predefined_samples) {
                        connectivity_tup_.emplace_back(0, iteration, 1);
                        // OMPL_INFORM("SSC[0,%ld] := 1", iteration);
                        if (pds_flags_.rows() > 0 && (pds_flags_(iteration) & PDS_FLAG_TERMINATE)) {
                            OMPL_INFORM("Early Termination at Iteration %ld: connected to open space", iteration);
                            solution = motion;
                            sat = true;
                        }
                    }
                    break;
                }
	    }
	}
	iteration++;
	if (has_bloom_limit && nn_->size() > bloom_limit) {
	    OMPL_INFORM("%s: bloom_limit reached, leaving", getName().c_str());
	    sat = true;
	}
	if (has_ec_limit && si_->getCheckedMotionCount() > ec_limit) {
	    OMPL_INFORM("%s: ec_limit %ld reached after %ld connections, leaving",
			getName().c_str(),
			ec_limit,
			si_->getCheckedMotionCount());
	    sat = true;
	}
    }

    bool solved = false;
    bool approximate = false;
    if (solution == nullptr)
    {
        solution = approxsol;
        approximate = true;
    }

    if (solution != nullptr)
    {
        lastGoalMotion_ = solution;

        /* construct the solution path */
        std::vector<Motion*> mpath;
        while (solution != nullptr)
        {
            mpath.push_back(solution);
            solution = solution->parent;
        }

        /* set the solution path */
        PathGeometric *path = new PathGeometric(si_);
        for (int i = mpath.size() - 1 ; i >= 0 ; --i)
            path->append(mpath[i]->state);
        pdef_->addSolutionPath(base::PathPtr(path), approximate, approxdif, getName());
        solved = true;
    }

    si_->freeState(xstate);
    si_->freeState(restate);

    if (rmotion->state)
        si_->freeState(rmotion->state);
    delete rmotion;

    OMPL_INFORM("%s: Created %u states. Total edge connections %lu", getName().c_str(), nn_->size(),
                si_->getCheckedMotionCount());

    return base::PlannerStatus(solved, approximate);
}

void ompl::geometric::ReRRT::setStateInjection(size_t start_from, std::vector<std::vector<double>> samples)
{
    sample_injection_ = start_from;
    samples_to_inject_ = std::move(samples);
}

void ompl::geometric::ReRRT::getPlannerData(base::PlannerData &data) const
{
    Planner::getPlannerData(data);

    std::vector<Motion*> motions;
    if (nn_)
        nn_->list(motions);

    if (lastGoalMotion_)
        data.addGoalVertex(base::PlannerDataVertex(lastGoalMotion_->state));

    for (unsigned int i = 0 ; i < motions.size() ; ++i)
    {
        if (motions[i]->parent == nullptr)
            data.addStartVertex(base::PlannerDataVertex(motions[i]->state));
        else
            data.addEdge(base::PlannerDataVertex(motions[i]->parent->state),
                         base::PlannerDataVertex(motions[i]->state));
    }
}

void ompl::geometric::ReRRT::addWholePdsTree(ssize_t dc_pds_index,
                                             ompl::geometric::ReRRT::Motion* anchor_pds)
{
    int current_pds_tree_index = 0;
    const auto& Q = predefined_samples_;
    const auto& QB = predefined_sample_tree_bases_;
    const auto& QE = predefined_sample_edges_;
    const auto& QEB = predefined_sample_edge_bases_;
    int pds_tree_number = QB.rows();
    while (current_pds_tree_index < pds_tree_number) {
        int min_range, max_range;
        min_range = QB(current_pds_tree_index);
        if (current_pds_tree_index + 1 < QB.rows())
            max_range = QB(current_pds_tree_index + 1);
        else
            max_range = -1;
        if (min_range <= dc_pds_index && (max_range < 0 || dc_pds_index < max_range))
            break;
        current_pds_tree_index += 1;
    }
    // OMPL_INFORM("%s: add Whole PDS tree from %ld, tree id %d", getName().c_str(),
    //            dc_pds_index, current_pds_tree_index);
    int vbase = QB(current_pds_tree_index);
    int ebase = QEB(current_pds_tree_index);
    int erange_min = ebase;
    int erange_max = current_pds_tree_index + 1 < pds_tree_number ?
                     QEB(current_pds_tree_index + 1) :
                     QE.rows();
    std::unordered_map<int, std::vector<int>> tree_edges(std::min(5, erange_max - erange_min / 4));
    for (int ei = erange_min; ei < erange_max; ei++) {
        tree_edges[QE(ei, 0)].emplace_back(QE(ei, 1));
        tree_edges[QE(ei, 1)].emplace_back(QE(ei, 0));
    }
    std::function<void(ssize_t, Motion*)> dfs;
    dfs = [&](ssize_t pds_index, Motion* anchor) {
        for (const auto& ni : tree_edges[pds_index]) {
            if (pds_flags_[ni] & _PDS_FLAG_ALREADY_IN_TREE)
                continue;
            Motion *motion = new Motion(si_);
	    si_->getStateSpace()->copyFromEigen3(motion->state, Q.row(ni));
            motion->parent = anchor;
            motion->is_nouveau = false;
            motion->motion_type = Motion::MOTION_OF_SAMPLE;
            motion->motion_index = ni;
            nn_->add(motion);
            pds_flags_[ni] |= _PDS_FLAG_ALREADY_IN_TREE;
            connectivity_tup_.emplace_back(0, ni, 1);
            // OMPL_INFORM("\t%s: add Whole PDS ID %d", getName().c_str(), ni);
            dfs(ni, motion);
        }
    };
    dfs(dc_pds_index, anchor_pds);
}

// vim: sw=4 expandtab
