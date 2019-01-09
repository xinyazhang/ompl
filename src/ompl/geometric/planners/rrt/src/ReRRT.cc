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

    OMPL_INFORM("%s: Starting planning with %u states already in datastructure", getName().c_str(), nn_->size());
    OMPL_INFORM("%s: Max Distance %f", getName().c_str(), maxDistance_);

    Motion *solution  = nullptr;
    Motion *approxsol = nullptr;
    double  approxdif = std::numeric_limits<double>::infinity();
    Motion *rmotion   = new Motion(si_);
    base::State *rstate = rmotion->state;
    base::State *xstate = si_->allocState();
    base::State *restate = si_->allocState();
    ssize_t nsample_counter = 0; // Count the number of normal samples
    ssize_t sample_injection_off = 0; // Track the sample to inject with offset
    bool sat = false;

    while (ptc == false && !sat)
    {
	bool oor = (sample_injection_off >= samples_to_inject_.size()); // out of range
	bool injecting = (nsample_counter >= std::abs(sample_injection_) && !oor);
	if (injecting) {
	    si_->getStateSpace()->copyFromReals(rstate, samples_to_inject_[sample_injection_off]);
#if 0
	    std::cout << getName() << ": Injecting state ";
	    print_state(std::cout, rstate, si_) << std::endl;
#endif
	} else {
	    /* sample random state (with goal biasing) */
	    if (goal_s && rng_.uniform01() < goalBias_ && goal_s->canSample())
		goal_s->sampleGoal(rstate);
	    else
		sampler_->sampleUniform(rstate);
	}
	// Reset nsample_counter if periodicity is required
	if (sample_injection_ < 0 && oor) {
	    nsample_counter = 0;
	    sample_injection_off = 0;
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
	bool directly_connected = false;
	for (auto nmotion: nmotions) {
	    base::State *dstate = rstate;

	    /* find state to add */
	    double d = si_->distance(nmotion->state, rstate);
	    if (d > maxDistance_ && !injecting)
	    {
		si_->getStateSpace()->interpolate(nmotion->state, rstate, maxDistance_ / d, xstate);
		dstate = xstate;
	    }

	    bool to_create_motion = true;
	    std::pair<base::State*, double> lastValid(restate, 0.0);
	    directly_connected = si_->checkMotion(nmotion->state, dstate, lastValid);
	    if (!directly_connected) {
		if (lastValid.first != nullptr && lastValid.second < 1.0 - 1e-4 && lastValid.second > 1e-4) {
		    si_->copyState(dstate, lastValid.first);
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
#if 0
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
		nsample_counter++;
	    }
	    if (directly_connected)
		break;
	}
	// Special processing for sample injection case
	if (injecting) {
	    if (!directly_connected) {
		// Proceed to the next
		sample_injection_off += 1;
	    } else {
		// Remove current sample
		// Note: do NOT add the offset
		samples_to_inject_.erase(samples_to_inject_.begin() + sample_injection_off);
	    }
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

    OMPL_INFORM("%s: Created %u states", getName().c_str(), nn_->size());

    return base::PlannerStatus(solved, approximate);
}

void ompl::geometric::ReRRT::setStateInjection(ssize_t start_from, std::vector<std::vector<double>> samples)
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

// vim: sw=4
