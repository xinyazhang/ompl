/*********************************************************************
* Software License Agreement (BSD License)
*
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

/* Author: Ioan Sucan */

#include "ompl/geometric/planners/rrt/RRTForest.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/tools/config/SelfConfig.h"

#include <queue>

namespace ompl {
namespace geometric {

struct RRTForest::Private {
    using SpaceInformationPtr = base::SpaceInformationPtr;

    enum {
        START_TREE = 0,
        GOAL_TREE = 1,
        TOTAL_INITIAL_TREE
    };
    SpaceInformationPtr si_;

    Private() = delete;
    Private(const SpaceInformationPtr& si)
    {
        si_ = si;
    }

    /** \brief collection of trees */
    std::vector<TreeData> forest_;

    /** \brief Tree Adjacency Matrix (direct, signed) **/
    Eigen::MatrixXi treeAdj_;

    /** \brief data for disjoint-set algorithm **/
    Eigen::VectorXi treeDisjointSetId_;

    struct InterTreeConnection {
        TreeData tree1, tree2;
        Motion *motion1, *motion2;

        InterTreeConnection(TreeData t1, TreeData t2, Motion *m1, Motion *m2)
            : tree1(t1), tree2(t2), motion1(m1), motion2(m2)
        {
            assert(!!t1);
            assert(!!t2);
            assert(!!m1);
            assert(!!m2);
        }
    };

    /** \brief Connections between trees **/
    std::vector<InterTreeConnection> treeConnections_;

    /** \brief Initialize the forest.

        roots: samples per column.
        distance_function: distance function for NearestNeighbors

            Side effects:
                forest_ consists of roots.cols() + 2 TreeData objects
                forest_[0] and [1] are reserved for initial tree and goal tree respectively. **/
    void initForest(const Eigen::MatrixXd& roots,
                    RRTForest* rrt_forest,
                    std::function<double(Motion*, Motion*)> distance_function);

    void unionTrees(const InterTreeConnection&);

    int findTreeSetId(TreeData tree);

    void freeForest();

    void clearForest(); // CAVEAT: this does not free the indirect references.

    std::shared_ptr<PathGeometric>
    getSolutionPath();

    // Add the segment from the root of the tree to the next inter-tree
    // connection
    void addToPath(std::shared_ptr<PathGeometric> path,
                   TreeData tree,
                   int itc_next);

    // Add the segment between two inter-tree connections to path
    void addToPath(std::shared_ptr<PathGeometric> path,
                   int itc_prev,
                   int itc_next);

    // Add the segment from the next inter-tree
    // connection to the root of the tree
    void addToPath(std::shared_ptr<PathGeometric> path,
                   int itc_prev,
                   TreeData tree);

    // Add the segment between two nodes with lowest common ancestor
    // algorithm. All three variants addToPath use this for their
    // implementations
    //
    // Note: this is only called ONCE per tree, so no need to use Tarjan's
    // algorithm
    void addPathLCA(std::shared_ptr<PathGeometric> path,
                    Motion *node_from,
                    Motion *node_to);
};

void RRTForest::Private::initForest(const Eigen::MatrixXd& roots,
                             RRTForest *rrt_forest,
                             std::function<double(Motion*, Motion*)> distance_function)
{
    // Note: roots are col-major and each column stores one root
    size_t N = roots.cols();
    size_t W = roots.rows(); // W: Width of the state vector
    treeAdj_.setZero(N + TOTAL_INITIAL_TREE, N + TOTAL_INITIAL_TREE);
    treeDisjointSetId_.resize(N + TOTAL_INITIAL_TREE);

    forest_.clear();
    forest_.reserve(N + TOTAL_INITIAL_TREE);
    for (size_t i = 0; i < N + TOTAL_INITIAL_TREE; i++) {
        treeDisjointSetId_(i) = i;
        auto tree = std::make_shared<Tree>(int(i));
        tree->nn.reset(tools::SelfConfig::getDefaultNearestNeighbors<Motion *>(rrt_forest));
        tree->nn->setDistanceFunction(distance_function);
        forest_.emplace_back(std::move(tree));
    }

    // Setup roots for sample trees
    std::vector<double> sample(W);
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < W; j++)
            sample[j] = roots(j, i);
        auto *motion = new Motion(si_);
        si_->getStateSpace()->copyFromReals(motion->state, sample);
        motion->root = motion->state;

        auto& tree = forest_[i+TOTAL_INITIAL_TREE]; // Set forest_ from 2
        tree->nn->add(motion);
    }

    // TODO: union initial roots

    treeConnections_.clear();
}

void RRTForest::Private::unionTrees(const InterTreeConnection& itc)
{
    auto root1 = findTreeSetId(itc.tree1);
    auto root2 = findTreeSetId(itc.tree2);
    if (root1 != root2)
        treeDisjointSetId_(root1) = root2;

    treeConnections_.emplace_back(itc);
    int itc_id = treeConnections_.size();

    treeAdj_(root1, root2) = itc_id;
    treeAdj_(root2, root1) = -itc_id;
}

int RRTForest::Private::findTreeSetId(TreeData tree)
{
    if (tree->id != treeDisjointSetId_(tree->id)) {
        int parent_id = treeDisjointSetId_(tree->id);
        treeDisjointSetId_(tree->id) = findTreeSetId(forest_[parent_id]);
    }
    return treeDisjointSetId_(tree->id);
}

void RRTForest::Private::freeForest()
{
    std::vector<Motion *> motions;
    for (auto &tree : forest_) {
        tree->nn->list(motions);
        for (auto &motion : motions) {
            if (motion->state != nullptr)
                si_->freeState(motion->state);
            delete motion;
        }
    }
}

void RRTForest::Private::clearForest()
{
    for (auto &tree : forest_)
        tree->nn->clear();
}


std::shared_ptr<PathGeometric>
RRTForest::Private::getSolutionPath()
{
    // Path of Trees
    Eigen::VectorXi distance_from_start = Eigen::VectorXi::Constant(treeDisjointSetId_.size(), -1);
    distance_from_start(START_TREE) = START_TREE;

    Eigen::VectorXi parent_tree = Eigen::VectorXi::Constant(treeDisjointSetId_.size(), -1);
    Eigen::VectorXi parent_itc = Eigen::VectorXi::Zero(treeDisjointSetId_.size());
    parent_tree(START_TREE) = START_TREE;

    std::queue<int> bfs_queue;
    bfs_queue.push(START_TREE);
    while (distance_from_start(GOAL_TREE) < 0) {
        int now = bfs_queue.front();
        bfs_queue.pop();
        int distance = distance_from_start(now);
        for (size_t i = 0; i < forest_.size(); i++) {
            if (treeAdj_(now, i) == 0)
                continue;
            if (distance_from_start(i) >= 0)
                continue;
            distance_from_start(i) = distance + 1;
            parent_tree(i) = now;
            parent_itc(i) = treeAdj_(now, i);
            bfs_queue.push(i);
        }
    }
    std::vector<int> path_tree_level;
    for (int now = GOAL_TREE; now != START_TREE; now = parent_tree(now))
        path_tree_level.push_back(now);
    std::reverse(path_tree_level.begin(), path_tree_level.end());

    auto path = std::make_shared<PathGeometric>(si_);
    // Note: START_TREE is missing in path_tree_level and needs special
    //       handling
    addToPath(path, forest_[START_TREE], parent_itc(path_tree_level[0]));
    for (size_t itc_now = 0; itc_now + 1 < path_tree_level.size(); itc_now++) {
        auto itc_next = itc_now + 1;
        addToPath(path, itc_now, itc_next);
    }
    addToPath(path, parent_itc(path_tree_level[GOAL_TREE]), forest_[GOAL_TREE]);

    return path;
}

void
RRTForest::Private::addToPath(std::shared_ptr<PathGeometric> path,
                              TreeData tree,
                              int itc_next)
{
    Motion *from, *to;
    if (itc_next > 0) {
        int index = itc_next - 1;
        assert(treeConnections_[index].tree1 == tree);
        from = to = treeConnections_[index].motion1;
    } else {
        int index = (-itc_next) - 1;
        assert(treeConnections_[index].tree2 == tree);
        from = to = treeConnections_[index].motion2;
    }
    while (from->parent != nullptr)
        from = from->parent;
    addPathLCA(path, from, to);
}

void
RRTForest::Private::addToPath(std::shared_ptr<PathGeometric> path,
                              int itc_prev,
                              int itc_next)
{
    Motion *from, *to;
    TreeData sc1, sc2;
    if (itc_prev > 0) {
        int index = itc_prev - 1;
        sc1 = treeConnections_[index].tree2;
        from = treeConnections_[index].motion2;
    } else {
        int index = (-itc_prev) - 1;
        sc1 = treeConnections_[index].tree1;
        from = treeConnections_[index].motion1;
    }

    if (itc_next > 0) {
        int index = itc_next - 1;
        sc2 = treeConnections_[index].tree1;
        to = treeConnections_[index].motion1;
    } else {
        int index = (-itc_next) - 1;
        sc2 = treeConnections_[index].tree2;
        to = treeConnections_[index].motion2;
    }

    assert(sc1 == sc2);

    addPathLCA(path, from, to);
}

void
RRTForest::Private::addToPath(std::shared_ptr<PathGeometric> path,
                              int itc_prev,
                              TreeData tree)
{
    Motion *from, *to;
    if (itc_prev > 0) {
        int index = itc_prev - 1;
        assert(treeConnections_[index].tree2 == tree);
        from = to = treeConnections_[index].motion2;
    } else {
        int index = (-itc_prev) - 1;
        assert(treeConnections_[index].tree1 == tree);
        from = to = treeConnections_[index].motion1;
    }
    while (to->parent != nullptr)
        to = to->parent;
    addPathLCA(path, from, to);
}

void
RRTForest::Private::addPathLCA(std::shared_ptr<PathGeometric> path,
                               RRTForest::Motion *node_from,
                               RRTForest::Motion *node_to)
{
    std::vector<Motion *> from_path;
    std::unordered_set<Motion *> from_set;
    for (auto now = node_from; now != nullptr; now = now->parent) {
        from_path.emplace_back(now);
        from_set.insert(now);
    }

    std::vector<Motion *> to_path;
    for (auto now = node_to; now != nullptr; now = now->parent) {
        to_path.emplace_back(now);
        if (from_set.find(now) != from_set.end()) // i.e. from_set.contains(now)
            break;
    }
    // Note: the common Motion is IN the to_path
    for (auto now : from_path) {
        if (now == to_path.back())
            break;
        path->append(now->state);
    }
    for (auto iter = to_path.rbegin(); iter != to_path.rend(); iter++) {
        path->append((*iter)->state);
    }
}

RRTForest::RRTForest(const base::SpaceInformationPtr &si) : base::Planner(si, "RRTForest")
{
    specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
    specs_.directed = true;

    Planner::declareParam<double>("range", this, &RRTForest::setRange, &RRTForest::getRange, "0.:1.:10000.");
    connectionPoint_ = std::make_pair<base::State *, base::State *>(nullptr, nullptr);

    d_.reset(new Private(si));
}

RRTForest::~RRTForest()
{
    freeMemory();
}

void RRTForest::setup()
{
    Planner::setup();
    tools::SelfConfig sc(si_, getName());
    sc.configurePlannerRange(maxDistance_);

#if 0
    if (!tStart_)
        tStart_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Motion *>(this));
    if (!tGoal_)
        tGoal_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Motion *>(this));
    tStart_->setDistanceFunction([this](const Motion *a, const Motion *b)
                                 {
                                     return distanceFunction(a, b);
                                 });
    tGoal_->setDistanceFunction([this](const Motion *a, const Motion *b)
                                {
                                    return distanceFunction(a, b);
                                });
#endif
    auto distance_function = [this](const Motion *a, const Motion *b)
                                 {
                                     return distanceFunction(a, b);
                                 };
    d_->initForest(additional_samples_, this, distance_function);
    tStart_ = d_->forest_[Private::START_TREE];
    tGoal_ = d_->forest_[Private::GOAL_TREE];
}

void RRTForest::setSamples(const Eigen::MatrixXd& samples)
{
    additional_samples_ = samples;
}

void RRTForest::freeMemory()
{
    d_->freeForest();
#if 0
    std::vector<Motion *> motions;

    if (tStart_)
    {
        tStart_->list(motions);
        for (auto &motion : motions)
        {
            if (motion->state != nullptr)
                si_->freeState(motion->state);
            delete motion;
        }
    }

    if (tGoal_)
    {
        tGoal_->list(motions);
        for (auto &motion : motions)
        {
            if (motion->state != nullptr)
                si_->freeState(motion->state);
            delete motion;
        }
    }
#endif
}

void RRTForest::clear()
{
    Planner::clear();
    sampler_.reset();
    freeMemory();
#if 0
    if (tStart_)
        tStart_->clear();
    if (tGoal_)
        tGoal_->clear();
#endif
    d_->clearForest();
    connectionPoint_ = std::make_pair<base::State *, base::State *>(nullptr, nullptr);
}

//
// Note: tgi.xmotion is the newly created motion added to the tree. Hence it
// is output, rather than the input as indicated by its position.
RRTForest::GrowState RRTForest::growTree(TreeData &tree_package,
                                         TreeGrowingInfo &tgi,
                                         Motion *rmotion)
{
    auto& tree = tree_package->nn;
    /* find closest state in the tree */
    Motion *nmotion = tree->nearest(rmotion);

    /* assume we can reach the state we go towards */
    bool reach = true;

    /* find state to add */
    base::State *dstate = rmotion->state;
    double d = si_->distance(nmotion->state, rmotion->state);
    if (d > maxDistance_)
    {
        si_->getStateSpace()->interpolate(nmotion->state, rmotion->state, maxDistance_ / d, tgi.xstate);
        dstate = tgi.xstate;
        reach = false;
    }
    // if we are in the start tree, we just check the motion like we normally do;
    // if we are in the goal tree, we need to check the motion in reverse, but checkMotion() assumes the first state it
    // receives as argument is valid,
    // so we check that one first
    //
    // TODO: Check if we really need to handle this special case in the paper.
    bool validMotion = tgi.start ?
                           si_->checkMotion(nmotion->state, dstate) :
                           si_->getStateValidityChecker()->isValid(dstate) && si_->checkMotion(dstate, nmotion->state);

    if (validMotion)
    {
        /* create a motion */
        auto *motion = new Motion(si_);
        si_->copyState(motion->state, dstate);
        motion->parent = nmotion;
        motion->root = nmotion->root;
        tgi.xmotion = motion;

        tree->add(motion);
        return reach ? REACHED : ADVANCED;
    }
    else
        return TRAPPED;
}

ompl::base::PlannerStatus
RRTForest::solve(const base::PlannerTerminationCondition &ptc)
{
    checkValidity();
    auto *goal = dynamic_cast<base::GoalSampleableRegion *>(pdef_->getGoal().get());

    if (goal == nullptr)
    {
        OMPL_ERROR("%s: Unknown type of goal", getName().c_str());
        return base::PlannerStatus::UNRECOGNIZED_GOAL_TYPE;
    }

    while (const base::State *st = pis_.nextStart())
    {
        auto *motion = new Motion(si_);
        si_->copyState(motion->state, st);
        motion->root = motion->state;
        tStart_->nn->add(motion);
    }

    if (tStart_->nn->size() == 0)
    {
        OMPL_ERROR("%s: Motion planning start tree could not be initialized!", getName().c_str());
        return base::PlannerStatus::INVALID_START;
    }

    if (!goal->couldSample())
    {
        OMPL_ERROR("%s: Insufficient states in sampleable goal region", getName().c_str());
        return base::PlannerStatus::INVALID_GOAL;
    }

    if (!sampler_)
        sampler_ = si_->allocStateSampler();

    OMPL_INFORM("%s: Starting planning with %d states already in datastructure", getName().c_str(),
                (int)(tStart_->nn->size() + tGoal_->nn->size()));

    std::vector<TreeGrowingInfo> tgis(d_->forest_.size());
    for (auto& tgi : tgis)
        tgi.xstate = si_->allocState();
    std::vector<GrowState> gss(d_->forest_.size());

    auto *rmotion = new Motion(si_);
    base::State *rstate = rmotion->state;
    std::cerr << "rmotion: " << rmotion << std::endl;
    std::cerr << "rstate: " << rstate << std::endl;
    bool startTree = true;
    bool solved = false;
    TreeGrowingInfo temp_tgi;
    temp_tgi.start = true;
    temp_tgi.xstate = si_->allocState();

    while (!ptc)
    {
#if 0
        // "Flip", we don't need this though
        TreeData &tree = startTree ? tStart_ : tGoal_;
        tgi.start = startTree;
        startTree = !startTree;
        TreeData &otherTree = startTree ? tStart_ : tGoal_;
#endif

        if (tGoal_->nn->size() == 0 || pis_.getSampledGoalsCount() < tGoal_->nn->size() / 2)
        {
                const base::State *st = tGoal_->nn->size() == 0 ? pis_.nextGoal(ptc) : pis_.nextGoal();
                if (st != nullptr)
                {
                        auto *motion = new Motion(si_);
                        si_->copyState(motion->state, st);
                        motion->root = motion->state;
                        tGoal_->nn->add(motion);
                }

                if (tGoal_->nn->size() == 0)
                {
                        OMPL_ERROR("%s: Unable to sample any valid states for goal tree", getName().c_str());
                        break;
                }
        }

        /* sample random state */
        sampler_->sampleUniform(rstate);

        /* grow every tree in the forest */
        for (size_t i = 0; i < d_->forest_.size(); i++) {
                auto &tree = d_->forest_[i];
                auto &tgi = tgis[i];
                tgi.start = true; // We do not treat goal tree specially.
                auto &gs = gss[i];
                gs = growTree(tree, tgi, rmotion);
        }
        std::cerr << __LINE__ << " rstate: " << rstate << std::endl;

        /* attempt to connect trees */

        bool merged = false;
        /* 1. Collect trees that connected to the new sample */
        for (size_t to = 0; to < d_->forest_.size(); to++) {
            auto gs = gss[to];
            if (gs == TRAPPED)
                continue; // No new sample, skipped
            // Setup the new added milestone as the target
            Motion *addedMotion = tgis[to].xmotion;
            if (gs != REACHED)
                si_->copyState(rstate, tgis[to].xstate);
            // Iterate through all other trees
            for (size_t from = 0; from < d_->forest_.size(); from++) {
                if (to == from)
                    continue;
                GrowState gsc = ADVANCED;
                // Try to connect the other tree to the new milestone in "to" tree
                auto otherTree = d_->forest_[from];
                while (gsc == ADVANCED)
                    gsc = growTree(otherTree, temp_tgi, rmotion);
                if (gsc != REACHED)
                    continue;
                Private::InterTreeConnection itc(d_->forest_[from], d_->forest_[to],
                                                 temp_tgi.xmotion, addedMotion);
                d_->unionTrees(itc);
                merged = true;
            }
        }

#if 0
        /* remember which motion was just added */
        Motion *addedMotion = tgi.xmotion;

        /* if reached, it means we used rstate directly, no need top copy again */
        if (gs != REACHED)
                si_->copyState(rstate, tgi.xstate);

        GrowState gsc = ADVANCED;
        tgi.start = startTree;
        while (gsc == ADVANCED)
                gsc = growTree(otherTree, tgi, rmotion);

        Motion *startMotion = startTree ? tgi.xmotion : addedMotion;
        Motion *goalMotion = startTree ? addedMotion : tgi.xmotion;

        /* if we connected the trees in a valid way (start and goal pair is valid)*/
#endif
        // Note: start and goal may connect indirectly. Hence we need to
        //       query the disjoint set each time when union occurs.
        if (merged && d_->findTreeSetId(tStart_) == d_->findTreeSetId(tGoal_)) {
            // TODO:
            // Build the solution tree is tedious, offload this to another
            // function.
            pdef_->addSolutionPath(d_->getSolutionPath(), false, 0.0, getName());
            solved = true;
            break;
#if 0
            // it must be the case that either the start tree or the goal tree has made some progress
            // so one of the parents is not nullptr. We go one step 'back' to avoid having a duplicate state
            // on the solution path
            if (startMotion->parent != nullptr)
                startMotion = startMotion->parent;
            else
                goalMotion = goalMotion->parent;

            connectionPoint_ = std::make_pair(startMotion->state, goalMotion->state);

            /* construct the solution path */
            Motion *solution = startMotion;
            std::vector<Motion *> mpath1;
            while (solution != nullptr)
            {
                mpath1.push_back(solution);
                solution = solution->parent;
            }

            solution = goalMotion;
            std::vector<Motion *> mpath2;
            while (solution != nullptr)
            {
                mpath2.push_back(solution);
                solution = solution->parent;
            }

            auto path(std::make_shared<PathGeometric>(si_));
            path->getStates().reserve(mpath1.size() + mpath2.size());
            for (int i = mpath1.size() - 1; i >= 0; --i)
                path->append(mpath1[i]->state);
            for (auto &i : mpath2)
                path->append(i->state);

            pdef_->addSolutionPath(path, false, 0.0, getName());
            solved = true;
            break;
#endif
        }
    }

    for (auto& tgi : tgis)
        si_->freeState(tgi.xstate);
    si_->freeState(rstate);
    si_->freeState(temp_tgi.xstate);
    delete rmotion;

    OMPL_INFORM("%s: Created %u states (%u start + %u goal)",
                getName().c_str(), tStart_->nn->size() + tGoal_->nn->size(),
                tStart_->nn->size(), tGoal_->nn->size());

    return solved ? base::PlannerStatus::EXACT_SOLUTION : base::PlannerStatus::TIMEOUT;
}

void RRTForest::getPlannerData(base::PlannerData &data) const
{
    Planner::getPlannerData(data);

    std::vector<Motion *> motions;
    for (auto &tree : d_->forest_) {
        tree->nn->list(motions);
        for (auto &motion : motions) {
            if (motion->parent == nullptr)
                data.addStartVertex(base::PlannerDataVertex(motion->state, 1));
            else
                data.addEdge(base::PlannerDataVertex(motion->parent->state, 1), base::PlannerDataVertex(motion->state, 1));
        }
        motions.clear();
    }
    for (const auto& itc: d_->treeConnections_) {
        data.addEdge(data.vertexIndex(itc.motion1->state),
                     data.vertexIndex(itc.motion2->state));
    }
#if 0
    if (tStart_)
        tStart_->list(motions);

    for (auto &motion : motions)
    {
        if (motion->parent == nullptr)
            data.addStartVertex(base::PlannerDataVertex(motion->state, 1));
        else
        {
            data.addEdge(base::PlannerDataVertex(motion->parent->state, 1), base::PlannerDataVertex(motion->state, 1));
        }
    }

    motions.clear();
    if (tGoal_)
        tGoal_->list(motions);

    for (auto &motion : motions)
    {
        if (motion->parent == nullptr)
            data.addGoalVertex(base::PlannerDataVertex(motion->state, 2));
        else
        {
            // The edges in the goal tree are reversed to be consistent with start tree
            data.addEdge(base::PlannerDataVertex(motion->state, 2), base::PlannerDataVertex(motion->parent->state, 2));
        }
    }

    // Add the edge connecting the two trees
    data.addEdge(data.vertexIndex(connectionPoint_.first), data.vertexIndex(connectionPoint_.second));
#endif
}

}
}

// vim: tabstop=4:shiftwidth=4:softtabstop=4:expandtab
