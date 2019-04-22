/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2010, Rice University
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
*   * Neither the name of the Rice University nor the names of its
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

#include "ompl/base/Planner.h"
#include "ompl/util/Exception.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include <sstream>
#include <thread>
#include <utility>

ompl::base::Planner::Planner(SpaceInformationPtr si, std::string name)
  : si_(std::move(si)), pis_(this), name_(std::move(name)), setup_(false)
{
    if (!si_)
        throw Exception(name_, "Invalid space information instance for planner");
}

const ompl::base::PlannerSpecs &ompl::base::Planner::getSpecs() const
{
    return specs_;
}

const std::string &ompl::base::Planner::getName() const
{
    return name_;
}

void ompl::base::Planner::setName(const std::string &name)
{
    name_ = name;
}

const ompl::base::SpaceInformationPtr &ompl::base::Planner::getSpaceInformation() const
{
    return si_;
}

const ompl::base::ProblemDefinitionPtr &ompl::base::Planner::getProblemDefinition() const
{
    return pdef_;
}

void ompl::base::Planner::setProblemDefinition(const ProblemDefinitionPtr &pdef)
{
    pdef_ = pdef;
    pis_.update();
}

const ompl::base::PlannerInputStates &ompl::base::Planner::getPlannerInputStates() const
{
    return pis_;
}

void ompl::base::Planner::setup()
{
    if (!si_->isSetup())
    {
        OMPL_INFORM("%s: Space information setup was not yet called. Calling now.", getName().c_str());
        si_->setup();
    }

    if (setup_)
        OMPL_WARN("%s: Planner setup called multiple times", getName().c_str());
    else
        setup_ = true;
}

void ompl::base::Planner::checkValidity()
{
    if (!isSetup())
        setup();
    pis_.checkValidity();
}

bool ompl::base::Planner::isSetup() const
{
    return setup_;
}

void ompl::base::Planner::clear()
{
    pis_.clear();
    pis_.update();
}

void ompl::base::Planner::getPlannerData(PlannerData &data) const
{
    for (const auto &plannerProgressProperty : plannerProgressProperties_)
        data.properties[plannerProgressProperty.first] = plannerProgressProperty.second();
}

void ompl::base::Planner::addGraph(const Eigen::Ref<const Eigen::MatrixXd> V, const Eigen::Ref<const Eigen::SparseMatrix<uint8_t>> E)
{
    throw std::runtime_error("addGraph not implemented in Planner " + getName());
}

ompl::base::PlannerStatus ompl::base::Planner::solve(const PlannerTerminationConditionFn &ptc, double checkInterval)
{
    return solve(PlannerTerminationCondition(ptc, checkInterval));
}

ompl::base::PlannerStatus ompl::base::Planner::solve(double solveTime)
{
    if (solveTime < 1.0)
        return solve(timedPlannerTerminationCondition(solveTime));
    return solve(timedPlannerTerminationCondition(solveTime, std::min(solveTime / 100.0, 0.1)));
}

void ompl::base::Planner::printProperties(std::ostream &out) const
{
    out << "Planner " + getName() + " specs:" << std::endl;
    out << "Multithreaded:                 " << (getSpecs().multithreaded ? "Yes" : "No") << std::endl;
    out << "Reports approximate solutions: " << (getSpecs().approximateSolutions ? "Yes" : "No") << std::endl;
    out << "Can optimize solutions:        " << (getSpecs().optimizingPaths ? "Yes" : "No") << std::endl;
    out << "Aware of the following parameters:";
    std::vector<std::string> params;
    params_.getParamNames(params);
    for (auto &param : params)
        out << " " << param;
    out << std::endl;
}

void ompl::base::Planner::printSettings(std::ostream &out) const
{
    out << "Declared parameters for planner " << getName() << ":" << std::endl;
    params_.print(out);
}

void ompl::base::Planner::setSampleSet(const Eigen::Ref<const Eigen::MatrixXd> )
{
    throw std::runtime_error("setSampleSet not implemented in Planner " + getName());
}

void ompl::base::Planner::setSampleSetFlags(const Eigen::Ref<const Eigen::Matrix<uint32_t, -1, 1>> QF)
{
    throw std::runtime_error("setSampleSetFlags not implemented in Planner " + getName());
}


void ompl::base::Planner::getSampleSetConnectivity(Eigen::SparseMatrix<int>& )
{
    throw std::runtime_error("getSampleSetConnectivity not implemented in Planner " + getName());
}

void ompl::base::Planner::getCompactGraph(Eigen::Matrix<int64_t, -1, 1>& nouveau_vertex_id,
                                          Eigen::MatrixXd& nouveau_vertices,
                                          Eigen::Matrix<int64_t, -1, 2>& edges) const
{
    throw std::runtime_error("getCompactGraph not implemented in Planner " + getName());
}

ssize_t ompl::base::Planner::peekPlannerDataSize() const
{
    throw std::runtime_error("peekPlannerDataSize not implemented in Planner " + getName());
}

namespace {
	// Note: we MUST use double dash instead of single dash, otherwise
 	// negative real values will be considered as an option key
	const std::string key_prefix = "--";

	bool is_option_key(const std::string& key)
	{
		return key.compare(0, key_prefix.size(), key_prefix) == 0;

	}

	std::string sanitize_option_key(const std::string& key)
	{
		if (!is_option_key(key))
			return key;
		size_t off = key_prefix.size();
		while (off < key.size() and key[off] == '-')
			off += 1;
		return key.substr(off, key.size() - off);
	}
}

/*
 * Option vector: variable length arguments similary to long options
 */
void ompl::base::Planner::setOptionVector(const std::vector<std::string>& ovec)
{
	option_dict_.clear();
	OMPL_INFORM("Is '--retract' option key? %d", is_option_key("--retract"));
	for (const auto& s : ovec) {
		OMPL_INFORM("[OptionVector]: %s", s.c_str());
	}
	for (size_t i = 0; i < ovec.size(); i++) {
		// Search for option key
		while (i < ovec.size() && !is_option_key(ovec[i]))
			i++;
		if (i >= ovec.size())
			break;
		OMPL_INFORM("Current key location %d", i);
		std::string current_key = sanitize_option_key(ovec[i]);
		auto current_key_index = i;
		i++; // find next key
		while (i < ovec.size() && !is_option_key(ovec[i]))
			i++;
		OMPL_INFORM("Next key location %d", i);
		std::vector<std::string> current_args;
		for (auto j = current_key_index + 1; j < i; j++)
			current_args.emplace_back(ovec[j]);
		option_dict_[current_key] = current_args;
	}
	OMPL_INFORM("[OptionDict]:");
	for (const auto& kv : option_dict_) {
		std::stringstream ss;
		ss << "\t[" << kv.first << "]: ";
		for (const auto& s: kv.second) {
			ss << s;
		}
		OMPL_INFORM("%s", ss.str().c_str());
	}
}

bool ompl::base::Planner::getOptionBool(const std::string& key) const
{
	// Can be option_dict_.contains(key) in C++ 20
	return option_dict_.find(key) != option_dict_.end();
}

bool ompl::base::Planner::getOptionStr(const std::string& key, std::string& value) const
{
	auto iter = option_dict_.find(key);
	if (iter == option_dict_.end())
	    return false;
	if (iter->second.empty())
	    return false;
	value = iter->second[0];
	return true;
}

bool ompl::base::Planner::getOptionInt(const std::string& key, long& value) const
{
	auto iter = option_dict_.find(key);
	if (iter == option_dict_.end())
	    return false;
	if (iter->second.empty())
	    return false;
	value = std::stol(iter->second[0]);
	return true;
}

bool ompl::base::Planner::getOptionReal(const std::string& key, double& value) const
{
	auto iter = option_dict_.find(key);
	if (iter == option_dict_.end())
	    return false;
	if (iter->second.empty())
	    return false;
	value = std::stod(iter->second[0]);
	return true;
}


void ompl::base::PlannerInputStates::clear()
{
    if (tempState_ != nullptr)
    {
        si_->freeState(tempState_);
        tempState_ = nullptr;
    }
    addedStartStates_ = 0;
    sampledGoalsCount_ = 0;
    pdef_ = nullptr;
    si_ = nullptr;
}

void ompl::base::PlannerInputStates::restart()
{
    addedStartStates_ = 0;
    sampledGoalsCount_ = 0;
}

bool ompl::base::PlannerInputStates::update()
{
    if (planner_ == nullptr)
        throw Exception("No planner set for PlannerInputStates");
    return use(planner_->getProblemDefinition());
}

void ompl::base::PlannerInputStates::checkValidity() const
{
    std::string error;

    if (pdef_ == nullptr)
        error = "Problem definition not specified";
    else
    {
        if (pdef_->getStartStateCount() <= 0)
            error = "No start states specified";
        else if (!pdef_->getGoal())
            error = "No goal specified";
    }

    if (!error.empty())
    {
        if (planner_ != nullptr)
            throw Exception(planner_->getName(), error);
        else
            throw Exception(error);
    }
}

bool ompl::base::PlannerInputStates::use(const ProblemDefinitionPtr &pdef)
{
    if (pdef)
        return use(pdef.get());

    clear();
    return true;
}

bool ompl::base::PlannerInputStates::use(const ProblemDefinition *pdef)
{
    if (pdef_ != pdef)
    {
        clear();
        pdef_ = pdef;
        si_ = pdef->getSpaceInformation().get();
        return true;
    }
    return false;
}

const ompl::base::State *ompl::base::PlannerInputStates::nextStart()
{
    if (pdef_ == nullptr || si_ == nullptr)
    {
        std::string error = "Missing space information or problem definition";
        if (planner_ != nullptr)
            throw Exception(planner_->getName(), error);
        else
            throw Exception(error);
    }

    while (addedStartStates_ < pdef_->getStartStateCount())
    {
        const base::State *st = pdef_->getStartState(addedStartStates_);
        addedStartStates_++;
        bool bounds = si_->satisfiesBounds(st);
        bool valid = bounds ? si_->isValid(st) : false;
        if (bounds && valid)
            return st;

        OMPL_WARN("%s: Skipping invalid start state (invalid %s)",
                  planner_ ? planner_->getName().c_str() : "PlannerInputStates", bounds ? "state" : "bounds");
        std::stringstream ss;
        si_->printState(st, ss);
        OMPL_DEBUG("%s: Discarded start state %s", planner_ ? planner_->getName().c_str() : "PlannerInputStates",
                   ss.str().c_str());
    }
    return nullptr;
}

const ompl::base::State *ompl::base::PlannerInputStates::nextGoal()
{
    // This initialization is safe since we are in a non-const function anyway.
    static PlannerTerminationCondition ptc = plannerAlwaysTerminatingCondition();
    return nextGoal(ptc);
}

const ompl::base::State *ompl::base::PlannerInputStates::nextGoal(const PlannerTerminationCondition &ptc)
{
    if (pdef_ == nullptr || si_ == nullptr)
    {
        std::string error = "Missing space information or problem definition";
        if (planner_ != nullptr)
            throw Exception(planner_->getName(), error);
        else
            throw Exception(error);
    }

    const GoalSampleableRegion *goal =
        pdef_->getGoal()->hasType(GOAL_SAMPLEABLE_REGION) ? pdef_->getGoal()->as<GoalSampleableRegion>() : nullptr;

    if (goal != nullptr)
    {
        time::point start_wait;
        bool first = true;
        bool attempt = true;
        while (attempt)
        {
            attempt = false;

            if (sampledGoalsCount_ < goal->maxSampleCount() && goal->canSample())
            {
                if (tempState_ == nullptr)
                    tempState_ = si_->allocState();
                do
                {
                    goal->sampleGoal(tempState_);
                    sampledGoalsCount_++;
                    bool bounds = si_->satisfiesBounds(tempState_);
                    bool valid = bounds ? si_->isValid(tempState_) : false;
                    if (bounds && valid)
                    {
                        if (!first)  // if we waited, show how long
                        {
                            OMPL_DEBUG("%s: Waited %lf seconds for the first goal sample.",
                                       planner_ ? planner_->getName().c_str() : "PlannerInputStates",
                                       time::seconds(time::now() - start_wait));
                        }
                        return tempState_;
                    }

                    OMPL_WARN("%s: Skipping invalid goal state (invalid %s)",
                              planner_ ? planner_->getName().c_str() : "PlannerInputStates",
                              bounds ? "state" : "bounds");
                    std::stringstream ss;
                    si_->printState(tempState_, ss);
                    OMPL_DEBUG("%s: Discarded goal state:\n%s",
                               planner_ ? planner_->getName().c_str() : "PlannerInputStates", ss.str().c_str());
                } while (!ptc && sampledGoalsCount_ < goal->maxSampleCount() && goal->canSample());
            }
            if (goal->couldSample() && !ptc)
            {
                if (first)
                {
                    first = false;
                    start_wait = time::now();
                    OMPL_DEBUG("%s: Waiting for goal region samples ...",
                               planner_ ? planner_->getName().c_str() : "PlannerInputStates");
                }
                std::this_thread::sleep_for(time::seconds(0.01));
                attempt = !ptc;
            }
        }
    }

    return nullptr;
}

bool ompl::base::PlannerInputStates::haveMoreStartStates() const
{
    if (pdef_ != nullptr)
        return addedStartStates_ < pdef_->getStartStateCount();
    return false;
}

bool ompl::base::PlannerInputStates::haveMoreGoalStates() const
{
    if ((pdef_ != nullptr) && pdef_->getGoal())
        if (pdef_->getGoal()->hasType(GOAL_SAMPLEABLE_REGION))
            return sampledGoalsCount_ < pdef_->getGoal()->as<GoalSampleableRegion>()->maxSampleCount();
    return false;
}
