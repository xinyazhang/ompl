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

#ifndef OMPL_GEOMETRIC_PLANNERS_RERRT_RRT_
#define OMPL_GEOMETRIC_PLANNERS_RERRT_RRT_

#include "ompl/geometric/planners/PlannerIncludes.h"
#include "ompl/datastructures/NearestNeighbors.h"

namespace ompl
{

    namespace geometric
    {

        /**
           @anchor ReRRT
           @par Short description
	   Retraction based RRT moves generated samples q in the C-Obs onto
	   the contact space q', and pull q' to q'' on contact space, where
	   \|q''-q\| is local mimimum.
	   This is a simplified version of the Retraction-based RRT.
           @par External documentation
           [[PDF]](http://gamma.cs.unc.edu/RRRT/RRRT.pdf)
           [[more]](http://gamma.cs.unc.edu/RRRT/)
        */

        /** \brief Retraction-based RRT */
        class ReRRT : public base::Planner
        {
        public:

            /** \brief Constructor */
            ReRRT(const base::SpaceInformationPtr &si);

            virtual ~ReRRT();

            virtual void getPlannerData(base::PlannerData &data) const;

            virtual base::PlannerStatus solve(const base::PlannerTerminationCondition &ptc);

            virtual void clear();

            /** \brief Set the goal bias

                In the process of randomly selecting states in
                the state space to attempt to go towards, the
                algorithm may in fact choose the actual goal state, if
                it knows it, with some probability. This probability
                is a real number between 0.0 and 1.0; its value should
                usually be around 0.05 and should not be too large. It
                is probably a good idea to use the default value. */
            void setGoalBias(double goalBias)
            {
                goalBias_ = goalBias;
            }

            /** \brief Get the goal bias the planner is using */
            double getGoalBias() const
            {
                return goalBias_;
            }

            /** \brief Set the range the planner is supposed to use.

                This parameter greatly influences the runtime of the
                algorithm. It represents the maximum length of a
                motion to be added in the tree of motions. */
            void setRange(double distance)
            {
                maxDistance_ = distance;
            }

            /** \brief Get the range the planner is using */
            double getRange() const
            {
                return maxDistance_;
            }

            /** \brief Set a different nearest neighbors datastructure */
            template<template<typename T> class NN>
            void setNearestNeighbors()
            {
                nn_.reset(new NN<Motion*>());
            }

            virtual void setup();

	    void setStateInjection(size_t start_from, std::vector<std::vector<double>> samples);
	    void setKNearest(int K) { knearest_ = K; }

	    virtual void setSampleSet(const Eigen::Ref<const Eigen::MatrixXd> Q) override;
	    virtual void setSampleSetFlags(const Eigen::Ref<const Eigen::Matrix<uint32_t, -1, 1>> QF) override;
	    virtual void getSampleSetConnectivity(Eigen::SparseMatrix<int>& ) override;
	    virtual void getCompactGraph(Eigen::Matrix<int64_t, -1, 1>& nouveau_vertex_id,
	                                 Eigen::MatrixXd& nouveau_vertices,
	                                 Eigen::Matrix<int64_t, -1, 2>& edges) const;

            /** \brief Representation of a motion

                This only contains pointers to parent motions as we
                only need to go backwards in the tree. */
            class Motion
            {
            public:

                Motion() : state(nullptr), parent(nullptr)
                {
                }

                /** \brief Constructor that allocates memory for the state */
                Motion(const base::SpaceInformationPtr &si) : state(si->allocState()), parent(nullptr)
                {
                }

                ~Motion()
                {
                }

                /** \brief The state contained by the motion */
                base::State       *state;

                /** \brief The parent motion in the exploration tree */
                Motion            *parent;

		int64_t motion_index = -65535;
		int64_t forest_index = 0;
		bool is_nouveau = true;
		enum {
			MOTION_OF_START,
			MOTION_OF_GOAL,
			MOTION_OF_SAMPLE,
		} motion_type;
            };

            //
            // Return a writable pointer to nn_
            // This may destroy the planner
            // USE IT WITH CAUTION.
            //
            std::shared_ptr< NearestNeighbors<Motion*> > _accessNearestNeighbors()
            {
                return nn_;
            }

            double distanceFunction(const Motion *a, const Motion *b) const
            {
                return si_->distance(a->state, b->state);
            }

        protected:

            /** \brief Free the memory allocated by this planner */
            void freeMemory();

            /** \brief Compute distance between motions (actually distance between contained states) */
            /** \brief State sampler */
            base::StateSamplerPtr                          sampler_;

            /** \brief A nearest-neighbors datastructure containing the tree of motions */
            std::shared_ptr< NearestNeighbors<Motion*> > nn_;

            /** \brief The fraction of time the goal is picked as the state to expand towards (if such a state is available) */
            double                                         goalBias_;

            /** \brief The maximum length of a motion to be added to a tree */
            double                                         maxDistance_;

            /** \brief The random number generator */
            RNG                                            rng_;

            /** \brief The most recent goal motion.  Used for PlannerData computation */
            Motion                                         *lastGoalMotion_;

	    size_t sample_injection_;
	    std::vector<std::vector<double>> samples_to_inject_;
	    int knearest_;

	    Eigen::MatrixXd predefined_samples_;
	    std::vector<Eigen::Triplet<int>> connectivity_tup_;

	    size_t motion_of_start_size_;
	    Eigen::Matrix<uint32_t, -1, 1> pds_flags_;
        };

    }
}

#endif
