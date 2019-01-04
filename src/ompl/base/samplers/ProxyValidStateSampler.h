/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2019, University of Texas at Austin
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

/* Author: Xinya Zhang. Modified based upon UniformValidStateSampler.h */

#ifndef OMPL_BASE_SAMPLERS_PROXY_VALID_STATE_SAMPLER_
#define OMPL_BASE_SAMPLERS_PROXY_VALID_STATE_SAMPLER_

#include "ompl/base/ValidStateSampler.h"
#include "ompl/base/StateSampler.h"
#include <vector>

namespace ompl
{
    namespace base
    {
        /** \brief A proxy state sampler that only samples valid states. */
        class ProxyValidStateSampler : public ValidStateSampler
        {
        public:
            /** \brief Constructor */
            ProxyValidStateSampler(const SpaceInformation *si,
                                   base::ValidStateSamplerPtr real_sampler);

            ~ProxyValidStateSampler() override;

            void cacheState(size_t begin,
                            std::vector<std::vector<double>> states); // Call with std::move to save object copy
            bool sample(State *state) override;
            bool sampleNear(State *state, const State *near, double distance) override;

        protected:
            /** \brief The sampler to proxy */
            base::ValidStateSamplerPtr real_sampler_;

            /** \brief The cache */
            struct ProxyCache;
            std::shared_ptr<ProxyCache> cache_;
        };
    }
}

#endif
