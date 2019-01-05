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

/* Author: Xinya Zhang. */

#include "ompl/base/ProxyStateSampler.h"
#include "ompl/base/SpaceInformation.h"

struct ProxyCache;
namespace ompl {
namespace base {

struct ProxyStateSampler::ProxyCache {
	using StateArray = std::vector<std::vector<double>>;
	std::vector<size_t> begins_;
	std::vector<StateArray> arrays_;
	size_t index_ = 0;
	const StateSpace *space_;

	ProxyCache(const StateSpace* space)
		: space_(space)
	{
	}

	void cache(size_t begin,
		   StateArray&& states)
	{
		begins_.emplace_back(begin);
		arrays_.emplace_back(states);
	}

	bool supply_sample(State *state)
	{
		bool ret = false;
		for (size_t i = 0; i < begins_.size(); i++) {
			// std::cerr << "Looking up the cache: " << index_ << " in [" << begins_[i] << ", " << begins_[i] + arrays_[i].size() << "]\n";
			if (index_ >= begins_[i] && index_ < begins_[i] + arrays_[i].size()) {
				size_t off = index_ - begins_[i];
				// std::cerr << "Supply sample from cached block " << i << " at offset " << off << std::endl;
				space_->copyFromReals(state, arrays_[i][off]);
				ret = true;
				break;
			}
		}
		index_++;
		return ret;
	}
};

ProxyStateSampler::ProxyStateSampler(const StateSpace *space,
                                     base::StateSamplerPtr real_sampler)
	: StateSampler(space), real_sampler_(real_sampler)
{
	cache_.reset(new ProxyCache(space));
}

// We need to put the destructor here because ProxyCache is not defined in the
// header.
ProxyStateSampler::~ProxyStateSampler()
{
}

void ProxyStateSampler::cacheState(size_t begin,
                                   std::vector<std::vector<double>> states)
{
	cache_->cache(begin, std::move(states));
}

void ProxyStateSampler::sampleUniform(State *state)
{
	if (cache_->supply_sample(state))
		return ;
	return real_sampler_->sampleUniform(state);
}

void ProxyStateSampler::sampleUniformNear(State *state, const State *near, double distance)
{
	if (cache_->supply_sample(state))
		return ;
	return real_sampler_->sampleUniformNear(state, near, distance);
}

void ProxyStateSampler::sampleGaussian(State *state, const State *mean, double stdDev)
{
	if (cache_->supply_sample(state))
		return ;
	return real_sampler_->sampleGaussian(state, mean, stdDev);
}

} // namespace base
} // namespace ompl
