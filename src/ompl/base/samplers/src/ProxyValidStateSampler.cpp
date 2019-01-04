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

#include "ompl/base/samplers/ProxyValidStateSampler.h"
#include "ompl/base/SpaceInformation.h"

struct ProxyCache;
namespace ompl {
namespace base {

struct ProxyValidStateSampler::ProxyCache {
	using StateArray = std::vector<std::vector<double>>;
	std::vector<size_t> begins_;
	std::vector<StateArray> arrays_;
	size_t index_ = 0;
	StateSpacePtr sspace_;

	ProxyCache(const SpaceInformation* si)
	{
		sspace_ = si->getStateSpace();
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
			if (index_ >= begins_[i] && index_ < begins_[i] + arrays_[i].size()) {
				size_t off = index_ - begins_[i];
				sspace_->copyFromReals(state, arrays_[i][off]);
				ret = true;
				break;
			}
		}
		index_++;
		return ret;
	}
};

ProxyValidStateSampler::ProxyValidStateSampler(const SpaceInformation *si,
                                               base::ValidStateSamplerPtr real_sampler)
	: ValidStateSampler(si), real_sampler_(real_sampler)
{
	name_ = "proxy_unifom";
	cache_.reset(new ProxyCache(si));
}

// We need to put the destructor here because ProxyCache is not defined in the
// header.
ProxyValidStateSampler::~ProxyValidStateSampler()
{
}

void ProxyValidStateSampler::cacheState(size_t begin,
                                        std::vector<std::vector<double>> states)
{
	cache_->cache(begin, std::move(states));
}

bool ProxyValidStateSampler::sample(State *state)
{
	if (cache_->supply_sample(state))
		return true;
	return real_sampler_->sample(state);
}

bool ProxyValidStateSampler::sampleNear(State *state, const State *near, const double distance)
{
	return real_sampler_->sampleNear(state, near, distance);
}

} // namespace base
} // namespace ompl
