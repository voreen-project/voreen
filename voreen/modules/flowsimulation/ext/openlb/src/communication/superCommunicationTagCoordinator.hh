/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef SUPER_COMMUNICATION_TAG_COORDINATOR_HH
#define SUPER_COMMUNICATION_TAG_COORDINATOR_HH

#include "superCommunicationTagCoordinator.h"

#include <vector>
#include <algorithm>

namespace olb {

#ifdef PARALLEL_MODE_MPI

template <typename T>
class SuperCommunicationTagCoordinator<T>::ChannelId {
private:
  const int _iC;
  const int _jC;

public:
  ChannelId(int iC, int jC):
    _iC(util::min(iC,jC)),
    _jC(util::max(iC,jC)) { }

  bool operator==(const ChannelId& rhs) const
  {
    return _iC == rhs._iC
        && _jC == rhs._jC;
  }

  bool operator<(const ChannelId& rhs) const
  {
    return  _iC  < rhs._iC
        || (_iC == rhs._iC && _jC < rhs._jC);
  }

};

template <typename T>
SuperCommunicationTagCoordinator<T>::SuperCommunicationTagCoordinator(LoadBalancer<T>& loadBalancer):
  _loadBalancer(loadBalancer) { }

template <typename T>
template <unsigned D>
void SuperCommunicationTagCoordinator<T>::coordinate(
  std::vector<BlockCommunicationNeighborhood<T,D>>& neighborhood)
{
  for (int iC = 0; iC < _loadBalancer.size(); ++iC) {
    neighborhood[iC].forRemoteNeighbors([&](int jC) {
      _tags[_loadBalancer.rank(jC)][{_loadBalancer.glob(iC),jC}] = -1;
    });
  }

  for (auto& [rank, tags] : _tags) {
    int i=0;
    for (auto tag=tags.begin(); tag != tags.end(); ++tag, ++i) {
      (*tag).second = i;
    }
  }
}

template <typename T>
int SuperCommunicationTagCoordinator<T>::get(int iC, int jC, int iGroup)
{
  auto& tags = _tags[_loadBalancer.rank(iC)];
  return iGroup*tags.size() + tags[{iC,jC}];
}

#endif // PARALLEL_MODE_MPI

}

#endif
