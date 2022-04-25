/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Jan E. Marquardt
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

#include <utilities/omath.h>

#ifndef ONGRIDCONTACTDETECTION_H
#define ONGRIDCONTACTDETECTION_H

namespace olb {

// Helps storing of 4 IDs (1 to 65534) in one uint64_t

namespace contact {
// Necessary to shift particle IDs, since 0 should be the default value.
void storePid(uint64_t& previous, uint16_t value)
{
  previous = (previous << 16) + value + 1;
}

// Index starts at 0
int getPid(const uint64_t& ids, uint16_t index)
{
  return ((ids >> (index * 16)) & 0xffff) - 1;
}
}

} // namespace olb

#endif /*ONGRIDCONTACTDETECTION_H*/
