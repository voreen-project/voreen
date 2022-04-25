/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2011, 2014 Mathias J. Krause, Simon Zimny
 *                2021 Clara Schragmann, Adrian Kummerlaender
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

/** \file
 * Representation of a block geometry -- generic implementation.
 */

#ifndef BLOCK_GEOMETRY_HH
#define BLOCK_GEOMETRY_HH

#include "geometry/blockGeometry.h"

#include "functors/lattice/indicator/blockIndicatorBaseF2D.h"
#include "functors/lattice/indicator/blockIndicatorBaseF3D.h"

#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"

namespace olb {


template<typename T, unsigned D>
BlockGeometry<T,D>::BlockGeometry(Cuboid<T,D>& cuboid, int padding, int iCglob)
  : BlockStructureD<D>(cuboid.getExtent(), padding),
    _data(this->getNcells()),
    _communicatable(_data),
    _cuboid(cuboid),
    _iCglob(iCglob),
    _statistics(this),
    clout(std::cout, ("BlockGeometry" + std::to_string(D) + "D"))
{
  addToStatisticsList(&_statistics.getStatisticsStatus());
}

template<typename T, unsigned D>
BlockGeometryStatistics<T,D>& BlockGeometry<T,D>::getStatistics(bool verbose)
{
  return _statistics;
}

template<typename T, unsigned D>
BlockGeometryStatistics<T,D> const& BlockGeometry<T,D>::getStatistics(bool verbose) const
{
  return _statistics;
}

template<typename T, unsigned D>
Vector<T,D> BlockGeometry<T,D>::getOrigin() const
{
  return _cuboid.getOrigin();
}

template<typename T, unsigned D>
T BlockGeometry<T,D>::getDeltaR() const
{
  return _cuboid.getDeltaR();
}

template<typename T, unsigned D>
int& BlockGeometry<T,D>::get(LatticeR<D>latticeR)
{
  resetStatistics();
  return _data[0][this->getCellId(latticeR)];
}

template<typename T, unsigned D>
int& BlockGeometry<T,D>::get(const int latticeR[D])
{
  return get(latticeR);
}

template<typename T, unsigned D>
int& BlockGeometry<T,D>::get(std::size_t iCell)
{
  resetStatistics();
  return _data[0][iCell];
}

template<typename T, unsigned D>
int BlockGeometry<T,D>::get(LatticeR<D> latticeR) const
{
  return _data[0][this->getCellId(latticeR)];
}

template<typename T, unsigned D>
int BlockGeometry<T,D>::get(const int latticeR[D]) const
{
  return get(latticeR);
}

template<typename T, unsigned D>
int BlockGeometry<T,D>::get(std::size_t iCell) const
{
  return _data[0][iCell];
}

template<typename T, unsigned D>
int BlockGeometry<T,D>::getMaterial(LatticeR<D> latticeR) const
{
  if (this->isInside(latticeR)) {
    return _data[0][this->getCellId(latticeR)];
  }
  else {
    return 0;
  }
}

template<typename T, unsigned D>
int BlockGeometry<T,D>::getMaterial(const int latticeR[D]) const
{
  LatticeR<D> loc{latticeR};
  if (this->isInside(loc)) {
    return _data[0][this->getCellId(loc)];
  }
  else {
    return 0;
  }
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::getPhysR(T physR[D], const int latticeR[D]) const
{
  LatticeR<D> loc{latticeR};
  getPhysR(physR, loc);
  return;
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::getPhysR(T physR[D], LatticeR<D> loc) const
{
  _cuboid.getPhysR(physR, loc);
  return;
}

template<typename T, unsigned D>
int const& BlockGeometry<T,D>::getIcGlob() const
{
  return _iCglob;
}

template<typename T, unsigned D>
Vector<int,D> BlockGeometry<T,D>::getExtent() const
{
  if constexpr (D == 3) {
    return Vector<int,3>(this->getNx(), this->getNy(), this->getNz());
  } else {
    return Vector<int,2>(this->getNx(), this->getNy());
  }
  __builtin_unreachable();
}

template<typename T, unsigned D>
int BlockGeometry<T,D>::clean(bool verbose, std::vector<int> bulkMaterials)
{
  using DESCRIPTOR = std::conditional_t<D==2,descriptors::D2Q9<>,descriptors::D3Q27<>>;
  int counter=0;
  bool toClean = true;
  this->forCoreSpatialLocations([&](LatticeR<D> latticeR) {
    // material not 0 and not in bulkMaterials
    if (get(latticeR) != 0 && (! util::isContained(bulkMaterials, get(latticeR))) ) {
      toClean = true;
      for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++){
        if ( util::isContained(bulkMaterials, get(latticeR + (descriptors::c<DESCRIPTOR>(iPop)))) ){
          toClean = false;
        }
      }
      if (toClean){
        get(latticeR) = 0;
        counter++;
      }
    }
  });
  if (verbose) {
    clout << "cleaned "<< counter << " outer boundary voxel(s)" << std::endl;
  }
  return counter;
}

template<typename T, unsigned D>
int BlockGeometry<T,D>::outerClean(bool verbose, std::vector<int> bulkMaterials)
{
  int counter=0;
  using DESCRIPTOR = std::conditional_t<D==2,descriptors::D2Q9<>,descriptors::D3Q27<>>;
  bool toClean = false;
  this->forCoreSpatialLocations([&](LatticeR<D> latticeR) {
    if (util::isContained(bulkMaterials, get(latticeR))) {
      toClean = false;
      for(int iPop = 1; iPop < DESCRIPTOR::q; iPop++){
        if(getMaterial(latticeR + (descriptors::c<DESCRIPTOR>(iPop))) == 0){
          toClean = true;
        }
      }
      if(toClean){
        get(latticeR) = 0;
        counter++;
      }
    }
  });
  if (verbose) {
    clout << "cleaned "<< counter << " outer fluid voxel(s)" << std::endl;
  }
  return counter;
}

template<typename T, unsigned D>
int BlockGeometry<T,D>::innerClean(bool verbose)
{
  int count = 0;
  int count2 = 0;
  using DESCRIPTOR = std::conditional_t<D==2,descriptors::D2Q5<>,descriptors::D3Q7<>>;

  this->forCoreSpatialLocations([&](LatticeR<D> latticeR) {
    if (get(latticeR) != 1 && get(latticeR) != 0) {
      count++;
      if constexpr (D==3){
        bool var[7] = {false};
        for(int iPop = 1; iPop < DESCRIPTOR::q; iPop++){
          if(getMaterial(latticeR + (descriptors::c<DESCRIPTOR>(iPop)))==1){
            var[iPop] = true;
          }
        }
        int comb[12][3]={{1,4,2},{1,4,3},{1,4,5},{1,4,6},{2,5,1},{2,5,3},{2,5,4},{2,5,6},{3,6,1},{3,6,2},{3,6,4},{3,6,5}};
        for(int i = 0; i < 12; i++){
          if (var[(comb[i][0])] == true
              && var[(comb[i][1])] == true
              && var[(comb[i][2])] == true){
            get(latticeR) = 1;
            count2++;//count2 is the same value as before, count2 gets increased even if this cell has already been cleaned
          }
        }
     } else {
        int var = 0;
        for(int iPop = 1; iPop < DESCRIPTOR::q; iPop++){
          if(getMaterial(latticeR + (descriptors::c<DESCRIPTOR>(iPop)))==1){
            var = var+1;
          }
        }
        if(var >= 3){
          get(latticeR) = 1;
          count2++;//count2 differs from original count2
        }
      }
    }
  });

  if (verbose) {
    clout << "cleaned "<< count2 << " inner boundary voxel(s)" << std::endl;
  }
  return count2;
}

template<typename T, unsigned D>
int BlockGeometry<T,D>::innerClean(int fromM, bool verbose)
{
  int count = 0;
  int count2 = 0;
  using DESCRIPTOR = std::conditional_t<D==2,descriptors::D2Q5<>,descriptors::D3Q7<>>;

  this->forCoreSpatialLocations([&](LatticeR<D> latticeR) {
    if (get(latticeR) != 1 && get(latticeR) != 0 && get(latticeR) == fromM) {
      count++;
      if constexpr (D==3){
        bool var[7] = {false};
        for(int iPop = 1; iPop < DESCRIPTOR::q; iPop++){
          if(getMaterial(latticeR + (descriptors::c<DESCRIPTOR>(iPop)))==1){
            var[iPop] = true;
          }
        }
        int comb[12][3]={{1,4,2},{1,4,3},{1,4,5},{1,4,6},{2,5,1},{2,5,3},{2,5,4},{2,5,6},{3,6,1},{3,6,2},{3,6,4},{3,6,5}};
        for(int i = 0; i < 12; i++){
          if (var[(comb[i][0])] == true
              && var[(comb[i][1])] == true
              && var[(comb[i][2])] == true){
            get(latticeR) = 1;
            count2++;//count2 is the same value as before, count2 gets increased even if this cell has already been cleaned
          }
        }
     } else {
        int var = 0;
        for(int iPop = 1; iPop < DESCRIPTOR::q; iPop++){
          if(getMaterial(latticeR + (descriptors::c<DESCRIPTOR>(iPop)))==1){
            var = var+1;
          }
        }
        if(var >= 3){
          get(latticeR) = 1;
          count2++;//count2 differs from original count2
        }
      }
    }
  });

  if (verbose){
    clout << "cleaned "<< count2
    << " inner boundary voxel(s) of Type " << fromM << std::endl;
  }
  return count2;
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::reset(IndicatorF<T,D>& domain)
{
  this->forCoreSpatialLocations([&](LatticeR<D> latticeR) {
    T physR[D] { };
    bool output{};
    getPhysR(physR, latticeR);
    domain(&output, physR);
    if (output) {
      get(latticeR) = 0;
    }
  });
}

template<typename T, unsigned D>
bool BlockGeometry<T,D>::find(int material, std::vector<unsigned> offset,
                                std::vector<int> var)
{
  bool found = false;
  for (var[0] = 0; var[0] < this->getNx(); var[0]++) {
    for (var[1] = 0; var[1] < this->getNy(); var[1]++) {
      if constexpr(D==3){
        for (var[2] = 0; var[2] < this->getNz(); var[2]++) {
          found = check(material, var, offset);
        }
      } else {
        found = check(material, var, offset);
      }
      if (found) {
        return found;
      }
    }
  }
  return found;
}

template<typename T, unsigned D>
bool BlockGeometry<T,D>::check(int material, std::vector<int> var,
                                 std::vector<unsigned> offset)
{
  bool found = true;
  for (int iOffsetX = -offset[0]; iOffsetX <= (int) offset[0]; ++iOffsetX) {
    for (int iOffsetY = -offset[1]; iOffsetY <= (int) offset[1]; ++iOffsetY) {
      if constexpr (D==3){
        for (int iOffsetZ = -offset[2]; iOffsetZ <= (int) offset[2]; ++iOffsetZ) {
          if (getMaterial({var[0] + iOffsetX, var[1] + iOffsetY, var[2] + iOffsetZ}) != material) {
            found = false;
          }
        }
      } else {
        if (getMaterial({var[0] + iOffsetX, var[1] + iOffsetY}) != material) {
          found = false;
        }
      }
    }
  }
  return found;
}

template<typename T, unsigned D>
bool BlockGeometry<T,D>::checkForErrors(bool verbose) const
{
  bool error = false;
  using DESCRIPTOR = std::conditional_t<D==2,descriptors::D2Q9<>,descriptors::D3Q27<>>;
  bool errorFound = false;
  this->forSpatialLocations([&](LatticeR<D> latticeR)
  {
    if (get(latticeR) == 0) {
          errorFound = false;
          for(int iPop = 1; iPop < DESCRIPTOR::q; iPop++){
            if(getMaterial(latticeR + (descriptors::c<DESCRIPTOR>(iPop))) == 1){
              errorFound = true;
            }
          }
      if(errorFound){
        error = true;
      }
    }
  });

  if (verbose) {
    if (error) {
      clout << "error!" << std::endl;
    }
    else {
      clout << "the model is correct!" << std::endl;
    }
  }
  return error;
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::rename(int fromM, int toM)
{
  this->forCoreSpatialLocations([&](LatticeR<D> latticeR) {
    if (get(latticeR) == fromM) {
      get(latticeR) = toM;
    }
  });
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::rename(int fromM, int toM, IndicatorF<T,D>& condition)
{
  T physR[D];
  this->forSpatialLocations([&](LatticeR<D> latticeR) {
    if (get(latticeR) == fromM) {
      getPhysR(physR, latticeR);
      bool inside[1];
      condition(inside, physR);
      if (inside[0]) {
        get(latticeR) = toM;
      }
    }
  });
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::rename(int fromM, int toM, LatticeR<D> offset)
{
  this->forCoreSpatialLocations([&](LatticeR<D> latticeR) {
    if (get(latticeR) == fromM) {
      bool found = true;
      for (int iOffsetX = -offset[0]; iOffsetX <= (int) offset[0]; ++iOffsetX) {
        for (int iOffsetY = -offset[1]; iOffsetY <= (int) offset[1]; ++iOffsetY) {
          if constexpr (D == 3) {
            for (int iOffsetZ = -offset[2]; iOffsetZ <= (int) offset[2]; ++iOffsetZ) {
              if (getMaterial({latticeR[0] + iOffsetX, latticeR[1] + iOffsetY, latticeR[2] + iOffsetZ}) != fromM) {
                if (getMaterial({latticeR[0] + iOffsetX, latticeR[1] + iOffsetY, latticeR[2] + iOffsetZ}) != 1245) {
                  found = false;
                }
              }
            }
          } else {
            if (getMaterial({latticeR[0] + iOffsetX, latticeR[1] + iOffsetY}) != fromM) {
              if (getMaterial({latticeR[0] + iOffsetX, latticeR[1] + iOffsetY}) != 1245) {
                found = false;
              }
            }
          }
        }
      }
      if (found) {
        get(latticeR) = 1245;
      }
    }
  });
  rename(1245,toM);
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::rename(int fromM, int toM, int testM,
    std::vector<int> testDirection)
{
  this->forCoreSpatialLocations([&](LatticeR<D> latticeR) {
    if (get(latticeR) == fromM) {
      // flag that indicates the renaming of the current voxel, valid voxels are not renamed
      bool isValid = true;
      for (int iOffsetX = util::min(testDirection[0],0); iOffsetX <= util::max(testDirection[0],0); ++iOffsetX) {
        for (int iOffsetY = util::min(testDirection[1],0); iOffsetY <= util::max(testDirection[1],0); ++iOffsetY) {
          if constexpr (D == 3){
            for (int iOffsetZ = util::min(testDirection[2],0); iOffsetZ <= util::max(testDirection[2],0); ++iOffsetZ) {
              if (iOffsetX!=0 || iOffsetY!=0 || iOffsetZ!=0) {
                if (getMaterial({latticeR[0] + iOffsetX, latticeR[1] + iOffsetY, latticeR[2] + iOffsetZ}) != testM) {
                  isValid = false;
                }
              }
            }
          } else {
            if (iOffsetX!=0 || iOffsetY!=0) {
              if (getMaterial({latticeR[0] + iOffsetX, latticeR[1] + iOffsetY}) != testM) {
                isValid = false;
              }
            }
          }
        }
      }
      if (!isValid) {
        get(latticeR) = toM;
      }
    }
  });
}


template<typename T, unsigned D>
void BlockGeometry<T,D>::rename(int fromM, int toM, int fluidM,
    IndicatorF<T,D>& condition, Vector<int,D> discreteNormal)
{
  rename(fromM, toM, condition);
  Vector<int,D> testDirection(discreteNormal);
  T physR[D];
  this->forCoreSpatialLocations([&](LatticeR<D> latticeR) {
    if (get(latticeR) == toM) {
      getPhysR(physR, latticeR);
      bool inside[1];
      condition(inside, physR);
      if (inside[0]) {
        if constexpr (D==3){
          if (getMaterial({latticeR[0]+testDirection[0],latticeR[1]+testDirection[1],latticeR[2]+testDirection[2]})!=fluidM ||
              getMaterial({latticeR[0]+2*testDirection[0],latticeR[1]+2*testDirection[1],latticeR[2]+2*testDirection[2]})!=fluidM ||
              getMaterial({latticeR[0]-testDirection[0],latticeR[1]-testDirection[1],latticeR[2]-testDirection[2]})!=0 ) {
            get(latticeR) = fromM;
          }
        } else {
          if (getMaterial({latticeR[0]+testDirection[0],latticeR[1]+testDirection[1]}) != fluidM ||
              getMaterial({latticeR[0]+2*testDirection[0],latticeR[1]+2*testDirection[1]}) != fluidM ||
              getMaterial({latticeR[0]-testDirection[0],latticeR[1]-testDirection[1]}) != 0) {
            get(latticeR) = fromM;
          }
        }
      }
    }
  });
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::rename(int fromM, int toM, int fluidM,
    IndicatorF<T,D>& condition)
{
  rename(fromM, toM, condition);
  std::vector<int> testDirection = getStatistics().computeDiscreteNormal(toM);
  T physR[D];
  // values that have been incorrectly changed from "fromM to "toM" get assigned back to "fromM"
  this->forCoreSpatialLocations([&](LatticeR<D> latticeR) {
    if (get(latticeR) == toM) {
      getPhysR(physR, latticeR);
      bool inside[1];
      condition(inside, physR);
      if (inside[0]) {
        if constexpr(D==3){
          if (getMaterial({latticeR[0]+testDirection[0],latticeR[1]+testDirection[1],latticeR[2]+testDirection[2]})!=fluidM ||
              getMaterial({latticeR[0]+2*testDirection[0],latticeR[1]+2*testDirection[1],latticeR[2]+2*testDirection[2]})!=fluidM ||
              getMaterial({latticeR[0]-testDirection[0],latticeR[1]-testDirection[1],latticeR[2]-testDirection[2]})!=0 ) {
            get(latticeR) = fromM;
          }
        } else {
          if (getMaterial({latticeR[0]+testDirection[0],latticeR[1]+testDirection[1]})!=fluidM ||
              getMaterial({latticeR[0]+2*testDirection[0],latticeR[1]+2*testDirection[1]})!=fluidM ||
              getMaterial({latticeR[0]-testDirection[0],latticeR[1]-testDirection[1]})!=0 ) {
            get(latticeR) = fromM;
          }
        }
      }
    }
  });
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::copyMaterialLayer(IndicatorF3D<T>& condition, int discreteNormal[D], int numberOfLayers)
{
  T physR[D];
  this->forCoreSpatialLocations([&](LatticeR<D> latticeR) {
    getPhysR(physR, latticeR);
    bool inside[1];
    condition(inside, physR);
    if (inside[0]) {
      for (int i = 0; i < numberOfLayers; i++) {
        if (0 <= latticeR[0] + i * discreteNormal[0] && latticeR[0] + i * discreteNormal[0] < this->getNx() &&
            0 <= latticeR[1] + i * discreteNormal[1] && latticeR[1] + i * discreteNormal[1] < this->getNy()){
          if constexpr(D==3){
            if(0 <= latticeR[2] + i * discreteNormal[2] && latticeR[2] + i * discreteNormal[2] < this->getNz()){
              get({latticeR[0] + i * discreteNormal[0], latticeR[1] + i * discreteNormal[1], latticeR[2] + i * discreteNormal[2]}) = get(latticeR);
            }
          } else {
            get({latticeR[0] + i * discreteNormal[0], latticeR[1] + i * discreteNormal[1]}) = get(latticeR);
          }
        }
      }
    }
  });
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::regionGrowing(int fromM, int toM, LatticeR<D> seed,
                                       std::vector<int> offset, std::map<std::vector<int>, int>* tmp)
{
  std::map<std::vector<int>, int> tmp2;
  bool firstCall = false;
  if (tmp == nullptr) {
    tmp = &tmp2;
    firstCall = true;
  }

  if (getMaterial(seed) == fromM) {
    std::vector<int> found;
    found.push_back(seed[0]);
    found.push_back(seed[1]);
    if constexpr(D==3){
      found.push_back(seed[2]);
      if (tmp->count(found) == 0) {
        (*tmp)[found] = 2;
        if (offset[0] != 0) {
          regionGrowing(fromM, toM, {seed[0] + 1, seed[1], seed[2]}, offset, tmp);
          regionGrowing(fromM, toM, {seed[0] - 1, seed[1], seed[2]}, offset, tmp);
        }
        if (offset[1] != 0) {
          regionGrowing(fromM, toM, {seed[0], seed[1] + 1, seed[2]}, offset, tmp);
          regionGrowing(fromM, toM, {seed[0], seed[1] - 1, seed[2]}, offset, tmp);
        }
        if (offset[2] != 0) {
          regionGrowing(fromM, toM, {seed[0], seed[1], seed[2] + 1}, offset, tmp);
          regionGrowing(fromM, toM, {seed[0], seed[1], seed[2] - 1}, offset, tmp);
        }
      }
    } else {
      if (tmp->count(found) == 0) {
        (*tmp)[found] = 2;
        if (offset[0] != 0) {
          regionGrowing(fromM, toM, {seed[0] + 1, seed[1]}, offset, tmp);
          regionGrowing(fromM, toM, {seed[0] - 1, seed[1]}, offset, tmp);
        }
        if (offset[1] != 0) {
          regionGrowing(fromM, toM, {seed[0], seed[1] + 1}, offset, tmp);
          regionGrowing(fromM, toM, {seed[0], seed[1] - 1}, offset, tmp);
        }
      }
    }
  }
  if (firstCall) {
    std::map<std::vector<int>, int>::iterator iter;
    for (iter = tmp->begin(); iter != tmp->end(); iter++) {
      if constexpr(D==3){
        get((iter->first)[0],(iter->first)[1],(iter->first)[2]) = toM;
      } else {
        get((iter->first)[0],(iter->first)[1]) = toM;
      }
    }
  }
  return;
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::addToStatisticsList(bool* statisticStatus)
{
  _statisticsUpdateNeeded.push_back(statisticStatus);
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::removeFromStatisticsList(bool* statisticStatus)
{
  _statisticsUpdateNeeded.remove(statisticStatus);
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::printLayer(std::vector<int> min, std::vector<int> max, bool linenumber)//each call of this function must be carefully changed(previous:(xo,x1,y0,y1),now:({x0,yo},{x1,z1}))
{
  for (int x = min[0]; x <= max[0]; x++) {
    if (linenumber) {
      clout << x << ": ";
    }
    for (int y = min[1]; y <= max[1]; y++) {
      if constexpr (D==3){
        for (int z = min[2]; z <= max[2]; z++) {
          clout << getMaterial({x, y, z}) << " ";
        }
        if (max[1] - min[1] != 0 && max[2] - min[2] != 0) {
          clout << std::endl;
        }
      } else {
        clout << getMaterial({x, y}) << " ";
      }
    }
    if (max[0] - min[0] != 0) {
      clout << std::endl;
    }
  }
  clout << std::endl;
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::printLayer(int direction, int layer, bool linenumber)
{
  assert(direction >= 0 && direction <= 2);
  if constexpr(D==3){
    switch (direction) {
    case 0:
      printLayer({layer, 0, 0},{layer, this->getNy() - 1, this->getNz() - 1}, linenumber);
      break;
    case 1:
      printLayer({0, layer, 0}, {this->getNx() - 1, layer, this->getNz() - 1}, linenumber);
      break;
    case 2:
      printLayer({0, 0, layer}, {this->getNx() - 1, this->getNy() - 1, layer}, linenumber);
      break;
    }
  } else {
    switch (direction) {
    case 0:
      printLayer({layer, 0}, {layer, this->getNy() - 1}, linenumber);
      break;
    case 1:
      printLayer({0, layer}, {this->getNx() - 1, layer}, linenumber);
      break;
    }
  }
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::printNode(std::vector<int> loc)
{
  for (int x = loc[0] - 1; x <= loc[0] + 1; x++) {
    clout << "x=" << x << std::endl;
    for (int y = loc[1] - 1; y <= loc[1] + 1; y++) {
      if constexpr (D==3){
        for (int z = loc[2] - 1; z <= loc[2] + 1; z++) {
          clout << getMaterial({x, y, z}) << " ";
        }
        clout << std::endl;
      } else {
        clout << getMaterial({x, y}) << " ";
      }
    }
    clout << std::endl;
  }
  clout << std::endl;
}

template<typename T, unsigned D>
void BlockGeometry<T,D>::resetStatistics()
{
  for (bool* update : _statisticsUpdateNeeded) {
    *update = true;
  }
}

} // namespace olb

#endif
