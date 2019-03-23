/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007-2014 Mathias J. Krause
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
 * The description of a 2D super lattice -- header file.
 */

#ifndef SUPER_LATTICE_2D_H
#define SUPER_LATTICE_2D_H

#include <vector>

#include "blockLattice2D.h"
#include "blockLatticeView2D.h"
#include "communication/communicator2D.h"
#include "postProcessing.h"
#include "serializer.h"
#include "communication/superStructure2D.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T> class Communicator2D;
template<typename T> class CuboidGeometry2D;
template<typename T, template<typename U> class Lattice> class BlockLattice2D;
template<typename T, template<typename U> class Lattice> class BlockLatticeView2D;
template<typename T> class LoadBalancer;
template<typename T> class SuperGeometry2D;
template<typename T, template<typename U> class Lattice> class SuperLatticeF2D;
template<typename T> class SuperStructure2D;
template<typename T> class SuperIndicatorF2D;


/// A super lattice combines a number of block lattices that are ordered
/// in a cuboid geometry.
/** The communication between the block lattices is done by two
 * communicators. One (_commStream) is responible to provide the data for
 * the streaming the other (_commBC) for the non-local boundary conditions.
 * To simplify the code structure ghost cells in an overlap of size
 * (_overlap) is indrocuced. It depends on the non-locality of the
 * boundary conditions but is at least one because of the streaming
 *
 * The algorithm is parallelized with mpi. The load balancer (_load)
 * distributes the block lattices to processes.
 *
 * WARNING: For unstructured grids there is an iterpolation needed
 * for the method buffer_outData in coboidNeighbourhood which is not
 * yet implemented! Moreover this class needs to be chanced
 * that the number of times steps for the collision and streaming is
 * is dependent of the refinement level.
 *
 * This class is not intended to be derived from.
 */
template<typename T, template<typename U> class Lattice>
class SuperLattice2D : public SuperStructure2D<T>, public BufferSerializable {

private:
  /// Lattices with ghost cell layer of size overlap
  std::vector<BlockLattice2D<T,Lattice> >     _extendedBlockLattices;
  /// View of the lattices without overlap
  std::vector<BlockLatticeView2D<T,Lattice> > _blockLattices;

  /// Size of the refinement overlap
  int                                         _overlapRefinement;
  /// This communicator handels the communication for the streaming
  Communicator2D<T>                           _commStream;
  /// This communicator handels the communication for the postprocessors
  Communicator2D<T>                           _commBC;
  /// Specifies if there is communication for non local boundary conditions
  /// needed. It is automatically swichted on if overlapBC >= 1 by the
  /// calling the constructer. (default =false)
  bool                                        _commBC_on;

  /// Statistic of the super structure
  LatticeStatistics<T>                        *_statistics;
  /// Specifies if there is statistic calculated. It is always
  /// needed for the ConstRhoBGK dynamics. (default =true)
  bool                                        _statistics_on;
public:
  /// Construction of a super lattice
  SuperLattice2D(SuperGeometry2D<T>& superGeometry, int overlapRefinement=0);
  SuperLattice2D(const SuperLattice2D&) = delete;
  ~SuperLattice2D();
  /// Read and write access to a block lattice
  BlockLattice2D<T,Lattice>& getExtendedBlockLattice(int locIC)
  {
    return _extendedBlockLattices[locIC];
  };
  /// Read only access to a block lattice
  BlockLattice2D<T,Lattice> const& getExtendedBlockLattice(int locIC) const
  {
    return _extendedBlockLattices[locIC];
  };
  /// Read and write access to a lattice (block lattice view, one
  /// without overlap).
  BlockLatticeView2D<T,Lattice>& getBlockLattice(int locIC)
  {
    return _blockLattices[locIC];
  };
  /// Read only access to a lattice
  BlockLatticeView2D<T,Lattice> const& getBlockLattice(int locIC) const
  {
    return _blockLattices[locIC];
  };

  /// Read and write access to the streaming communicator
  Communicator2D<T>& get_commStream()
  {
    return _commStream;
  };
  /// Read only access to the streaming communicator
  Communicator2D<T> const& get_commStream() const
  {
    return _commStream;
  };
  /// Read and write access to the boundary communicator
  Communicator2D<T>& get_commBC()
  {
    return _commBC;
  };
  /// Read only access to the boundary communicator
  Communicator2D<T> const& get_commBC() const
  {
    return _commBC;
  };

  /// Return a handle to the LatticeStatistics object
  LatticeStatistics<T>& getStatistics();
  /// Return a constant handle to the LatticeStatistics object
  LatticeStatistics<T> const& getStatistics() const;

  /// Write access to lattice cells that returns false if
  /// iX/iY is not in any of the cuboids
  bool set(T iX, T iY, Cell<T,Lattice> const& cell);
  /// Read only access to lattice cells that returns false if
  /// iX/iY is not in any of the cuboids
  bool get(T iX, T iY, Cell<T,Lattice>& cell) const;
  /// Read only access to lattice cells over the cuboid number
  /// and local coordinates   WARNING!!! NO ERROR HANDLING IMPLEMENTED!!!
  Cell<T,Lattice> get(int iC, int locX, int locY) const;

  /// Write access to the memory of the data of the super structure
  bool* operator() (int iCloc, int iX, int iY, int iData) override
  {
    return (bool*)&getExtendedBlockLattice(iCloc).get(iX+this->_overlap, iY+this->_overlap)[iData];
  };
  /// Read only access to the dim of the data of the super structure
  int getDataSize() const override
  {
    return Lattice<T>::q;
  };
  /// Read only access to the data type dim of the data of the super structure
  int getDataTypeSize() const override
  {
    return sizeof(T);
  };
  /// Initialize all lattice cells to become ready for simulation
  void initialize();

  /// Defines the dynamics by material
  void defineDynamics(SuperGeometry2D<T>& superGeometry, int material, Dynamics<T,Lattice>* dynamics);
  /// Defines the dynamics on a domain described by an indicator reference
  void defineDynamics(SuperIndicatorF2D<T>& indicator, Dynamics<T,Lattice>* dynamics);
  /// Defines the dynamics on a domain described by an indicator pointer
  void defineDynamics(std::unique_ptr<SuperIndicatorF2D<T>> const& indicator, Dynamics<T, Lattice>* dynamics);
  /// Defines rho on a rectangular domain
  void defineRhoU (T x0, T x1, T y0, T y1, T rho, const T u[Lattice<T>::d] );
  /// Defines rho and u on a domain with a particular material number
  void defineRhoU(SuperGeometry2D<T>& sGeometry, int material,
                  AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u);
  /// Defines rho on a rectangular domain
  void defineRho (T x0, T x1, T y0, T y1, T rho );
  /// Defines rho on a domain with a particular material number
  void defineRho(SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& rho);
  /// Defines u on a rectangular domain
  void defineU (T x0, T x1, T y0, T y1, const T u[Lattice<T>::d] );
  /// Defines u on a domain with a particular material number
  void defineU(SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& u);
  // Defines a Population on a domain with a particular material number
  void definePopulations(SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& Pop);
  /// Defines an external field on a rectangular domain
  void defineExternalField (T x0, T x1, T y0, T y1, int fieldBeginsAt, int sizeOfField, T* field );
  /// Defines an external field on a domain with a particular material number
  void defineExternalField(SuperGeometry2D<T>& sGeometry, int material,
                           int fieldBeginsAt, int sizeOfField, AnalyticalF2D<T,T>& field);
  /// Defines an external field on a domain with a particular indicator
  void defineExternalField(SuperGeometry2D<T>& sGeometry, IndicatorF2D<T>& indicator,
                           int fieldBeginsAt, int sizeOfField, AnalyticalF2D<T,T>& field);
  /// Resets an external Particle Field
  void resetExternalParticleField(SuperGeometry2D<T>& sGeometry, IndicatorF2D<T>& indicator);
  void setExternalParticleField(SuperGeometry2D<T>& sGeometry, AnalyticalF2D<T,T>& velocity,
                                ParticleIndicatorF2D<T,T>& sIndicator);
  void addExternalField(SuperGeometry2D<T>& sGeometry, IndicatorF2D<T>& indicator,
                        int fieldBeginsAt, int sizeOfField, AnalyticalF2D<T,T>& field);
  void addExternalField(SuperGeometry2D<T>& sGeometry, IndicatorF2D<T>& indicator,
                        int fieldBeginsAt, int sizeOfField, AnalyticalF2D<T,T>& field,
                        AnalyticalF2D<T,T>& porous);
  void multiplyExternalField(SuperGeometry2D<T>& sGeometry, IndicatorF2D<T>& indicator,
                             int fieldBeginsAt, int sizeOfField, AnalyticalF2D<T,T>& field);
  /// Defines an external field on a domain with a particular material number
  void defineExternalField(SuperGeometry2D<T>& sGeometry, int material, int fieldBeginsAt,
                           int sizeOfField, SuperLatticeF2D<T,Lattice>& field);

  /// Initializes the equilibrium
  void iniEquilibrium (T x0, T x1, T y0, T y1, T rho, const T u[Lattice<T>::d] );
  /// Initializes the equilibrium on a domain with a particular material number
  void iniEquilibrium(SuperGeometry2D<T>& sGeometry, int material,
                      AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u);

  /// Apply collision step to a rectangular domain
  void collide(T x0, T x1, T y0, T y1);
  /// Apply collision step to the whole domain
  void collide();
  /// TO BE DONE: Apply collision step to a rectangular domain,
  /// with fixed velocity
  // void staticCollide(T x0, T x1, T y0, T y1,
  //                  TensorField2D<T,2> const& u);
  /// TO BE DONE: Apply collision step to the whole domain,
  /// with fixed velocity
  // void staticCollide(TensorField2D<T,2> const& u);
  /// Apply streaming step to a rectangular domain
  void stream(T x0, T x1, T y0, T y1);
  /// Apply streaming step to the whole domain
  void stream();
  /// TO BE DONE: Apply first collision, then streaming step
  /// to a rectangular domain
  // void collideAndStream(T x0, T x1, T y0, T y1);
  /// Apply first collision, then streaming step
  /// to the whole domain
  void collideAndStream();
  /// Subtract a constant offset from the density within the whole domain
  void stripeOffDensityOffset (int x0_, int x1_, int y0_, int y1_, T offset );
  /// Subtract a constant offset from the density within a rect. domain
  void stripeOffDensityOffset(T offset);
  /// Switches Statistics on (default on)
  void statisticsOn()
  {
    _statistics_on = true;
  };
  /// Switches Statistics off (default on). That speeds up
  /// the execution time.
  void statisticsOff()
  {
    _statistics_on = false;
  };

  /// Adds a coupling generator for one partner superLattice
  template<template<typename U> class Slattice>
  void addLatticeCoupling(SuperGeometry2D<T>& sGeometry, int material,
                          LatticeCouplingGenerator2D<T, Lattice> const& lcGen,
                          SuperLattice2D<T,Slattice>& partnerLattice );
  /// Executes coupling generatur for one partner superLattice
  void executeCoupling();

  //void communicate(bool verbose=true);

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;
private:
  /// Resets and reduce the statistics
  void reset_statistics();
};

} // namespace olb

#endif
