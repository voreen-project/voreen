/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias J. Krause, Jonathan Jeppener-Haltenhoff
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
 * The description of a 3D super lattice -- header file.
 */

#ifndef SUPER_LATTICE_3D_H
#define SUPER_LATTICE_3D_H

#include <vector>
#include <map>

#include "blockLattice3D.h"
#include "cellD.h"
#include "blockLatticeView3D.h"
#include "superData3D.h"
#include "communication/communicator3D.h"
#include "communication/superPropagation.h"
#include "postProcessing.h"
#include "serializer.h"
#include "communication/superStructure3D.hh"
#include "utilities/functorPtr.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T, typename DESCRIPTOR> class SuperLatticeF3D;
template<typename T> class SuperGeometry3D;
template<typename T> class SuperIndicatorF3D;
template<typename T> class LoadBalancer;
template<typename T> class CuboidGeometry3D;
template <typename T> class IndicatorSphere3D;
template<typename T, typename DESCRIPTOR> class BlockLatticeView3D;


/// Super class maintaining block lattices for a cuboid decomposition
template<typename T, typename DESCRIPTOR>
class SuperLattice3D : public SuperStructure3D<T>, public BufferSerializable {
private:
  /// Lattices with ghost cell layer of size overlap
  std::vector<BlockLattice3D<T, DESCRIPTOR> > _extendedBlockLattices;
  /// View of the lattices without overlap
  std::vector<BlockLatticeView3D<T, DESCRIPTOR> > _blockLattices;

#ifdef NEW_INTERIM_BLOCK_PROPAGATION
  /// Communicator for propagation of populations between blocks
  SuperPropagationCommunicator<T,DESCRIPTOR> _commStream;
#else
  /// Communicator for propagation of populations between blocks
  Communicator3D<T> _commStream;
#endif

  /// Communicator for non-local boundary conditions
  Communicator3D<T> _commBC;
  /// Specifies if there is communication for non local boundary conditions
  /// needed. It is automatically swichted on if overlapBC >= 1 by the
  /// calling the constructer. (default = false)
  bool _commBC_on;

  /// Statistics of the super structure
  LatticeStatistics<T> _statistics;
  /// Specifies if statistics are to be calculated
  /**
   * Always needed for the ConstRhoBGK dynamics. (default = true)
   **/
  bool _statistics_on;

public:
  /// Construct lattice for the cuboid decomposition of superGeometry 
  SuperLattice3D(SuperGeometry3D<T>& superGeometry);

  SuperLattice3D(const SuperLattice3D&) = delete;
  ~SuperLattice3D() = default;

  /// Read and write access to a block lattice
  BlockLattice3D<T, DESCRIPTOR>& getExtendedBlockLattice(int locIC)
  {
    return _extendedBlockLattices[locIC];
  }
  ;
  /// Read only access to a block lattice
  BlockLattice3D<T, DESCRIPTOR> const& getExtendedBlockLattice(int locIC) const
  {
    return _extendedBlockLattices[locIC];
  }
  ;
  /// Read and write access to a lattice (block lattice view, one
  /// without overlap).
  BlockLatticeView3D<T, DESCRIPTOR>& getBlockLattice(int locIC)
  {
    return _blockLattices[locIC];
  }
  ;
  /// Read only access to a lattice
  BlockLatticeView3D<T, DESCRIPTOR> const& getBlockLattice(int locIC) const
  {
    return _blockLattices[locIC];
  }
  ;

  /// Read and write access to the boundary communicator
  Communicator3D<T>& get_commBC()
  {
    return _commBC;
  }
  ;
  /// Read only access to the boundary communicator
  Communicator3D<T> const& get_commBC() const
  {
    return _commBC;
  }
  ;

  /// Return a handle to the LatticeStatistics object
  LatticeStatistics<T>& getStatistics();
  /// Return a constant handle to the LatticeStatistics object
  LatticeStatistics<T> const& getStatistics() const;

  /// Return local Cell
  Cell<T,DESCRIPTOR> get(int iC, int iX, int iY, int iZ);
  Cell<T,DESCRIPTOR> get(int latticeR[DESCRIPTOR::d+1])
  {
    return get(latticeR[0], latticeR[1], latticeR[2], latticeR[3]);
  }
  Cell<T,DESCRIPTOR> get(const std::vector<int>& latticeR)
  {
    return get(latticeR[0], latticeR[1], latticeR[2], latticeR[3]);
  }

  /// Return Cell copy from anywhere
  CellD<T,DESCRIPTOR> getCopy(int iC, int iX, int iY, int iZ) const;
  CellD<T,DESCRIPTOR> getCopy(Vector<int,4> pos) const;

  CellD<T,DESCRIPTOR> getCopy(T locX, T locY, T locZ) const;
  CellD<T,DESCRIPTOR> getCopy(Vector<T,3> pos) const;

  void setCopy(int iC, int iX, int iY, int iZ, ConstCell<T,DESCRIPTOR>& cell);
  void setCopy(Vector<int,4> pos, ConstCell<T,DESCRIPTOR>& cell);

  bool setCopy(T locX, T locY, T locZ, ConstCell<T,DESCRIPTOR>& cell);
  bool setCopy(Vector<T,3> pos, ConstCell<T,DESCRIPTOR>& cell);

  /// Write access to the memory of the data of the super structure
  std::uint8_t* operator() (int iCloc, int iX, int iY, int iZ, int iPop) override
  {
    OLB_ASSERT(iPop < DESCRIPTOR::q, "iPop out of bounds");
    return reinterpret_cast<std::uint8_t*>(&getBlockLattice(iCloc).getPop(iX,iY,iZ,iPop));
  };
  std::uint8_t* operator() (int iCloc, std::size_t iCell, int iPop) override
  {
    OLB_ASSERT(iPop < DESCRIPTOR::q, "iPop out of bounds");
    return reinterpret_cast<std::uint8_t*>(&getExtendedBlockLattice(iCloc).getPop(iCell,iPop));
  };

  /// Read only access to the dim of the data of the super structure
  int getDataSize() const override
  {
    return DESCRIPTOR::q;
  };
  /// Read only access to the data type dim of the data of the super structure
  int getDataTypeSize() const override
  {
    return sizeof(T);
  };

  /// Initialize all lattice cells to become ready for simulation
  void initialize();

  /// Define the dynamics on a domain described by an indicator
  void defineDynamics(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, Dynamics<T, DESCRIPTOR>* dynamics);
  /// Define the dynamics on a domain with a particular material number
  void defineDynamics(SuperGeometry3D<T>& sGeometry, int material,
                      Dynamics<T, DESCRIPTOR>* dynamics);

  /// Define rho on a domain described by an indicator
  /**
   * \param indicator Indicator describing the target domain
   * \param rho       Analytical functor
   **/
  void defineRho(FunctorPtr<SuperIndicatorF3D<T>>&&, AnalyticalF<3,T,T>& rho);
  /// Define rho on a domain with a particular material number
  void defineRho(SuperGeometry3D<T>& sGeometry, int material, AnalyticalF<3,T,T>& rho);

  /// Define u on a domain described by an indicator
  /**
   * \param indicator Indicator describing the target domain
   * \param u         Analytical functor
   **/
  void defineU(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, AnalyticalF<3,T,T>& u);
  /// Define u on a domain with a particular material number
  void defineU(SuperGeometry3D<T>& sGeometry, int material, AnalyticalF<3,T,T>& u);

  /// Define rho and u on a domain described by an indicator
  /**
   * \param indicator Indicator describing the target domain
   * \param rho       Analytical functor
   * \param u         Analytical functor
   **/
  void defineRhoU(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                  AnalyticalF<3,T,T>& rho, AnalyticalF<3,T,T>& u);
  /// Define rho and u on a domain with a particular material number
  void defineRhoU(SuperGeometry3D<T>& sGeometry, int material,
                  AnalyticalF<3,T,T>& rho, AnalyticalF<3,T,T>& u);

  /// Define a population on a domain described by an indicator
  /**
   * \param indicator Indicator describing the target domain
   * \param Pop       Analytical functor, target dimension DESCRIPTOR::q
   **/
  void definePopulations(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                         AnalyticalF<3,T,T>& Pop);
  /// Define a population on a domain with a particular material number
  void definePopulations(SuperGeometry3D<T>& sGeometry, int material,
                         AnalyticalF<3,T,T>& Pop);
  /**
   * \param indicator Indicator describing the target domain
   * \param Pop       Super functor, target dimension DESCRIPTOR::q
   **/
  void definePopulations(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                         SuperF3D<T,T>& Pop);
  /// Define a population on a domain with a particular material number
  void definePopulations(SuperGeometry3D<T>& sGeometry, int material,
                         SuperF3D<T,T>& Pop);

  /// Define an external field on a domain with a particular material number
  template <typename FIELD>
  void defineField(SuperGeometry3D<T>& sGeometry, int material,
                   SuperF3D<T,T>& field);

  /// Define an external field on an indicated domain with a particular material number
  template <typename FIELD>
  void defineField(SuperGeometry3D<T>& sGeometry, IndicatorF3D<T>& indicator,
                   int material, SuperF3D<T,T>& field);

  /// Defines a field on a domain described by an indicator
  template <typename FIELD>
  void defineField(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                   AnalyticalF<3,T,T>& field)
  {
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      _extendedBlockLattices[iC].template defineField<FIELD>(
        indicator->getExtendedBlockIndicatorF(iC), field);
    }
  }
  /// Defines a field on a domain with a particular material number
  template <typename FIELD>
  void defineField(SuperGeometry3D<T>& sGeometry, int material,
                   AnalyticalF<3,T,T>& field)

  {
    defineField<FIELD>(sGeometry.getMaterialIndicator(material), field);
  }
  /// Defines a field on a indicated domain
  template <typename FIELD>
  void defineField(SuperGeometry3D<T>& sGeometry, IndicatorF3D<T>& indicator,
                   AnalyticalF<3,T,T>& field);
  /// Initialize by equilibrium on a domain described by an indicator
  /**
   * \param indicator Indicator describing the target domain
   * \param rho       Analytical functor (global)
   * \param u         Analytical functor (global)
   **/
  void iniEquilibrium(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                      AnalyticalF<3,T,T>& rho, AnalyticalF<3,T,T>& u);
  /// Initialize by equilibrium on a domain with a particular material number
  void iniEquilibrium(SuperGeometry3D<T>& sGeometry, int material,
                      AnalyticalF<3,T,T>& rho, AnalyticalF<3,T,T>& u);
  /// Initialize by non- and equilibrium on a domain described by an indicator
  /**
   * \param indicator Indicator describing the target domain
   * \param rho       Analytical functor (global)
   * \param u         Analytical functor (global)
   * \param pi        Analytical functor (global)
   **/
  void iniRegularized(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                      AnalyticalF<3,T,T>& rho, AnalyticalF<3,T,T>& u,
                      AnalyticalF<3,T,T>& pi);
  /// Initialize by non- and equilibrium on a domain with a particular material number
  void iniRegularized(SuperGeometry3D<T>& sGeometry, int material,
                      AnalyticalF<3,T,T>& rho, AnalyticalF<3,T,T>& u,
                      AnalyticalF<3,T,T>& pi);

  void collideAndStream();

  /// Subtract a constant offset from the density within the whole domain
  void stripeOffDensityOffset (
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T offset );
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
  template<typename Slattice>
  void addLatticeCoupling(LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
                          SuperLattice3D<T,Slattice>& partnerLattice );
  /// Adds a coupling generator for one partner superLattice
  template<typename Slattice>
  void addLatticeCoupling(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                          LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
                          SuperLattice3D<T,Slattice>& partnerLattice );
  /// Adds a coupling generator for one partner superLattice
  template<typename Slattice>
  void addLatticeCoupling(SuperGeometry3D<T>& sGeometry, int material,
                          LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
                          SuperLattice3D<T,Slattice>& partnerLattice );
  /// Adds a coupling generator for a multiple partner superLattice
  template<typename Slattice>
  void addLatticeCoupling(LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
                          std::vector<SuperLattice3D<T,Slattice>*> partnerLattices );
  /// Adds a coupling generator for a multiple partner superLattice
  template<typename Slattice>
  void addLatticeCoupling(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                          LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
                          std::vector<SuperLattice3D<T,Slattice>*> partnerLattices );
  /// Adds a coupling generator for a multiple partner superLattice
  template<typename Slattice>
  void addLatticeCoupling(SuperGeometry3D<T>& sGeometry, int material,
                          LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
                          std::vector<SuperLattice3D<T,Slattice>*> partnerLattices );
  /// Add a non-local post-processing step
  void addPostProcessor(PostProcessorGenerator3D<T, DESCRIPTOR> const& ppGen);
  void addPostProcessor(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                        PostProcessorGenerator3D<T, DESCRIPTOR> const& ppGen);
  void addPostProcessor(SuperGeometry3D<T>& sGeometry, int material,
                        PostProcessorGenerator3D<T, DESCRIPTOR> const& ppGen);
  /// Executes coupling generator for one partner superLattice
  void executeCoupling();

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

  /// Reinit population structure after deserialization
  void postLoad();

private:
  /// Resets and reduce the statistics
  void reset_statistics();
};

////////// FREE FUNCTIONS //////////

template<typename T, typename DESCRIPTOR>
void setSuperExternalParticleField( SuperGeometry3D<T>& sGeometry, AnalyticalF<3,T,T>& velocity,
                                    SmoothIndicatorF3D<T,T,true>&,
                                    SuperLattice3D<T, DESCRIPTOR>& sLattice );

template<typename T, typename DESCRIPTOR>
void setSuperExternalParticleField( SuperGeometry3D<T>& sGeometry, AnalyticalF<3,T,T>& velocity,
                                    SmoothIndicatorF3D<T,T,true>&,
                                    SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                    Vector<bool,3> periodicity );

} // namespace olb

#endif
