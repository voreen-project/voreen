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
#include "cellD.h"
#include "blockLatticeView2D.h"
#include "communication/communicator2D.h"
#include "communication/superPropagation.h"
#include "postProcessing.h"
#include "serializer.h"
#include "communication/superStructure2D.h"
#include "utilities/functorPtr.h"

#include "core/olbDebug.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T> class CuboidGeometry2D;
template<typename T, typename DESCRIPTOR> class BlockLattice2D;
template<typename T, typename DESCRIPTOR> class BlockLatticeView2D;
template<typename T> class LoadBalancer;
template<typename T> class SuperGeometry2D;
template<typename T, typename DESCRIPTOR> class SuperLatticeF2D;
template<typename T> class SuperStructure2D;
template<typename T> class SuperIndicatorF2D;


/// Super class maintaining block lattices for a cuboid decomposition
template<typename T, typename DESCRIPTOR>
class SuperLattice2D : public SuperStructure2D<T>, public BufferSerializable {
private:
  /// Lattices with ghost cell layer of size overlap
  std::vector<BlockLattice2D<T,DESCRIPTOR> >     _extendedBlockLattices;
  /// View of the lattices without overlap
  std::vector<BlockLatticeView2D<T,DESCRIPTOR> > _blockLattices;

#ifdef NEW_INTERIM_BLOCK_PROPAGATION
  /// Communicator for propagation of populations between blocks
  SuperPropagationCommunicator<T,DESCRIPTOR> _commStream;
#else
  /// Communicator for propagation of populations between blocks
  Communicator2D<T> _commStream;
#endif

  /// Communicator for non-local boundary conditions
  Communicator2D<T>                           _commBC;
  /// Specifies if there is communication for non local boundary conditions
  /// needed. It is automatically swichted on if overlapBC >= 1 by the
  /// calling the constructer. (default = false)
  bool                                        _commBC_on;

  /// Statistics of the super structure
  LatticeStatistics<T>                        _statistics;
  /// Specifies if statistics are to be calculated
  /**
   * Always needed for the ConstRhoBGK dynamics. (default = true)
   **/
  bool                                        _statistics_on;

public:
  /// Construct lattice for the cuboid decomposition of superGeometry 
  SuperLattice2D(SuperGeometry2D<T>& superGeometry);

  SuperLattice2D(const SuperLattice2D&) = delete;
  ~SuperLattice2D() = default;

  /// Read and write access to a block lattice
  BlockLattice2D<T,DESCRIPTOR>& getExtendedBlockLattice(int locIC)
  {
    return _extendedBlockLattices[locIC];
  };
  /// Read only access to a block lattice
  BlockLattice2D<T,DESCRIPTOR> const& getExtendedBlockLattice(int locIC) const
  {
    return _extendedBlockLattices[locIC];
  };

  /// Read and write access to a lattice (block lattice view, one
  /// without overlap).
  BlockLatticeView2D<T,DESCRIPTOR>& getBlockLattice(int locIC)
  {
    return _blockLattices[locIC];
  };
  /// Read only access to a lattice
  BlockLatticeView2D<T,DESCRIPTOR> const& getBlockLattice(int locIC) const
  {
    return _blockLattices[locIC];
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

  /// Return local Cell
  Cell<T,DESCRIPTOR> get(int iC, int iX, int iY);
  Cell<T,DESCRIPTOR> get(int latticeR[DESCRIPTOR::d+1])
  {
    return get(latticeR[0], latticeR[1], latticeR[2]);
  }
  Cell<T,DESCRIPTOR> get(const std::vector<int>& latticeR)
  {
    return get(latticeR[0], latticeR[1], latticeR[2]);
  }

  CellD<T,DESCRIPTOR> getCopy(int iC, int iX, int iY) const;
  CellD<T,DESCRIPTOR> getCopy(Vector<int,3> pos) const;

  CellD<T,DESCRIPTOR> getCopy(T locX, T locY) const;
  CellD<T,DESCRIPTOR> getCopy(Vector<T,2> pos) const;

  void setCopy(int iC, int iX, int iY, ConstCell<T,DESCRIPTOR>& cell);
  void setCopy(Vector<int,3> pos, ConstCell<T,DESCRIPTOR>& cell);

  bool setCopy(T locX, T locY, ConstCell<T,DESCRIPTOR>& cell);
  bool setCopy(Vector<T,2> pos, ConstCell<T,DESCRIPTOR>& cell);

  /// Write access to the memory of the data of the super structure
  std::uint8_t* operator() (int iCloc, int iX, int iY, int iPop) override
  {
    OLB_ASSERT(iPop < DESCRIPTOR::q, "iPop out of bounds");
    return reinterpret_cast<std::uint8_t*>(&getBlockLattice(iCloc).getPop(iX,iY,iPop));
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

  /// Defines the dynamics on a domain described by an indicator
  void defineDynamics(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, Dynamics<T, DESCRIPTOR>* dynamics);
  /// Defines the dynamics by material number
  void defineDynamics(SuperGeometry2D<T>& superGeometry, int material, Dynamics<T,DESCRIPTOR>* dynamics);

  /// Defines rho on a rectangular domain
  void defineRhoU (T x0, T x1, T y0, T y1, T rho, const T u[DESCRIPTOR::d] );
  /// Defines rho and u on a domain described by an indicator
  void defineRhoU(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                  AnalyticalF<2,T,T>& rho, AnalyticalF<2,T,T>& u);
  /// Defines rho and u on a domain with a particular material number
  void defineRhoU(SuperGeometry2D<T>& sGeometry, int material,
                  AnalyticalF<2,T,T>& rho, AnalyticalF<2,T,T>& u);

  /// Defines rho on a rectangular domain
  void defineRho (T x0, T x1, T y0, T y1, T rho );
  /// Defines rho on a domain described by an indicator
  void defineRho(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, AnalyticalF<2,T,T>& rho);
  /// Defines rho on a domain with a particular material number
  void defineRho(SuperGeometry2D<T>& sGeometry, int material, AnalyticalF<2,T,T>& rho);

  /// Defines u on a rectangular domain
  void defineU (T x0, T x1, T y0, T y1, const T u[DESCRIPTOR::d] );
  /// Defines u on a domain described by an indicator
  void defineU(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, AnalyticalF<2,T,T>& u);
  /// Defines u on a domain with a particular material number
  void defineU(SuperGeometry2D<T>& sGeometry, int material, AnalyticalF<2,T,T>& u);

  // Defines a population on a domain described by an indicator
  void definePopulations(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, AnalyticalF<2,T,T>& Pop);
  // Defines a population on a domain with a particular material number
  void definePopulations(SuperGeometry2D<T>& sGeometry, int material, AnalyticalF<2,T,T>& Pop);

  /// Defines a field on a rectangular domain
  template <typename FIELD>
  void defineField (T x0, T x1, T y0, T y1, T* field );
  /// Defines a field on a domain described by an indicator
  template <typename FIELD>
  void defineField(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                   SuperLatticeF2D<T,DESCRIPTOR>& field);
  /// Defines a field on a domain with a particular material number
  template <typename FIELD>
  void defineField(SuperGeometry2D<T>& sGeometry, int material,
                   SuperLatticeF2D<T,DESCRIPTOR>& field);

  /// Defines a field on a domain described by an indicator
  template <typename FIELD>
  void defineField(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                   AnalyticalF<2,T,T>& field)
  {
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      _extendedBlockLattices[iC].template defineField<FIELD>(
        indicator->getExtendedBlockIndicatorF(iC), field);
    }
  }

  /// Defines a field on a domain with a particular material number
  template <typename FIELD>
  void defineField(SuperGeometry2D<T>& sGeometry, int material,
                   AnalyticalF<2,T,T>& field)
  {
    defineField<FIELD>(sGeometry.getMaterialIndicator(material), field);
  }

  /// Defines a field on a indicated domain
  template <typename FIELD>
  void defineField(SuperGeometry2D<T>& sGeometry, IndicatorF2D<T>& indicator,
                   AnalyticalF<2,T,T>& field);


  template <typename FIELD>
  void addField(SuperGeometry2D<T>& sGeometry, IndicatorF2D<T>& indicator,
                AnalyticalF<2,T,T>& field);
  template <typename FIELD>
  void addField(SuperGeometry2D<T>& sGeometry, IndicatorF2D<T>& indicator,
                AnalyticalF<2,T,T>& field,
                AnalyticalF<2,T,T>& porous);
  template <typename FIELD>
  void multiplyField(SuperGeometry2D<T>& sGeometry, IndicatorF2D<T>& indicator,
                     AnalyticalF<2,T,T>& field);

  /// Initializes the equilibrium
  void iniEquilibrium (T x0, T x1, T y0, T y1, T rho, const T u[DESCRIPTOR::d] );
  /// Initializes the equilibrium on a domain described by an indicator
  void iniEquilibrium(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                      AnalyticalF<2,T,T>& rho, AnalyticalF<2,T,T>& u);
  /// Initializes the equilibrium on a domain with a particular material number
  void iniEquilibrium(SuperGeometry2D<T>& sGeometry, int material,
                      AnalyticalF<2,T,T>& rho, AnalyticalF<2,T,T>& u);

  /// Initializes populations with equilibrium and non-equilibrium
  void iniRegularized (T x0, T x1, T y0, T y1, T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR >::n] );
  /// Initializes populations on a domain described by an indicator
  void iniRegularized(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                      AnalyticalF<2,T,T>& rho, AnalyticalF<2,T,T>& u, AnalyticalF<2,T,T>& pi);
  /// Initializes populations on a domain with a particular material number
  void iniRegularized(SuperGeometry2D<T>& sGeometry, int material,
                      AnalyticalF<2,T,T>& rho, AnalyticalF<2,T,T>& u, AnalyticalF<2,T,T>& pi);

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
  template<typename Slattice>
  void addLatticeCoupling(LatticeCouplingGenerator2D<T, DESCRIPTOR> const& lcGen,
                          SuperLattice2D<T,Slattice>& partnerLattice );
  /// Adds a coupling generator for one partner superLattice
  template<typename Slattice>
  void addLatticeCoupling(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                          LatticeCouplingGenerator2D<T, DESCRIPTOR> const& lcGen,
                          SuperLattice2D<T,Slattice>& partnerLattice );
  /// Adds a coupling generator for one partner superLattice
  template<typename Slattice>
  void addLatticeCoupling(SuperGeometry2D<T>& sGeometry, int material,
                          LatticeCouplingGenerator2D<T, DESCRIPTOR> const& lcGen,
                          SuperLattice2D<T,Slattice>& partnerLattice );
  /// Adds a coupling generator for a multiple partner superLattice
  template<typename Slattice>
  void addLatticeCoupling(LatticeCouplingGenerator2D<T, DESCRIPTOR> const& lcGen,
                          std::vector<SuperLattice2D<T,Slattice>*> partnerLattices );
  /// Adds a coupling generator for a multiple partner superLattice
  template<typename Slattice>
  void addLatticeCoupling(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                          LatticeCouplingGenerator2D<T, DESCRIPTOR> const& lcGen,
                          std::vector<SuperLattice2D<T,Slattice>*> partnerLattices );
  /// Adds a coupling generator for a multiple partner superLattice
  template<typename Slattice>
  void addLatticeCoupling(SuperGeometry2D<T>& sGeometry, int material,
                          LatticeCouplingGenerator2D<T, DESCRIPTOR> const& lcGen,
                          std::vector<SuperLattice2D<T,Slattice>*> partnerLattices );
  /// Add a non-local post-processing step
  void addPostProcessor(PostProcessorGenerator2D<T, DESCRIPTOR> const& ppGen);
  void addPostProcessor(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                        PostProcessorGenerator2D<T, DESCRIPTOR> const& ppGen);
  void addPostProcessor(SuperGeometry2D<T>& sGeometry, int material,
                        PostProcessorGenerator2D<T, DESCRIPTOR> const& ppGen);
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
void setSuperExternalParticleField( SuperGeometry2D<T>& sGeometry, AnalyticalF<2,T,T>& velocity,
                                    SmoothIndicatorF2D<T,T,true>&,
                                    SuperLattice2D<T, DESCRIPTOR>& sLattice );

template<typename T, typename DESCRIPTOR>
void setSuperExternalParticleField( SuperGeometry2D<T>& sGeometry, AnalyticalF<2,T,T>& velocity,
                                    SmoothIndicatorF2D<T,T,true>&,
                                    SuperLattice2D<T, DESCRIPTOR>& sLattice,
                                    Vector<bool,2> periodicity );

//Geng2019
template<typename T, typename DESCRIPTOR>
void setSuperZetaParticleField( SuperGeometry2D<T>& sGeometry, AnalyticalF<2,T,T>& velocity,
                                SmoothIndicatorF2D<T,T,true>& sIndicator,
                                SuperLattice2D<T, DESCRIPTOR>& sLattice );

} // namespace olb

#endif
