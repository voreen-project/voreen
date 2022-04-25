/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Davide Dapelo
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

/* Drag force models for Lagrangian two-way coupling -- header file.
 */

#ifndef LB_DRAG_MODELS_H
#define LB_DRAG_MODELS_H

#include "functors/lattice/reductionF3D.h"
#include "ReynoldsNumbers3D.h"

namespace olb {

/** Abstact base class for DragModelBase.
  * Its raison d'etre consists of not being templetized in Lattice.
  */
template<typename T, template<typename V> class Particle>
class DragModel {
public:
  /// Returns the scalar drag coefficient to overload. globicFull = { globic, latticeRoundedP[0, ..., 2] }
  virtual T operator() (Particle<T>* p, T latticeVelF[], T latticeVelP[], int globicFull[], T continuousPhaseFraction=1.)=0;
  virtual std::string name()=0;
protected:
  /// Functional to compute particle Reynolds number
  std::shared_ptr<ParticleReynoldsNumber<T, Particle> > _reP;
};

/** Abstact class for all the drag models.
  * The virtual method computeDragCoeff returns the drag coefficient.
  * Input parameters in attice units.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class DragModelBase : public DragModel<T,Particle> {
public:
  /// Constructor
  DragModelBase(UnitConverter<T, Lattice>& converter);
  /// Returns the scalar drag coefficient to overload.
  //virtual T operator() (Particle<T>* p, T latticeVelF[], T latticeVelP[], int globic)=0;
protected:
  UnitConverter<T, Lattice>& _converter; // reference to a UnitConverter
};

/** Class to compute a drag coefficient Cd=1.83 for low-Re Stokes drag.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class StokesSimplifiedDragModel : public DragModelBase<T,Lattice,Particle> {
public:
  /// Constructor
  StokesSimplifiedDragModel(UnitConverter<T, Lattice>& converter);
  /// Returns the scalar drag coefficient. globicFull = { globic, latticeRoundedP[0, ..., 2] }
  virtual T operator() (Particle<T>* p, T latticeVelF[], T latticeVelP[], int globicFull[], T continuousPhaseFraction=1.) override;
  std::string name() override {return "StokesSimplifiedDragModel";}
};

/** Class to compute the standard drag coefficient
  * as in Morsi and Alexander (1972).
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class MorsiDragModel : public DragModelBase<T,Lattice,Particle> {
public:
  /// Constructor
  MorsiDragModel(UnitConverter<T, Lattice>& converter);
  /// Returns the scalar drag coefficient. globicFull = { globic, latticeRoundedP[0, ..., 2] }
  virtual T operator() (Particle<T>* p, T latticeVelF[], T latticeVelP[], int globicFull[], T continuousPhaseFraction=1.) override;
  std::string name() override {return "MorsiDragModel";}
};

/** Class to compute the standard drag coefficient
  * as in Morsi and Alexander (1972), in a power-law fluid.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class PowerLawMorsiDragModel : public MorsiDragModel<T,Lattice,Particle> {
public:
  /// Constructor
  PowerLawMorsiDragModel ( UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice );
  std::string name() override {return "PowerLawMorsiDragModel";}
};

/** Class to compute the standard drag coefficient
  * as in Schiller and Naumann (1935).
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class SchillerNaumannDragModel : public DragModelBase<T,Lattice,Particle> {
public:
  /// Constructor
  SchillerNaumannDragModel(UnitConverter<T, Lattice>& converter);
  /// Returns the scalar drag coefficient. globicFull = { globic, latticeRoundedP[0, ..., 2] }
  virtual T operator() (Particle<T>* p, T latticeVelF[], T latticeVelP[], int globicFull[], T continuousPhaseFraction=1.) override;
  std::string name() override {return "SchillerNaumannDragModel";}
};

/** Class to compute the standard drag coefficient
  * as in Schiller and Naumann (1935), in a power-law fluid.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class PowerLawSchillerNaumannDragModel : public SchillerNaumannDragModel<T,Lattice,Particle> {
public:
  /// Constructor
  PowerLawSchillerNaumannDragModel ( UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice );
  std::string name() override {return "PowerLawSchillerNaumannDragModel";}
};

/** Class to compute the drag coefficient for gas bubbles in a liquid fluid phase
  * as in Dewsbury et al. (1999).
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class DewsburyDragModel : public DragModelBase<T,Lattice,Particle> {
public:
  /// Constructor
  DewsburyDragModel(UnitConverter<T, Lattice>& converter);
  /// Returns the scalar drag coefficient. globicFull = { globic, latticeRoundedP[0, ..., 2] }
  virtual T operator() (Particle<T>* p, T latticeVelF[], T latticeVelP[], int globicFull[], T continuousPhaseFraction=1.) override;
  std::string name() override {return "DewsburyDragModel";}
};

/** Class to compute the drag coefficient for gas bubbles in a liquid fluid phase
  * as in Dewsbury et al. (1999), in a power-law fluid.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class PowerLawDewsburyDragModel : public DewsburyDragModel<T,Lattice,Particle> {
public:
  /// Constructor
  PowerLawDewsburyDragModel ( UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice );
  std::string name() override {return "PowerLawDewsburyDragModel";}
};

/** Class to compute the drag coefficient for gas bubbles in a liquid fluid phase
  * as in Sun, Guo, Wang et al. (2015).
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class SunDragModel : public DragModelBase<T,Lattice,Particle> {
public:
  /// Constructor
  SunDragModel(UnitConverter<T, Lattice>& converter);
  /// Returns the scalar drag coefficient. globicFull = { globic, latticeRoundedP[0, ..., 2] }
  virtual T operator() (Particle<T>* p, T latticeVelF[], T latticeVelP[], int globicFull[], T continuousPhaseFraction=1.) override;
  std::string name() override {return "SunDragModel";}
};

/** Class to compute the drag coefficient for gas bubbles in a liquid fluid phase
  * as in Sun, Guo, Wang et al. (2015), in a power-law fluid.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class PowerLawSunDragModel : public SunDragModel<T,Lattice,Particle> {
public:
  /// Constructor
  PowerLawSunDragModel ( UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice );
  std::string name() override {return "PowerLawSunDragModel";}
};


}

#endif
