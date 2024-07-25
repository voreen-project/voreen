/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause, Albert Mink
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

#ifndef ANALYTICAL_F_H
#define ANALYTICAL_F_H

#include <vector>
#include <random>

#include "analyticalBaseF.h"
#include "geometry/superGeometry.h"
#include "indicator/smoothIndicatorF2D.h"
#include "indicator/smoothIndicatorF3D.h"
#include "utilities/adHelpers.h"
#include "utilities/dimensionConverter.h"


/**
 *  The functor dimensions are given by F: S^m -> T^n  (S=source, T=target)
 *  and are implemented via GenericF(n,m).
 *  Don't get confused by the flipped order of source and target.
 */

namespace olb {

template<typename T, typename S, bool> class SmoothIndicatorSphere3D;
template<typename T,  typename DESCRIPTOR> class RadiativeUnitConverter;

////////////////////////////////////////////////////////////////////////////////
////////implementation of several 1d,2d,3d functors (analyticalFXD)/////////////
////////////////////////////////////////////////////////////////////////////////

template <unsigned D, typename T, typename S>
class AnalyticalComposed final : public AnalyticalF<D,T,S> {
private:
  std::vector<std::reference_wrapper<AnalyticalF<D,T,S>>> _f;
public:
  template <unsigned otherD=D, typename = typename std::enable_if_t<otherD==2>>
  AnalyticalComposed(AnalyticalF<D,T,S>& f0, AnalyticalF<D,T,S>& f1)
    : AnalyticalF<D,T,S>(2), _f{f0, f1}
  {
    this->getName() = "composed";
  }
  template <unsigned otherD=D, typename = typename std::enable_if_t<otherD==3>>
  AnalyticalComposed(AnalyticalF<D,T,S>& f0, AnalyticalF<D,T,S>& f1, AnalyticalF<D,T,S>& f2)
    : AnalyticalF<D,T,S>(3), _f{f0, f1, f2}
  {
    this->getName() = "composed";
  }
  AnalyticalComposed(std::vector<AnalyticalF<D,T,S>>& f);
  bool operator() (T output[], const S x[]) override;
};


/// AnalyticalConst: DD -> XD, where XD is defined by value.size()
template <unsigned D, typename T, typename S>
class AnalyticalConst final: public AnalyticalF<D,T,S> {
private:
  // is constant return value of operator()
  std::vector<T> _c;
public:
  AnalyticalConst(T value);
  AnalyticalConst(T value0, T value1);
  AnalyticalConst(T value0, T value1, T value2);
  AnalyticalConst(const Vector<T,2>& value);
  AnalyticalConst(const Vector<T,3>& value);
  AnalyticalConst(const std::vector<T>& value);
  bool operator() (T output[], const S x[]) override;

  template<typename V, typename U>
  using exchange_type = AnalyticalConst<D,V,U>;

  template<typename V, typename U>
  auto copyAs() const {
    return exchange_type<V,U>(util::copyAs<V,T>(_c));
  }
};

/// AnalyticalNormal: DD -> XD, where XD is defined by value.size()
template <unsigned D, typename T, typename S>
class AnalyticalNormal final: public AnalyticalF<D,T,S> {
private:
  // is constant return value of operator()
  std::vector<T> _mean;
  T _stdDev;
public:
  AnalyticalNormal(std::vector<T> mean, T stdDev);
  bool operator() (T output[], const S x[]) override;
};

/// AnalyticalRandomBase: virtual base class for all the random functionals
template <unsigned D, typename T, typename S>
class AnalyticalRandomBase : public AnalyticalF<D,T,S> {
protected:
  AnalyticalRandomBase();
  std::random_device rd;
  std::mt19937 gen;
};

/// AnalyticalRandomUniform: DD -> 1D with random image in (0,1)
template <unsigned D, typename T, typename S>
class AnalyticalRandomUniform : public AnalyticalRandomBase<D,T,S> {
public:
  AnalyticalRandomUniform(T minVal=0., T maxVal=1.);
  bool operator() (T output[], const S x[]) override;
private:
  std::uniform_real_distribution<T> distro;
};

/// AnalyticalRandomNormal: DD -> 1D with random image in (0,1)
template <unsigned D, typename T, typename S>
class AnalyticalRandomNormal : public AnalyticalRandomBase<D,T,S> {
public:
  AnalyticalRandomNormal(T mean=0., T stdDev=1.);
  bool operator() (T output[], const S x[]) override;
protected:
  std::normal_distribution<T> distro;
};

/// AnalyticalRandomNormal: DD -> 1D with random image in (0,1)
/// Normal distribution cut off outside [mean-n*stdDev, mean+n*stdDev]
template <unsigned D, typename T, typename S>
class AnalyticalRandomTruncatedNormal : public AnalyticalRandomNormal<D,T,S> {
public:
  AnalyticalRandomTruncatedNormal(T mean=0., T stdDev=1., T n=3.);
  bool operator() (T output[], const S x[]) override;
private:
  T _min;
  T _max;
};


/// AnalyticalRandomOld: DD -> 1D with random image in (0,1)
template <unsigned D, typename T, typename S>
class AnalyticalRandomOld : public AnalyticalF<D,T,S> {
public:
  AnalyticalRandomOld();
  bool operator() (T output[], const S x[]) override;
};


/** Computes resulting velocity of an object from translational and rotational velocity.
 * \param position the rotation center
 * \param velocity translational velocity of the object - expected in lattice units
 * \param angularVelocity rotational velocity of the object - expected in lattice units
 */
template <typename T, typename S, typename DESCRIPTOR>
class EccentricVelocityField : public AnalyticalF<DESCRIPTOR::d,T,S> {
protected:
  Vector<T,DESCRIPTOR::d> _position;
  Vector<T,DESCRIPTOR::d> _velocity;
  Vector<T,utilities::dimensions::convert<DESCRIPTOR::d>::rotation> _angularVelocity;
public:
  EccentricVelocityField( Vector<T,DESCRIPTOR::d> position,
                          Vector<T,DESCRIPTOR::d> velocity,
                          Vector<T,utilities::dimensions::convert<DESCRIPTOR::d>::rotation> angularVelocity );
  bool operator()(T output[], const S input[]) override;
};

/** Computes resulting lattice velocity of an object from translational and rotational velocity.
 * \param position the rotation center
 * \param velocity translational velocity of the object - expected in lattice units
 * \param angularVelocity rotational velocity of the object - expected in lattice units
 * \param converter unit converter to convert to lattice velocity
 */
template <typename T, typename S, typename DESCRIPTOR>
class EccentricLatticeVelocityField : public AnalyticalF<DESCRIPTOR::d,T,S> {
protected:
  Vector<T,DESCRIPTOR::d> _position;
  Vector<T,DESCRIPTOR::d> _velocity;
  Vector<T,utilities::dimensions::convert<DESCRIPTOR::d>::rotation> _angularVelocity;
  UnitConverter<T,DESCRIPTOR> const& _converter;
public:
  EccentricLatticeVelocityField( Vector<T,DESCRIPTOR::d> position,
                          Vector<T,DESCRIPTOR::d> velocity,
                          Vector<T,utilities::dimensions::convert<DESCRIPTOR::d>::rotation> angularVelocity,
                          UnitConverter<T,DESCRIPTOR> const& converter );
  bool operator()(T output[], const S input[]) override;
};



/// Square wave with given period length, amplitude, difference
/// (= length of positive time / length of period)
template <unsigned D, typename T, typename S>
class AnalyticalSquareWave : public AnalyticalF<D,T,S> {
public:
  AnalyticalSquareWave(T period=1, T amplitude=1, T difference=0.5);
  bool operator() (T output[], const S x[]) override;
protected:
  T _period;
  T _amplitude;
  T _difference;
};


/// Smoothed square wave. epsilon = width of the mollified interval
template <unsigned D, typename T, typename S>
class AnalyticalSmoothedSquareWave : public AnalyticalSquareWave<D,T,S> {
public:
  AnalyticalSmoothedSquareWave(T period=1, T amplitude=1, T difference=0.5, T epsilon=1.e-3);
  bool operator() (T output[], const S x[]) override;
protected:
  T _epsilon;
};


/** Concatenate an analytical functor with any other function.
 * \param _f: an analytical functor S^D -> T^n
 * \param _g: a function T^n -> U^k or T -> U
 * g can be e.g. a lambda expression or a free function.
 * If g: T -> U, then g is applied to each component of the result of f.
 * If g: T* -> U*, then g is applied to the complete vector. The user is
 * responsible for ensuring that the dimensions of f and g fit together.
 *
 * \example: Apply cosine to linearly transformed values:
 * AnalyticalLinear1D<double,double> f(1.1, 1.5);
 * AnalyticalConcatenation conc(f, std::cos);
 */
template <unsigned D, typename U, typename T, typename S,
  bool ComponentWise, bool ReturnArray>
class AnalyticalConcatenation : public AnalyticalF<D,U,S> {
protected:
  using return_type_g = std::conditional_t<ReturnArray,U*,U>;
  using function_t = std::conditional_t<ComponentWise,
    std::function<return_type_g(T)>, std::function<return_type_g(T*)>>;

  AnalyticalF<D,T,S>& _f;
  function_t _g;

public:
  /** concatenate functor f and some lambda expression g.
   * targetDim needs to be specified if g is vector (array)-valued
   */
  template <typename G>
  AnalyticalConcatenation(AnalyticalF<D,T,S>& f, G g, unsigned targetDim=1)
  : AnalyticalF<D,U,S>((ComponentWise) ? f.getTargetDim() : targetDim),
    _f(f), _g(g) { }

  /** concatenate functor f and some function g.
   * g maps value to value (no pointer) and is applied component-wise
   */
  AnalyticalConcatenation(AnalyticalF<D,T,S>& f, U (*g)(T))
  : AnalyticalF<D,U,S>(f.getTargetDim()), _f(f), _g(g) {
    static_assert(ComponentWise);
    static_assert(! std::is_pointer_v<U>);
  }

  /** concatenate functor f and some function g.
   * targetDim needs to be specified if g is vector (array)-valued
   */
  // wrapped_U equals U or a pointer to U
  template<typename wrapped_U>
  AnalyticalConcatenation(
    AnalyticalF<D,T,S>& f, wrapped_U (*g)(T*), unsigned targetDim=1)
  : AnalyticalF<D,U,S>(targetDim),
    _f(f), _g(g) {
    static_assert(! ComponentWise);
  }

  /** concatenate functor f and some function g.
   * targetDim needs to be specified if g is vector (array)-valued
   */
  // wrapped_U equals U or a pointer to U
  template<typename wrapped_U>
  AnalyticalConcatenation(
    AnalyticalF<D,T,S>& f, wrapped_U (*g)(const T*), unsigned targetDim=1)
  : AnalyticalF<D,U,S>(targetDim),
    _f(f), _g(g) {
    static_assert(! ComponentWise);
  }

  bool operator() (U output[], const S x[]) override {
    T outputTmp[_f.getTargetDim()];
    _f(outputTmp, x);
    if constexpr (ComponentWise) {
      for (int i = 0; i < _f.getTargetDim(); ++i) {
        output[i] = _g(outputTmp[i]);
      }
    } else {  // g works on the vector
      if constexpr (ReturnArray) {
        for (int i = 0; i < this->getTargetDim(); ++i) {
          const auto* outputTmp2 = _g(outputTmp);
          output[i] = outputTmp2[i];
        }
      } else {  // g returns value
        output[0] = _g(outputTmp);
      }
    }
    return true;
  }
};

template <unsigned D, typename T, typename S, typename G>
AnalyticalConcatenation(AnalyticalF<D,T,S>&, G g, unsigned _=1)
 -> AnalyticalConcatenation<D,
      std::remove_pointer_t<decltype(
        g(std::conditional_t<std::is_invocable_v<G,T>,T,T*>{}))>,T,S,
      std::is_invocable_v<G,T>,
      std::is_pointer_v<decltype(
        g(std::conditional_t<std::is_invocable_v<G,T>,T,T*>{}))>>;

// default for U seems to be necessary for multiply overloaded functions
template <unsigned D, typename T, typename S, typename U=T>
AnalyticalConcatenation(AnalyticalF<D,T,S>&, U (*g)(T))
 -> AnalyticalConcatenation<D,U,T,S,true,false>;

template <unsigned D, typename wrapped_U, typename T, typename S>
AnalyticalConcatenation(AnalyticalF<D,T,S>&, wrapped_U(T*), unsigned)
 -> AnalyticalConcatenation<D,std::remove_pointer_t<wrapped_U>,T,S,
      false,std::is_pointer_v<wrapped_U>>;

template <unsigned D, typename wrapped_U, typename T, typename S>
AnalyticalConcatenation(AnalyticalF<D,T,S>&, wrapped_U(const T*), unsigned)
 -> AnalyticalConcatenation<D,std::remove_pointer_t<wrapped_U>,T,S,
      false,std::is_pointer_v<wrapped_U>>;

////////////// CONVERSION FROM NEW TO OLD IMPLEMENTATION //////////////////////////

template <typename T, typename S>
using AnalyticalConst1D = AnalyticalConst<1,T,S>;
template <typename T, typename S>
using AnalyticalConst2D = AnalyticalConst<2,T,S>;
template <typename T, typename S>
using AnalyticalConst3D = AnalyticalConst<3,T,S>;

template <typename T, typename S>
using AnalyticalComposed2D = AnalyticalComposed<2,T,S>;
template <typename T, typename S>
using AnalyticalComposed3D = AnalyticalComposed<3,T,S>;

template <typename T, typename S>
using AnalyticalRandom1D = AnalyticalRandomOld<1,T,S>;
template <typename T, typename S>
using AnalyticalRandom2D = AnalyticalRandomOld<2,T,S>;
template <typename T, typename S>
using AnalyticalRandom3D = AnalyticalRandomOld<3,T,S>;






////////////// OLD IMPLEMENTATION //////////////////////////


//////////////////////////////////1D////////////////////////////////////////////

/// AnalyticalLinear1D: 1D -> 1D troughout given points (x0,v0) and (x1,v1)
//  Punktsteigungsform
template <typename T, typename S>
class AnalyticalLinear1D : public AnalyticalF1D<T,S> {
private:
  T _a;
  T _b;
public:
  AnalyticalLinear1D(T a, T b);
  AnalyticalLinear1D(S x0, T v0, S x1, T v1);
  bool operator() (T output[], const S x[]) override; ///< returns line _a*x + _b

  template<typename V, typename U>
  using exchange_type = AnalyticalLinear1D<V,U>;

  template<typename V, typename U>
  auto copyAs() const {
    return exchange_type<V,U>(_a, _b);
  }
};


/// represents an inverse parabola profile like it is used in Poiseuille inflow
/// note: output depends only on first parameter, maps 1D,2D,3D->1D
template <typename T, typename S>
class AnalyticalSquare1D : public AnalyticalF1D<T,S> {
private:
  S _cp;
  S _r;
  T _maxi;
public:
  AnalyticalSquare1D(S cp, S r, T maxi);
  bool operator() (T output[], const S x[]) override;
};


/// SinusStartScale: 1D -> 1D a start curve based on sinus for a continuous transition at 0 and 1
template <typename T, typename S>
class SinusStartScale : public AnalyticalF1D<T,S> {
protected:
  S _numTimeSteps;
  T _maxValue;
public:
  SinusStartScale(int numTimeSteps=1, T maxValue=1);
  bool operator() (T output[], const S x[]) override;
};


/// PolynomialStartScale: 1D -> 1D a start curve based on a polynomial fifth order for a continuous transition at 0 and 1: maxValue*(6*y^5-15*y^4+10*y^3)
template <typename T, typename S>
class PolynomialStartScale : public AnalyticalF1D<T,S> {
protected:
  S _numTimeSteps;
  T _maxValue;
public:
  PolynomialStartScale(S numTimeSteps=S(1), T maxValue=T(1));
  bool operator() (T output[], const S x[]) override;
};

/// Sinus: Sinus with period and amplitude
template <typename T, typename S>
class Sinus : public AnalyticalF1D<T,S> {
protected:
  T _period;
  T _amplitude;
public:
  Sinus (T period=1, T amplitude=1);
  bool operator() (T output[], const S x[]) override;
};

/// Cosinus: Cosinus with period and amplitude
template <typename T, typename S>
class Cosinus : public AnalyticalF1D<T,S> {
protected:
  T _period;
  T _amplitude;
public:
  Cosinus(T period=1, T amplitude=1);
  bool operator() (T output[], const S x[]) override;

  template<typename V, typename U>
  using exchange_type = Cosinus<V,U>;

  template<typename V, typename U>
  auto copyAs() const {
    return exchange_type<V,U>(_period, _amplitude);
  }
};

/// CosinusComposite: Composition of two Cosinus to shift the low point within a period - difference denotes the share of the period in which the low point is located. Calculated with case discrimination (x%period < d or d <= x%period)
template <typename T, typename S>
class CosinusComposite : public AnalyticalF1D<T,S> {
protected:
  T _period;
  T _difference;
  T _amplitude;
public:
  CosinusComposite(T period=1, T amplitude=1, T difference = 1);
  bool operator() (T output[], const S x[]) override;
};



//////////////////////////////////2D////////////////////////////////////////////

/// AnalyticalLinear2D: 2D -> 1D troughout given points (x0,y0,v0), (x1,y1,v1), (x2,y2,v2)
template <typename T, typename S>
class AnalyticalLinear2D final : public AnalyticalF2D<T,S> {
protected:
  T _a;
  T _b;
  T _c;
public:
  AnalyticalLinear2D(T a, T b, T c);
  AnalyticalLinear2D(S x0, S y0, T v0, S x1, S y1, T v1, S x2, S y2, T v2);
  bool operator() (T output[], const S x[]) override;

  template<typename V, typename U>
  using exchange_type = AnalyticalLinear2D<V,U>;

  template<typename V, typename U>
  auto copyAs() const {
    return exchange_type<V,U>(_a, _b, _c);
  }
};

/// AnalyticalRandom2D: 2D -> 1D with maxValue in the center decreasing linearly with the distrance to the center to zero at the radius and zero outside
template <typename T, typename S>
class AnalyticalParticleAdsorptionLinear2D final : public AnalyticalF2D<T,S> {
protected:
  T _center[2];
  T _radius;
  T _maxValue;
public:
  AnalyticalParticleAdsorptionLinear2D(T center[], T radius, T maxValue);
  bool operator() (T output[], const S x[]);
};

//TODO: to be removed
/** Computes resulting velocity of an object from translational and rotational velocity.
 * \param indicator Class defining the object (needs to be a SmoothIndicatorF2D<T,T,true>)
 * \param u translational velocity of the object - expected in lattice units
 * \param omega rotational velocity of the object - expected in lattice units
 */
template <typename T, typename S, typename DESCRIPTOR>
class ParticleU2D : public AnalyticalF2D<T,S> {
protected:
  SmoothIndicatorF2D<T,T,true>& _indicator;
  UnitConverter<T,DESCRIPTOR> const& _converter;
public:
  ParticleU2D(SmoothIndicatorF2D<T,T,true>& indicator, UnitConverter<T,DESCRIPTOR> const& converter);
  bool operator()(T output[], const S input[]) override;
};


//////////////////////////////////3D////////////////////////////////////////////
/// AnalyticalLinear3D: 3D -> 1D troughout given points (x0,y0,z0,v0), (x1,y1,z1,v1), (x2,y2,z2,v2), (x3,y3,z3,v3)
template <typename T, typename S>
class AnalyticalLinear3D final : public AnalyticalF3D<T,S> {
protected:
  T _a;
  T _b;
  T _c;
  T _d;
public:
  AnalyticalLinear3D(T a, T b, T c, T d);
  AnalyticalLinear3D(S x0, S y0, S z0, T v0, S x1, S y1, S z1, T v1, S x2, S y2,
                     S z2, T v2, S x3, S y3, S z3, T v3);
  bool operator() (T output[], const S x[]) override;
};

/// AnalyticalScaled3D: 3D -> Image(AnalyticalF) scales AnalyticalF by _scale
template <typename T, typename S>
class AnalyticalScaled3D final : public AnalyticalF3D<T,S> {
private:
  AnalyticalF3D<T,S>& _f;
  T _scale;
public:
  AnalyticalScaled3D(AnalyticalF3D<T,S>& f, T scale);
  bool operator() (T output[], const S x[]) override;
};

/// see Mink et al. 2016 in Sec.3.1.
template <typename T, typename S, typename DESCRIPTOR>
class PLSsolution3D : public AnalyticalF3D<T,S> {
private:
  T _physSigmaEff;
  T _physDiffusionCoefficient;
public:
  PLSsolution3D(RadiativeUnitConverter<T,DESCRIPTOR> const& converter);
  bool operator()(T output[1], const S x[3]) override;
};

/// light source as a cylinder along z-axis
template <typename T, typename S, typename DESCRIPTOR>
class LightSourceCylindrical3D : public AnalyticalF3D<T,S> {
private:
  T _physSigmaEff;
  T _physDiffusionCoefficient;
  Vector<T,3> _center;
public:
  LightSourceCylindrical3D(RadiativeUnitConverter<T,DESCRIPTOR> const& converter, Vector<T,3> center = {T(0), T(0), T(0)});
  bool operator()(T output[1], const S x[3]) override;
};

/**
* \param _position of light source
* \param _orientation direction of light source (normalized)
* \param _falloff is power of the cosine
*/
template <typename T, typename S>
class Spotlight : public AnalyticalF3D<T,S> {
private:
  Vector<T,3> const _position;
  Vector<T,3> const _orientation;
  T const _falloff;
public:
  Spotlight(Vector<T,3> position, Vector<T,3> direction, T falloff);
  bool operator()(T output[1], const S x[3]) override;
};

/// 8.6.1 Gauss Hill inital values
template <typename T, typename S>
class GaussianHill2D : public AnalyticalF2D<T,S> {
private:
  T _sigma;
  Vector<T,2> _x0;
  T _c0;
public:
  GaussianHill2D(T sigma, Vector<T,2> x0, T c0);
  bool operator()(T output[1], const S x[2]) override;
};

/// 8.6.1 Gauss Hill time evolution
template <typename T, typename S>
class GaussianHillTimeEvolution2D : public AnalyticalF2D<T,S> {
private:
  T _sigma02;
  T _D;
  T _t;
  Vector<T,2> _x0;
  Vector<T,2> _u;
  T _c0;
public:
  GaussianHillTimeEvolution2D(T sigma0, T D, T t, Vector<T,2> x0, Vector<T,2> u, T c0);
  bool operator()(T output[1], const S x[2]) override;
};


//TODO: to be removed
/** Computes resulting velocity of an object from translational and rotational velocity.
 * \param indicator Class defining the object (needs to be a SmoothIndicatorF3D)
 * \param u translational velocity of the object - expected in lattice units
 * \param omega rotational velocity of the object - expected in lattice units
 */
template <typename T, typename S, typename DESCRIPTOR>
class ParticleU3D : public AnalyticalF3D<T,S> {
protected:
  SmoothIndicatorF3D<T, T, true>& _indicator;
  UnitConverter<T,DESCRIPTOR> const& _converter;
public:
  ParticleU3D(SmoothIndicatorF3D<T, T, true>& indicator, UnitConverter<T,DESCRIPTOR> const& converter);
  bool operator()(T output[], const S input[]) override;
};

/** Returns a constant value on every cuboids.
 * The cuboid decomposition is independent of the simulation geometry.
 * SuperStructure etc. are only needed because our functor stucture requires
 * this.
 * Scaling of this functor is problematic because of a primitive underlying
 * search algorithm. Hence, this functor should not be applied in any
 * performance critical context.
 */
template <unsigned D, typename T, typename S>
class AnalyticalCuboidwiseConst final : public AnalyticalF<D,T,S>{
private:
  const unsigned _numberOfCuboids;

  std::shared_ptr<const CuboidGeometry<T,D>>  _cuboids;
  const std::vector<T>                        _values;

public:
  AnalyticalCuboidwiseConst(SuperGeometry<T,D>& sGeometry,
    const std::vector<T>& values, unsigned targetDim=D)
  : AnalyticalF<D,T,S>(targetDim),
    _numberOfCuboids(values.size()),
    _cuboids(std::make_shared<CuboidGeometry<T,D>>(
      sGeometry.getCuboidGeometry().getMotherCuboid(),
      _numberOfCuboids)),
    _values(values)
  {
    this->getName() = "cuboidwiseConst";
  }

  bool operator() (T output[], const S input[]) override {
    const auto index = _cuboids->get_iC(input[0], input[1], input[2]);
    for (int i = 0; i < this->getTargetDim(); ++i) {
      output[i] = _values[index];
    }
    return true;
  }
};

} // end namespace olb
#endif
