/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 202021 Max Sagebaum, Nicolas Gauger
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

#ifndef A_DIFF_TAPE_H
#define A_DIFF_TAPE_H

#include <cstdint>
#include <vector>

namespace olb {

namespace util {

namespace ADtape {

template<typename T_Tape>
struct ActiveType {
  using Tape = T_Tape;
  using Real = typename Tape::Real;
  using Gradient = Real;

private:
  static Tape globalTape;

  Real value_;
  int id_;

public:
  ActiveType() = default;

  ActiveType(ActiveType const& value) = default; // We do not need to store copy operations

  ActiveType(Real const& value) : value_(value), id_(0) {}

  ActiveType& operator=(Real const& value) {
    value_ = value;
    id_ = 0;

    return *this;
  }

  ActiveType& operator =(ActiveType const& value) = default; // We do not need to store copy operations

#define OPERATOR(OP, rhs) \
  ActiveType& OP(ActiveType const& value) { \
    *this = rhs; \
    return *this; \
  }

  OPERATOR(operator +=, *this + value)
  OPERATOR(operator -=, *this - value)
  OPERATOR(operator *=, *this * value)
  OPERATOR(operator /=, *this / value)

#undef OPERATOR

  Real& value() {
    return value_;
  }

  Real const& value() const {
    return value_;
  }

  int& getIdentifier() {
    return id_;
  }

  int const& getIdentifier() const {
    return id_;
  }

  Real& gradient() {
    return getTape().gradient(id_);
  }

  Real const& gradient() const {
    return getTape().gradient(id_);
  }

  static Tape& getTape() {
    return globalTape;
  }
};

template<typename Tape>
Tape ActiveType<Tape>::globalTape;

#define BINARY_OPERATOR(OP, VAL, JAC1, JAC2) \
template<typename Tape> \
ActiveType<Tape> OP(ActiveType<Tape> const& arg1, ActiveType<Tape> const& arg2) { \
  using Real = typename Tape::Real; \
  \
  Real val1 = arg1.value(); \
  Real val2 = arg2.value(); \
  ActiveType<Tape> lhs; \
  ActiveType<Tape>::getTape().store(lhs, VAL, arg1, JAC1, arg2, JAC2); \
  \
  return lhs; \
} \
\
template<typename Tape> \
ActiveType<Tape> OP(typename Tape::Real const& val1, ActiveType<Tape> const& arg2) { \
  using Real = typename Tape::Real; \
  \
  Real val2 = arg2.value(); \
  ActiveType<Tape> lhs; \
  ActiveType<Tape>::getTape().store(lhs, VAL, arg2, JAC2); \
  \
  return lhs; \
} \
\
template<typename Tape> \
ActiveType<Tape> OP(ActiveType<Tape> const& arg1, typename Tape::Real const& val2) { \
  using Real = typename Tape::Real; \
  \
  Real val1 = arg1.value(); \
  ActiveType<Tape> lhs; \
  ActiveType<Tape>::getTape().store(lhs, VAL, arg1, JAC1); \
  \
  return lhs; \
}

using std::atan2;
using std::hypot;
using std::max;
using std::min;
using std::pow;

BINARY_OPERATOR(operator+, val1 + val2, 1.0, 1.0)
BINARY_OPERATOR(operator-, val1 - val2, 1.0, -1.0)
BINARY_OPERATOR(operator*, val1 * val2, val2, val1)
BINARY_OPERATOR(operator/, val1 / val2, 1.0 / val2, -val1 / val2 / val2)
BINARY_OPERATOR(atan2, atan2(val1, val2), val2 / (val1 * val1 + val2 * val2), -val1 / (val1 * val1 + val2 * val2))
BINARY_OPERATOR(hypot, hypot(val1, val2), val1 / hypot(val1, val2), val2 / hypot(val1, val2))
BINARY_OPERATOR(max, max(val1, val2), val1 > val2 ? 1.0 : 0.0, val1 > val2 ? 0.0 : 1.0)
BINARY_OPERATOR(min, min(val1, val2), val1 < val2 ? 1.0 : 0.0, val1 < val2 ? 0.0 : 1.0)
BINARY_OPERATOR(pow, pow(val1, val2), val2 * pow(val1, val2 - 1.0), log(val1) * pow(val1, val2))
#undef BINARY_OPERATOR

#define UNARY_OPERATOR(OP, VAL, JAC) \
template<typename Tape> \
ActiveType<Tape> OP(ActiveType<Tape> const& arg) { \
  using Real = typename Tape::Real; \
  \
  Real val = arg.value(); \
  ActiveType<Tape> lhs; \
  ActiveType<Tape>::getTape().store(lhs, VAL, arg, JAC); \
  \
  return lhs; \
}

using std::abs;
using std::acos;
using std::asin;
using std::atan;
using std::atanh;
using std::cbrt;
using std::cos;
using std::cosh;
using std::exp;
using std::log;
using std::log10;
using std::sin;
using std::sinh;
using std::sqrt;
using std::tan;
using std::tanh;

UNARY_OPERATOR(operator+, val, 1.0)
UNARY_OPERATOR(operator-, -val, -1.0)
UNARY_OPERATOR(abs, abs(val), val == 0.0 ? 0.0 : (val > 0.0 ? 1.0 : -1.0))
UNARY_OPERATOR(acos, acos(val), -1.0 / sqrt(1.0 - val * val))
UNARY_OPERATOR(asin, asin(val), 1.0 / sqrt(1.0 - val * val))
UNARY_OPERATOR(atan, atan(val), 1.0 / (1.0 + val * val))
UNARY_OPERATOR(atanh, atanh(val), 1.0 / (1.0 - val * val))
UNARY_OPERATOR(cbrt, cbrt(val), 1.0 / (3.0 * cbrt(val) * cbrt(val)))
UNARY_OPERATOR(cos, cos(val), -sin(val))
UNARY_OPERATOR(cosh, cosh(val), sinh(val))
UNARY_OPERATOR(exp, exp(val), exp(val))
UNARY_OPERATOR(log, log(val), 1.0 / val)
UNARY_OPERATOR(log10, log10(val), 0.434294481903252 / val)
UNARY_OPERATOR(sin, sin(val), cos(val))
UNARY_OPERATOR(sinh, sinh(val), cosh(val))
UNARY_OPERATOR(sqrt, sqrt(val), 0.5 / sqrt(val))
UNARY_OPERATOR(tan, tan(val), 1.0 / cos(val) / cos(val))
UNARY_OPERATOR(tanh, tanh(val), 1.0 - tanh(val) * tanh(val))
#undef UNARY_OPERATOR

#define CONDITION_OPERATOR(OP) \
template<typename Tape> \
bool operator OP(ActiveType<Tape> const& arg1, ActiveType<Tape> const& arg2) { \
  return arg1.value() OP arg2.value(); \
}

CONDITION_OPERATOR(==)
CONDITION_OPERATOR(!=)
CONDITION_OPERATOR(>)
CONDITION_OPERATOR(<)
CONDITION_OPERATOR(>=)
CONDITION_OPERATOR(<=)
CONDITION_OPERATOR(&&)
CONDITION_OPERATOR(||)

#undef CONDITION_OPERATOR

#define CONDITION_OPERATOR(OP) \
template<typename Tape> \
bool operator OP(ActiveType<Tape> const& arg) { \
  return OP arg.value(); \
}

using std::swap;

template<typename Tape>
void swap(ActiveType<Tape>& arg1, ActiveType<Tape>& arg2) {
  swap(arg1.value(), arg2.value());
  swap(arg1.getIdentifier(), arg2.getIdentifier());
}

CONDITION_OPERATOR(!)
#undef CONDITION_OPERATOR

template<typename T_Real>
struct Tape {
public:
  using Real = T_Real;
  using OpType = uint8_t;

  int maxIdentifier;

  std::vector<Real> jacobianVector;
  std::vector<int> idVector;
  std::vector<OpType> opVector;

  std::vector<Real> adjoints;

private:
  void pushJacobian(int const& id, Real const& jacobi) {
    jacobianVector.push_back(jacobi);
    idVector.push_back(id);
  }

  void pushStmt(OpType const& opType) {
    opVector.push_back(opType);
  }

  void assignValue(ActiveType<Tape>& lhs, Real const& lhsValue) {
    maxIdentifier += 1;

    lhs.value() = lhsValue;
    lhs.getIdentifier() = maxIdentifier;
  }

  void resizeAdjointVector() {
    if((int)adjoints.size() <= maxIdentifier) {
      adjoints.resize(maxIdentifier + 1);
    }
  }

public:
  void store(ActiveType<Tape>& lhs, Real const& lhsValue) {
    // Register input always store a statement
    pushStmt(0);

    assignValue(lhs, lhsValue);
  }

  void store(ActiveType<Tape>& lhs, Real const& lhsValue, ActiveType<Tape> const& rhs, Real const& rhsJac) {
    pushJacobian(rhs.getIdentifier(), rhsJac);
    pushStmt(1);

    assignValue(lhs, lhsValue);
  }

  void store(ActiveType<Tape>& lhs, Real const& lhsValue,
             ActiveType<Tape> const& rhs1, Real const& rhs1Jac,
             ActiveType<Tape> const& rhs2, Real const& rhs2Jac) {
    pushJacobian(rhs1.getIdentifier(), rhs1Jac);
    pushJacobian(rhs2.getIdentifier(), rhs2Jac);
    pushStmt(2);

    assignValue(lhs, lhsValue);
  }

  void evaluate() {

    resizeAdjointVector();

    int lhsId = maxIdentifier;

    int stmtPos = (size_t)opVector.size();
    int argPos = (size_t)jacobianVector.size();

    while(lhsId > 0) {

      stmtPos -= 1;

      OpType nArgs = opVector[stmtPos];
      Real adjoint = adjoints[lhsId];

      for(OpType i = 0; i < nArgs; i += 1) {
        argPos -= 1;
        adjoints[idVector[argPos]] += jacobianVector[argPos] * adjoint;
      }

      lhsId -= 1;
    }
  }

  void reset() {
    maxIdentifier = 0;

    jacobianVector.clear();
    idVector.clear();
    opVector.clear();

    std::fill(adjoints.begin(), adjoints.end(), Real());
  }

  void registerInput(ActiveType<Tape>& value) {
    store(value, value.value());
  }

  void registerOutput(ActiveType<Tape>& value) {
    store(value, value.value(), value, Real(1.0));
  }

  Real& gradient(int const& id) {

    resizeAdjointVector();

    if(id <= maxIdentifier) {
      return adjoints[id];
    } else {
      return adjoints[0];
    }
  }
};

using RealReverse = ActiveType<Tape<double>>;

}

}

}

#endif
