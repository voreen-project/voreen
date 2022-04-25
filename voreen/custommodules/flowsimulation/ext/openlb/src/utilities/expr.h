/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef UTILITIES_EXPR_H
#define UTILITIES_EXPR_H

#include <string>
#include <sstream>
#include <memory>
#include <variant>

#include "core/meta.h"
#include "core/cellD.h"

namespace olb {

/// Basic value-substitute enabling extraction of expression trees for code generation
class Expr : public ExprBase {
public:
  enum struct Op {
    Add, Mul, Sub, Div, Sqrt, Abs
  };

private:
  struct Symbol {
    std::string name;
    std::string describe() const {
      return name;
    }
  };

  struct Constant {
    double value;
    std::string describe() const {
      return std::to_string(value);
    }
  };

  class Binary {
  private:
    std::shared_ptr<Expr> lhs;
    Op op;
    std::shared_ptr<Expr> rhs;

  public:
    Binary(Expr lhs, Op op, Expr rhs):
      lhs{std::make_shared<Expr>(lhs)},
      op{op},
      rhs{std::make_shared<Expr>(rhs)} { }

    std::string describe() const {
      switch (op) {
      case Op::Add: return "(" + lhs->describe() + ") + (" + rhs->describe() + ")";
      case Op::Sub: return "(" + lhs->describe() + ") - (" + rhs->describe() + ")";
      case Op::Mul: return "(" + lhs->describe() + ") * (" + rhs->describe() + ")";
      case Op::Div: return "(" + lhs->describe() + ") / (" + rhs->describe() + ")";
      default: throw std::invalid_argument("Unsupported binary operation");
      }
    }
  };

  class Unary {
  private:
    Op op;
    std::shared_ptr<Expr> arg;

  public:
    Unary(Op op, Expr arg):
      op{op},
      arg{std::make_shared<Expr>(arg)} { }

    std::string describe() const {
      switch (op) {
      case Op::Sqrt: return "sqrt(" + arg->describe() + ")";
      case Op::Abs:  return "Abs(" + arg->describe() + ")";
      default: throw std::invalid_argument("Unsupported unary operation");
      }
    }
  };

  std::variant<Symbol,Constant,Binary,Unary> _payload;

public:
  Expr(const Expr& rhs):
    _payload(rhs._payload) { }

  Expr(double v):
    _payload(Constant{v}) { }

  Expr():
    Expr(double{}) { }

  Expr(std::string name):
    _payload(Symbol{name}) { }

  Expr(Expr lhs, Op op, Expr rhs):
    _payload(Binary(lhs, op, rhs)) { }

  Expr(Op op, Expr rhs):
    _payload(Unary(op, rhs)) { }

  std::string describe() const {
    std::string out;
    std::visit([&out](auto& x) {
      out = x.describe();
    }, _payload);
    return out;
  };

  Expr& operator+=(Expr rhs);
  Expr& operator-=(Expr rhs);
  Expr& operator*=(Expr rhs);
  Expr& operator/=(Expr rhs);
};

Expr& Expr::operator+=(Expr rhs) {
  _payload = Binary(*this, Op::Add, rhs);
  return *this;
}

Expr& Expr::operator-=(Expr rhs) {
  _payload = Binary(*this, Op::Sub, rhs);
  return *this;
}

Expr& Expr::operator*=(Expr rhs) {
  _payload = Binary(*this, Op::Mul, rhs);
  return *this;
}

Expr& Expr::operator/=(Expr rhs) {
  _payload = Binary(*this, Op::Div, rhs);
  return *this;
}

Expr operator+(Expr lhs, Expr rhs) {
  return Expr(lhs, Expr::Op::Add, rhs);
}

Expr operator-(Expr lhs, Expr rhs) {
  return Expr(lhs, Expr::Op::Sub, rhs);
}

Expr operator*(Expr lhs, Expr rhs) {
  return Expr(lhs, Expr::Op::Mul, rhs);
}

Expr operator/(Expr lhs, Expr rhs) {
  return Expr(lhs, Expr::Op::Div, rhs);
}

Expr operator-(Expr rhs) {
  return Expr(-1.) * rhs;
}

namespace util {

Expr sqrt(Expr x) {
  return Expr(Expr::Op::Sqrt, x);
}

Expr fabs(Expr x) {
  return Expr(Expr::Op::Abs, x);
}

}

}

#endif
