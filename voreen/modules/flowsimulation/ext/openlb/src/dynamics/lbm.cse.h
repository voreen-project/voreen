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

#ifndef DYNAMICS_LBM_CSE_H
#define DYNAMICS_LBM_CSE_H


#ifndef DISABLE_CSE

#include "lbm.h"
#include "latticeDescriptors.h"

namespace olb {


template <typename... FIELDS>
struct lbm<descriptors::D2Q5<FIELDS...>> {

template <typename CELL, typename V=typename CELL::value_t>
static auto computeRho(CELL& cell) any_platform
{

return cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
}

template <typename CELL, typename J, typename V=typename CELL::value_t>
static void computeJ(CELL& cell, J& j) any_platform
{
j[0] = -V{1}*cell[1] + V{1}*cell[3];
j[1] = -V{1}*cell[2] + V{1}*cell[4];

}

template <typename CELL, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeRhoU(CELL& cell, RHO& rho, U& u) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
auto x1 = V{1}/x0;
rho = x0;
u[0] = -x1*(cell[1] - cell[3]);
u[1] = -x1*(cell[2] - cell[4]);

}

template <typename CELL, typename RHO, typename J, typename V=typename CELL::value_t>
static void computeRhoJ(CELL& cell, RHO& rho, J& j) any_platform
{
rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
j[0] = -V{1}*cell[1] + V{1}*cell[3];
j[1] = -V{1}*cell[2] + V{1}*cell[4];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeStress(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = V{0.333333333333333} - V{0.333333333333333}*rho;
pi[0] = -rho*u[0]*u[0] + x0 + V{1}*cell[1] + V{1}*cell[3];
pi[1] = -rho*u[0]*u[1];
pi[2] = -rho*u[1]*u[1] + x0 + V{1}*cell[2] + V{1}*cell[4];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeAllMomenta(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
auto x1 = cell[1] - cell[3];
auto x2 = V{1}/x0;
auto x3 = x1*x2;
auto x4 = cell[2] - cell[4];
auto x5 = -V{0.333333333333333}*cell[0];
rho = x0;
u[0] = -x3;
u[1] = -x2*x4;
pi[0] = -x2*x1*x1 + x5 + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4];
pi[1] = -x3*x4;
pi[2] = -x2*x4*x4 + x5 - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4];

}

template <typename CELL, typename FEQ, typename V=typename CELL::value_t>
static void computeFeq(CELL& cell, FEQ& fEq) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4];
auto x1 = x0 + V{1};
auto x2 = x0 + V{1};
auto x3 = V{1} / ((x2)*(x2));
auto x4 = V{1.5}*x3;
auto x5 = -cell[2] + cell[4];
auto x6 = x5*x5;
auto x12 = -cell[1] + cell[3];
auto x13 = x4*(x12*x12) + V{-1};
auto x14 = V{1} / (x2);
auto x15 = x14*(V{3}*cell[1] - V{3}*cell[3]);
auto x16 = V{3}*x3;
auto x17 = cell[1] - cell[3];
auto x18 = x17*x17;
auto x19 = cell[2] - cell[4];
auto x20 = x19*x19;
auto x21 = x16*x18 - x20*x4 + V{1};
auto x22 = V{3}*cell[4];
auto x23 = V{3}*cell[2];
fEq[0] = -V{0.333333333333333}*x1*(x13 + x4*x6) + V{-0.333333333333333};
fEq[1] = V{0.166666666666667}*x1*(x15 + x21) + V{-0.166666666666667};
fEq[2] = -V{0.166666666666667}*x1*(x13 + x14*(x22 - x23) - x16*x6) + V{-0.166666666666667};
fEq[3] = V{0.166666666666667}*x1*(-x15 + x21) + V{-0.166666666666667};
fEq[4] = -V{0.166666666666667}*x1*(x14*(-x22 + x23) - x16*x20 + x18*x4 + V{-1}) + V{-0.166666666666667};

}

template <typename CELL, typename FNEQ, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq, RHO& rho, U& u) any_platform
{
auto x0 = u[1]*u[1];
auto x1 = V{1.5}*x0;
auto x2 = u[0]*u[0];
auto x3 = V{1.5}*x2;
auto x4 = x3 + V{-1};
auto x5 = V{0.166666666666667}*rho;
auto x6 = V{3}*u[0];
auto x12 = -x1 + V{3}*x2 + V{1};
auto x13 = V{3}*u[1];
auto x14 = V{3}*x0;
fNeq[0] = V{0.333333333333333}*rho*(x1 + x4) + cell[0] + V{0.333333333333333};
fNeq[1] = -x5*(x12 - x6) + cell[1] + V{0.166666666666667};
fNeq[2] = x5*(x13 - x14 + x4) + cell[2] + V{0.166666666666667};
fNeq[3] = -x5*(x12 + x6) + cell[3] + V{0.166666666666667};
fNeq[4] = -x5*(x13 + x14 - x3 + V{1}) + cell[4] + V{0.166666666666667};

}

template <typename CELL, typename FNEQ, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
auto x1 = V{1} / ((x0)*(x0));
auto x2 = V{1.5}*x1;
auto x3 = cell[2] - cell[4];
auto x4 = x3*x3;
auto x5 = x2*x4;
auto x6 = cell[1] - cell[3];
auto x12 = x6*x6;
auto x13 = x12*x2 + V{-1};
auto x14 = V{0.166666666666667}*cell[0] + V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667};
auto x15 = V{1} / (x0);
auto x16 = x15*(V{3}*cell[1] - V{3}*cell[3]);
auto x17 = V{3}*x1;
auto x18 = x12*x17 - x5 + V{1};
auto x19 = x15*(V{3}*cell[2] - V{3}*cell[4]);
auto x20 = x13 - x17*x4;
fNeq[0] = (x13 + x5)*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}) + cell[0] + V{0.333333333333333};
fNeq[1] = -x14*(x16 + x18) + cell[1] + V{0.166666666666667};
fNeq[2] = x14*(-x19 + x20) + cell[2] + V{0.166666666666667};
fNeq[3] = -x14*(-x16 + x18) + cell[3] + V{0.166666666666667};
fNeq[4] = x14*(x19 + x20) + cell[4] + V{0.166666666666667};

}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto bgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x5 = omega + V{-1};
auto x6 = u[1]*u[1];
auto x7 = V{1.5}*x6;
auto x8 = u[0]*u[0];
auto x9 = V{1.5}*x8;
auto x10 = x9 + V{-1};
auto x11 = V{0.166666666666667}*omega;
auto x12 = V{3}*u[0];
auto x13 = -x7 + V{3}*x8 + V{1};
auto x14 = V{3}*u[1];
auto x15 = V{3}*x6;
cell[0] = -V{0.333333333333333}*omega*(rho*(x10 + x7) + V{1}) - x5*cell[0];
cell[1] = x11*(rho*(-x12 + x13) + V{-1}) - x5*cell[1];
cell[2] = -x11*(rho*(x10 + x14 - x15) + V{1}) - x5*cell[2];
cell[3] = x11*(rho*(x12 + x13) + V{-1}) - x5*cell[3];
cell[4] = x11*(rho*(x14 + x15 - x9 + V{1}) + V{-1}) - x5*cell[4];
return x6 + x8;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto adeBgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega)
{
auto x5 = omega + V{-1};
auto x6 = V{3}*u[0];
auto x7 = V{0.166666666666667}*omega;
auto x8 = V{3}*u[1];
cell[0] = V{0.333333333333333}*omega*(rho + V{-1}) - x5*cell[0];
cell[1] = -x5*cell[1] - x7*(rho*(x6 + V{-1}) + V{1});
cell[2] = -x5*cell[2] - x7*(rho*(x8 + V{-1}) + V{1});
cell[3] = -x5*cell[3] + x7*(rho*(x6 + V{1}) + V{-1});
cell[4] = -x5*cell[4] + x7*(rho*(x8 + V{1}) + V{-1});
return u[0]*u[0] + u[1]*u[1];
}

template <typename CELL, typename PRESSURE, typename J, typename OMEGA, typename V=typename CELL::value_t>
static auto incBgkCollision(CELL& cell, PRESSURE& pressure, J& j, OMEGA& omega) any_platform
{
auto x5 = j[0]*j[0];
auto x6 = V{0.5}*x5;
auto x7 = j[1]*j[1];
auto x8 = V{0.5}*x7;
auto x9 = omega + V{-1};
auto x10 = V{0.25}*x7;
auto x11 = V{0.5}*j[0];
auto x12 = V{0.5}*pressure;
auto x13 = V{0.166666666666667} - x12;
auto x14 = V{0.25}*x5;
auto x15 = V{0.5}*j[1];
auto x16 = x12 + V{-0.166666666666667};
cell[0] = -omega*(-V{1}*pressure + x6 + x8 + V{0.333333333333333}) - x9*cell[0];
cell[1] = -omega*(x10 + x11 + x13 - x6) - x9*cell[1];
cell[2] = -omega*(x13 + x14 + x15 - x8) - x9*cell[2];
cell[3] = omega*(-x10 + x11 + x16 + x6) - x9*cell[3];
cell[4] = omega*(-x14 + x15 + x16 + x8) - x9*cell[4];
return x5 + x7;
}

template <typename CELL, typename RHO, typename U, typename RATIORHO, typename OMEGA, typename V=typename CELL::value_t>
static auto constRhoBgkCollision(CELL& cell, RHO& rho, U& u, RATIORHO& ratioRho, OMEGA& omega) any_platform
{
auto x5 = omega + V{-1};
auto x6 = V{0.333333333333333}*rho;
auto x7 = u[1]*u[1];
auto x8 = V{1.5}*x7;
auto x9 = u[0]*u[0];
auto x10 = V{1.5}*x9;
auto x11 = x10 + V{-1};
auto x12 = x11 + x8;
auto x13 = V{0.166666666666667}*rho;
auto x14 = V{3}*u[0];
auto x15 = -x8 + V{3}*x9 + V{1};
auto x16 = -x14 + x15;
auto x17 = ratioRho*x13;
auto x18 = V{3}*u[1];
auto x19 = V{3}*x7;
auto x20 = x11 + x18 - x19;
auto x21 = x14 + x15;
auto x22 = -x10 + x18 + x19 + V{1};
cell[0] = -ratioRho*x12*x6 - x5*(x12*x6 + cell[0] + V{0.333333333333333}) + V{-0.333333333333333};
cell[1] = x16*x17 - x5*(-x13*x16 + cell[1] + V{0.166666666666667}) + V{-0.166666666666667};
cell[2] = -x17*x20 - x5*(x13*x20 + cell[2] + V{0.166666666666667}) + V{-0.166666666666667};
cell[3] = x17*x21 - x5*(-x13*x21 + cell[3] + V{0.166666666666667}) + V{-0.166666666666667};
cell[4] = x17*x22 - x5*(-x13*x22 + cell[4] + V{0.166666666666667}) + V{-0.166666666666667};
return x7 + x9;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x5 = omega + V{-1};
auto x6 = V{3}*u[0];
auto x7 = V{0.0833333333333333}*rho;
auto x8 = x6 + V{1};
auto x9 = -x7*x8 - V{0.5}*cell[1] + V{0.5}*cell[3];
auto x10 = x6 + V{-1};
auto x11 = V{0.166666666666667}*rho;
auto x12 = V{3}*u[1];
auto x13 = x12 + V{1};
auto x14 = -x13*x7 - V{0.5}*cell[2] + V{0.5}*cell[4];
auto x15 = x12 + V{-1};
cell[0] = V{0.333333333333333}*rho + V{-0.333333333333333};
cell[1] = -x10*x11 + x5*(x7*(V{1} - x6) + x9) + V{-0.166666666666667};
cell[2] = -x11*x15 + x5*(x14 + x7*(V{1} - x12)) + V{-0.166666666666667};
cell[3] = x11*x8 - x5*(-x10*x7 + x9) + V{-0.166666666666667};
cell[4] = x11*x13 - x5*(x14 - x15*x7) + V{-0.166666666666667};
return u[0]*u[0] + u[1]*u[1];
}

template <typename CELL, typename RHO, typename U, typename PI, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, PI& pi, OMEGA& omega) any_platform
{
auto x5 = omega + V{-1};
auto x6 = u[1]*u[1];
auto x7 = V{1.5}*x6;
auto x8 = u[0]*u[0];
auto x9 = V{1.5}*x8;
auto x10 = x9 + V{-1};
auto x11 = V{0.166666666666667}*rho;
auto x12 = V{3}*u[0];
auto x13 = -x7 + V{3}*x8 + V{1};
auto x14 = -x5*(V{0.5}*pi[0] - V{0.25}*pi[2]) + V{-0.166666666666667};
auto x15 = V{3}*u[1];
auto x16 = V{3}*x6;
auto x17 = x5*(V{0.25}*pi[0] - V{0.5}*pi[2]) + V{-0.166666666666667};
cell[0] = -V{0.333333333333333}*rho*(x10 + x7) + V{0.5}*x5*(pi[0] + pi[2]) + V{-0.333333333333333};
cell[1] = x11*(-x12 + x13) + x14;
cell[2] = -x11*(x10 + x15 - x16) + x17;
cell[3] = x11*(x12 + x13) + x14;
cell[4] = x11*(x15 + x16 - x9 + V{1}) + x17;
return x6 + x8;
}

template <typename CELL, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineEqFirstOrder(CELL& cell, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x5 = V{3}*newU[0];
auto x6 = V{3}*newU[1];
cell[0] = V{0.333333333333333}*newRho + V{-0.333333333333333};
cell[1] = -V{0.166666666666667}*newRho*(x5 + V{-1}) + V{-0.166666666666667};
cell[2] = -V{0.166666666666667}*newRho*(x6 + V{-1}) + V{-0.166666666666667};
cell[3] = V{0.166666666666667}*newRho*(x5 + V{1}) + V{-0.166666666666667};
cell[4] = V{0.166666666666667}*newRho*(x6 + V{1}) + V{-0.166666666666667};

}

template <typename CELL, typename OLDRHO, typename OLDU, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineNEq(CELL& cell, OLDRHO& oldRho, OLDU& oldU, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x5 = oldU[1]*oldU[1];
auto x6 = V{1.5}*x5;
auto x7 = oldU[0]*oldU[0];
auto x8 = V{1.5}*x7;
auto x9 = x8 + V{-1};
auto x10 = newU[1]*newU[1];
auto x11 = V{1.5}*x10;
auto x12 = newU[0]*newU[0];
auto x13 = V{1.5}*x12;
auto x14 = x13 + V{-1};
auto x15 = V{0.166666666666667}*newRho;
auto x16 = V{3}*newU[0];
auto x17 = -x11 + V{3}*x12 + V{1};
auto x18 = V{0.166666666666667}*oldRho;
auto x19 = V{3}*oldU[0];
auto x20 = -x6 + V{3}*x7 + V{1};
auto x21 = V{3}*oldU[1];
auto x22 = V{3}*x5;
auto x23 = V{3}*newU[1];
auto x24 = V{3}*x10;
cell[0] = -V{0.333333333333333}*newRho*(x11 + x14) + V{0.333333333333333}*oldRho*(x6 + x9) + cell[0];
cell[1] = x15*(-x16 + x17) - x18*(-x19 + x20) + cell[1];
cell[2] = -x15*(x14 + x23 - x24) + x18*(x21 - x22 + x9) + cell[2];
cell[3] = x15*(x16 + x17) - x18*(x19 + x20) + cell[3];
cell[4] = x15*(-x13 + x23 + x24 + V{1}) - x18*(x21 + x22 - x8 + V{1}) + cell[4];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void defineNEqFromPi(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x5 = V{0.5}*pi[0];
auto x6 = V{0.5}*pi[2];
auto x7 = u[1]*u[1];
auto x8 = V{1.5}*x7;
auto x9 = u[0]*u[0];
auto x10 = V{1.5}*x9;
auto x11 = x10 + V{-1};
auto x12 = V{0.166666666666667}*rho;
auto x13 = V{3}*u[0];
auto x14 = -x8 + V{3}*x9 + V{1};
auto x15 = x5 - V{0.25}*pi[2] + V{-0.166666666666667};
auto x16 = V{3}*u[1];
auto x17 = V{3}*x7;
auto x18 = x6 - V{0.25}*pi[0] + V{-0.166666666666667};
cell[0] = -V{0.333333333333333}*rho*(x11 + x8) - x5 - x6 + V{-0.333333333333333};
cell[1] = x12*(-x13 + x14) + x15;
cell[2] = -x12*(x11 + x16 - x17) + x18;
cell[3] = x12*(x13 + x14) + x15;
cell[4] = x12*(-x10 + x16 + x17 + V{1}) + x18;

}

template <typename CELL, typename FORCE, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell, FORCE& force) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4];
auto x1 = x0 + V{1};
auto x2 = cell[1] - cell[3];
auto x3 = cell[2] - cell[4];
auto x4 = x0 + V{1};
auto x5 = x4*(x2*force[1] + x3*force[0]);
auto x6 = x2*x3;
auto x7 = V{0.333333333333333}*cell[0];
auto x8 = V{1} / (x1);
auto x9 = x4*x8;
auto x10 = x3*x9*force[1] + x7 + x8*(x3*x3) + V{0.333333333333333}*cell[1] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[3] - V{0.666666666666667}*cell[4];
auto x11 = x2*x9*force[0] + x7 + x8*(x2*x2) - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[2] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4];
return x10*x10 + x11*x11 + (V{0.5}*x5 + V{1}*x6)*(V{1}*x5 + V{2}*x6)/((x1)*(x1));
}

template <typename CELL, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
auto x1 = cell[1] - cell[3];
auto x2 = x1*x1;
auto x3 = cell[2] - cell[4];
auto x4 = x3*x3;
auto x5 = V{0.333333333333333}*cell[0];
auto x6 = V{1}/x0;
auto x7 = x4*x6 + x5 + V{0.333333333333333}*cell[1] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[3] - V{0.666666666666667}*cell[4];
auto x8 = x2*x6 + x5 - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[2] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4];
return x7*x7 + x8*x8 + V{2}*x2*x4/((x0)*(x0));
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addExternalForce(CELL& cell, RHO& rho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x5 = force[0]*u[0];
auto x6 = force[1]*u[1];
auto x7 = rho*(V{0.5}*omega + V{-1});
auto x8 = V{6}*u[0];
auto x9 = V{0.166666666666667}*force[0];
auto x10 = -V{0.5}*x6;
auto x11 = V{6}*u[1];
auto x12 = V{0.166666666666667}*force[1];
auto x13 = -V{0.5}*x5;
cell[0] = V{1}*x7*(x5 + x6) + cell[0];
cell[1] = -x7*(x10 + x9*(x8 + V{-3})) + cell[1];
cell[2] = -x7*(x12*(x11 + V{-3}) + x13) + cell[2];
cell[3] = -x7*(x10 + x9*(x8 + V{3})) + cell[3];
cell[4] = -x7*(x12*(x11 + V{3}) + x13) + cell[4];

}

};

template <typename... FIELDS>
struct lbm<descriptors::D2Q9<FIELDS...>> {

template <typename CELL, typename V=typename CELL::value_t>
static auto computeRho(CELL& cell) any_platform
{

return cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
}

template <typename CELL, typename J, typename V=typename CELL::value_t>
static void computeJ(CELL& cell, J& j) any_platform
{
auto x0 = -cell[3] + cell[7];
j[0] = V{1}*x0 - V{1}*cell[1] - V{1}*cell[2] + V{1}*cell[5] + V{1}*cell[6];
j[1] = V{1}*x0 + V{1}*cell[1] - V{1}*cell[4] - V{1}*cell[5] + V{1}*cell[8];

}

template <typename CELL, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeRhoU(CELL& cell, RHO& rho, U& u) any_platform
{
auto x0 = cell[2] + cell[3];
auto x1 = cell[7] + cell[8];
auto x2 = x0 + x1 + cell[0] + cell[1] + cell[4] + cell[5] + cell[6] + V{1};
auto x3 = cell[1] - cell[5];
auto x4 = V{1}/x2;
rho = x2;
u[0] = -x4*(x0 + x3 - cell[6] - cell[7]);
u[1] = x4*(x1 + x3 - cell[3] - cell[4]);

}

template <typename CELL, typename RHO, typename J, typename V=typename CELL::value_t>
static void computeRhoJ(CELL& cell, RHO& rho, J& j) any_platform
{
auto x0 = cell[5] + cell[6];
auto x1 = cell[1] + cell[8];
auto x2 = -cell[3] + cell[7];
rho = x0 + x1 + cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + V{1};
j[0] = V{1}*x0 + V{1}*x2 - V{1}*cell[1] - V{1}*cell[2];
j[1] = V{1}*x1 + V{1}*x2 - V{1}*cell[4] - V{1}*cell[5];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeStress(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = V{1}*cell[1];
auto x1 = V{1}*cell[5];
auto x2 = V{1}*cell[3] + V{1}*cell[7];
auto x3 = -V{0.333333333333333}*rho + x0 + x1 + x2 + V{0.333333333333333};
pi[0] = -rho*u[0]*u[0] + x3 + V{1}*cell[2] + V{1}*cell[6];
pi[1] = -rho*u[0]*u[1] - x0 - x1 + x2;
pi[2] = -rho*u[1]*u[1] + x3 + V{1}*cell[4] + V{1}*cell[8];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeAllMomenta(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = cell[1] + cell[2];
auto x1 = cell[7] + cell[8];
auto x2 = x0 + x1 + cell[0] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
auto x3 = -cell[5];
auto x4 = x3 + cell[3];
auto x5 = x0 + x4 - cell[6] - cell[7];
auto x6 = V{1} / (x2);
auto x7 = V{1}*x6;
auto x8 = x1 + x3 + cell[1] - cell[3] - cell[4];
auto x9 = -V{0.333333333333333}*cell[0] + V{0.666666666666667}*cell[1] + V{0.666666666666667}*cell[3] + V{0.666666666666667}*cell[5] + V{0.666666666666667}*cell[7];
rho = x2;
u[0] = -x5*x7;
u[1] = x7*x8;
pi[0] = -x7*x5*x5 + x9 + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[4] + V{0.666666666666667}*cell[6] - V{0.333333333333333}*cell[8];
pi[1] = V{1}*x4 + V{1}*x5*x6*x8 - V{1}*cell[1] + V{1}*cell[7];
pi[2] = -x7*x8*x8 + x9 - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[4] - V{0.333333333333333}*cell[6] + V{0.666666666666667}*cell[8];

}

template <typename CELL, typename FEQ, typename V=typename CELL::value_t>
static void computeFeq(CELL& cell, FEQ& fEq) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x1 = x0 + V{1};
auto x2 = x0 + V{1};
auto x3 = V{1} / ((x2)*(x2));
auto x4 = V{1.5}*x3;
auto x5 = -cell[3] + cell[7];
auto x6 = -cell[4] + cell[8];
auto x7 = cell[1] - cell[5];
auto x8 = x5 + x6 + x7;
auto x9 = x8*x8;
auto x10 = x4*x9;
auto x20 = -cell[2] + cell[6];
auto x21 = x20 + x5 - cell[1] + cell[5];
auto x22 = x4*(x21*x21) + V{-1};
auto x23 = x10 + x22;
auto x24 = V{1} / (x2);
auto x25 = V{3}*cell[2];
auto x26 = V{3}*cell[3];
auto x27 = V{3}*cell[6];
auto x28 = V{3}*cell[7];
auto x29 = V{3}*cell[1];
auto x30 = V{3}*cell[5];
auto x31 = x29 - x30;
auto x32 = x24*(x25 + x26 - x27 - x28 + x31);
auto x33 = V{4.5}*x3;
auto x34 = cell[2] - cell[6];
auto x35 = x34 + x6 + V{2}*cell[1] - V{2}*cell[5];
auto x36 = x33*(x35*x35);
auto x37 = V{1} - x10;
auto x38 = -x26 + x28;
auto x39 = x24*(x31 + x38 - V{3}*cell[4] + V{3}*cell[8]);
auto x40 = x34 + x7 + cell[3] - cell[7];
auto x41 = x40*x40;
auto x42 = x4*x41;
auto x43 = x39 - x42;
auto x44 = x37 + x43;
auto x45 = V{3}*x3;
auto x46 = x37 + x41*x45;
auto x47 = V{2}*cell[3];
auto x48 = V{2}*cell[7];
auto x49 = x20 - x47 + x48 + x6;
auto x50 = x45*x9;
auto x51 = -x32;
auto x52 = x34 + x47 - x48 + cell[4] - cell[8];
fEq[0] = -V{0.444444444444444}*x1*x23 + V{-0.444444444444444};
fEq[1] = V{0.0277777777777778}*x1*(x32 + x36 + x44) + V{-0.0277777777777778};
fEq[2] = V{0.111111111111111}*x1*(x32 + x46) + V{-0.111111111111111};
fEq[3] = -V{0.0277777777777778}*(x1*(x23 + x24*(-x25 + x27 - x29 + x30 + x38) - x33*x49*x49 + x39) + V{1});
fEq[4] = -V{0.111111111111111}*x1*(x22 + x39 - x50) + V{-0.111111111111111};
fEq[5] = -V{0.0277777777777778}*x1*(x10 + x32 - x36 + x39 + x42 + V{-1}) + V{-0.0277777777777778};
fEq[6] = V{0.111111111111111}*x1*(x46 + x51) + V{-0.111111111111111};
fEq[7] = V{0.0277777777777778}*(x1*(x33*(x52*x52) + x44 + x51) + V{-1});
fEq[8] = V{0.111111111111111}*x1*(x43 + x50 + V{1}) + V{-0.111111111111111};

}

template <typename CELL, typename FNEQ, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq, RHO& rho, U& u) any_platform
{
auto x0 = u[1]*u[1];
auto x1 = V{1.5}*x0;
auto x2 = u[0]*u[0];
auto x3 = V{1.5}*x2;
auto x4 = x3 + V{-1};
auto x5 = x1 + x4;
auto x6 = V{0.0277777777777778}*rho;
auto x7 = V{3}*u[0];
auto x8 = -x7;
auto x9 = V{3}*u[1];
auto x10 = -x3;
auto x20 = x10 + x9;
auto x21 = u[0] - u[1];
auto x22 = V{1} - x1;
auto x23 = x22 + V{4.5}*(x21*x21);
auto x24 = V{0.111111111111111}*rho;
auto x25 = V{3}*x2 + x22;
auto x26 = u[0] + u[1];
auto x27 = V{4.5}*(x26*x26);
auto x28 = V{3}*x0;
fNeq[0] = V{0.444444444444444}*rho*x5 + cell[0] + V{0.444444444444444};
fNeq[1] = -x6*(x20 + x23 + x8) + cell[1] + V{0.0277777777777778};
fNeq[2] = -x24*(x25 + x8) + cell[2] + V{0.111111111111111};
fNeq[3] = x6*(-x27 + x5 + x7 + x9) + cell[3] + V{0.0277777777777778};
fNeq[4] = x24*(-x28 + x4 + x9) + cell[4] + V{0.111111111111111};
fNeq[5] = -x6*(x10 + x23 + x7 - x9) + cell[5] + V{0.0277777777777778};
fNeq[6] = -x24*(x25 + x7) + cell[6] + V{0.111111111111111};
fNeq[7] = -x6*(x20 + x22 + x27 + x7) + cell[7] + V{0.0277777777777778};
fNeq[8] = -x24*(x20 + x28 + V{1}) + cell[8] + V{0.111111111111111};

}

template <typename CELL, typename FNEQ, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
auto x1 = V{1} / ((x0)*(x0));
auto x2 = V{1.5}*x1;
auto x3 = cell[1] - cell[5];
auto x4 = -cell[4] + cell[8];
auto x5 = x3 + x4 - cell[3] + cell[7];
auto x6 = x5*x5;
auto x7 = x2*x6;
auto x8 = cell[2] - cell[6];
auto x9 = x3 + x8 + cell[3] - cell[7];
auto x10 = x9*x9;
auto x20 = x10*x2;
auto x21 = x20 + V{-1};
auto x22 = x21 + x7;
auto x23 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778};
auto x24 = V{1} / (x0);
auto x25 = V{3}*cell[3];
auto x26 = V{3}*cell[7];
auto x27 = V{3}*cell[1] - V{3}*cell[5];
auto x28 = x24*(x25 - x26 + x27 + V{3}*cell[2] - V{3}*cell[6]);
auto x29 = V{4.5}*x1;
auto x30 = x4 + x8 + V{2}*cell[1] - V{2}*cell[5];
auto x31 = x29*(x30*x30);
auto x32 = V{1} - x7;
auto x33 = x24*(-x25 + x26 + x27 - V{3}*cell[4] + V{3}*cell[8]);
auto x34 = -x20 + x33;
auto x35 = x32 + x34;
auto x36 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111};
auto x37 = V{3}*x1;
auto x38 = x10*x37 + x32;
auto x39 = -x28;
auto x40 = V{2}*cell[3];
auto x41 = V{2}*cell[7];
auto x42 = x4 - x40 + x41 - cell[2] + cell[6];
auto x43 = x22 + x33;
auto x44 = x37*x6;
auto x45 = x40 - x41 + x8 + cell[4] - cell[8];
fNeq[0] = x22*(V{0.444444444444444}*cell[0] + V{0.444444444444444}*cell[1] + V{0.444444444444444}*cell[2] + V{0.444444444444444}*cell[3] + V{0.444444444444444}*cell[4] + V{0.444444444444444}*cell[5] + V{0.444444444444444}*cell[6] + V{0.444444444444444}*cell[7] + V{0.444444444444444}*cell[8] + V{0.444444444444444}) + cell[0] + V{0.444444444444444};
fNeq[1] = -x23*(x28 + x31 + x35) + cell[1] + V{0.0277777777777778};
fNeq[2] = -x36*(x28 + x38) + cell[2] + V{0.111111111111111};
fNeq[3] = x23*(-x29*x42*x42 + x39 + x43) + cell[3] + V{0.0277777777777778};
fNeq[4] = x36*(x21 + x33 - x44) + cell[4] + V{0.111111111111111};
fNeq[5] = x23*(x28 - x31 + x43) + cell[5] + V{0.0277777777777778};
fNeq[6] = -x36*(x38 + x39) + cell[6] + V{0.111111111111111};
fNeq[7] = -x23*(x29*(x45*x45) + x35 + x39) + cell[7] + V{0.0277777777777778};
fNeq[8] = -x36*(x34 + x44 + V{1}) + cell[8] + V{0.111111111111111};

}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto bgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x9 = omega + V{-1};
auto x10 = u[1]*u[1];
auto x11 = V{1.5}*x10;
auto x12 = u[0]*u[0];
auto x13 = V{1.5}*x12;
auto x14 = x13 + V{-1};
auto x15 = x11 + x14;
auto x16 = V{0.0277777777777778}*omega;
auto x17 = V{3}*u[0];
auto x18 = -x17;
auto x19 = V{3}*u[1];
auto x20 = -x13;
auto x21 = x19 + x20;
auto x22 = u[0] - u[1];
auto x23 = V{1} - x11;
auto x24 = x23 + V{4.5}*(x22*x22);
auto x25 = V{0.111111111111111}*omega;
auto x26 = V{3}*x12 + x23;
auto x27 = u[0] + u[1];
auto x28 = V{4.5}*(x27*x27);
auto x29 = V{3}*x10;
cell[0] = -V{0.444444444444444}*omega*(rho*x15 + V{1}) - x9*cell[0];
cell[1] = x16*(rho*(x18 + x21 + x24) + V{-1}) - x9*cell[1];
cell[2] = x25*(rho*(x18 + x26) + V{-1}) - x9*cell[2];
cell[3] = -x16*(rho*(x15 + x17 + x19 - x28) + V{1}) - x9*cell[3];
cell[4] = -x25*(rho*(x14 + x19 - x29) + V{1}) - x9*cell[4];
cell[5] = x16*(rho*(x17 - x19 + x20 + x24) + V{-1}) - x9*cell[5];
cell[6] = x25*(rho*(x17 + x26) + V{-1}) - x9*cell[6];
cell[7] = x16*(rho*(x17 + x21 + x23 + x28) + V{-1}) - x9*cell[7];
cell[8] = x25*(rho*(x21 + x29 + V{1}) + V{-1}) - x9*cell[8];
return x10 + x12;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto adeBgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega)
{
auto x9 = omega + V{-1};
auto x10 = V{3}*u[0];
auto x11 = V{3}*u[1];
auto x12 = x11 + V{1};
auto x13 = V{0.0277777777777778}*omega;
auto x14 = x10 + V{-1};
auto x15 = V{0.111111111111111}*omega;
auto x16 = x10 + V{1};
cell[0] = V{0.444444444444444}*omega*(rho + V{-1}) - x9*cell[0];
cell[1] = x13*(rho*(-x10 + x12) + V{-1}) - x9*cell[1];
cell[2] = -x15*(rho*x14 + V{1}) - x9*cell[2];
cell[3] = -x13*(rho*(x11 + x14) + V{1}) - x9*cell[3];
cell[4] = -x15*(rho*(x11 + V{-1}) + V{1}) - x9*cell[4];
cell[5] = x13*(rho*(-x11 + x16) + V{-1}) - x9*cell[5];
cell[6] = x15*(rho*x16 + V{-1}) - x9*cell[6];
cell[7] = x13*(rho*(x11 + x16) + V{-1}) - x9*cell[7];
cell[8] = x15*(rho*x12 + V{-1}) - x9*cell[8];
return u[0]*u[0] + u[1]*u[1];
}

template <typename CELL, typename PRESSURE, typename J, typename OMEGA, typename V=typename CELL::value_t>
static auto incBgkCollision(CELL& cell, PRESSURE& pressure, J& j, OMEGA& omega) any_platform
{
auto x9 = j[0]*j[0];
auto x10 = j[1]*j[1];
auto x11 = omega + V{-1};
auto x12 = -j[0] + j[1];
auto x13 = V{0.0833333333333333}*j[1];
auto x14 = V{0.0833333333333333}*j[0];
auto x15 = V{0.0416666666666667}*x9;
auto x16 = V{0.0416666666666667}*x10;
auto x17 = V{0.0833333333333333}*pressure;
auto x18 = x15 + x16 - x17 + V{0.0277777777777778};
auto x19 = x14 + x18;
auto x20 = V{0.166666666666667}*x10;
auto x21 = V{0.333333333333333}*j[0];
auto x22 = V{0.333333333333333}*x9;
auto x23 = V{0.333333333333333}*pressure;
auto x24 = V{0.111111111111111} - x23;
auto x25 = j[0] + j[1];
auto x26 = V{0.125}*(x25*x25);
auto x27 = V{0.166666666666667}*x9;
auto x28 = V{0.333333333333333}*j[1];
auto x29 = V{0.333333333333333}*x10;
auto x30 = j[0] - j[1];
auto x31 = x23 + V{-0.111111111111111};
cell[0] = -omega*(-V{1.33333333333333}*pressure + V{0.666666666666667}*x10 + V{0.666666666666667}*x9 + V{0.444444444444444}) - x11*cell[0];
cell[1] = -omega*(-x13 + x19 - V{0.125}*x12*x12) - x11*cell[1];
cell[2] = -omega*(x20 + x21 - x22 + x24) - x11*cell[2];
cell[3] = -omega*(x13 + x19 - x26) - x11*cell[3];
cell[4] = -omega*(x24 + x27 + x28 - x29) - x11*cell[4];
cell[5] = -omega*(x13 - x14 + x18 - V{0.125}*x30*x30) - x11*cell[5];
cell[6] = omega*(-x20 + x21 + x22 + x31) - x11*cell[6];
cell[7] = omega*(x13 + x14 - x15 - x16 + x17 + x26 + V{-0.0277777777777778}) - x11*cell[7];
cell[8] = omega*(-x27 + x28 + x29 + x31) - x11*cell[8];
return x10 + x9;
}

template <typename CELL, typename RHO, typename U, typename RATIORHO, typename OMEGA, typename V=typename CELL::value_t>
static auto constRhoBgkCollision(CELL& cell, RHO& rho, U& u, RATIORHO& ratioRho, OMEGA& omega) any_platform
{
auto x9 = omega + V{-1};
auto x10 = V{0.444444444444444}*rho;
auto x11 = u[1]*u[1];
auto x12 = V{1.5}*x11;
auto x13 = u[0]*u[0];
auto x14 = V{1.5}*x13;
auto x15 = x14 + V{-1};
auto x16 = x12 + x15;
auto x17 = V{0.0277777777777778}*rho;
auto x18 = V{3}*u[0];
auto x19 = -x18;
auto x20 = V{3}*u[1];
auto x21 = -x14;
auto x22 = x20 + x21;
auto x23 = u[0] - u[1];
auto x24 = V{1} - x12;
auto x25 = x24 + V{4.5}*(x23*x23);
auto x26 = x19 + x22 + x25;
auto x27 = ratioRho*x17;
auto x28 = V{0.111111111111111}*rho;
auto x29 = V{3}*x13 + x24;
auto x30 = x19 + x29;
auto x31 = ratioRho*x28;
auto x32 = u[0] + u[1];
auto x33 = V{4.5}*(x32*x32);
auto x34 = x16 + x18 + x20 - x33;
auto x35 = V{3}*x11;
auto x36 = x15 + x20 - x35;
auto x37 = x18 - x20 + x21 + x25;
auto x38 = x18 + x29;
auto x39 = x18 + x22 + x24 + x33;
auto x40 = x22 + x35 + V{1};
cell[0] = -ratioRho*x10*x16 - x9*(x10*x16 + cell[0] + V{0.444444444444444}) + V{-0.444444444444444};
cell[1] = x26*x27 - x9*(-x17*x26 + cell[1] + V{0.0277777777777778}) + V{-0.0277777777777778};
cell[2] = x30*x31 - x9*(-x28*x30 + cell[2] + V{0.111111111111111}) + V{-0.111111111111111};
cell[3] = -x27*x34 - x9*(x17*x34 + cell[3] + V{0.0277777777777778}) + V{-0.0277777777777778};
cell[4] = -x31*x36 - x9*(x28*x36 + cell[4] + V{0.111111111111111}) + V{-0.111111111111111};
cell[5] = x27*x37 - x9*(-x17*x37 + cell[5] + V{0.0277777777777778}) + V{-0.0277777777777778};
cell[6] = x31*x38 - x9*(-x28*x38 + cell[6] + V{0.111111111111111}) + V{-0.111111111111111};
cell[7] = x27*x39 - x9*(-x17*x39 + cell[7] + V{0.0277777777777778}) + V{-0.0277777777777778};
cell[8] = x31*x40 - x9*(-x28*x40 + cell[8] + V{0.111111111111111}) + V{-0.111111111111111};
return x11 + x13;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x9 = omega + V{-1};
auto x10 = V{0.0833333333333333}*cell[6];
auto x11 = V{0.0833333333333333}*cell[2];
auto x12 = V{3}*u[0];
auto x13 = V{3}*u[1];
auto x14 = x13 + V{1};
auto x15 = -x12 + x14;
auto x16 = V{0.00462962962962963}*rho;
auto x17 = x12 + V{1};
auto x18 = -x13 + x17;
auto x19 = V{0.00925925925925926}*rho;
auto x20 = x17*x19;
auto x21 = x12 + V{-1};
auto x22 = x19*x21;
auto x23 = x13 + V{-1};
auto x24 = x14*x19 + x19*x23 + V{0.0833333333333333}*cell[4] - V{0.0833333333333333}*cell[8];
auto x25 = x9*(x10 - x11 + x15*x16 - x16*x18 - x20 - x22 + x24 - V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[5]);
auto x26 = V{0.0277777777777778}*rho;
auto x27 = V{0.333333333333333}*cell[1];
auto x28 = V{0.333333333333333}*cell[5];
auto x29 = x18*x19;
auto x30 = V{0.037037037037037}*rho;
auto x31 = x15*x19;
auto x32 = x12 + x14;
auto x33 = x12 + x23;
auto x34 = x19*x32 + x19*x33 + V{0.333333333333333}*cell[3] - V{0.333333333333333}*cell[7] + V{4.62592926927149e-18};
auto x35 = x9*(x17*x30 + x21*x30 + x27 - x28 + x29 - x31 + x34 + V{0.333333333333333}*cell[2] - V{0.333333333333333}*cell[6]);
auto x36 = V{0.111111111111111}*rho;
auto x37 = x9*(-x10 + x11 + x16*x32 + x16*x33 + x20 + x22 + x24 + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[7] + V{2.31296463463574e-18});
auto x38 = x9*(x14*x30 + x23*x30 - x27 + x28 - x29 + x31 + x34 + V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[8]);
cell[0] = V{0.444444444444444}*rho + V{-0.444444444444444};
cell[1] = x15*x26 + x25 + V{-0.0277777777777778};
cell[2] = -x21*x36 - x35 + V{-0.111111111111111};
cell[3] = -x26*x33 - x37 + V{-0.0277777777777778};
cell[4] = -x23*x36 - x38 + V{-0.111111111111111};
cell[5] = x18*x26 - x25 + V{-0.0277777777777778};
cell[6] = x17*x36 + x35 + V{-0.111111111111111};
cell[7] = x26*x32 + x37 + V{-0.0277777777777778};
cell[8] = x14*x36 + x38 + V{-0.111111111111111};
return u[0]*u[0] + u[1]*u[1];
}

template <typename CELL, typename RHO, typename U, typename PI, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, PI& pi, OMEGA& omega) any_platform
{
auto x9 = omega + V{-1};
auto x10 = u[1]*u[1];
auto x11 = V{1.5}*x10;
auto x12 = u[0]*u[0];
auto x13 = V{1.5}*x12;
auto x14 = x13 + V{-1};
auto x15 = x11 + x14;
auto x16 = V{0.0277777777777778}*rho;
auto x17 = V{3}*u[0];
auto x18 = -x17;
auto x19 = V{3}*u[1];
auto x20 = -x13;
auto x21 = x19 + x20;
auto x22 = u[0] - u[1];
auto x23 = V{1} - x11;
auto x24 = x23 + V{4.5}*(x22*x22);
auto x25 = V{0.25}*pi[1];
auto x26 = V{0.0833333333333333}*pi[0] + V{0.0833333333333333}*pi[2];
auto x27 = -x9*(-x25 + x26) + V{-0.0277777777777778};
auto x28 = V{0.111111111111111}*rho;
auto x29 = V{3}*x12 + x23;
auto x30 = -x9*(V{0.333333333333333}*pi[0] - V{0.166666666666667}*pi[2]) + V{-0.111111111111111};
auto x31 = x9*(x25 + x26);
auto x32 = u[0] + u[1];
auto x33 = V{4.5}*(x32*x32);
auto x34 = V{3}*x10;
auto x35 = x9*(V{0.166666666666667}*pi[0] - V{0.333333333333333}*pi[2]) + V{-0.111111111111111};
cell[0] = -V{0.444444444444444}*rho*x15 + V{0.666666666666667}*x9*(pi[0] + pi[2]) + V{-0.444444444444444};
cell[1] = x16*(x18 + x21 + x24) + x27;
cell[2] = x28*(x18 + x29) + x30;
cell[3] = -x16*(x15 + x17 + x19 - x33) - x31 + V{-0.0277777777777778};
cell[4] = -x28*(x14 + x19 - x34) + x35;
cell[5] = x16*(x17 - x19 + x20 + x24) + x27;
cell[6] = x28*(x17 + x29) + x30;
cell[7] = x16*(x17 + x21 + x23 + x33) - x31 + V{-0.0277777777777778};
cell[8] = x28*(x21 + x34 + V{1}) + x35;
return x10 + x12;
}

template <typename CELL, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineEqFirstOrder(CELL& cell, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x9 = V{3}*newU[0];
auto x10 = V{3}*newU[1];
auto x11 = x10 + V{1};
auto x12 = x9 + V{-1};
auto x13 = x9 + V{1};
cell[0] = V{0.444444444444444}*newRho + V{-0.444444444444444};
cell[1] = V{0.0277777777777778}*newRho*(x11 - x9) + V{-0.0277777777777778};
cell[2] = -V{0.111111111111111}*newRho*x12 + V{-0.111111111111111};
cell[3] = -V{0.0277777777777778}*newRho*(x10 + x12) + V{-0.0277777777777778};
cell[4] = -V{0.111111111111111}*newRho*(x10 + V{-1}) + V{-0.111111111111111};
cell[5] = V{0.0277777777777778}*newRho*(-x10 + x13) + V{-0.0277777777777778};
cell[6] = V{0.111111111111111}*newRho*x13 + V{-0.111111111111111};
cell[7] = V{0.0277777777777778}*newRho*(x10 + x13) + V{-0.0277777777777778};
cell[8] = V{0.111111111111111}*newRho*x11 + V{-0.111111111111111};

}

template <typename CELL, typename OLDRHO, typename OLDU, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineNEq(CELL& cell, OLDRHO& oldRho, OLDU& oldU, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x9 = oldU[1]*oldU[1];
auto x10 = V{1.5}*x9;
auto x11 = oldU[0]*oldU[0];
auto x12 = V{1.5}*x11;
auto x13 = x12 + V{-1};
auto x14 = x10 + x13;
auto x15 = newU[1]*newU[1];
auto x16 = V{1.5}*x15;
auto x17 = newU[0]*newU[0];
auto x18 = V{1.5}*x17;
auto x19 = x18 + V{-1};
auto x20 = x16 + x19;
auto x21 = V{0.0277777777777778}*newRho;
auto x22 = V{3}*newU[0];
auto x23 = -x22;
auto x24 = V{3}*newU[1];
auto x25 = -x18;
auto x26 = x24 + x25;
auto x27 = newU[0] - newU[1];
auto x28 = V{1} - x16;
auto x29 = x28 + V{4.5}*(x27*x27);
auto x30 = V{0.0277777777777778}*oldRho;
auto x31 = V{3}*oldU[0];
auto x32 = -x31;
auto x33 = V{3}*oldU[1];
auto x34 = -x12;
auto x35 = x33 + x34;
auto x36 = oldU[0] - oldU[1];
auto x37 = V{1} - x10;
auto x38 = x37 + V{4.5}*(x36*x36);
auto x39 = V{0.111111111111111}*newRho;
auto x40 = V{3}*x17 + x28;
auto x41 = V{0.111111111111111}*oldRho;
auto x42 = V{3}*x11 + x37;
auto x43 = oldU[0] + oldU[1];
auto x44 = V{4.5}*(x43*x43);
auto x45 = newU[0] + newU[1];
auto x46 = V{4.5}*(x45*x45);
auto x47 = V{3}*x9;
auto x48 = V{3}*x15;
cell[0] = -V{0.444444444444444}*newRho*x20 + V{0.444444444444444}*oldRho*x14 + cell[0];
cell[1] = x21*(x23 + x26 + x29) - x30*(x32 + x35 + x38) + cell[1];
cell[2] = x39*(x23 + x40) - x41*(x32 + x42) + cell[2];
cell[3] = -x21*(x20 + x22 + x24 - x46) + x30*(x14 + x31 + x33 - x44) + cell[3];
cell[4] = -x39*(x19 + x24 - x48) + x41*(x13 + x33 - x47) + cell[4];
cell[5] = x21*(x22 - x24 + x25 + x29) - x30*(x31 - x33 + x34 + x38) + cell[5];
cell[6] = x39*(x22 + x40) - x41*(x31 + x42) + cell[6];
cell[7] = x21*(x22 + x26 + x28 + x46) - x30*(x31 + x35 + x37 + x44) + cell[7];
cell[8] = x39*(x26 + x48 + V{1}) - x41*(x35 + x47 + V{1}) + cell[8];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void defineNEqFromPi(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x9 = u[1]*u[1];
auto x10 = V{1.5}*x9;
auto x11 = u[0]*u[0];
auto x12 = V{1.5}*x11;
auto x13 = x12 + V{-1};
auto x14 = x10 + x13;
auto x15 = V{0.0277777777777778}*rho;
auto x16 = V{3}*u[0];
auto x17 = -x16;
auto x18 = V{3}*u[1];
auto x19 = -x12;
auto x20 = x18 + x19;
auto x21 = u[0] - u[1];
auto x22 = V{1} - x10;
auto x23 = x22 + V{4.5}*(x21*x21);
auto x24 = V{0.25}*pi[1];
auto x25 = V{0.0833333333333333}*pi[0] + V{0.0833333333333333}*pi[2] + V{-0.0277777777777778};
auto x26 = -x24 + x25;
auto x27 = V{0.111111111111111}*rho;
auto x28 = V{3}*x11 + x22;
auto x29 = V{0.333333333333333}*pi[0] - V{0.166666666666667}*pi[2] + V{-0.111111111111111};
auto x30 = u[0] + u[1];
auto x31 = V{4.5}*(x30*x30);
auto x32 = x24 + x25;
auto x33 = V{3}*x9;
auto x34 = -V{0.166666666666667}*pi[0] + V{0.333333333333333}*pi[2] + V{-0.111111111111111};
cell[0] = -V{0.444444444444444}*rho*x14 - V{0.666666666666667}*pi[0] - V{0.666666666666667}*pi[2] + V{-0.444444444444444};
cell[1] = x15*(x17 + x20 + x23) + x26;
cell[2] = x27*(x17 + x28) + x29;
cell[3] = -x15*(x14 + x16 + x18 - x31) + x32;
cell[4] = -x27*(x13 + x18 - x33) + x34;
cell[5] = x15*(x16 - x18 + x19 + x23) + x26;
cell[6] = x27*(x16 + x28) + x29;
cell[7] = x15*(x16 + x20 + x22 + x31) + x32;
cell[8] = x27*(x20 + x33 + V{1}) + x34;

}

template <typename CELL, typename FORCE, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell, FORCE& force) any_platform
{
auto x0 = cell[7] + cell[8];
auto x1 = cell[2] + cell[3];
auto x2 = x0 + x1 + cell[0] + cell[1] + cell[4] + cell[5] + cell[6];
auto x3 = V{1} / (x2 + V{1});
auto x4 = V{1}*x3;
auto x5 = cell[1] - cell[5];
auto x6 = x1 + x5 - cell[6] - cell[7];
auto x7 = x0 + x5 - cell[3] - cell[4];
auto x8 = x6*x7;
auto x9 = x2 + V{1};
auto x10 = x9*(x6*force[1] - x7*force[0]);
auto x11 = x3*x9;
auto x12 = -V{0.333333333333333}*cell[0] + V{0.666666666666667}*cell[1] + V{0.666666666666667}*cell[3] + V{0.666666666666667}*cell[5] + V{0.666666666666667}*cell[7];
auto x13 = -x11*x6*force[0] + x12 - x3*x6*x6 + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[4] + V{0.666666666666667}*cell[6] - V{0.333333333333333}*cell[8];
auto x14 = x11*x7*force[1] + x12 - x3*x7*x7 - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[4] - V{0.333333333333333}*cell[6] + V{0.666666666666667}*cell[8];
return (-V{0.5}*x10*x3 + x4*x8 - V{1}*cell[1] + V{1}*cell[3] - V{1}*cell[5] + V{1}*cell[7])*(-x10*x4 + V{2}*x3*x8 - V{2}*cell[1] + V{2}*cell[3] - V{2}*cell[5] + V{2}*cell[7]) + x13*x13 + x14*x14;
}

template <typename CELL, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell) any_platform
{
auto x0 = -cell[7];
auto x1 = cell[7] + cell[8];
auto x2 = cell[1] + cell[2] + cell[3];
auto x3 = V{1} / (x1 + x2 + cell[0] + cell[4] + cell[5] + cell[6] + V{1});
auto x4 = -cell[5];
auto x5 = cell[1] - cell[3];
auto x6 = x1 + x4 + x5 - cell[4];
auto x7 = x0 + x2 + x4 - cell[6];
auto x8 = x0 - x3*x6*x7 + x5 + cell[5];
auto x9 = V{1}*x3;
auto x10 = -V{0.333333333333333}*cell[0] + V{0.666666666666667}*cell[1] + V{0.666666666666667}*cell[3] + V{0.666666666666667}*cell[5] + V{0.666666666666667}*cell[7];
auto x11 = x10 - x9*x7*x7 + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[4] + V{0.666666666666667}*cell[6] - V{0.333333333333333}*cell[8];
auto x12 = x10 - x9*x6*x6 - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[4] - V{0.333333333333333}*cell[6] + V{0.666666666666667}*cell[8];
return x11*x11 + x12*x12 + V{2}*(x8*x8);
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addExternalForce(CELL& cell, RHO& rho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x9 = force[0]*u[0];
auto x10 = force[1]*u[1];
auto x11 = rho*(V{0.5}*omega + V{-1});
auto x12 = V{9}*u[0];
auto x13 = V{6}*u[1];
auto x14 = x13 + V{3};
auto x15 = V{9}*u[1];
auto x16 = V{6}*u[0];
auto x17 = V{0.0277777777777778}*x11;
auto x18 = x16 + V{-3};
auto x19 = V{0.111111111111111}*force[0];
auto x20 = -V{0.333333333333333}*x10;
auto x21 = x13 + V{-3};
auto x22 = V{0.111111111111111}*force[1];
auto x23 = -V{0.333333333333333}*x9;
auto x24 = x16 + V{3};
cell[0] = V{1.33333333333333}*x11*(x10 + x9) + cell[0];
cell[1] = -x17*((-x12 + x14)*force[1] - (x15 - x16 + V{3})*force[0]) + cell[1];
cell[2] = -x11*(x18*x19 + x20) + cell[2];
cell[3] = -x17*((x12 + x21)*force[1] + (x15 + x18)*force[0]) + cell[3];
cell[4] = -x11*(x21*x22 + x23) + cell[4];
cell[5] = -x17*((-x15 + x24)*force[0] - (x12 - x13 + V{3})*force[1]) + cell[5];
cell[6] = -x11*(x19*x24 + x20) + cell[6];
cell[7] = -x17*((x12 + x14)*force[1] + (x15 + x24)*force[0]) + cell[7];
cell[8] = -x11*(x14*x22 + x23) + cell[8];

}

};

template <typename... FIELDS>
struct lbm<descriptors::D3Q7<FIELDS...>> {

template <typename CELL, typename V=typename CELL::value_t>
static auto computeRho(CELL& cell) any_platform
{

return cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
}

template <typename CELL, typename J, typename V=typename CELL::value_t>
static void computeJ(CELL& cell, J& j) any_platform
{
j[0] = -V{1}*cell[1] + V{1}*cell[4];
j[1] = -V{1}*cell[2] + V{1}*cell[5];
j[2] = -V{1}*cell[3] + V{1}*cell[6];

}

template <typename CELL, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeRhoU(CELL& cell, RHO& rho, U& u) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
auto x1 = V{1}/x0;
rho = x0;
u[0] = -x1*(cell[1] - cell[4]);
u[1] = -x1*(cell[2] - cell[5]);
u[2] = -x1*(cell[3] - cell[6]);

}

template <typename CELL, typename RHO, typename J, typename V=typename CELL::value_t>
static void computeRhoJ(CELL& cell, RHO& rho, J& j) any_platform
{
rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
j[0] = -V{1}*cell[1] + V{1}*cell[4];
j[1] = -V{1}*cell[2] + V{1}*cell[5];
j[2] = -V{1}*cell[3] + V{1}*cell[6];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeStress(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = V{0.25} - V{0.25}*rho;
auto x1 = rho*u[0];
pi[0] = -rho*u[0]*u[0] + x0 + V{1}*cell[1] + V{1}*cell[4];
pi[1] = -x1*u[1];
pi[2] = -x1*u[2];
pi[3] = -rho*u[1]*u[1] + x0 + V{1}*cell[2] + V{1}*cell[5];
pi[4] = -rho*u[1]*u[2];
pi[5] = -rho*u[2]*u[2] + x0 + V{1}*cell[3] + V{1}*cell[6];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeAllMomenta(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
auto x1 = cell[1] - cell[4];
auto x2 = V{1}/x0;
auto x3 = x1*x2;
auto x4 = cell[2] - cell[5];
auto x5 = x2*x4;
auto x6 = cell[3] - cell[6];
auto x7 = -V{0.25}*cell[0];
auto x8 = x7 - V{0.25}*cell[3] - V{0.25}*cell[6];
auto x9 = -V{0.25}*cell[2] - V{0.25}*cell[5];
auto x20 = -V{0.25}*cell[1] - V{0.25}*cell[4];
rho = x0;
u[0] = -x3;
u[1] = -x5;
u[2] = -x2*x6;
pi[0] = -x2*x1*x1 + x8 + x9 + V{0.75}*cell[1] + V{0.75}*cell[4];
pi[1] = -x3*x4;
pi[2] = -x3*x6;
pi[3] = -x2*x4*x4 + x20 + x8 + V{0.75}*cell[2] + V{0.75}*cell[5];
pi[4] = -x5*x6;
pi[5] = -x2*x6*x6 + x20 + x7 + x9 + V{0.75}*cell[3] + V{0.75}*cell[6];

}

template <typename CELL, typename FEQ, typename V=typename CELL::value_t>
static void computeFeq(CELL& cell, FEQ& fEq) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6];
auto x1 = x0 + V{1};
auto x2 = x0 + V{1};
auto x3 = V{1} / ((x2)*(x2));
auto x4 = V{2}*x3;
auto x5 = -cell[1] + cell[4];
auto x6 = x5*x5;
auto x7 = x4*x6;
auto x8 = -cell[2] + cell[5];
auto x9 = x8*x8;
auto x17 = x4*x9;
auto x18 = -cell[3] + cell[6];
auto x19 = x18*x18;
auto x20 = x19*x4;
auto x21 = x17 + x20 + V{-1};
auto x22 = V{1} / (x2);
auto x23 = V{4}*cell[4];
auto x24 = V{4}*cell[1];
auto x25 = V{6}*x3;
auto x26 = V{4}*cell[5];
auto x27 = V{4}*cell[2];
auto x28 = x7 + V{-1};
auto x29 = V{4}*cell[6];
auto x30 = V{4}*cell[3];
auto x31 = cell[2] - cell[5];
auto x32 = x31*x31;
auto x33 = x32*x4;
auto x34 = cell[1] - cell[4];
auto x35 = x34*x34;
auto x36 = cell[3] - cell[6];
auto x37 = x36*x36;
auto x38 = x37*x4 + V{-1};
auto x39 = x35*x4;
fEq[0] = -V{0.25}*x1*(x21 + x7) + V{-0.25};
fEq[1] = -V{0.125}*x1*(x21 + x22*(x23 - x24) - x25*x6) + V{-0.125};
fEq[2] = -V{0.125}*x1*(x20 + x22*(x26 - x27) - x25*x9 + x28) + V{-0.125};
fEq[3] = -V{0.125}*x1*(x17 - x19*x25 + x22*(x29 - x30) + x28) + V{-0.125};
fEq[4] = -V{0.125}*x1*(x22*(-x23 + x24) - x25*x35 + x33 + x38) + V{-0.125};
fEq[5] = -V{0.125}*x1*(x22*(-x26 + x27) - x25*x32 + x38 + x39) + V{-0.125};
fEq[6] = -V{0.125}*x1*(x22*(-x29 + x30) - x25*x37 + x33 + x39 + V{-1}) + V{-0.125};

}

template <typename CELL, typename FNEQ, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq, RHO& rho, U& u) any_platform
{
auto x0 = u[0]*u[0];
auto x1 = V{2}*x0;
auto x2 = u[1]*u[1];
auto x3 = V{2}*x2;
auto x4 = u[2]*u[2];
auto x5 = V{2}*x4;
auto x6 = x3 + x5 + V{-1};
auto x7 = V{0.125}*rho;
auto x8 = V{4}*u[0];
auto x9 = V{6}*x0;
auto x17 = V{4}*u[1];
auto x18 = V{6}*x2;
auto x19 = x1 + V{-1};
auto x20 = V{4}*u[2];
auto x21 = V{6}*x4;
auto x22 = -x3;
auto x23 = V{1} - x5;
auto x24 = -x1;
fNeq[0] = V{0.25}*rho*(x1 + x6) + cell[0] + V{0.25};
fNeq[1] = x7*(x6 + x8 - x9) + cell[1] + V{0.125};
fNeq[2] = x7*(x17 - x18 + x19 + x5) + cell[2] + V{0.125};
fNeq[3] = x7*(x19 + x20 - x21 + x3) + cell[3] + V{0.125};
fNeq[4] = -x7*(x22 + x23 + x8 + x9) + cell[4] + V{0.125};
fNeq[5] = -x7*(x17 + x18 + x23 + x24) + cell[5] + V{0.125};
fNeq[6] = -x7*(x20 + x21 + x22 + x24 + V{1}) + cell[6] + V{0.125};

}

template <typename CELL, typename FNEQ, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
auto x1 = V{1} / ((x0)*(x0));
auto x2 = V{2}*x1;
auto x3 = cell[1] - cell[4];
auto x4 = x3*x3;
auto x5 = x2*x4;
auto x6 = cell[2] - cell[5];
auto x7 = x6*x6;
auto x8 = x2*x7;
auto x9 = cell[3] - cell[6];
auto x17 = x9*x9;
auto x18 = x17*x2;
auto x19 = x18 + x8 + V{-1};
auto x20 = V{0.125}*cell[0] + V{0.125}*cell[1] + V{0.125}*cell[2] + V{0.125}*cell[3] + V{0.125}*cell[4] + V{0.125}*cell[5] + V{0.125}*cell[6] + V{0.125};
auto x21 = V{1} / (x0);
auto x22 = x21*(V{4}*cell[1] - V{4}*cell[4]);
auto x23 = V{6}*x1;
auto x24 = x19 - x23*x4;
auto x25 = x21*(V{4}*cell[2] - V{4}*cell[5]);
auto x26 = x5 + V{-1};
auto x27 = x18 - x23*x7 + x26;
auto x28 = x21*(V{4}*cell[3] - V{4}*cell[6]);
auto x29 = -x17*x23 + x26 + x8;
fNeq[0] = (x19 + x5)*(V{0.25}*cell[0] + V{0.25}*cell[1] + V{0.25}*cell[2] + V{0.25}*cell[3] + V{0.25}*cell[4] + V{0.25}*cell[5] + V{0.25}*cell[6] + V{0.25}) + cell[0] + V{0.25};
fNeq[1] = x20*(-x22 + x24) + cell[1] + V{0.125};
fNeq[2] = x20*(-x25 + x27) + cell[2] + V{0.125};
fNeq[3] = x20*(-x28 + x29) + cell[3] + V{0.125};
fNeq[4] = x20*(x22 + x24) + cell[4] + V{0.125};
fNeq[5] = x20*(x25 + x27) + cell[5] + V{0.125};
fNeq[6] = x20*(x28 + x29) + cell[6] + V{0.125};

}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto bgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x7 = omega + V{-1};
auto x8 = u[0]*u[0];
auto x9 = V{2}*x8;
auto x10 = u[1]*u[1];
auto x11 = V{2}*x10;
auto x12 = u[2]*u[2];
auto x13 = V{2}*x12;
auto x14 = x11 + x13 + V{-1};
auto x15 = V{0.125}*omega;
auto x16 = V{4}*u[0];
auto x17 = V{6}*x8;
auto x18 = V{4}*u[1];
auto x19 = V{6}*x10;
auto x20 = x9 + V{-1};
auto x21 = V{4}*u[2];
auto x22 = V{6}*x12;
auto x23 = -x11;
auto x24 = V{1} - x13;
auto x25 = -x9;
cell[0] = -V{0.25}*omega*(rho*(x14 + x9) + V{1}) - x7*cell[0];
cell[1] = -x15*(rho*(x14 + x16 - x17) + V{1}) - x7*cell[1];
cell[2] = -x15*(rho*(x13 + x18 - x19 + x20) + V{1}) - x7*cell[2];
cell[3] = -x15*(rho*(x11 + x20 + x21 - x22) + V{1}) - x7*cell[3];
cell[4] = x15*(rho*(x16 + x17 + x23 + x24) + V{-1}) - x7*cell[4];
cell[5] = x15*(rho*(x18 + x19 + x24 + x25) + V{-1}) - x7*cell[5];
cell[6] = x15*(rho*(x21 + x22 + x23 + x25 + V{1}) + V{-1}) - x7*cell[6];
return x10 + x12 + x8;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto adeBgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega)
{
auto x7 = omega + V{-1};
auto x8 = V{4}*u[0];
auto x9 = V{0.125}*omega;
auto x10 = V{4}*u[1];
auto x11 = V{4}*u[2];
cell[0] = V{0.25}*omega*(rho + V{-1}) - x7*cell[0];
cell[1] = -x7*cell[1] - x9*(rho*(x8 + V{-1}) + V{1});
cell[2] = -x7*cell[2] - x9*(rho*(x10 + V{-1}) + V{1});
cell[3] = -x7*cell[3] - x9*(rho*(x11 + V{-1}) + V{1});
cell[4] = -x7*cell[4] + x9*(rho*(x8 + V{1}) + V{-1});
cell[5] = -x7*cell[5] + x9*(rho*(x10 + V{1}) + V{-1});
cell[6] = -x7*cell[6] + x9*(rho*(x11 + V{1}) + V{-1});
return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

template <typename CELL, typename PRESSURE, typename J, typename OMEGA, typename V=typename CELL::value_t>
static auto incBgkCollision(CELL& cell, PRESSURE& pressure, J& j, OMEGA& omega) any_platform
{
auto x7 = j[0]*j[0];
auto x8 = j[1]*j[1];
auto x9 = j[2]*j[2];
auto x10 = omega + V{-1};
auto x11 = V{0.25}*x8;
auto x12 = V{0.5}*j[0];
auto x13 = V{0.75}*x7;
auto x14 = V{0.25}*x9;
auto x15 = V{0.5}*pressure;
auto x16 = -x15;
auto x17 = x14 + x16 + V{0.125};
auto x18 = V{0.25}*x7;
auto x19 = V{0.5}*j[1];
auto x20 = V{0.75}*x8;
auto x21 = V{0.5}*j[2];
auto x22 = V{0.75}*x9;
auto x23 = -x11;
auto x24 = -x14 + x15 + V{-0.125};
auto x25 = -x18;
cell[0] = -omega*(-V{1}*pressure + V{0.5}*x7 + V{0.5}*x8 + V{0.5}*x9 + V{0.25}) - x10*cell[0];
cell[1] = -omega*(x11 + x12 - x13 + x17) - x10*cell[1];
cell[2] = -omega*(x17 + x18 + x19 - x20) - x10*cell[2];
cell[3] = -omega*(x11 + x16 + x18 + x21 - x22 + V{0.125}) - x10*cell[3];
cell[4] = omega*(x12 + x13 + x23 + x24) - x10*cell[4];
cell[5] = omega*(x19 + x20 + x24 + x25) - x10*cell[5];
cell[6] = omega*(x15 + x21 + x22 + x23 + x25 + V{-0.125}) - x10*cell[6];
return x7 + x8 + x9;
}

template <typename CELL, typename RHO, typename U, typename RATIORHO, typename OMEGA, typename V=typename CELL::value_t>
static auto constRhoBgkCollision(CELL& cell, RHO& rho, U& u, RATIORHO& ratioRho, OMEGA& omega) any_platform
{
auto x7 = omega + V{-1};
auto x8 = V{0.25}*rho;
auto x9 = u[0]*u[0];
auto x10 = V{2}*x9;
auto x11 = u[1]*u[1];
auto x12 = V{2}*x11;
auto x13 = u[2]*u[2];
auto x14 = V{2}*x13;
auto x15 = x12 + x14 + V{-1};
auto x16 = x10 + x15;
auto x17 = V{0.125}*rho;
auto x18 = V{4}*u[0];
auto x19 = V{6}*x9;
auto x20 = x15 + x18 - x19;
auto x21 = ratioRho*x17;
auto x22 = V{4}*u[1];
auto x23 = V{6}*x11;
auto x24 = x10 + V{-1};
auto x25 = x14 + x22 - x23 + x24;
auto x26 = V{4}*u[2];
auto x27 = V{6}*x13;
auto x28 = x12 + x24 + x26 - x27;
auto x29 = -x12;
auto x30 = V{1} - x14;
auto x31 = x18 + x19 + x29 + x30;
auto x32 = -x10;
auto x33 = x22 + x23 + x30 + x32;
auto x34 = x26 + x27 + x29 + x32 + V{1};
cell[0] = -ratioRho*x16*x8 - x7*(x16*x8 + cell[0] + V{0.25}) + V{-0.25};
cell[1] = -x20*x21 - x7*(x17*x20 + cell[1] + V{0.125}) + V{-0.125};
cell[2] = -x21*x25 - x7*(x17*x25 + cell[2] + V{0.125}) + V{-0.125};
cell[3] = -x21*x28 - x7*(x17*x28 + cell[3] + V{0.125}) + V{-0.125};
cell[4] = x21*x31 - x7*(-x17*x31 + cell[4] + V{0.125}) + V{-0.125};
cell[5] = x21*x33 - x7*(-x17*x33 + cell[5] + V{0.125}) + V{-0.125};
cell[6] = x21*x34 - x7*(-x17*x34 + cell[6] + V{0.125}) + V{-0.125};
return x11 + x13 + x9;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x7 = omega + V{-1};
auto x8 = V{4}*u[0];
auto x9 = V{0.0625}*rho;
auto x10 = x8 + V{1};
auto x11 = -x10*x9 - V{0.5}*cell[1] + V{0.5}*cell[4];
auto x12 = x8 + V{-1};
auto x13 = V{0.125}*rho;
auto x14 = V{4}*u[1];
auto x15 = x14 + V{1};
auto x16 = -x15*x9 - V{0.5}*cell[2] + V{0.5}*cell[5];
auto x17 = x14 + V{-1};
auto x18 = V{4}*u[2];
auto x19 = x18 + V{1};
auto x20 = -x19*x9 - V{0.5}*cell[3] + V{0.5}*cell[6];
auto x21 = x18 + V{-1};
cell[0] = V{0.25}*rho + V{-0.25};
cell[1] = -x12*x13 + x7*(x11 + x9*(V{1} - x8)) + V{-0.125};
cell[2] = -x13*x17 + x7*(x16 + x9*(V{1} - x14)) + V{-0.125};
cell[3] = -x13*x21 + x7*(x20 + x9*(V{1} - x18)) + V{-0.125};
cell[4] = x10*x13 - x7*(x11 - x12*x9) + V{-0.125};
cell[5] = x13*x15 - x7*(x16 - x17*x9) + V{-0.125};
cell[6] = x13*x19 - x7*(x20 - x21*x9) + V{-0.125};
return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

template <typename CELL, typename RHO, typename U, typename PI, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, PI& pi, OMEGA& omega) any_platform
{
auto x7 = omega + V{-1};
auto x8 = u[0]*u[0];
auto x9 = V{2}*x8;
auto x10 = u[1]*u[1];
auto x11 = V{2}*x10;
auto x12 = u[2]*u[2];
auto x13 = V{2}*x12;
auto x14 = x11 + x13 + V{-1};
auto x15 = V{0.125}*rho;
auto x16 = V{4}*u[0];
auto x17 = V{6}*x8;
auto x18 = V{0.25}*pi[3];
auto x19 = V{0.25}*pi[5];
auto x20 = x7*(x18 + x19 - V{0.75}*pi[0]) + V{-0.125};
auto x21 = V{4}*u[1];
auto x22 = V{6}*x10;
auto x23 = x9 + V{-1};
auto x24 = V{0.25}*pi[0];
auto x25 = x7*(x19 + x24 - V{0.75}*pi[3]) + V{-0.125};
auto x26 = V{4}*u[2];
auto x27 = V{6}*x12;
auto x28 = x7*(x18 + x24 - V{0.75}*pi[5]) + V{-0.125};
auto x29 = -x11;
auto x30 = V{1} - x13;
auto x31 = -x9;
cell[0] = -V{0.25}*rho*(x14 + x9) + V{0.5}*x7*(pi[0] + pi[3] + pi[5]) + V{-0.25};
cell[1] = -x15*(x14 + x16 - x17) + x20;
cell[2] = -x15*(x13 + x21 - x22 + x23) + x25;
cell[3] = -x15*(x11 + x23 + x26 - x27) + x28;
cell[4] = x15*(x16 + x17 + x29 + x30) + x20;
cell[5] = x15*(x21 + x22 + x30 + x31) + x25;
cell[6] = x15*(x26 + x27 + x29 + x31 + V{1}) + x28;
return x10 + x12 + x8;
}

template <typename CELL, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineEqFirstOrder(CELL& cell, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x7 = V{4}*newU[0];
auto x8 = V{4}*newU[1];
auto x9 = V{4}*newU[2];
cell[0] = V{0.25}*newRho + V{-0.25};
cell[1] = -V{0.125}*newRho*(x7 + V{-1}) + V{-0.125};
cell[2] = -V{0.125}*newRho*(x8 + V{-1}) + V{-0.125};
cell[3] = -V{0.125}*newRho*(x9 + V{-1}) + V{-0.125};
cell[4] = V{0.125}*newRho*(x7 + V{1}) + V{-0.125};
cell[5] = V{0.125}*newRho*(x8 + V{1}) + V{-0.125};
cell[6] = V{0.125}*newRho*(x9 + V{1}) + V{-0.125};

}

template <typename CELL, typename OLDRHO, typename OLDU, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineNEq(CELL& cell, OLDRHO& oldRho, OLDU& oldU, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x7 = oldU[0]*oldU[0];
auto x8 = V{2}*x7;
auto x9 = oldU[1]*oldU[1];
auto x10 = V{2}*x9;
auto x11 = oldU[2]*oldU[2];
auto x12 = V{2}*x11;
auto x13 = x10 + x12 + V{-1};
auto x14 = newU[0]*newU[0];
auto x15 = V{2}*x14;
auto x16 = newU[1]*newU[1];
auto x17 = V{2}*x16;
auto x18 = newU[2]*newU[2];
auto x19 = V{2}*x18;
auto x20 = x17 + x19 + V{-1};
auto x21 = V{0.125}*oldRho;
auto x22 = V{4}*oldU[0];
auto x23 = V{6}*x7;
auto x24 = V{0.125}*newRho;
auto x25 = V{4}*newU[0];
auto x26 = V{6}*x14;
auto x27 = V{4}*oldU[1];
auto x28 = V{6}*x9;
auto x29 = x8 + V{-1};
auto x30 = V{4}*newU[1];
auto x31 = V{6}*x16;
auto x32 = x15 + V{-1};
auto x33 = V{4}*oldU[2];
auto x34 = V{6}*x11;
auto x35 = V{4}*newU[2];
auto x36 = V{6}*x18;
auto x37 = -x17;
auto x38 = V{1} - x19;
auto x39 = -x10;
auto x40 = V{1} - x12;
auto x41 = -x15;
auto x42 = -x8;
cell[0] = -V{0.25}*newRho*(x15 + x20) + V{0.25}*oldRho*(x13 + x8) + cell[0];
cell[1] = x21*(x13 + x22 - x23) - x24*(x20 + x25 - x26) + cell[1];
cell[2] = x21*(x12 + x27 - x28 + x29) - x24*(x19 + x30 - x31 + x32) + cell[2];
cell[3] = x21*(x10 + x29 + x33 - x34) - x24*(x17 + x32 + x35 - x36) + cell[3];
cell[4] = -x21*(x22 + x23 + x39 + x40) + x24*(x25 + x26 + x37 + x38) + cell[4];
cell[5] = -x21*(x27 + x28 + x40 + x42) + x24*(x30 + x31 + x38 + x41) + cell[5];
cell[6] = -x21*(x33 + x34 + x39 + x42 + V{1}) + x24*(x35 + x36 + x37 + x41 + V{1}) + cell[6];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void defineNEqFromPi(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x7 = u[0]*u[0];
auto x8 = V{2}*x7;
auto x9 = u[1]*u[1];
auto x10 = V{2}*x9;
auto x11 = u[2]*u[2];
auto x12 = V{2}*x11;
auto x13 = x10 + x12 + V{-1};
auto x14 = V{0.125}*rho;
auto x15 = V{4}*u[0];
auto x16 = V{6}*x7;
auto x17 = -V{0.25}*pi[3];
auto x18 = -V{0.25}*pi[5] + V{-0.125};
auto x19 = x17 + x18 + V{0.75}*pi[0];
auto x20 = V{4}*u[1];
auto x21 = V{6}*x9;
auto x22 = x8 + V{-1};
auto x23 = -V{0.25}*pi[0];
auto x24 = x18 + x23 + V{0.75}*pi[3];
auto x25 = V{4}*u[2];
auto x26 = V{6}*x11;
auto x27 = x17 + x23 + V{0.75}*pi[5] + V{-0.125};
auto x28 = -x10;
auto x29 = V{1} - x12;
auto x30 = -x8;
cell[0] = -V{0.25}*rho*(x13 + x8) - V{0.5}*pi[0] - V{0.5}*pi[3] - V{0.5}*pi[5] + V{-0.25};
cell[1] = -x14*(x13 + x15 - x16) + x19;
cell[2] = -x14*(x12 + x20 - x21 + x22) + x24;
cell[3] = -x14*(x10 + x22 + x25 - x26) + x27;
cell[4] = x14*(x15 + x16 + x28 + x29) + x19;
cell[5] = x14*(x20 + x21 + x29 + x30) + x24;
cell[6] = x14*(x25 + x26 + x28 + x30 + V{1}) + x27;

}

template <typename CELL, typename FORCE, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell, FORCE& force) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6];
auto x1 = x0 + V{1};
auto x2 = V{1} / ((x1)*(x1));
auto x3 = cell[1] - cell[4];
auto x4 = cell[3] - cell[6];
auto x5 = x0 + V{1};
auto x6 = V{0.5}*x5;
auto x7 = V{1}*x3;
auto x8 = x4*x7 + x6*(x3*force[2] + x4*force[0]);
auto x9 = cell[2] - cell[5];
auto x10 = x3*force[1] + x9*force[0];
auto x11 = V{1}*x5;
auto x12 = V{2}*x9;
auto x13 = x4*force[1] + x9*force[2];
auto x14 = V{1} / (x1);
auto x15 = x14*x5;
auto x16 = V{0.25}*cell[0];
auto x17 = x16 + V{0.25}*cell[1] + V{0.25}*cell[4];
auto x18 = V{0.25}*cell[2] + V{0.25}*cell[5];
auto x19 = x14*(x4*x4) + x15*x4*force[2] + x17 + x18 - V{0.75}*cell[3] - V{0.75}*cell[6];
auto x20 = V{0.25}*cell[3] + V{0.25}*cell[6];
auto x21 = x14*(x9*x9) + x15*x9*force[1] + x17 + x20 - V{0.75}*cell[2] - V{0.75}*cell[5];
auto x22 = x14*(x3*x3) + x15*x3*force[0] + x16 + x18 + x20 - V{0.75}*cell[1] - V{0.75}*cell[4];
return x2*(x10*x11 + x12*x3)*(x10*x6 + x7*x9) + x2*(x11*x13 + x12*x4)*(x13*x6 + V{1}*x4*x9) + 2*x2*(x8*x8) + x19*x19 + x21*x21 + x22*x22;
}

template <typename CELL, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell) any_platform
{
auto x0 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
auto x1 = V{2}/((x0)*(x0));
auto x2 = cell[1] - cell[4];
auto x3 = x2*x2;
auto x4 = cell[2] - cell[5];
auto x5 = x4*x4;
auto x6 = cell[3] - cell[6];
auto x7 = x6*x6;
auto x8 = V{1}/x0;
auto x9 = V{0.25}*cell[0];
auto x10 = x9 + V{0.25}*cell[1] + V{0.25}*cell[4];
auto x11 = V{0.25}*cell[2] + V{0.25}*cell[5];
auto x12 = x10 + x11 + x7*x8 - V{0.75}*cell[3] - V{0.75}*cell[6];
auto x13 = V{0.25}*cell[3] + V{0.25}*cell[6];
auto x14 = x10 + x13 + x5*x8 - V{0.75}*cell[2] - V{0.75}*cell[5];
auto x15 = x11 + x13 + x3*x8 + x9 - V{0.75}*cell[1] - V{0.75}*cell[4];
return x1*x3*x5 + x1*x3*x7 + x1*x5*x7 + x12*x12 + x14*x14 + x15*x15;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addExternalForce(CELL& cell, RHO& rho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x7 = force[0]*u[0];
auto x8 = force[1]*u[1];
auto x9 = force[2]*u[2];
auto x10 = rho*(V{0.5}*omega + V{-1});
auto x11 = V{12}*u[0];
auto x12 = V{0.125}*force[0];
auto x13 = V{0.5}*x8;
auto x14 = V{0.5}*x9;
auto x15 = x13 + x14;
auto x16 = V{12}*u[1];
auto x17 = V{0.125}*force[1];
auto x18 = V{0.5}*x7;
auto x19 = x14 + x18;
auto x20 = V{12}*u[2];
auto x21 = V{0.125}*force[2];
auto x22 = x13 + x18;
cell[0] = V{1}*x10*(x7 + x8 + x9) + cell[0];
cell[1] = x10*(-x12*(x11 + V{-4}) + x15) + cell[1];
cell[2] = x10*(-x17*(x16 + V{-4}) + x19) + cell[2];
cell[3] = x10*(-x21*(x20 + V{-4}) + x22) + cell[3];
cell[4] = x10*(-x12*(x11 + V{4}) + x15) + cell[4];
cell[5] = x10*(-x17*(x16 + V{4}) + x19) + cell[5];
cell[6] = x10*(-x21*(x20 + V{4}) + x22) + cell[6];

}

};

template <typename... FIELDS>
struct lbm<descriptors::D3Q19<FIELDS...>> {

template <typename CELL, typename V=typename CELL::value_t>
static auto computeRho(CELL& cell) any_platform
{

return cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
}

template <typename CELL, typename J, typename V=typename CELL::value_t>
static void computeJ(CELL& cell, J& j) any_platform
{
auto x0 = cell[13] - cell[4];
auto x1 = cell[15] - cell[6];
auto x2 = cell[17] - cell[8];
j[0] = V{1}*x0 + V{1}*x1 + V{1}*cell[10] + V{1}*cell[14] + V{1}*cell[16] - V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7];
j[1] = V{1}*x0 + V{1}*x2 + V{1}*cell[11] - V{1}*cell[14] + V{1}*cell[18] - V{1}*cell[2] + V{1}*cell[5] - V{1}*cell[9];
j[2] = V{1}*x1 + V{1}*x2 + V{1}*cell[12] - V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + V{1}*cell[7] + V{1}*cell[9];

}

template <typename CELL, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeRhoU(CELL& cell, RHO& rho, U& u) any_platform
{
auto x0 = cell[10] + cell[14] + cell[16];
auto x1 = cell[11] + cell[18] + cell[5];
auto x2 = cell[12] + cell[7] + cell[9];
auto x3 = x0 + x1 + x2 + cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + V{1};
auto x4 = cell[13] - cell[4];
auto x5 = cell[15] - cell[6];
auto x6 = V{1}/x3;
auto x7 = cell[17] - cell[8];
rho = x3;
u[0] = x6*(x0 + x4 + x5 - cell[1] - cell[5] - cell[7]);
u[1] = x6*(x1 + x4 + x7 - cell[14] - cell[2] - cell[9]);
u[2] = x6*(x2 + x5 + x7 - cell[16] - cell[18] - cell[3]);

}

template <typename CELL, typename RHO, typename J, typename V=typename CELL::value_t>
static void computeRhoJ(CELL& cell, RHO& rho, J& j) any_platform
{
auto x0 = cell[10] + cell[14] + cell[16];
auto x1 = cell[11] + cell[18] + cell[5];
auto x2 = cell[12] + cell[7] + cell[9];
auto x3 = cell[13] - cell[4];
auto x4 = cell[15] - cell[6];
auto x5 = cell[17] - cell[8];
rho = x0 + x1 + x2 + cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + V{1};
j[0] = V{1}*x0 + V{1}*x3 + V{1}*x4 - V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7];
j[1] = V{1}*x1 + V{1}*x3 + V{1}*x5 - V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9];
j[2] = V{1}*x2 + V{1}*x4 + V{1}*x5 - V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeStress(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = V{1}*cell[5];
auto x1 = V{1}*cell[14];
auto x2 = -V{0.333333333333333}*rho;
auto x3 = V{1}*cell[13] + V{1}*cell[4];
auto x4 = x0 + x1 + x2 + x3 + V{0.333333333333333};
auto x5 = V{1}*cell[7];
auto x6 = V{1}*cell[16];
auto x7 = V{1}*cell[15] + V{1}*cell[6];
auto x8 = x5 + x6 + x7;
auto x9 = rho*u[0];
auto x10 = V{1}*cell[9];
auto x11 = V{1}*cell[18];
auto x12 = V{1}*cell[17] + V{1}*cell[8];
auto x13 = x10 + x11 + x12;
pi[0] = -rho*u[0]*u[0] + x4 + x8 + V{1}*cell[10] + V{1}*cell[1];
pi[1] = -x0 - x1 + x3 - x9*u[1];
pi[2] = -x5 - x6 + x7 - x9*u[2];
pi[3] = -rho*u[1]*u[1] + x13 + x4 + V{1}*cell[11] + V{1}*cell[2];
pi[4] = -rho*u[1]*u[2] - x10 - x11 + x12;
pi[5] = -rho*u[2]*u[2] + x13 + x2 + x8 + V{1}*cell[12] + V{1}*cell[3] + V{0.333333333333333};

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeAllMomenta(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = cell[10] + cell[14] + cell[16];
auto x1 = cell[11] + cell[13] + cell[18] + cell[5];
auto x2 = cell[12] + cell[15] + cell[17] + cell[7] + cell[9];
auto x3 = x0 + x1 + x2 + cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + V{1};
auto x4 = -cell[4];
auto x5 = -cell[6];
auto x6 = cell[13] - cell[5];
auto x7 = cell[15] - cell[7];
auto x8 = x0 + x4 + x5 + x6 + x7 - cell[1];
auto x9 = V{1} / (x3);
auto x10 = V{1}*x9;
auto x11 = -cell[14];
auto x12 = -cell[8];
auto x13 = cell[17] - cell[9];
auto x14 = x1 + x11 + x12 + x13 + x4 - cell[2];
auto x15 = -cell[16];
auto x16 = -cell[18];
auto x17 = x12 + x15 + x16 + x2 + x5 - cell[3];
auto x18 = -V{0.333333333333333}*cell[0];
auto x19 = x18 - V{0.333333333333333}*cell[12] + V{0.666666666666667}*cell[13] + V{0.666666666666667}*cell[14] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5];
auto x20 = -V{0.333333333333333}*cell[11] + V{0.666666666666667}*cell[15] + V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x21 = x8*x9;
auto x32 = -V{0.333333333333333}*cell[10] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
rho = x3;
u[0] = x10*x8;
u[1] = x10*x14;
u[2] = x10*x17;
pi[0] = -x10*x8*x8 + x19 + x20 + V{0.666666666666667}*cell[10] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9];
pi[1] = V{1}*x11 - V{1}*x14*x21 + V{1}*x6 + V{1}*cell[4];
pi[2] = V{1}*x15 - V{1}*x17*x21 + V{1}*x7 + V{1}*cell[6];
pi[3] = -x10*x14*x14 + x19 + x32 + V{0.666666666666667}*cell[11] - V{0.333333333333333}*cell[15] - V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7];
pi[4] = V{1}*x13 - V{1}*x14*x17*x9 + V{1}*x16 + V{1}*cell[8];
pi[5] = -x10*x17*x17 + x18 + x20 + x32 + V{0.666666666666667}*cell[12] - V{0.333333333333333}*cell[13] - V{0.333333333333333}*cell[14] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5];

}

template <typename CELL, typename FEQ, typename V=typename CELL::value_t>
static void computeFeq(CELL& cell, FEQ& fEq) any_platform
{
auto x0 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x1 = x0 + V{1};
auto x2 = x0 + V{1};
auto x3 = V{1} / ((x2)*(x2));
auto x4 = V{1.5}*x3;
auto x5 = cell[13] - cell[4];
auto x6 = cell[15] - cell[6];
auto x7 = x5 + x6;
auto x8 = cell[14] - cell[5];
auto x9 = -cell[1];
auto x10 = cell[16] - cell[7];
auto x11 = x10 + x9 + cell[10];
auto x12 = x11 + x7 + x8;
auto x13 = x12*x12;
auto x14 = x13*x4;
auto x15 = cell[17] - cell[8];
auto x16 = x15 + x5;
auto x17 = cell[18] - cell[9];
auto x18 = -cell[2];
auto x19 = -cell[14] + cell[5];
auto x20 = x18 + x19 + cell[11];
auto x21 = x16 + x17 + x20;
auto x41 = x21*x21;
auto x42 = x4*x41;
auto x43 = x15 + x6;
auto x44 = -cell[16] + cell[7];
auto x45 = -cell[3];
auto x46 = -cell[18] + cell[9];
auto x47 = x45 + x46 + cell[12];
auto x48 = x43 + x44 + x47;
auto x49 = x48*x48;
auto x50 = x4*x49;
auto x51 = x42 + x50 + V{-1};
auto x52 = x14 + x51;
auto x53 = V{1} / (x2);
auto x54 = V{3}*cell[14];
auto x55 = V{3}*cell[16];
auto x56 = V{3}*cell[5];
auto x57 = V{3}*cell[7];
auto x58 = V{3}*cell[13] - V{3}*cell[4];
auto x59 = V{3}*cell[15] - V{3}*cell[6];
auto x60 = x53*(x54 + x55 - x56 - x57 + x58 + x59 + V{3}*cell[10] - V{3}*cell[1]);
auto x61 = V{3}*x3;
auto x62 = x13*x61;
auto x63 = V{3}*cell[18];
auto x64 = V{3}*cell[9];
auto x65 = V{3}*cell[17] - V{3}*cell[8];
auto x66 = x53*(-x54 + x56 + x58 + x63 - x64 + x65 + V{3}*cell[11] - V{3}*cell[2]);
auto x67 = x41*x61;
auto x68 = x14 + V{-1};
auto x69 = x53*(-x55 + x57 + x59 - x63 + x64 + x65 + V{3}*cell[12] - V{3}*cell[3]);
auto x70 = x49*x61;
auto x71 = V{4.5}*x3;
auto x72 = x17 + x18 + cell[11];
auto x73 = x11 + x43 + x72 + V{2}*cell[13] - V{2}*cell[4];
auto x74 = x71*(x73*x73);
auto x75 = x52 + x60;
auto x76 = -x66;
auto x77 = V{2}*cell[14];
auto x78 = V{2}*cell[5];
auto x79 = -cell[15] + cell[6];
auto x80 = x15 - cell[10] + cell[1];
auto x81 = x44 + x72 - x77 + x78 + x79 + x80;
auto x82 = x8 + x9 + cell[10];
auto x83 = x16 + x47 + x82 + V{2}*cell[15] - V{2}*cell[6];
auto x84 = x71*(x83*x83);
auto x85 = -x69;
auto x86 = V{2}*cell[16];
auto x87 = V{2}*cell[7];
auto x88 = -cell[13] + cell[4];
auto x89 = x19 + x47 + x80 - x86 + x87 + x88;
auto x90 = x44 + x45 + cell[12];
auto x91 = x20 + x7 + x90 + V{2}*cell[17] - V{2}*cell[8];
auto x92 = x71*(x91*x91);
auto x93 = x52 + x66;
auto x94 = V{2}*cell[18];
auto x95 = V{2}*cell[9];
auto x96 = x6 - cell[11] + cell[2];
auto x97 = x8 + x88 + x90 - x94 + x95 + x96;
auto x98 = -x42;
auto x99 = V{1} - x50;
auto x100 = x98 + x99;
auto x101 = x100 + x60;
auto x102 = -x14;
auto x103 = x102 + x66;
auto x104 = x102 + x69;
auto x105 = -x60;
auto x106 = -cell[17] + cell[8];
auto x107 = x106 + x11 + x46 + x77 - x78 + x96;
auto x108 = x5 - cell[12] + cell[3];
auto x109 = x106 + x108 + x17 + x82 + x86 - x87;
auto x110 = x52 + x69;
auto x111 = x10 + x108 + x20 + x79 + x94 - x95;
fEq[0] = -V{0.333333333333333}*x1*x52 + V{-0.333333333333333};
fEq[1] = -V{0.0555555555555556}*x1*(x51 + x60 - x62) + V{-0.0555555555555556};
fEq[2] = -V{0.0555555555555556}*x1*(x50 + x66 - x67 + x68) + V{-0.0555555555555556};
fEq[3] = -V{0.0555555555555556}*x1*(x42 + x68 + x69 - x70) + V{-0.0555555555555556};
fEq[4] = -V{0.0277777777777778}*x1*(x66 - x74 + x75) + V{-0.0277777777777778};
fEq[5] = -V{0.0277777777777778}*(x1*(-x71*x81*x81 + x75 + x76) + V{1});
fEq[6] = -V{0.0277777777777778}*x1*(x69 + x75 - x84) + V{-0.0277777777777778};
fEq[7] = -V{0.0277777777777778}*(x1*(-x71*x89*x89 + x75 + x85) + V{1});
fEq[8] = -V{0.0277777777777778}*x1*(x69 - x92 + x93) + V{-0.0277777777777778};
fEq[9] = -V{0.0277777777777778}*(x1*(-x71*x97*x97 + x85 + x93) + V{1});
fEq[10] = V{0.0555555555555556}*x1*(x101 + x62) + V{-0.0555555555555556};
fEq[11] = V{0.0555555555555556}*x1*(x103 + x67 + x99) + V{-0.0555555555555556};
fEq[12] = V{0.0555555555555556}*x1*(x104 + x70 + x98 + V{1}) + V{-0.0555555555555556};
fEq[13] = V{0.0277777777777778}*x1*(x101 + x103 + x74) + V{-0.0277777777777778};
fEq[14] = -V{0.0277777777777778}*(x1*(x105 - x71*x107*x107 + x93) + V{1});
fEq[15] = V{0.0277777777777778}*x1*(x101 + x104 + x84) + V{-0.0277777777777778};
fEq[16] = -V{0.0277777777777778}*(x1*(x105 + x110 - x71*x109*x109) + V{1});
fEq[17] = V{0.0277777777777778}*x1*(x100 + x103 + x69 + x92) + V{-0.0277777777777778};
fEq[18] = -V{0.0277777777777778}*(x1*(x110 - x71*x111*x111 + x76) + V{1});

}

template <typename CELL, typename FNEQ, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq, RHO& rho, U& u) any_platform
{
auto x0 = u[0]*u[0];
auto x1 = V{1.5}*x0;
auto x2 = u[1]*u[1];
auto x3 = V{1.5}*x2;
auto x4 = u[2]*u[2];
auto x5 = V{1.5}*x4;
auto x6 = x3 + x5 + V{-1};
auto x7 = x1 + x6;
auto x8 = V{0.0555555555555556}*rho;
auto x9 = V{3}*u[0];
auto x10 = V{3}*x0;
auto x11 = V{3}*u[1];
auto x12 = V{3}*x2;
auto x13 = x1 + V{-1};
auto x14 = V{3}*u[2];
auto x15 = V{3}*x4;
auto x16 = V{0.0277777777777778}*rho;
auto x17 = u[0] + u[1];
auto x18 = V{4.5}*(x17*x17);
auto x19 = x7 + x9;
auto x20 = -x11;
auto x21 = u[0] - u[1];
auto x41 = -V{4.5}*x21*x21;
auto x42 = u[0] + u[2];
auto x43 = V{4.5}*(x42*x42);
auto x44 = -x14;
auto x45 = -u[2];
auto x46 = x45 + u[0];
auto x47 = -V{4.5}*x46*x46;
auto x48 = u[1] + u[2];
auto x49 = V{4.5}*(x48*x48);
auto x50 = x11 + x7;
auto x51 = x45 + u[1];
auto x52 = -V{4.5}*x51*x51;
auto x53 = -x3;
auto x54 = V{1} - x5;
auto x55 = x53 + x54;
auto x56 = x55 + x9;
auto x57 = -x1;
auto x58 = x11 + x57;
auto x59 = x14 + x57;
auto x60 = -x9;
auto x61 = x14 + x7;
fNeq[0] = V{0.333333333333333}*rho*x7 + cell[0] + V{0.333333333333333};
fNeq[1] = x8*(-x10 + x6 + x9) + cell[1] + V{0.0555555555555556};
fNeq[2] = x8*(x11 - x12 + x13 + x5) + cell[2] + V{0.0555555555555556};
fNeq[3] = x8*(x13 + x14 - x15 + x3) + cell[3] + V{0.0555555555555556};
fNeq[4] = x16*(x11 - x18 + x19) + cell[4] + V{0.0277777777777778};
fNeq[5] = x16*(x19 + x20 + x41) + cell[5] + V{0.0277777777777778};
fNeq[6] = x16*(x14 + x19 - x43) + cell[6] + V{0.0277777777777778};
fNeq[7] = x16*(x19 + x44 + x47) + cell[7] + V{0.0277777777777778};
fNeq[8] = x16*(x14 - x49 + x50) + cell[8] + V{0.0277777777777778};
fNeq[9] = x16*(x44 + x50 + x52) + cell[9] + V{0.0277777777777778};
fNeq[10] = -x8*(x10 + x56) + cell[10] + V{0.0555555555555556};
fNeq[11] = -x8*(x12 + x54 + x58) + cell[11] + V{0.0555555555555556};
fNeq[12] = -x8*(x15 + x53 + x59 + V{1}) + cell[12] + V{0.0555555555555556};
fNeq[13] = -x16*(x18 + x56 + x58) + cell[13] + V{0.0277777777777778};
fNeq[14] = x16*(x41 + x50 + x60) + cell[14] + V{0.0277777777777778};
fNeq[15] = -x16*(x43 + x56 + x59) + cell[15] + V{0.0277777777777778};
fNeq[16] = x16*(x47 + x60 + x61) + cell[16] + V{0.0277777777777778};
fNeq[17] = -x16*(x14 + x49 + x55 + x58) + cell[17] + V{0.0277777777777778};
fNeq[18] = x16*(x20 + x52 + x61) + cell[18] + V{0.0277777777777778};

}

template <typename CELL, typename FNEQ, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq) any_platform
{
auto x0 = cell[10] + cell[14];
auto x1 = cell[12] + cell[7];
auto x2 = x0 + x1 + cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + V{1};
auto x3 = V{1} / ((x2)*(x2));
auto x4 = V{1.5}*x3;
auto x5 = cell[13] - cell[4];
auto x6 = cell[15] - cell[6];
auto x7 = x5 + x6;
auto x8 = -cell[1];
auto x9 = cell[16] - cell[7];
auto x10 = x8 + x9;
auto x11 = x0 - cell[5];
auto x12 = x10 + x11 + x7;
auto x13 = x12*x12;
auto x14 = x13*x4;
auto x15 = cell[17] - cell[8];
auto x16 = x15 + x5;
auto x17 = cell[18] - cell[9];
auto x18 = -cell[2];
auto x19 = x18 + cell[11] - cell[14] + cell[5];
auto x20 = x16 + x17 + x19;
auto x21 = x20*x20;
auto x41 = x21*x4;
auto x42 = x15 + x6;
auto x43 = -cell[3];
auto x44 = -cell[18] + cell[9];
auto x45 = x43 + x44;
auto x46 = x1 - cell[16];
auto x47 = x42 + x45 + x46;
auto x48 = x47*x47;
auto x49 = x4*x48;
auto x50 = x41 + x49 + V{-1};
auto x51 = x14 + x50;
auto x52 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x53 = V{1} / (x2);
auto x54 = V{3}*cell[14];
auto x55 = V{3}*cell[16];
auto x56 = V{3}*cell[5];
auto x57 = V{3}*cell[7];
auto x58 = V{3}*cell[13] - V{3}*cell[4];
auto x59 = V{3}*cell[15] - V{3}*cell[6];
auto x60 = x53*(x54 + x55 - x56 - x57 + x58 + x59 + V{3}*cell[10] - V{3}*cell[1]);
auto x61 = V{3}*x3;
auto x62 = x13*x61;
auto x63 = V{3}*cell[18];
auto x64 = V{3}*cell[9];
auto x65 = V{3}*cell[17] - V{3}*cell[8];
auto x66 = x53*(-x54 + x56 + x58 + x63 - x64 + x65 + V{3}*cell[11] - V{3}*cell[2]);
auto x67 = x21*x61;
auto x68 = x14 + V{-1};
auto x69 = x53*(-x55 + x57 + x59 - x63 + x64 + x65 + V{3}*cell[12] - V{3}*cell[3]);
auto x70 = x48*x61;
auto x71 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x72 = V{4.5}*x3;
auto x73 = x10 + cell[10];
auto x74 = x17 + x18 + x42 + x73 + cell[11] + V{2}*cell[13] - V{2}*cell[4];
auto x75 = x72*(x74*x74);
auto x76 = x51 + x60;
auto x77 = -x66;
auto x78 = -cell[17] + cell[8];
auto x79 = x44 + x6 + x73 + x78 - cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5];
auto x80 = -x72*x79*x79;
auto x81 = x11 + x8;
auto x82 = x16 + x45 + x81 + cell[12] + V{2}*cell[15] - V{2}*cell[6];
auto x83 = x72*(x82*x82);
auto x84 = -x69;
auto x85 = x5 - cell[12] + cell[3];
auto x86 = x17 + x78 + x81 + x85 + V{2}*cell[16] - V{2}*cell[7];
auto x87 = -x72*x86*x86;
auto x88 = x19 + x43 + x46 + x7 + V{2}*cell[17] - V{2}*cell[8];
auto x89 = x72*(x88*x88);
auto x90 = x51 + x66;
auto x91 = x19 + x85 + x9 - cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9];
auto x92 = -x72*x91*x91;
auto x93 = -x41;
auto x94 = V{1} - x49;
auto x95 = x93 + x94;
auto x96 = x60 + x95;
auto x97 = -x14;
auto x98 = x66 + x97;
auto x99 = x69 + x97;
auto x100 = -x60;
auto x101 = x51 + x69;
fNeq[0] = x51*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.333333333333333}) + cell[0] + V{0.333333333333333};
fNeq[1] = x52*(x50 + x60 - x62) + cell[1] + V{0.0555555555555556};
fNeq[2] = x52*(x49 + x66 - x67 + x68) + cell[2] + V{0.0555555555555556};
fNeq[3] = x52*(x41 + x68 + x69 - x70) + cell[3] + V{0.0555555555555556};
fNeq[4] = x71*(x66 - x75 + x76) + cell[4] + V{0.0277777777777778};
fNeq[5] = x71*(x76 + x77 + x80) + cell[5] + V{0.0277777777777778};
fNeq[6] = x71*(x69 + x76 - x83) + cell[6] + V{0.0277777777777778};
fNeq[7] = x71*(x76 + x84 + x87) + cell[7] + V{0.0277777777777778};
fNeq[8] = x71*(x69 - x89 + x90) + cell[8] + V{0.0277777777777778};
fNeq[9] = x71*(x84 + x90 + x92) + cell[9] + V{0.0277777777777778};
fNeq[10] = -x52*(x62 + x96) + cell[10] + V{0.0555555555555556};
fNeq[11] = -x52*(x67 + x94 + x98) + cell[11] + V{0.0555555555555556};
fNeq[12] = -x52*(x70 + x93 + x99 + V{1}) + cell[12] + V{0.0555555555555556};
fNeq[13] = -x71*(x75 + x96 + x98) + cell[13] + V{0.0277777777777778};
fNeq[14] = x71*(x100 + x80 + x90) + cell[14] + V{0.0277777777777778};
fNeq[15] = -x71*(x83 + x96 + x99) + cell[15] + V{0.0277777777777778};
fNeq[16] = x71*(x100 + x101 + x87) + cell[16] + V{0.0277777777777778};
fNeq[17] = -x71*(x69 + x89 + x95 + x98) + cell[17] + V{0.0277777777777778};
fNeq[18] = x71*(x101 + x77 + x92) + cell[18] + V{0.0277777777777778};

}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto bgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x19 = omega + V{-1};
auto x20 = u[0]*u[0];
auto x21 = V{1.5}*x20;
auto x22 = u[1]*u[1];
auto x23 = V{1.5}*x22;
auto x24 = u[2]*u[2];
auto x25 = V{1.5}*x24;
auto x26 = x23 + x25 + V{-1};
auto x27 = x21 + x26;
auto x28 = V{0.0555555555555556}*omega;
auto x29 = V{3}*u[0];
auto x30 = V{3}*x20;
auto x31 = V{3}*u[1];
auto x32 = V{3}*x22;
auto x33 = x21 + V{-1};
auto x34 = V{3}*u[2];
auto x35 = V{3}*x24;
auto x36 = V{0.0277777777777778}*omega;
auto x37 = u[0] + u[1];
auto x38 = V{4.5}*(x37*x37);
auto x39 = x27 + x29;
auto x40 = -x31;
auto x41 = -u[0];
auto x42 = x41 + u[1];
auto x43 = u[0] + u[2];
auto x44 = V{4.5}*(x43*x43);
auto x45 = -x34;
auto x46 = x41 + u[2];
auto x47 = u[1] + u[2];
auto x48 = V{4.5}*(x47*x47);
auto x49 = x27 + x31;
auto x50 = -u[1];
auto x51 = x50 + u[2];
auto x52 = -x23;
auto x53 = V{1} - x25;
auto x54 = x52 + x53;
auto x55 = x29 + x54;
auto x56 = -x21;
auto x57 = x31 + x56;
auto x58 = x34 + x56;
auto x59 = -x29;
auto x60 = x50 + u[0];
auto x61 = -u[2];
auto x62 = x61 + u[0];
auto x63 = x27 + x34;
auto x64 = x61 + u[1];
cell[0] = -V{0.333333333333333}*omega*(rho*x27 + V{1}) - x19*cell[0];
cell[1] = -x19*cell[1] - x28*(rho*(x26 + x29 - x30) + V{1});
cell[2] = -x19*cell[2] - x28*(rho*(x25 + x31 - x32 + x33) + V{1});
cell[3] = -x19*cell[3] - x28*(rho*(x23 + x33 + x34 - x35) + V{1});
cell[4] = -x19*cell[4] - x36*(rho*(x31 - x38 + x39) + V{1});
cell[5] = -x19*cell[5] - x36*(rho*(x39 + x40 - V{4.5}*x42*x42) + V{1});
cell[6] = -x19*cell[6] - x36*(rho*(x34 + x39 - x44) + V{1});
cell[7] = -x19*cell[7] - x36*(rho*(x39 + x45 - V{4.5}*x46*x46) + V{1});
cell[8] = -x19*cell[8] - x36*(rho*(x34 - x48 + x49) + V{1});
cell[9] = -x19*cell[9] - x36*(rho*(x45 + x49 - V{4.5}*x51*x51) + V{1});
cell[10] = -x19*cell[10] + x28*(rho*(x30 + x55) + V{-1});
cell[11] = -x19*cell[11] + x28*(rho*(x32 + x53 + x57) + V{-1});
cell[12] = -x19*cell[12] + x28*(rho*(x35 + x52 + x58 + V{1}) + V{-1});
cell[13] = -x19*cell[13] + x36*(rho*(x38 + x55 + x57) + V{-1});
cell[14] = -x19*cell[14] - x36*(rho*(x49 + x59 - V{4.5}*x60*x60) + V{1});
cell[15] = -x19*cell[15] + x36*(rho*(x44 + x55 + x58) + V{-1});
cell[16] = -x19*cell[16] - x36*(rho*(x59 + x63 - V{4.5}*x62*x62) + V{1});
cell[17] = -x19*cell[17] + x36*(rho*(x34 + x48 + x54 + x57) + V{-1});
cell[18] = -x19*cell[18] - x36*(rho*(x40 + x63 - V{4.5}*x64*x64) + V{1});
return x20 + x22 + x24;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto adeBgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega)
{
auto x19 = omega + V{-1};
auto x20 = V{3}*u[0];
auto x21 = x20 + V{-1};
auto x22 = V{0.0555555555555556}*omega;
auto x23 = V{3}*u[1];
auto x24 = x23 + V{-1};
auto x25 = V{3}*u[2];
auto x26 = V{0.0277777777777778}*omega;
auto x27 = -x20;
auto x28 = x23 + V{1};
auto x29 = x25 + V{1};
auto x30 = -x23;
auto x31 = x20 + V{1};
auto x32 = -x25;
cell[0] = V{0.333333333333333}*omega*(rho + V{-1}) - x19*cell[0];
cell[1] = -x19*cell[1] - x22*(rho*x21 + V{1});
cell[2] = -x19*cell[2] - x22*(rho*x24 + V{1});
cell[3] = -x19*cell[3] - x22*(rho*(x25 + V{-1}) + V{1});
cell[4] = -x19*cell[4] - x26*(rho*(x21 + x23) + V{1});
cell[5] = -x19*cell[5] + x26*(rho*(x27 + x28) + V{-1});
cell[6] = -x19*cell[6] - x26*(rho*(x21 + x25) + V{1});
cell[7] = -x19*cell[7] + x26*(rho*(x27 + x29) + V{-1});
cell[8] = -x19*cell[8] - x26*(rho*(x24 + x25) + V{1});
cell[9] = -x19*cell[9] + x26*(rho*(x29 + x30) + V{-1});
cell[10] = -x19*cell[10] + x22*(rho*x31 + V{-1});
cell[11] = -x19*cell[11] + x22*(rho*x28 + V{-1});
cell[12] = -x19*cell[12] + x22*(rho*x29 + V{-1});
cell[13] = -x19*cell[13] + x26*(rho*(x23 + x31) + V{-1});
cell[14] = -x19*cell[14] + x26*(rho*(x30 + x31) + V{-1});
cell[15] = -x19*cell[15] + x26*(rho*(x25 + x31) + V{-1});
cell[16] = -x19*cell[16] + x26*(rho*(x31 + x32) + V{-1});
cell[17] = -x19*cell[17] + x26*(rho*(x25 + x28) + V{-1});
cell[18] = -x19*cell[18] + x26*(rho*(x28 + x32) + V{-1});
return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

template <typename CELL, typename PRESSURE, typename J, typename OMEGA, typename V=typename CELL::value_t>
static auto incBgkCollision(CELL& cell, PRESSURE& pressure, J& j, OMEGA& omega) any_platform
{
auto x19 = j[0]*j[0];
auto x20 = j[1]*j[1];
auto x21 = j[2]*j[2];
auto x22 = omega + V{-1};
auto x23 = V{0.0833333333333333}*x20;
auto x24 = V{0.166666666666667}*j[0];
auto x25 = V{0.166666666666667}*x19;
auto x26 = V{0.0833333333333333}*x21;
auto x27 = V{0.166666666666667}*pressure;
auto x28 = -x27;
auto x29 = x26 + x28 + V{0.0555555555555556};
auto x30 = V{0.0833333333333333}*x19;
auto x31 = V{0.166666666666667}*j[1];
auto x32 = V{0.166666666666667}*x20;
auto x33 = V{0.166666666666667}*j[2];
auto x34 = V{0.166666666666667}*x21;
auto x35 = j[0] + j[1];
auto x36 = V{0.125}*(x35*x35);
auto x37 = V{0.0833333333333333}*j[0];
auto x38 = V{0.0833333333333333}*j[1];
auto x39 = x37 + x38;
auto x40 = V{0.0416666666666667}*x19;
auto x41 = V{0.0416666666666667}*x20;
auto x42 = V{0.0416666666666667}*x21;
auto x43 = V{0.0833333333333333}*pressure;
auto x44 = x40 + x41 + x42 - x43 + V{0.0277777777777778};
auto x45 = -j[0];
auto x46 = x45 + j[1];
auto x47 = -x38;
auto x48 = x37 + x44;
auto x49 = V{0.0833333333333333}*j[2];
auto x50 = j[0] + j[2];
auto x51 = V{0.125}*(x50*x50);
auto x52 = x45 + j[2];
auto x53 = -x49;
auto x54 = j[1] + j[2];
auto x55 = V{0.125}*(x54*x54);
auto x56 = x38 + x44;
auto x57 = -j[1];
auto x58 = x57 + j[2];
auto x59 = -x23;
auto x60 = -x26 + x27 + V{-0.0555555555555556};
auto x61 = -x30;
auto x62 = -x40 - x41 - x42 + x43 + V{-0.0277777777777778};
auto x63 = x57 + j[0];
auto x64 = -x37;
auto x65 = x49 + x62;
auto x66 = -j[2];
auto x67 = x66 + j[0];
auto x68 = x44 + x49;
auto x69 = x66 + j[1];
cell[0] = -omega*(-V{1}*pressure + V{0.5}*x19 + V{0.5}*x20 + V{0.5}*x21 + V{0.333333333333333}) - x22*cell[0];
cell[1] = -omega*(x23 + x24 - x25 + x29) - x22*cell[1];
cell[2] = -omega*(x29 + x30 + x31 - x32) - x22*cell[2];
cell[3] = -omega*(x23 + x28 + x30 + x33 - x34 + V{0.0555555555555556}) - x22*cell[3];
cell[4] = -omega*(-x36 + x39 + x44) - x22*cell[4];
cell[5] = -omega*(x47 + x48 - V{0.125}*x46*x46) - x22*cell[5];
cell[6] = -omega*(x48 + x49 - x51) - x22*cell[6];
cell[7] = -omega*(x48 + x53 - V{0.125}*x52*x52) - x22*cell[7];
cell[8] = -omega*(x49 - x55 + x56) - x22*cell[8];
cell[9] = -omega*(x53 + x56 - V{0.125}*x58*x58) - x22*cell[9];
cell[10] = omega*(x24 + x25 + x59 + x60) - x22*cell[10];
cell[11] = omega*(x31 + x32 + x60 + x61) - x22*cell[11];
cell[12] = omega*(x27 + x33 + x34 + x59 + x61 + V{-0.0555555555555556}) - x22*cell[12];
cell[13] = omega*(x36 + x39 + x62) - x22*cell[13];
cell[14] = -omega*(x56 + x64 - V{0.125}*x63*x63) - x22*cell[14];
cell[15] = omega*(x37 + x51 + x65) - x22*cell[15];
cell[16] = -omega*(x64 + x68 - V{0.125}*x67*x67) - x22*cell[16];
cell[17] = omega*(x38 + x55 + x65) - x22*cell[17];
cell[18] = -omega*(x47 + x68 - V{0.125}*x69*x69) - x22*cell[18];
return x19 + x20 + x21;
}

template <typename CELL, typename RHO, typename U, typename RATIORHO, typename OMEGA, typename V=typename CELL::value_t>
static auto constRhoBgkCollision(CELL& cell, RHO& rho, U& u, RATIORHO& ratioRho, OMEGA& omega) any_platform
{
auto x19 = omega + V{-1};
auto x20 = V{0.333333333333333}*rho;
auto x21 = u[0]*u[0];
auto x22 = V{1.5}*x21;
auto x23 = u[1]*u[1];
auto x24 = V{1.5}*x23;
auto x25 = u[2]*u[2];
auto x26 = V{1.5}*x25;
auto x27 = x24 + x26 + V{-1};
auto x28 = x22 + x27;
auto x29 = V{0.0555555555555556}*rho;
auto x30 = V{3}*u[0];
auto x31 = V{3}*x21;
auto x32 = x27 + x30 - x31;
auto x33 = ratioRho*x29;
auto x34 = V{3}*u[1];
auto x35 = V{3}*x23;
auto x36 = x22 + V{-1};
auto x37 = x26 + x34 - x35 + x36;
auto x38 = V{3}*u[2];
auto x39 = V{3}*x25;
auto x40 = x24 + x36 + x38 - x39;
auto x41 = V{0.0277777777777778}*rho;
auto x42 = u[0] + u[1];
auto x43 = V{4.5}*(x42*x42);
auto x44 = x28 + x30;
auto x45 = x34 - x43 + x44;
auto x46 = ratioRho*x41;
auto x47 = -x34;
auto x48 = -u[0];
auto x49 = x48 + u[1];
auto x50 = x44 + x47 - V{4.5}*x49*x49;
auto x51 = u[0] + u[2];
auto x52 = V{4.5}*(x51*x51);
auto x53 = x38 + x44 - x52;
auto x54 = -x38;
auto x55 = x48 + u[2];
auto x56 = x44 + x54 - V{4.5}*x55*x55;
auto x57 = u[1] + u[2];
auto x58 = V{4.5}*(x57*x57);
auto x59 = x28 + x34;
auto x60 = x38 - x58 + x59;
auto x61 = -u[1];
auto x62 = x61 + u[2];
auto x63 = x54 + x59 - V{4.5}*x62*x62;
auto x64 = -x24;
auto x65 = V{1} - x26;
auto x66 = x64 + x65;
auto x67 = x30 + x66;
auto x68 = x31 + x67;
auto x69 = -x22;
auto x70 = x34 + x69;
auto x71 = x35 + x65 + x70;
auto x72 = x38 + x69;
auto x73 = x39 + x64 + x72 + V{1};
auto x74 = x43 + x67 + x70;
auto x75 = -x30;
auto x76 = x61 + u[0];
auto x77 = x59 + x75 - V{4.5}*x76*x76;
auto x78 = x52 + x67 + x72;
auto x79 = -u[2];
auto x80 = x79 + u[0];
auto x81 = x28 + x38;
auto x82 = x75 + x81 - V{4.5}*x80*x80;
auto x83 = x38 + x58 + x66 + x70;
auto x84 = x79 + u[1];
auto x85 = x47 + x81 - V{4.5}*x84*x84;
cell[0] = -ratioRho*x20*x28 - x19*(x20*x28 + cell[0] + V{0.333333333333333}) + V{-0.333333333333333};
cell[1] = -x19*(x29*x32 + cell[1] + V{0.0555555555555556}) - x32*x33 + V{-0.0555555555555556};
cell[2] = -x19*(x29*x37 + cell[2] + V{0.0555555555555556}) - x33*x37 + V{-0.0555555555555556};
cell[3] = -x19*(x29*x40 + cell[3] + V{0.0555555555555556}) - x33*x40 + V{-0.0555555555555556};
cell[4] = -x19*(x41*x45 + cell[4] + V{0.0277777777777778}) - x45*x46 + V{-0.0277777777777778};
cell[5] = -x19*(x41*x50 + cell[5] + V{0.0277777777777778}) - x46*x50 + V{-0.0277777777777778};
cell[6] = -x19*(x41*x53 + cell[6] + V{0.0277777777777778}) - x46*x53 + V{-0.0277777777777778};
cell[7] = -x19*(x41*x56 + cell[7] + V{0.0277777777777778}) - x46*x56 + V{-0.0277777777777778};
cell[8] = -x19*(x41*x60 + cell[8] + V{0.0277777777777778}) - x46*x60 + V{-0.0277777777777778};
cell[9] = -x19*(x41*x63 + cell[9] + V{0.0277777777777778}) - x46*x63 + V{-0.0277777777777778};
cell[10] = -x19*(-x29*x68 + cell[10] + V{0.0555555555555556}) + x33*x68 + V{-0.0555555555555556};
cell[11] = -x19*(-x29*x71 + cell[11] + V{0.0555555555555556}) + x33*x71 + V{-0.0555555555555556};
cell[12] = -x19*(-x29*x73 + cell[12] + V{0.0555555555555556}) + x33*x73 + V{-0.0555555555555556};
cell[13] = -x19*(-x41*x74 + cell[13] + V{0.0277777777777778}) + x46*x74 + V{-0.0277777777777778};
cell[14] = -x19*(x41*x77 + cell[14] + V{0.0277777777777778}) - x46*x77 + V{-0.0277777777777778};
cell[15] = -x19*(-x41*x78 + cell[15] + V{0.0277777777777778}) + x46*x78 + V{-0.0277777777777778};
cell[16] = -x19*(x41*x82 + cell[16] + V{0.0277777777777778}) - x46*x82 + V{-0.0277777777777778};
cell[17] = -x19*(-x41*x83 + cell[17] + V{0.0277777777777778}) + x46*x83 + V{-0.0277777777777778};
cell[18] = -x19*(x41*x85 + cell[18] + V{0.0277777777777778}) - x46*x85 + V{-0.0277777777777778};
return x21 + x23 + x25;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x19 = omega + V{-1};
auto x20 = V{3}*u[0];
auto x21 = x20 + V{1};
auto x22 = V{0.00925925925925926}*rho;
auto x23 = x20 + V{-1};
auto x24 = V{3}*u[1];
auto x25 = -x24;
auto x26 = x21 + x25;
auto x27 = V{0.00462962962962963}*rho;
auto x28 = -x20;
auto x29 = x24 + V{1};
auto x30 = x28 + x29;
auto x31 = x26*x27 - x27*x30 - V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[5];
auto x32 = V{3}*u[2];
auto x33 = -x32;
auto x34 = x21 + x33;
auto x35 = x32 + V{1};
auto x36 = x28 + x35;
auto x37 = x27*x34 - x27*x36 - V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[7];
auto x38 = V{0.166666666666667}*cell[4];
auto x39 = V{0.166666666666667}*cell[13];
auto x40 = x21 + x24;
auto x41 = x27*x40;
auto x42 = x23 + x24;
auto x43 = x27*x42;
auto x44 = x38 - x39 + x41 + x43;
auto x45 = V{0.166666666666667}*cell[6];
auto x46 = V{0.166666666666667}*cell[15];
auto x47 = x21 + x32;
auto x48 = x27*x47;
auto x49 = x23 + x32;
auto x50 = x27*x49;
auto x51 = x45 - x46 + x48 + x50;
auto x52 = x19*(x21*x22 + x22*x23 + x31 + x37 + x44 + x51 - V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] + V{4.62592926927149e-18});
auto x53 = V{0.0555555555555556}*rho;
auto x54 = x25 + V{1};
auto x55 = x28 + x54;
auto x56 = x29 + x32;
auto x57 = -x27*x56 + V{0.166666666666667}*cell[17] - V{0.166666666666667}*cell[8];
auto x58 = x27*(x33 + x54) + x57;
auto x59 = V{0.166666666666667}*cell[18];
auto x60 = V{0.166666666666667}*cell[9];
auto x61 = x32 + x54;
auto x62 = x27*x61;
auto x63 = x29 + x33;
auto x64 = x27*x63;
auto x65 = x59 - x60 + x62 - x64;
auto x66 = -x22*x29 + x31 - x38 + x39 - x41 + x65 + V{0.166666666666667}*cell[11] - V{0.166666666666667}*cell[2];
auto x67 = x24 + V{-1};
auto x68 = x33 + V{1};
auto x69 = x28 + x68;
auto x70 = -x22*x35 + x37 - x45 + x46 - x48 - x59 + x60 - x62 + x64 + V{0.166666666666667}*cell[12] - V{0.166666666666667}*cell[3];
auto x71 = x32 + V{-1};
auto x72 = V{0.0833333333333333}*cell[2];
auto x73 = V{0.0833333333333333}*cell[11];
auto x74 = x27*x29;
auto x75 = x27*x67;
auto x76 = V{0.0833333333333333}*cell[7];
auto x77 = V{0.0833333333333333}*cell[16];
auto x78 = V{0.00231481481481481}*rho;
auto x79 = x34*x78;
auto x80 = x36*x78;
auto x81 = x76 - x77 + x79 - x80;
auto x82 = V{0.0833333333333333}*cell[1];
auto x83 = -V{0.0833333333333333}*cell[10];
auto x84 = x21*x27;
auto x85 = x23*x27;
auto x86 = V{0.0833333333333333}*cell[6];
auto x87 = V{0.0833333333333333}*cell[15];
auto x88 = x47*x78;
auto x89 = x49*x78;
auto x90 = x86 - x87 + x88 + x89;
auto x91 = x82 + x83 + x84 + x85 + x90 + V{2.31296463463574e-18};
auto x92 = V{0.0833333333333333}*cell[8];
auto x93 = V{0.0833333333333333}*cell[17];
auto x94 = x56*x78;
auto x95 = x32 + x67;
auto x96 = x78*x95;
auto x97 = x92 - x93 + x94 + x96;
auto x98 = V{0.0833333333333333}*cell[9];
auto x99 = V{0.0833333333333333}*cell[18];
auto x100 = x63*x78;
auto x101 = x61*x78;
auto x102 = x100 - x101 + x98 - x99;
auto x103 = x19*(x102 + x44 + x72 - x73 + x74 + x75 + x81 + x91 + x97);
auto x104 = V{0.0277777777777778}*rho;
auto x105 = -x75;
auto x106 = -x72;
auto x107 = -x74;
auto x108 = x106 + x107 + x73 + x81;
auto x109 = -x100 + x101 - x98 + x99;
auto x110 = -x92 + x93 - x94 - x96;
auto x111 = x19*(x105 + x108 + x109 + x110 + x31 + x91);
auto x112 = V{0.0833333333333333}*cell[3];
auto x113 = V{0.0833333333333333}*cell[12];
auto x114 = x27*x35;
auto x115 = x27*x71;
auto x116 = x26*x78 - x30*x78 - V{0.0833333333333333}*cell[14] + V{0.0833333333333333}*cell[5];
auto x117 = x112 - x113 + x114 + x115 + x116;
auto x118 = V{0.0833333333333333}*cell[4];
auto x119 = V{0.0833333333333333}*cell[13];
auto x120 = x40*x78;
auto x121 = x42*x78;
auto x122 = x118 - x119 + x120 + x121 + x82 + x83 + x84 + x85 + V{2.31296463463574e-18};
auto x123 = x19*(x109 + x117 + x122 + x51 + x97);
auto x124 = -x115;
auto x125 = -x112 + x113 - x114 + x116;
auto x126 = x19*(x102 + x110 + x122 + x124 + x125 + x37);
auto x127 = -x118 + x119 - x120;
auto x128 = x108 + x125 + x127 - x86 + x87 - x88;
auto x129 = x105 - x121;
auto x130 = x19*(x106 + x107 + x117 + x127 + x129 + x65 + x73 - x76 + x77 - x79 + x80 + x90);
auto x131 = -x27*x95 + x57;
cell[0] = V{0.333333333333333}*rho + V{-0.333333333333333};
cell[1] = -x23*x53 - x52 + V{-0.0555555555555556};
cell[2] = x19*(x22*x54 + x27*x55 + x58 + x66) - x53*x67 + V{-0.0555555555555556};
cell[3] = x19*(x22*x68 + x27*x69 + x58 + x70) - x53*x71 + V{-0.0555555555555556};
cell[4] = -x103 - x104*x42 + V{-0.0277777777777778};
cell[5] = x104*x30 - x111 + V{-0.0277777777777778};
cell[6] = -x104*x49 - x123 + V{-0.0277777777777778};
cell[7] = x104*x36 - x126 + V{-0.0277777777777778};
cell[8] = -x104*x95 + x19*(x128 + x27*x54 + x27*x68 + x55*x78 + x58 + x69*x78) + V{-0.0277777777777778};
cell[9] = x104*x61 + x130 + V{-0.0277777777777778};
cell[10] = x21*x53 + x52 + V{-0.0555555555555556};
cell[11] = -x19*(x131 - x22*x67 - x43 + x66) + x29*x53 + V{-0.0555555555555556};
cell[12] = -x19*(x131 - x22*x71 - x50 + x70) + x35*x53 + V{-0.0555555555555556};
cell[13] = x103 + x104*x40 + V{-0.0277777777777778};
cell[14] = x104*x26 + x111 + V{-0.0277777777777778};
cell[15] = x104*x47 + x123 + V{-0.0277777777777778};
cell[16] = x104*x34 + x126 + V{-0.0277777777777778};
cell[17] = x104*x56 - x19*(x124 + x128 + x129 + x131 - x89) + V{-0.0277777777777778};
cell[18] = x104*x63 - x130 + V{-0.0277777777777778};
return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

template <typename CELL, typename RHO, typename U, typename PI, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, PI& pi, OMEGA& omega) any_platform
{
auto x19 = omega + V{-1};
auto x20 = u[0]*u[0];
auto x21 = V{1.5}*x20;
auto x22 = u[1]*u[1];
auto x23 = V{1.5}*x22;
auto x24 = u[2]*u[2];
auto x25 = V{1.5}*x24;
auto x26 = x23 + x25 + V{-1};
auto x27 = x21 + x26;
auto x28 = V{0.0555555555555556}*rho;
auto x29 = V{3}*u[0];
auto x30 = V{3}*x20;
auto x31 = V{0.0833333333333333}*pi[3];
auto x32 = V{0.0833333333333333}*pi[5];
auto x33 = x19*(x31 + x32 - V{0.166666666666667}*pi[0]) + V{-0.0555555555555556};
auto x34 = V{3}*u[1];
auto x35 = V{3}*x22;
auto x36 = x21 + V{-1};
auto x37 = V{0.0833333333333333}*pi[0];
auto x38 = x19*(x32 + x37 - V{0.166666666666667}*pi[3]) + V{-0.0555555555555556};
auto x39 = V{3}*u[2];
auto x40 = V{3}*x24;
auto x41 = x19*(x31 + x37 - V{0.166666666666667}*pi[5]) + V{-0.0555555555555556};
auto x42 = V{0.25}*pi[1];
auto x43 = V{0.0833333333333333}*pi[0];
auto x44 = V{0.0833333333333333}*pi[3];
auto x45 = x43 + x44 - V{0.0416666666666667}*pi[5];
auto x46 = x19*(x42 + x45);
auto x47 = V{0.0277777777777778}*rho;
auto x48 = u[0] + u[1];
auto x49 = V{4.5}*(x48*x48);
auto x50 = x27 + x29;
auto x51 = -x34;
auto x52 = -u[0];
auto x53 = x52 + u[1];
auto x54 = x19*(-x42 + x45) + V{0.0277777777777778};
auto x55 = V{0.25}*pi[2];
auto x56 = V{0.0833333333333333}*pi[5];
auto x57 = x43 + x56 - V{0.0416666666666667}*pi[3];
auto x58 = x19*(x55 + x57);
auto x59 = u[0] + u[2];
auto x60 = V{4.5}*(x59*x59);
auto x61 = -x39;
auto x62 = x52 + u[2];
auto x63 = x19*(-x55 + x57) + V{0.0277777777777778};
auto x64 = V{0.25}*pi[4];
auto x65 = V{0.0416666666666667}*pi[0];
auto x66 = x19*(x44 + x56 + x64 - x65);
auto x67 = u[1] + u[2];
auto x68 = V{4.5}*(x67*x67);
auto x69 = x27 + x34;
auto x70 = -u[1];
auto x71 = x70 + u[2];
auto x72 = x19*(-x44 - x56 + x64 + x65) + V{-0.0277777777777778};
auto x73 = -x23;
auto x74 = V{1} - x25;
auto x75 = x73 + x74;
auto x76 = x29 + x75;
auto x77 = -x21;
auto x78 = x34 + x77;
auto x79 = x39 + x77;
auto x80 = -x29;
auto x81 = x70 + u[0];
auto x82 = -u[2];
auto x83 = x82 + u[0];
auto x84 = x27 + x39;
auto x85 = x82 + u[1];
cell[0] = -V{0.333333333333333}*rho*x27 + V{0.5}*x19*(pi[0] + pi[3] + pi[5]) + V{-0.333333333333333};
cell[1] = -x28*(x26 + x29 - x30) + x33;
cell[2] = -x28*(x25 + x34 - x35 + x36) + x38;
cell[3] = -x28*(x23 + x36 + x39 - x40) + x41;
cell[4] = -x46 - x47*(x34 - x49 + x50) + V{-0.0277777777777778};
cell[5] = -x47*(x50 + x51 - V{4.5}*x53*x53) - x54;
cell[6] = -x47*(x39 + x50 - x60) - x58 + V{-0.0277777777777778};
cell[7] = -x47*(x50 + x61 - V{4.5}*x62*x62) - x63;
cell[8] = -x47*(x39 - x68 + x69) - x66 + V{-0.0277777777777778};
cell[9] = -x47*(x61 + x69 - V{4.5}*x71*x71) + x72;
cell[10] = x28*(x30 + x76) + x33;
cell[11] = x28*(x35 + x74 + x78) + x38;
cell[12] = x28*(x40 + x73 + x79 + V{1}) + x41;
cell[13] = -x46 + x47*(x49 + x76 + x78) + V{-0.0277777777777778};
cell[14] = -x47*(x69 + x80 - V{4.5}*x81*x81) - x54;
cell[15] = x47*(x60 + x76 + x79) - x58 + V{-0.0277777777777778};
cell[16] = -x47*(x80 + x84 - V{4.5}*x83*x83) - x63;
cell[17] = x47*(x39 + x68 + x75 + x78) - x66 + V{-0.0277777777777778};
cell[18] = -x47*(x51 + x84 - V{4.5}*x85*x85) + x72;
return x20 + x22 + x24;
}

template <typename CELL, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineEqFirstOrder(CELL& cell, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x19 = V{3}*newU[0];
auto x20 = x19 + V{-1};
auto x21 = V{3}*newU[1];
auto x22 = x21 + V{-1};
auto x23 = V{3}*newU[2];
auto x24 = -x19;
auto x25 = x21 + V{1};
auto x26 = x23 + V{1};
auto x27 = -x21;
auto x28 = x19 + V{1};
auto x29 = -x23;
cell[0] = V{0.333333333333333}*newRho + V{-0.333333333333333};
cell[1] = -V{0.0555555555555556}*newRho*x20 + V{-0.0555555555555556};
cell[2] = -V{0.0555555555555556}*newRho*x22 + V{-0.0555555555555556};
cell[3] = -V{0.0555555555555556}*newRho*(x23 + V{-1}) + V{-0.0555555555555556};
cell[4] = -V{0.0277777777777778}*newRho*(x20 + x21) + V{-0.0277777777777778};
cell[5] = V{0.0277777777777778}*newRho*(x24 + x25) + V{-0.0277777777777778};
cell[6] = -V{0.0277777777777778}*newRho*(x20 + x23) + V{-0.0277777777777778};
cell[7] = V{0.0277777777777778}*newRho*(x24 + x26) + V{-0.0277777777777778};
cell[8] = -V{0.0277777777777778}*newRho*(x22 + x23) + V{-0.0277777777777778};
cell[9] = V{0.0277777777777778}*newRho*(x26 + x27) + V{-0.0277777777777778};
cell[10] = V{0.0555555555555556}*newRho*x28 + V{-0.0555555555555556};
cell[11] = V{0.0555555555555556}*newRho*x25 + V{-0.0555555555555556};
cell[12] = V{0.0555555555555556}*newRho*x26 + V{-0.0555555555555556};
cell[13] = V{0.0277777777777778}*newRho*(x21 + x28) + V{-0.0277777777777778};
cell[14] = V{0.0277777777777778}*newRho*(x27 + x28) + V{-0.0277777777777778};
cell[15] = V{0.0277777777777778}*newRho*(x23 + x28) + V{-0.0277777777777778};
cell[16] = V{0.0277777777777778}*newRho*(x28 + x29) + V{-0.0277777777777778};
cell[17] = V{0.0277777777777778}*newRho*(x23 + x25) + V{-0.0277777777777778};
cell[18] = V{0.0277777777777778}*newRho*(x25 + x29) + V{-0.0277777777777778};

}

template <typename CELL, typename OLDRHO, typename OLDU, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineNEq(CELL& cell, OLDRHO& oldRho, OLDU& oldU, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x19 = oldU[0]*oldU[0];
auto x20 = V{1.5}*x19;
auto x21 = oldU[1]*oldU[1];
auto x22 = V{1.5}*x21;
auto x23 = oldU[2]*oldU[2];
auto x24 = V{1.5}*x23;
auto x25 = x22 + x24 + V{-1};
auto x26 = x20 + x25;
auto x27 = newU[0]*newU[0];
auto x28 = V{1.5}*x27;
auto x29 = newU[1]*newU[1];
auto x30 = V{1.5}*x29;
auto x31 = newU[2]*newU[2];
auto x32 = V{1.5}*x31;
auto x33 = x30 + x32 + V{-1};
auto x34 = x28 + x33;
auto x35 = V{0.0555555555555556}*oldRho;
auto x36 = V{3}*oldU[0];
auto x37 = V{3}*x19;
auto x38 = V{0.0555555555555556}*newRho;
auto x39 = V{3}*newU[0];
auto x40 = V{3}*x27;
auto x41 = V{3}*oldU[1];
auto x42 = V{3}*x21;
auto x43 = x20 + V{-1};
auto x44 = V{3}*newU[1];
auto x45 = V{3}*x29;
auto x46 = x28 + V{-1};
auto x47 = V{3}*oldU[2];
auto x48 = V{3}*x23;
auto x49 = V{3}*newU[2];
auto x50 = V{3}*x31;
auto x51 = V{0.0277777777777778}*oldRho;
auto x52 = oldU[0] + oldU[1];
auto x53 = V{4.5}*(x52*x52);
auto x54 = x26 + x36;
auto x55 = V{0.0277777777777778}*newRho;
auto x56 = newU[0] + newU[1];
auto x57 = V{4.5}*(x56*x56);
auto x58 = x34 + x39;
auto x59 = -x41;
auto x60 = oldU[0] - oldU[1];
auto x61 = -V{4.5}*x60*x60;
auto x62 = -x44;
auto x63 = newU[0] - newU[1];
auto x64 = -V{4.5}*x63*x63;
auto x65 = oldU[0] + oldU[2];
auto x66 = V{4.5}*(x65*x65);
auto x67 = newU[0] + newU[2];
auto x68 = V{4.5}*(x67*x67);
auto x69 = -x47;
auto x70 = -oldU[2];
auto x71 = x70 + oldU[0];
auto x72 = -V{4.5}*x71*x71;
auto x73 = -x49;
auto x74 = -newU[2];
auto x75 = x74 + newU[0];
auto x76 = -V{4.5}*x75*x75;
auto x77 = oldU[1] + oldU[2];
auto x78 = V{4.5}*(x77*x77);
auto x79 = x26 + x41;
auto x80 = newU[1] + newU[2];
auto x81 = V{4.5}*(x80*x80);
auto x82 = x34 + x44;
auto x83 = x70 + oldU[1];
auto x84 = -V{4.5}*x83*x83;
auto x85 = x74 + newU[1];
auto x86 = -V{4.5}*x85*x85;
auto x87 = -x30;
auto x88 = V{1} - x32;
auto x89 = x87 + x88;
auto x90 = x39 + x89;
auto x91 = -x22;
auto x92 = V{1} - x24;
auto x93 = x91 + x92;
auto x94 = x36 + x93;
auto x95 = -x28;
auto x96 = x44 + x95;
auto x97 = -x20;
auto x98 = x41 + x97;
auto x99 = x49 + x95;
auto x100 = x47 + x97;
auto x101 = -x36;
auto x102 = -x39;
auto x103 = x26 + x47;
auto x104 = x34 + x49;
cell[0] = -V{0.333333333333333}*newRho*x34 + V{0.333333333333333}*oldRho*x26 + cell[0];
cell[1] = x35*(x25 + x36 - x37) - x38*(x33 + x39 - x40) + cell[1];
cell[2] = x35*(x24 + x41 - x42 + x43) - x38*(x32 + x44 - x45 + x46) + cell[2];
cell[3] = x35*(x22 + x43 + x47 - x48) - x38*(x30 + x46 + x49 - x50) + cell[3];
cell[4] = x51*(x41 - x53 + x54) - x55*(x44 - x57 + x58) + cell[4];
cell[5] = x51*(x54 + x59 + x61) - x55*(x58 + x62 + x64) + cell[5];
cell[6] = x51*(x47 + x54 - x66) - x55*(x49 + x58 - x68) + cell[6];
cell[7] = x51*(x54 + x69 + x72) - x55*(x58 + x73 + x76) + cell[7];
cell[8] = x51*(x47 - x78 + x79) - x55*(x49 - x81 + x82) + cell[8];
cell[9] = x51*(x69 + x79 + x84) - x55*(x73 + x82 + x86) + cell[9];
cell[10] = -x35*(x37 + x94) + x38*(x40 + x90) + cell[10];
cell[11] = -x35*(x42 + x92 + x98) + x38*(x45 + x88 + x96) + cell[11];
cell[12] = -x35*(x100 + x48 + x91 + V{1}) + x38*(x50 + x87 + x99 + V{1}) + cell[12];
cell[13] = -x51*(x53 + x94 + x98) + x55*(x57 + x90 + x96) + cell[13];
cell[14] = x51*(x101 + x61 + x79) - x55*(x102 + x64 + x82) + cell[14];
cell[15] = -x51*(x100 + x66 + x94) + x55*(x68 + x90 + x99) + cell[15];
cell[16] = x51*(x101 + x103 + x72) - x55*(x102 + x104 + x76) + cell[16];
cell[17] = -x51*(x47 + x78 + x93 + x98) + x55*(x49 + x81 + x89 + x96) + cell[17];
cell[18] = x51*(x103 + x59 + x84) - x55*(x104 + x62 + x86) + cell[18];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void defineNEqFromPi(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x19 = u[0]*u[0];
auto x20 = V{1.5}*x19;
auto x21 = u[1]*u[1];
auto x22 = V{1.5}*x21;
auto x23 = u[2]*u[2];
auto x24 = V{1.5}*x23;
auto x25 = x22 + x24 + V{-1};
auto x26 = x20 + x25;
auto x27 = V{0.0555555555555556}*rho;
auto x28 = V{3}*u[0];
auto x29 = V{3}*x19;
auto x30 = -V{0.0833333333333333}*pi[3];
auto x31 = -V{0.0833333333333333}*pi[5] + V{-0.0555555555555556};
auto x32 = x30 + x31 + V{0.166666666666667}*pi[0];
auto x33 = V{3}*u[1];
auto x34 = V{3}*x21;
auto x35 = x20 + V{-1};
auto x36 = -V{0.0833333333333333}*pi[0];
auto x37 = x31 + x36 + V{0.166666666666667}*pi[3];
auto x38 = V{3}*u[2];
auto x39 = V{3}*x23;
auto x40 = x30 + x36 + V{0.166666666666667}*pi[5] + V{-0.0555555555555556};
auto x41 = V{0.0277777777777778}*rho;
auto x42 = u[0] + u[1];
auto x43 = V{4.5}*(x42*x42);
auto x44 = x26 + x28;
auto x45 = V{0.25}*pi[1];
auto x46 = V{0.0833333333333333}*pi[3];
auto x47 = V{0.0833333333333333}*pi[0] + V{-0.0277777777777778};
auto x48 = x46 + x47 - V{0.0416666666666667}*pi[5];
auto x49 = x45 + x48;
auto x50 = -x33;
auto x51 = -u[0];
auto x52 = x51 + u[1];
auto x53 = -x45 + x48;
auto x54 = u[0] + u[2];
auto x55 = V{4.5}*(x54*x54);
auto x56 = V{0.25}*pi[2];
auto x57 = V{0.0833333333333333}*pi[5];
auto x58 = x47 + x57 - V{0.0416666666666667}*pi[3];
auto x59 = x56 + x58;
auto x60 = -x38;
auto x61 = x51 + u[2];
auto x62 = -x56 + x58;
auto x63 = u[1] + u[2];
auto x64 = V{4.5}*(x63*x63);
auto x65 = x26 + x33;
auto x66 = V{0.25}*pi[4];
auto x67 = x46 + x57 - V{0.0416666666666667}*pi[0] + V{-0.0277777777777778};
auto x68 = x66 + x67;
auto x69 = -u[1];
auto x70 = x69 + u[2];
auto x71 = -x66 + x67;
auto x72 = -x22;
auto x73 = V{1} - x24;
auto x74 = x72 + x73;
auto x75 = x28 + x74;
auto x76 = -x20;
auto x77 = x33 + x76;
auto x78 = x38 + x76;
auto x79 = -x28;
auto x80 = x69 + u[0];
auto x81 = -u[2];
auto x82 = x81 + u[0];
auto x83 = x26 + x38;
auto x84 = x81 + u[1];
cell[0] = -V{0.333333333333333}*rho*x26 - V{0.5}*pi[0] - V{0.5}*pi[3] - V{0.5}*pi[5] + V{-0.333333333333333};
cell[1] = -x27*(x25 + x28 - x29) + x32;
cell[2] = -x27*(x24 + x33 - x34 + x35) + x37;
cell[3] = -x27*(x22 + x35 + x38 - x39) + x40;
cell[4] = -x41*(x33 - x43 + x44) + x49;
cell[5] = -x41*(x44 + x50 - V{4.5}*x52*x52) + x53;
cell[6] = -x41*(x38 + x44 - x55) + x59;
cell[7] = -x41*(x44 + x60 - V{4.5}*x61*x61) + x62;
cell[8] = -x41*(x38 - x64 + x65) + x68;
cell[9] = -x41*(x60 + x65 - V{4.5}*x70*x70) + x71;
cell[10] = x27*(x29 + x75) + x32;
cell[11] = x27*(x34 + x73 + x77) + x37;
cell[12] = x27*(x39 + x72 + x78 + V{1}) + x40;
cell[13] = x41*(x43 + x75 + x77) + x49;
cell[14] = -x41*(x65 + x79 - V{4.5}*x80*x80) + x53;
cell[15] = x41*(x55 + x75 + x78) + x59;
cell[16] = -x41*(x79 + x83 - V{4.5}*x82*x82) + x62;
cell[17] = x41*(x38 + x64 + x74 + x77) + x68;
cell[18] = -x41*(x50 + x83 - V{4.5}*x84*x84) + x71;

}

template <typename CELL, typename FORCE, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell, FORCE& force) any_platform
{
auto x0 = cell[15] - cell[6];
auto x1 = cell[13] - cell[4];
auto x2 = cell[10] + cell[14] + cell[16];
auto x3 = x0 + x1 + x2 - cell[1] - cell[5] - cell[7];
auto x4 = cell[17] - cell[8];
auto x5 = cell[12] + cell[7] + cell[9];
auto x6 = x0 + x4 + x5 - cell[16] - cell[18] - cell[3];
auto x7 = cell[11] + cell[18] + cell[5];
auto x8 = x2 + x5 + x7 + cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8];
auto x9 = V{1} / (x8 + V{1});
auto x10 = x9*(x8 + V{1});
auto x11 = V{0.5}*x10;
auto x12 = x3*x9;
auto x13 = V{1}*x12;
auto x14 = x11*(x3*force[2] + x6*force[0]) - x13*x6 + V{1}*cell[15] - V{1}*cell[16] + V{1}*cell[6] - V{1}*cell[7];
auto x15 = x1 + x4 + x7 - cell[14] - cell[2] - cell[9];
auto x16 = x15*force[0] + x3*force[1];
auto x17 = V{1}*x10;
auto x18 = V{2}*x15;
auto x19 = x15*force[2] + x6*force[1];
auto x20 = x6*x9;
auto x21 = -V{0.333333333333333}*cell[0];
auto x22 = x21 - V{0.333333333333333}*cell[12] + V{0.666666666666667}*cell[13] + V{0.666666666666667}*cell[14] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5];
auto x23 = -V{0.333333333333333}*cell[11] + V{0.666666666666667}*cell[15] + V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x24 = x10*x3*force[0] + x22 + x23 - x9*x3*x3 + V{0.666666666666667}*cell[10] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9];
auto x25 = -V{0.333333333333333}*cell[10] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
auto x26 = x10*x15*force[1] + x22 + x25 - x9*x15*x15 + V{0.666666666666667}*cell[11] - V{0.333333333333333}*cell[15] - V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7];
auto x27 = x10*x6*force[2] + x21 + x23 + x25 - x9*x6*x6 + V{0.666666666666667}*cell[12] - V{0.333333333333333}*cell[13] - V{0.333333333333333}*cell[14] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5];
return (x11*x16 - x13*x15 + V{1}*cell[13] - V{1}*cell[14] + V{1}*cell[4] - V{1}*cell[5])*(-x12*x18 + x16*x17 + V{2}*cell[13] - V{2}*cell[14] + V{2}*cell[4] - V{2}*cell[5]) + (x11*x19 - V{1}*x15*x20 + V{1}*cell[17] - V{1}*cell[18] + V{1}*cell[8] - V{1}*cell[9])*(x17*x19 - x18*x20 + V{2}*cell[17] - V{2}*cell[18] + V{2}*cell[8] - V{2}*cell[9]) + 2*(x14*x14) + x24*x24 + x26*x26 + x27*x27;
}

template <typename CELL, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell) any_platform
{
auto x0 = -cell[4];
auto x1 = -cell[8];
auto x2 = x1 + cell[18];
auto x3 = cell[11] + cell[17] + cell[5];
auto x4 = x0 + x2 + x3 + cell[13] - cell[14] - cell[2] - cell[9];
auto x5 = cell[10] + cell[13] + cell[15];
auto x6 = cell[12] + cell[7] + cell[9];
auto x7 = V{1} / (x3 + x5 + x6 + cell[0] + cell[14] + cell[16] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + V{1});
auto x8 = x0 + cell[14];
auto x9 = -cell[6];
auto x10 = x9 + cell[16];
auto x11 = x10 + x5 + x8 - cell[1] - cell[5] - cell[7];
auto x12 = x11*x7;
auto x13 = x12*x4 + x8 - cell[13] + cell[5];
auto x14 = x1 + x6 + x9 + cell[15] - cell[16] + cell[17] - cell[18] - cell[3];
auto x15 = x10 + x12*x14 - cell[15] + cell[7];
auto x16 = x14*x4*x7 + x2 - cell[17] + cell[9];
auto x17 = V{1}*x7;
auto x18 = V{0.333333333333333}*cell[0];
auto x19 = x18 + V{0.333333333333333}*cell[10] - V{0.666666666666667}*cell[17] - V{0.666666666666667}*cell[18] + V{0.333333333333333}*cell[1] - V{0.666666666666667}*cell[8] - V{0.666666666666667}*cell[9];
auto x20 = V{0.333333333333333}*cell[11] - V{0.666666666666667}*cell[15] - V{0.666666666666667}*cell[16] + V{0.333333333333333}*cell[2] - V{0.666666666666667}*cell[6] - V{0.666666666666667}*cell[7];
auto x21 = x17*(x14*x14) + x19 + x20 - V{0.666666666666667}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5];
auto x22 = V{0.333333333333333}*cell[12] - V{0.666666666666667}*cell[13] - V{0.666666666666667}*cell[14] + V{0.333333333333333}*cell[3] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[5];
auto x23 = x17*(x4*x4) + x19 + x22 - V{0.666666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7];
auto x24 = x17*(x11*x11) + x18 + x20 + x22 - V{0.666666666666667}*cell[10] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9];
return V{2}*(x13*x13) + V{2}*(x15*x15) + V{2}*(x16*x16) + x21*x21 + x23*x23 + x24*x24;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addExternalForce(CELL& cell, RHO& rho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x19 = force[0]*u[0];
auto x20 = force[1]*u[1];
auto x21 = force[2]*u[2];
auto x22 = rho*(V{0.5}*omega + V{-1});
auto x23 = V{6}*u[0];
auto x24 = x23 + V{-3};
auto x25 = V{0.0555555555555556}*force[0];
auto x26 = V{0.166666666666667}*x20;
auto x27 = V{0.166666666666667}*x21;
auto x28 = x26 + x27;
auto x29 = V{6}*u[1];
auto x30 = x29 + V{-3};
auto x31 = V{0.0555555555555556}*force[1];
auto x32 = V{0.166666666666667}*x19;
auto x33 = x27 + x32;
auto x34 = V{6}*u[2];
auto x35 = x34 + V{-3};
auto x36 = V{0.0555555555555556}*force[2];
auto x37 = x26 + x32;
auto x38 = V{9}*u[1];
auto x39 = V{0.0277777777777778}*force[0];
auto x40 = V{9}*u[0];
auto x41 = V{0.0277777777777778}*force[1];
auto x42 = -V{0.0833333333333333}*x21;
auto x43 = -x40;
auto x44 = x29 + V{3};
auto x45 = V{3} - x23;
auto x46 = V{9}*u[2];
auto x47 = V{0.0277777777777778}*force[2];
auto x48 = -V{0.0833333333333333}*x20;
auto x49 = x34 + V{3};
auto x50 = -V{0.0833333333333333}*x19;
auto x51 = -x38;
auto x52 = V{3} - x29;
auto x53 = x23 + V{3};
auto x54 = -x46;
auto x55 = V{3} - x34;
cell[0] = V{1}*x22*(x19 + x20 + x21) + cell[0];
cell[1] = x22*(-x24*x25 + x28) + cell[1];
cell[2] = x22*(-x30*x31 + x33) + cell[2];
cell[3] = x22*(-x35*x36 + x37) + cell[3];
cell[4] = -x22*(x39*(x24 + x38) + x41*(x30 + x40) + x42) + cell[4];
cell[5] = -x22*(-x39*(x38 + x45) + x41*(x43 + x44) + x42) + cell[5];
cell[6] = -x22*(x39*(x24 + x46) + x47*(x35 + x40) + x48) + cell[6];
cell[7] = -x22*(-x39*(x45 + x46) + x47*(x43 + x49) + x48) + cell[7];
cell[8] = -x22*(x41*(x30 + x46) + x47*(x35 + x38) + x50) + cell[8];
cell[9] = -x22*(-x41*(x46 + x52) + x47*(x49 + x51) + x50) + cell[9];
cell[10] = x22*(-x25*x53 + x28) + cell[10];
cell[11] = x22*(-x31*x44 + x33) + cell[11];
cell[12] = x22*(-x36*x49 + x37) + cell[12];
cell[13] = -x22*(x39*(x38 + x53) + x41*(x40 + x44) + x42) + cell[13];
cell[14] = -x22*(x39*(x51 + x53) - x41*(x40 + x52) + x42) + cell[14];
cell[15] = -x22*(x39*(x46 + x53) + x47*(x40 + x49) + x48) + cell[15];
cell[16] = -x22*(x39*(x53 + x54) - x47*(x40 + x55) + x48) + cell[16];
cell[17] = -x22*(x41*(x44 + x46) + x47*(x38 + x49) + x50) + cell[17];
cell[18] = -x22*(x41*(x44 + x54) - x47*(x38 + x55) + x50) + cell[18];

}

};

template <typename... FIELDS>
struct lbm<descriptors::D3Q27<FIELDS...>> {

template <typename CELL, typename V=typename CELL::value_t>
static auto computeRho(CELL& cell) any_platform
{

return cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[19] + cell[1] + cell[20] + cell[21] + cell[22] + cell[23] + cell[24] + cell[25] + cell[26] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
}

template <typename CELL, typename J, typename V=typename CELL::value_t>
static void computeJ(CELL& cell, J& j) any_platform
{
auto x0 = -cell[10];
auto x1 = x0 - cell[11] + cell[17] + cell[23] + cell[24] - cell[4];
auto x2 = -cell[12] + cell[19] + cell[25] - cell[6];
auto x3 = cell[13] + cell[21] - cell[26] - cell[8];
j[0] = V{1}*x1 + V{1}*x2 - V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[18] - V{1}*cell[1] + V{1}*cell[20] + V{1}*cell[26] - V{1}*cell[5] - V{1}*cell[7];
j[1] = V{1}*x1 + V{1}*x3 + V{1}*cell[12] + V{1}*cell[15] - V{1}*cell[18] + V{1}*cell[22] - V{1}*cell[25] - V{1}*cell[2] + V{1}*cell[5] - V{1}*cell[9];
j[2] = V{1}*x0 + V{1}*x2 + V{1}*x3 + V{1}*cell[11] + V{1}*cell[16] - V{1}*cell[20] - V{1}*cell[22] + V{1}*cell[23] - V{1}*cell[24] - V{1}*cell[3] + V{1}*cell[7] + V{1}*cell[9];

}

template <typename CELL, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeRhoU(CELL& cell, RHO& rho, U& u) any_platform
{
auto x0 = cell[13] + cell[1] + cell[5] + cell[7];
auto x1 = cell[18] + cell[25] + cell[2] + cell[9];
auto x2 = cell[10] + cell[20] + cell[22] + cell[24] + cell[3];
auto x3 = x0 + x1 + x2 + cell[0] + cell[11] + cell[12] + cell[14] + cell[15] + cell[16] + cell[17] + cell[19] + cell[21] + cell[23] + cell[26] + cell[4] + cell[6] + cell[8] + V{1};
auto x4 = -cell[23];
auto x5 = x4 + cell[10] + cell[11] - cell[17] - cell[24] + cell[4];
auto x6 = cell[12] - cell[19] - cell[25] + cell[6];
auto x7 = V{1}/x3;
auto x8 = -cell[13] - cell[21] + cell[26] + cell[8];
rho = x3;
u[0] = -x7*(x0 + x5 + x6 - cell[14] - cell[18] - cell[20] - cell[26]);
u[1] = -x7*(x1 + x5 + x8 - cell[12] - cell[15] - cell[22] - cell[5]);
u[2] = -x7*(x2 + x4 + x6 + x8 - cell[11] - cell[16] - cell[7] - cell[9]);

}

template <typename CELL, typename RHO, typename J, typename V=typename CELL::value_t>
static void computeRhoJ(CELL& cell, RHO& rho, J& j) any_platform
{
auto x0 = cell[14] + cell[18] + cell[20] + cell[26];
auto x1 = cell[12] + cell[15] + cell[22] + cell[5];
auto x2 = cell[11] + cell[16] + cell[23] + cell[7] + cell[9];
auto x3 = -cell[10];
auto x4 = x3 - cell[11] + cell[17] + cell[23] + cell[24] - cell[4];
auto x5 = -cell[12] + cell[19] + cell[25] - cell[6];
auto x6 = cell[13] + cell[21] - cell[26] - cell[8];
rho = x0 + x1 + x2 + cell[0] + cell[10] + cell[13] + cell[17] + cell[19] + cell[1] + cell[21] + cell[24] + cell[25] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + V{1};
j[0] = V{1}*x0 + V{1}*x4 + V{1}*x5 - V{1}*cell[13] - V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7];
j[1] = V{1}*x1 + V{1}*x4 + V{1}*x6 - V{1}*cell[18] - V{1}*cell[25] - V{1}*cell[2] - V{1}*cell[9];
j[2] = V{1}*x2 + V{1}*x3 + V{1}*x5 + V{1}*x6 - V{1}*cell[20] - V{1}*cell[22] - V{1}*cell[24] - V{1}*cell[3];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeStress(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = V{1}*cell[10];
auto x1 = V{1}*cell[23];
auto x2 = V{1}*cell[13];
auto x3 = V{1}*cell[26];
auto x4 = x0 + x1 + x2 + x3;
auto x5 = V{1}*cell[5];
auto x6 = V{1}*cell[18];
auto x7 = -V{0.333333333333333}*rho;
auto x8 = V{1}*cell[11];
auto x9 = V{1}*cell[24];
auto x10 = x8 + x9;
auto x11 = x10 + V{1}*cell[17] + V{1}*cell[4];
auto x12 = x11 + x5 + x6 + x7 + V{0.333333333333333};
auto x13 = V{1}*cell[7];
auto x14 = V{1}*cell[20];
auto x15 = V{1}*cell[12];
auto x16 = V{1}*cell[25];
auto x17 = x15 + x16;
auto x18 = x17 + V{1}*cell[19] + V{1}*cell[6];
auto x19 = x13 + x14 + x18;
auto x20 = rho*u[0];
auto x21 = x0 + x1 - x2 - x3;
auto x22 = -x15 - x16;
auto x23 = -x8 - x9;
auto x24 = V{1}*cell[9];
auto x25 = V{1}*cell[22];
auto x26 = x4 + V{1}*cell[21] + V{1}*cell[8];
auto x27 = x24 + x25 + x26;
pi[0] = -rho*u[0]*u[0] + x12 + x19 + x4 + V{1}*cell[14] + V{1}*cell[1];
pi[1] = x11 - x20*u[1] + x21 + x22 - x5 - x6;
pi[2] = -x13 - x14 + x18 - x20*u[2] + x21 + x23;
pi[3] = -rho*u[1]*u[1] + x12 + x17 + x27 + V{1}*cell[15] + V{1}*cell[2];
pi[4] = -rho*u[1]*u[2] + x22 + x23 - x24 - x25 + x26;
pi[5] = -rho*u[2]*u[2] + x10 + x19 + x27 + x7 + V{1}*cell[16] + V{1}*cell[3] + V{0.333333333333333};

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void computeAllMomenta(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x0 = cell[17] + cell[24];
auto x1 = cell[19] + cell[25];
auto x2 = cell[10] + cell[13];
auto x3 = x2 + cell[21] + cell[23];
auto x4 = cell[1] + cell[5] + cell[7];
auto x5 = cell[11] + cell[18] + cell[2] + cell[4] + cell[9];
auto x6 = cell[12] + cell[20] + cell[22] + cell[26] + cell[3] + cell[6] + cell[8];
auto x7 = x0 + x1 + x3 + x4 + x5 + x6 + cell[0] + cell[14] + cell[15] + cell[16] + V{1};
auto x8 = -cell[26];
auto x9 = cell[11] - cell[18] + cell[4];
auto x10 = cell[12] - cell[20] + cell[6];
auto x11 = -cell[24];
auto x12 = -cell[23];
auto x13 = x11 + x12 - cell[17];
auto x14 = -cell[25];
auto x15 = x14 - cell[19];
auto x16 = x10 + x13 + x15 + x2 + x4 + x8 + x9 - cell[14];
auto x17 = V{1} / (x7);
auto x18 = V{1}*x17;
auto x19 = -cell[5];
auto x20 = -cell[12];
auto x21 = -cell[22] + cell[26] + cell[8];
auto x22 = -cell[13];
auto x23 = x22 + cell[10] - cell[21];
auto x24 = x13 + x19 + x20 + x21 + x23 + x5 - cell[15] + cell[25];
auto x25 = -cell[7];
auto x26 = -cell[11];
auto x27 = -cell[9];
auto x28 = x12 + x15 + x23 + x25 + x26 + x27 + x6 - cell[16] + cell[24];
auto x29 = V{0.666666666666667}*cell[10];
auto x40 = V{0.666666666666667}*cell[11];
auto x41 = V{0.666666666666667}*cell[12];
auto x42 = V{0.666666666666667}*cell[13];
auto x43 = V{0.666666666666667}*cell[23];
auto x44 = V{0.666666666666667}*cell[24];
auto x45 = V{0.666666666666667}*cell[25];
auto x46 = V{0.666666666666667}*cell[26];
auto x47 = -V{0.333333333333333}*cell[0];
auto x48 = x29 + x40 + x41 + x42 + x43 + x44 + x45 + x46 + x47 - V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5];
auto x49 = -V{0.333333333333333}*cell[15] + V{0.666666666666667}*cell[19] + V{0.666666666666667}*cell[20] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x50 = x16*x17;
auto x51 = x22 + x8 + cell[10] + cell[23];
auto x52 = x14 + x20;
auto x53 = x11 + x26;
auto x54 = -V{0.333333333333333}*cell[14] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[21] + V{0.666666666666667}*cell[22] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
rho = x7;
u[0] = -x16*x18;
u[1] = -x18*x24;
u[2] = -x18*x28;
pi[0] = -x18*x16*x16 + x48 + x49 + V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[21] - V{0.333333333333333}*cell[22] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9];
pi[1] = V{1}*x0 + V{1}*x19 - V{1}*x24*x50 + V{1}*x51 + V{1}*x52 + V{1}*x9;
pi[2] = V{1}*x1 + V{1}*x10 + V{1}*x25 - V{1}*x28*x50 + V{1}*x51 + V{1}*x53;
pi[3] = -x18*x24*x24 + x48 + x54 + V{0.666666666666667}*cell[15] - V{0.333333333333333}*cell[19] - V{0.333333333333333}*cell[20] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7];
pi[4] = -V{1}*x17*x24*x28 + V{1}*x21 + V{1}*x27 + V{1}*x3 + V{1}*x52 + V{1}*x53;
pi[5] = -x18*x28*x28 + x29 + x40 + x41 + x42 + x43 + x44 + x45 + x46 + x47 + x49 + x54 + V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5];

}

template <typename CELL, typename FEQ, typename V=typename CELL::value_t>
static void computeFeq(CELL& cell, FEQ& fEq) any_platform
{
auto x0 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[19] + cell[1] + cell[20] + cell[21] + cell[22] + cell[23] + cell[24] + cell[25] + cell[26] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x1 = x0 + V{1};
auto x2 = x0 + V{1};
auto x3 = V{1} / ((x2)*(x2));
auto x4 = V{1.5}*x3;
auto x5 = -cell[4];
auto x6 = -cell[6];
auto x7 = x5 + x6 + cell[17] + cell[19];
auto x8 = -cell[7];
auto x9 = cell[18] - cell[5];
auto x10 = x8 + x9 + cell[20];
auto x11 = -cell[11] + cell[24];
auto x12 = -cell[10] + cell[23];
auto x13 = x11 + x12;
auto x14 = cell[14] - cell[1];
auto x15 = -cell[12] + cell[25];
auto x16 = -cell[13] + cell[26];
auto x17 = x15 + x16;
auto x18 = x10 + x13 + x14 + x17 + x7;
auto x19 = x18*x18;
auto x20 = x19*x4;
auto x21 = -cell[8];
auto x22 = cell[13] - cell[26];
auto x23 = x21 + x22 + cell[21];
auto x24 = -cell[9];
auto x25 = cell[12] - cell[25];
auto x26 = x24 + x25 + cell[22];
auto x27 = x5 + cell[17];
auto x28 = -cell[2];
auto x29 = -cell[18] + cell[5];
auto x57 = x28 + x29 + cell[15];
auto x58 = x27 + x57;
auto x59 = x13 + x23 + x26 + x58;
auto x60 = x59*x59;
auto x61 = x4*x60;
auto x62 = -cell[3];
auto x63 = x15 + x62 + cell[16];
auto x64 = cell[11] - cell[24];
auto x65 = -cell[22] + cell[9];
auto x66 = x64 + x65;
auto x67 = x6 + cell[19];
auto x68 = -cell[20] + cell[7];
auto x69 = x67 + x68;
auto x70 = x12 + x23 + x63 + x66 + x69;
auto x71 = x70*x70;
auto x72 = x4*x71;
auto x73 = x61 + x72 + V{-1};
auto x74 = x20 + x73;
auto x75 = V{1} / (x2);
auto x76 = V{3}*cell[14];
auto x77 = V{3}*cell[1];
auto x78 = V{3}*cell[17];
auto x79 = V{3}*cell[23];
auto x80 = V{3}*cell[24];
auto x81 = V{3}*cell[4];
auto x82 = V{3}*cell[10];
auto x83 = -x82;
auto x84 = V{3}*cell[11];
auto x85 = -x84;
auto x86 = x78 + x79 + x80 - x81 + x83 + x85;
auto x87 = V{3}*cell[19];
auto x88 = V{3}*cell[25];
auto x89 = V{3}*cell[6];
auto x90 = V{3}*cell[12];
auto x91 = -x90;
auto x92 = x87 + x88 - x89 + x91;
auto x93 = V{3}*cell[18];
auto x94 = V{3}*cell[26];
auto x95 = V{3}*cell[5];
auto x96 = V{3}*cell[13];
auto x97 = -x96;
auto x98 = x93 + x94 - x95 + x97;
auto x99 = V{3}*cell[20];
auto x100 = V{3}*cell[7];
auto x101 = -x100 + x99;
auto x102 = x75*(x101 + x76 - x77 + x86 + x92 + x98);
auto x103 = V{3}*x3;
auto x104 = V{3}*cell[15];
auto x105 = V{3}*cell[2];
auto x106 = V{3}*cell[21];
auto x107 = V{3}*cell[8];
auto x108 = -x94;
auto x109 = x106 - x107 + x108 + x96;
auto x110 = V{3}*cell[22];
auto x111 = V{3}*cell[9];
auto x112 = -x88;
auto x113 = x110 - x111 + x112 + x90;
auto x114 = -x93 + x95;
auto x115 = x75*(x104 - x105 + x109 + x113 + x114 + x86);
auto x116 = x20 + V{-1};
auto x117 = V{3}*cell[16];
auto x118 = V{3}*cell[3];
auto x119 = -x80;
auto x120 = -x110 + x111 + x119 + x84;
auto x121 = x100 - x99;
auto x122 = x75*(x109 + x117 - x118 + x120 + x121 + x79 + x83 + x92);
auto x123 = V{4.5}*x3;
auto x124 = x8 + cell[20];
auto x125 = x124 + x14 + x67;
auto x126 = x21 + cell[21];
auto x127 = V{2}*cell[10];
auto x128 = V{2}*cell[23];
auto x129 = -x127 + x128;
auto x130 = x126 + x129;
auto x131 = x24 + cell[22];
auto x132 = V{2}*cell[11];
auto x133 = V{2}*cell[24];
auto x134 = -x132 + x133;
auto x135 = x131 + x134;
auto x136 = V{2}*cell[4];
auto x137 = V{2}*cell[17];
auto x138 = x28 + cell[15];
auto x139 = -x136 + x137 + x138;
auto x140 = x125 + x130 + x135 + x139;
auto x141 = x102 + x74;
auto x142 = x115 + x141;
auto x143 = V{2}*cell[26];
auto x144 = -x143;
auto x145 = V{2}*cell[13];
auto x146 = x144 + x145 + x68;
auto x147 = V{2}*cell[25];
auto x148 = -x147;
auto x149 = V{2}*cell[12];
auto x150 = -cell[19] + cell[6];
auto x151 = x148 + x149 + x150;
auto x152 = -cell[14] + cell[1];
auto x153 = x126 + x152;
auto x154 = V{2}*cell[18];
auto x155 = V{2}*cell[5];
auto x156 = x138 - x154 + x155;
auto x157 = x131 + x146 + x151 + x153 + x156;
auto x158 = -x115;
auto x159 = x141 + x158;
auto x160 = x62 + cell[16];
auto x161 = x160 + x65;
auto x162 = x147 - x149;
auto x163 = x14 + x27 + x9;
auto x164 = V{2}*cell[6];
auto x165 = V{2}*cell[19];
auto x166 = -x164 + x165;
auto x167 = x130 + x161 + x162 + x163 + x166;
auto x168 = -x122;
auto x169 = x132 - x133;
auto x170 = V{2}*cell[20];
auto x171 = V{2}*cell[7];
auto x172 = -x170 + x171;
auto x173 = -cell[17] + cell[4];
auto x174 = x173 + x29;
auto x175 = x144 + x145 + x153 + x161 + x169 + x172 + x174;
auto x176 = V{2}*cell[8];
auto x177 = V{2}*cell[21];
auto x178 = -x176 + x177;
auto x179 = x160 + x178;
auto x180 = x129 + x146 + x179 + x57 + x7;
auto x181 = x115 + x74;
auto x182 = x122 + x181;
auto x183 = -cell[15] + cell[2];
auto x184 = x169 + x183;
auto x185 = x173 + x9;
auto x186 = V{2}*cell[22];
auto x187 = V{2}*cell[9];
auto x188 = -x186 + x187;
auto x189 = x160 + x188;
auto x190 = x162 + x184 + x185 + x189 + x69;
auto x191 = V{3}*cell[10];
auto x192 = V{3}*cell[23];
auto x193 = x22 + x63;
auto x194 = x139 + x14;
auto x195 = x11 + x166 + x178 - x191 + x192 + x193 + x194;
auto x196 = V{3}*cell[24];
auto x197 = V{3}*cell[11];
auto x198 = cell[10] - cell[23];
auto x199 = x136 - x137;
auto x200 = x152 + x172;
auto x201 = x183 + x188 + x193 - x196 + x197 + x198 + x199 + x200;
auto x202 = V{3}*cell[25];
auto x203 = V{3}*cell[12];
auto x204 = x152 + x164 - x165;
auto x205 = -cell[16];
auto x206 = x205 + cell[3];
auto x207 = x186 - x187 + x206;
auto x208 = x11 + x198;
auto x209 = x156 - x202 + x203 + x204 + x207 + x208 + x22;
auto x210 = -x79;
auto x211 = -x78 + x81;
auto x212 = -x87 + x89;
auto x213 = x75*(x108 + x112 + x114 + x119 + x121 + x210 + x211 + x212 - x76 + x77 + x82 + x84 + x90 + x96);
auto x214 = V{3}*cell[26];
auto x215 = V{3}*cell[13];
auto x216 = x25 + x64;
auto x217 = x12 + x156 + x179 + x200 - x214 + x215 + x216;
auto x218 = -x106 + x107 + x210 + x82;
auto x219 = x75*(x101 + x113 - x117 + x118 + x212 + x218 + x80 + x85 + x94 + x97);
auto x220 = -cell[21] + cell[8];
auto x221 = x183 + x220;
auto x222 = x17 + x185 + x198 + x221 + x66;
auto x223 = x222*x222;
auto x224 = x223*x4;
auto x225 = x150 + x198;
auto x226 = x205 + x220 + cell[3];
auto x227 = x11 + x124 + x16 + x225 + x226 + x26;
auto x228 = x227*x227;
auto x229 = x228*x4 + V{-1};
auto x230 = x224 + x229;
auto x231 = x75*(-x104 + x105 + x120 + x211 + x218 + x88 + x91 + x98);
auto x232 = x152 + x68;
auto x233 = x174 + x216 + x22 + x225 + x232;
auto x234 = x233*x233;
auto x235 = x234*x4;
auto x236 = x231 + x235;
auto x237 = x219 + x230 + x236;
auto x238 = x213 + x230;
auto x239 = x219 + x235;
auto x240 = x127 - x128;
auto x241 = x150 + x240;
auto x242 = x184 + x199 + x220 + x232 + x241 + x65;
auto x243 = x236 + x238;
auto x244 = -x102;
auto x245 = x143 - x145;
auto x246 = x154 - x155;
auto x247 = x125 + x162 + x221 + x245 + x246 + x65;
auto x248 = x131 + x148 + x149 + x174 + x204 + x226 + x240;
auto x249 = x238 + x239;
auto x250 = x170 - x171;
auto x251 = x135 + x163 + x226 + x245 + x250;
auto x252 = x122 + x74;
auto x253 = x176 - x177 + x183 + x206;
auto x254 = x10 + x173 + x241 + x245 + x253;
auto x255 = x124 + x134 + x151 + x207 + x58;
auto x256 = x16 + x191 - x192 + x199 + x204 + x216 + x253;
auto x257 = x12 + x16;
auto x258 = x194 + x196 - x197 + x207 + x25 + x250 + x257;
auto x259 = x14 + x246;
auto x260 = x166 + x183 + x189 + x202 - x203 + x257 + x259 + x64;
auto x261 = x15 + x208 + x214 - x215 + x250 + x253 + x259;
fEq[0] = -V{0.296296296296296}*x1*x74 + V{-0.296296296296296};
fEq[1] = -V{0.0740740740740741}*x1*(x102 - x103*x19 + x73) + V{-0.0740740740740741};
fEq[2] = -V{0.0740740740740741}*x1*(-x103*x60 + x115 + x116 + x72) + V{-0.0740740740740741};
fEq[3] = -V{0.0740740740740741}*x1*(-x103*x71 + x116 + x122 + x61) + V{-0.0740740740740741};
fEq[4] = -V{0.0185185185185185}*(x1*(-x123*x140*x140 + x142) + V{1});
fEq[5] = -V{0.0185185185185185}*(x1*(-x123*x157*x157 + x159) + V{1});
fEq[6] = -V{0.0185185185185185}*(x1*(x122 - x123*x167*x167 + x141) + V{1});
fEq[7] = -V{0.0185185185185185}*(x1*(-x123*x175*x175 + x141 + x168) + V{1});
fEq[8] = -V{0.0185185185185185}*(x1*(-x123*x180*x180 + x182) + V{1});
fEq[9] = -V{0.0185185185185185}*(x1*(-x123*x190*x190 + x168 + x181) + V{1});
fEq[10] = -V{0.00462962962962963}*(x1*(x122 - x123*x195*x195 + x142) + V{1});
fEq[11] = -V{0.00462962962962963}*(x1*(-x123*x201*x201 + x142 + x168) + V{1});
fEq[12] = -V{0.00462962962962963}*(x1*(x122 - x123*x209*x209 + x159) + V{1});
fEq[13] = -V{0.00462962962962963}*(x1*(-x123*x217*x217 - x213 + x237) + V{1});
fEq[14] = -V{0.0740740740740741}*x1*(-x103*x234 + x238) + V{-0.0740740740740741};
fEq[15] = -V{0.0740740740740741}*x1*(-x103*x223 + x229 + x236) + V{-0.0740740740740741};
fEq[16] = -V{0.0740740740740741}*x1*(-x103*x228 + x224 + x239 + V{-1}) + V{-0.0740740740740741};
fEq[17] = -V{0.0185185185185185}*(x1*(-x123*x242*x242 + x243) + V{1});
fEq[18] = -V{0.0185185185185185}*(x1*(-x123*x247*x247 + x181 + x244) + V{1});
fEq[19] = -V{0.0185185185185185}*(x1*(-x123*x248*x248 + x249) + V{1});
fEq[20] = -V{0.0185185185185185}*(x1*(-x123*x251*x251 + x244 + x252) + V{1});
fEq[21] = -V{0.0185185185185185}*(x1*(-x123*x254*x254 + x237) + V{1});
fEq[22] = -V{0.0185185185185185}*(x1*(-x123*x255*x255 + x158 + x252) + V{1});
fEq[23] = -V{0.00462962962962963}*(x1*(-x123*x256*x256 + x219 + x243) + V{1});
fEq[24] = -V{0.00462962962962963}*(x1*(-x123*x258*x258 - x219 + x243) + V{1});
fEq[25] = -V{0.00462962962962963}*(x1*(-x123*x260*x260 - x231 + x249) + V{1});
fEq[26] = -V{0.00462962962962963}*(x1*(-x123*x261*x261 + x182 + x244) + V{1});

}

template <typename CELL, typename FNEQ, typename RHO, typename U, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq, RHO& rho, U& u) any_platform
{
auto x0 = u[0]*u[0];
auto x1 = V{1.5}*x0;
auto x2 = u[1]*u[1];
auto x3 = V{1.5}*x2;
auto x4 = u[2]*u[2];
auto x5 = V{1.5}*x4;
auto x6 = x3 + x5 + V{-1};
auto x7 = x1 + x6;
auto x8 = V{0.0740740740740741}*rho;
auto x9 = V{3}*u[0];
auto x10 = V{3}*x0;
auto x11 = V{3}*u[1];
auto x12 = V{3}*x2;
auto x13 = x1 + V{-1};
auto x14 = V{3}*u[2];
auto x15 = V{3}*x4;
auto x16 = V{0.0185185185185185}*rho;
auto x17 = u[0] + u[1];
auto x18 = V{4.5}*(x17*x17);
auto x19 = x7 + x9;
auto x20 = x11 + x19;
auto x21 = u[0] - u[1];
auto x22 = -V{4.5}*x21*x21;
auto x23 = -x11;
auto x24 = x19 + x23;
auto x25 = u[0] + u[2];
auto x26 = V{4.5}*(x25*x25);
auto x27 = -x14;
auto x28 = -u[2];
auto x29 = x28 + u[0];
auto x57 = -V{4.5}*x29*x29;
auto x58 = u[1] + u[2];
auto x59 = V{4.5}*(x58*x58);
auto x60 = x11 + x7;
auto x61 = x14 + x60;
auto x62 = x28 + u[1];
auto x63 = -V{4.5}*x62*x62;
auto x64 = V{0.00462962962962963}*rho;
auto x65 = x17 + u[2];
auto x66 = V{4.5}*(x65*x65);
auto x67 = x17 + x28;
auto x68 = V{4.5}*(x67*x67);
auto x69 = x21 + u[2];
auto x70 = V{4.5}*(x69*x69);
auto x71 = -x9;
auto x72 = x58 - u[0];
auto x73 = V{4.5}*(x72*x72);
auto x74 = -x3;
auto x75 = V{1} - x5;
auto x76 = x74 + x75;
auto x77 = -x1;
auto x78 = x11 + x77;
auto x79 = x14 + x76 + x78;
auto x80 = x76 + x9;
auto x81 = x14 + x77;
auto x82 = x78 + x80;
auto x83 = x80 + x81;
auto x84 = x14 + x7;
fNeq[0] = V{0.296296296296296}*rho*x7 + cell[0] + V{0.296296296296296};
fNeq[1] = x8*(-x10 + x6 + x9) + cell[1] + V{0.0740740740740741};
fNeq[2] = x8*(x11 - x12 + x13 + x5) + cell[2] + V{0.0740740740740741};
fNeq[3] = x8*(x13 + x14 - x15 + x3) + cell[3] + V{0.0740740740740741};
fNeq[4] = x16*(-x18 + x20) + cell[4] + V{0.0185185185185185};
fNeq[5] = x16*(x22 + x24) + cell[5] + V{0.0185185185185185};
fNeq[6] = x16*(x14 + x19 - x26) + cell[6] + V{0.0185185185185185};
fNeq[7] = x16*(x19 + x27 + x57) + cell[7] + V{0.0185185185185185};
fNeq[8] = x16*(-x59 + x61) + cell[8] + V{0.0185185185185185};
fNeq[9] = x16*(x27 + x60 + x63) + cell[9] + V{0.0185185185185185};
fNeq[10] = x64*(x14 + x20 - x66) + cell[10] + V{0.00462962962962963};
fNeq[11] = x64*(x20 + x27 - x68) + cell[11] + V{0.00462962962962963};
fNeq[12] = x64*(x14 + x24 - x70) + cell[12] + V{0.00462962962962963};
fNeq[13] = -x64*(x71 + x73 + x79) + cell[13] + V{0.00462962962962963};
fNeq[14] = -x8*(x10 + x80) + cell[14] + V{0.0740740740740741};
fNeq[15] = -x8*(x12 + x75 + x78) + cell[15] + V{0.0740740740740741};
fNeq[16] = -x8*(x15 + x74 + x81 + V{1}) + cell[16] + V{0.0740740740740741};
fNeq[17] = -x16*(x18 + x82) + cell[17] + V{0.0185185185185185};
fNeq[18] = x16*(x22 + x60 + x71) + cell[18] + V{0.0185185185185185};
fNeq[19] = -x16*(x26 + x83) + cell[19] + V{0.0185185185185185};
fNeq[20] = x16*(x57 + x71 + x84) + cell[20] + V{0.0185185185185185};
fNeq[21] = -x16*(x59 + x79) + cell[21] + V{0.0185185185185185};
fNeq[22] = x16*(x23 + x63 + x84) + cell[22] + V{0.0185185185185185};
fNeq[23] = -x64*(x14 + x66 + x82) + cell[23] + V{0.00462962962962963};
fNeq[24] = -x64*(x27 + x68 + x82) + cell[24] + V{0.00462962962962963};
fNeq[25] = -x64*(x23 + x70 + x83) + cell[25] + V{0.00462962962962963};
fNeq[26] = x64*(x61 + x71 - x73) + cell[26] + V{0.00462962962962963};

}

template <typename CELL, typename FNEQ, typename V=typename CELL::value_t>
static void computeFneq(CELL& cell, FNEQ& fNeq) any_platform
{
auto x0 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[19] + cell[1] + cell[20] + cell[21] + cell[22] + cell[23] + cell[24] + cell[25] + cell[26] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
auto x1 = V{1} / ((x0)*(x0));
auto x2 = V{1.5}*x1;
auto x3 = -cell[23];
auto x4 = cell[12] - cell[25];
auto x5 = -cell[19] + cell[6];
auto x6 = x4 + x5;
auto x7 = -cell[17];
auto x8 = x7 + cell[4];
auto x9 = cell[13] - cell[26];
auto x10 = cell[11] - cell[24];
auto x11 = -cell[14] + cell[1];
auto x12 = x10 + x11;
auto x13 = -cell[20];
auto x14 = -cell[18] + cell[5];
auto x15 = x13 + x14 + cell[7];
auto x16 = x12 + x15 + x3 + x6 + x8 + x9 + cell[10];
auto x17 = x16*x16;
auto x18 = x17*x2;
auto x19 = -cell[12] + cell[25];
auto x20 = -cell[5];
auto x21 = -cell[22];
auto x22 = x20 + x21 + cell[18] + cell[9];
auto x23 = -cell[21] + cell[8];
auto x24 = x23 + x7 + cell[4];
auto x25 = -cell[13] + cell[26];
auto x26 = x25 + x3 + cell[10];
auto x27 = -cell[15] + cell[2];
auto x28 = x10 + x27;
auto x29 = x19 + x22 + x24 + x26 + x28;
auto x57 = x29*x29;
auto x58 = x2*x57;
auto x59 = -cell[11] + cell[24];
auto x60 = -cell[7];
auto x61 = x23 + x60 + cell[20];
auto x62 = -cell[16];
auto x63 = cell[22] - cell[9];
auto x64 = x62 + x63 + cell[3];
auto x65 = x26 + x59 + x6 + x61 + x64;
auto x66 = x65*x65;
auto x67 = x2*x66;
auto x68 = x58 + x67 + V{-1};
auto x69 = x18 + x68;
auto x70 = V{0.0740740740740741}*cell[0] + V{0.0740740740740741}*cell[10] + V{0.0740740740740741}*cell[11] + V{0.0740740740740741}*cell[12] + V{0.0740740740740741}*cell[13] + V{0.0740740740740741}*cell[14] + V{0.0740740740740741}*cell[15] + V{0.0740740740740741}*cell[16] + V{0.0740740740740741}*cell[17] + V{0.0740740740740741}*cell[18] + V{0.0740740740740741}*cell[19] + V{0.0740740740740741}*cell[1] + V{0.0740740740740741}*cell[20] + V{0.0740740740740741}*cell[21] + V{0.0740740740740741}*cell[22] + V{0.0740740740740741}*cell[23] + V{0.0740740740740741}*cell[24] + V{0.0740740740740741}*cell[25] + V{0.0740740740740741}*cell[26] + V{0.0740740740740741}*cell[2] + V{0.0740740740740741}*cell[3] + V{0.0740740740740741}*cell[4] + V{0.0740740740740741}*cell[5] + V{0.0740740740740741}*cell[6] + V{0.0740740740740741}*cell[7] + V{0.0740740740740741}*cell[8] + V{0.0740740740740741}*cell[9] + V{0.0740740740740741};
auto x71 = V{1} / (x0);
auto x72 = V{3}*cell[5];
auto x73 = V{3}*cell[7];
auto x74 = V{3}*cell[13];
auto x75 = V{3}*cell[18];
auto x76 = V{3}*cell[20];
auto x77 = V{3}*cell[26];
auto x78 = V{3}*cell[10];
auto x79 = V{3}*cell[11];
auto x80 = -V{3}*cell[23];
auto x81 = V{3}*cell[24];
auto x82 = x78 + x79 + x80 - x81 - V{3}*cell[17] + V{3}*cell[4];
auto x83 = V{3}*cell[12];
auto x84 = V{3}*cell[25];
auto x85 = x83 - x84 - V{3}*cell[19] + V{3}*cell[6];
auto x86 = x71*(x72 + x73 + x74 - x75 - x76 - x77 + x82 + x85 - V{3}*cell[14] + V{3}*cell[1]);
auto x87 = -x86;
auto x88 = V{3}*x1;
auto x89 = -x17*x88 + x68;
auto x90 = V{3}*cell[9];
auto x91 = V{3}*cell[22];
auto x92 = -x74 + x77 - V{3}*cell[21] + V{3}*cell[8];
auto x93 = x71*(-x72 + x75 + x82 - x83 + x84 + x90 - x91 + x92 - V{3}*cell[15] + V{3}*cell[2]);
auto x94 = -x93;
auto x95 = x18 + V{-1};
auto x96 = -x57*x88 + x67 + x95;
auto x97 = x71*(-x73 + x76 + x78 - x79 + x80 + x81 + x85 - x90 + x91 + x92 - V{3}*cell[16] + V{3}*cell[3]);
auto x98 = -x97;
auto x99 = x58 - x66*x88 + x95;
auto x100 = V{0.0185185185185185}*cell[0] + V{0.0185185185185185}*cell[10] + V{0.0185185185185185}*cell[11] + V{0.0185185185185185}*cell[12] + V{0.0185185185185185}*cell[13] + V{0.0185185185185185}*cell[14] + V{0.0185185185185185}*cell[15] + V{0.0185185185185185}*cell[16] + V{0.0185185185185185}*cell[17] + V{0.0185185185185185}*cell[18] + V{0.0185185185185185}*cell[19] + V{0.0185185185185185}*cell[1] + V{0.0185185185185185}*cell[20] + V{0.0185185185185185}*cell[21] + V{0.0185185185185185}*cell[22] + V{0.0185185185185185}*cell[23] + V{0.0185185185185185}*cell[24] + V{0.0185185185185185}*cell[25] + V{0.0185185185185185}*cell[26] + V{0.0185185185185185}*cell[2] + V{0.0185185185185185}*cell[3] + V{0.0185185185185185}*cell[4] + V{0.0185185185185185}*cell[5] + V{0.0185185185185185}*cell[6] + V{0.0185185185185185}*cell[7] + V{0.0185185185185185}*cell[8] + V{0.0185185185185185}*cell[9] + V{0.0185185185185185};
auto x101 = V{4.5}*x1;
auto x102 = -cell[6];
auto x103 = cell[14] - cell[1];
auto x104 = x102 + x103 + cell[19];
auto x105 = V{2}*cell[10];
auto x106 = -x105;
auto x107 = V{2}*cell[23];
auto x108 = x106 + x107 + cell[21] - cell[8];
auto x109 = x60 + cell[20];
auto x110 = cell[15] - cell[2];
auto x111 = V{2}*cell[11];
auto x112 = V{2}*cell[24];
auto x113 = -x111 + x112;
auto x114 = x109 + x110 + x113;
auto x115 = V{2}*cell[4];
auto x116 = V{2}*cell[17];
auto x117 = -x115 + x116;
auto x118 = x104 + x108 + x114 + x117 + x63;
auto x119 = x69 + x87;
auto x120 = x119 + x94;
auto x121 = V{2}*cell[12];
auto x122 = V{2}*cell[25];
auto x123 = -x121 + x122;
auto x124 = V{2}*cell[13];
auto x125 = V{2}*cell[26];
auto x126 = -x124 + x125;
auto x127 = V{2}*cell[5];
auto x128 = V{2}*cell[18];
auto x129 = -x127 + x128;
auto x130 = x21 + x27 + cell[9];
auto x131 = x104 + x123 + x126 + x129 + x130 + x61;
auto x132 = -x101*x131*x131;
auto x133 = x119 + x93;
auto x134 = cell[17] - cell[4];
auto x135 = x103 + x134;
auto x136 = V{2}*cell[6];
auto x137 = V{2}*cell[19];
auto x138 = cell[16] - cell[3];
auto x139 = -x136 + x137 + x138;
auto x140 = x108 + x123 + x135 + x139 + x22;
auto x141 = V{2}*cell[7];
auto x142 = V{2}*cell[20];
auto x143 = -x141 + x142;
auto x144 = x126 + x20 + cell[18];
auto x145 = x113 + x135 + x143 + x144 + x23 + x64;
auto x146 = -x101*x145*x145;
auto x147 = V{2}*cell[8];
auto x148 = V{2}*cell[21];
auto x149 = x110 - x147 + x148;
auto x150 = x138 + x149;
auto x151 = x102 + x106 + x107 + x124 - x125 + x134 + x15 + x150 + cell[19];
auto x152 = x69 + x94;
auto x153 = x152 + x98;
auto x154 = x121 - x122 + x14;
auto x155 = x62 + cell[3];
auto x156 = x155 + x5;
auto x157 = V{2}*cell[9];
auto x158 = V{2}*cell[22];
auto x159 = -x157 + x158;
auto x160 = x114 + x134 + x154 + x156 + x159;
auto x161 = -x101*x160*x160;
auto x162 = x152 + x97;
auto x163 = V{0.00462962962962963}*cell[0] + V{0.00462962962962963}*cell[10] + V{0.00462962962962963}*cell[11] + V{0.00462962962962963}*cell[12] + V{0.00462962962962963}*cell[13] + V{0.00462962962962963}*cell[14] + V{0.00462962962962963}*cell[15] + V{0.00462962962962963}*cell[16] + V{0.00462962962962963}*cell[17] + V{0.00462962962962963}*cell[18] + V{0.00462962962962963}*cell[19] + V{0.00462962962962963}*cell[1] + V{0.00462962962962963}*cell[20] + V{0.00462962962962963}*cell[21] + V{0.00462962962962963}*cell[22] + V{0.00462962962962963}*cell[23] + V{0.00462962962962963}*cell[24] + V{0.00462962962962963}*cell[25] + V{0.00462962962962963}*cell[26] + V{0.00462962962962963}*cell[2] + V{0.00462962962962963}*cell[3] + V{0.00462962962962963}*cell[4] + V{0.00462962962962963}*cell[5] + V{0.00462962962962963}*cell[6] + V{0.00462962962962963}*cell[7] + V{0.00462962962962963}*cell[8] + V{0.00462962962962963}*cell[9] + V{0.00462962962962963};
auto x164 = V{3}*cell[10];
auto x165 = V{3}*cell[23];
auto x166 = x103 + x117;
auto x167 = x139 + x149 - x164 + x165 + x166 + x19 + x59 + x9;
auto x168 = -cell[10];
auto x169 = x168 + x25 + cell[23];
auto x170 = x155 + x4;
auto x171 = x110 + x143 + x159 + x166 + x169 + x170 - V{3}*cell[11] + V{3}*cell[24];
auto x172 = -x101*x171*x171;
auto x173 = x103 + x129 + x139 + x157 - x158 + x169 + x28 - V{3}*cell[12] + V{3}*cell[25];
auto x174 = -x101*x173*x173;
auto x175 = x12 + x127 - x128 + x141 - x142 + x150 + x168 + x4 + V{3}*cell[13] + cell[23] - V{3}*cell[26];
auto x176 = -x101*x175*x175;
auto x177 = x105 - x107;
auto x178 = x11 + x177;
auto x179 = x115 - x116;
auto x180 = x111 - x112 + x13 + x130 + x178 + x179 + x23 + x5 + cell[7];
auto x181 = x69 + x86;
auto x182 = x181 + x93;
auto x183 = x136 - x137;
auto x184 = x154 + x178 + x183 + x24 + x64;
auto x185 = x147 - x148;
auto x186 = x109 + x144 + x156 + x177 + x185 + x27 + x8;
auto x187 = x69 + x93;
auto x188 = x11 + x164 - x165 + x170 + x179 + x183 + x185 + x25 + x28;
fNeq[0] = x69*(V{0.296296296296296}*cell[0] + V{0.296296296296296}*cell[10] + V{0.296296296296296}*cell[11] + V{0.296296296296296}*cell[12] + V{0.296296296296296}*cell[13] + V{0.296296296296296}*cell[14] + V{0.296296296296296}*cell[15] + V{0.296296296296296}*cell[16] + V{0.296296296296296}*cell[17] + V{0.296296296296296}*cell[18] + V{0.296296296296296}*cell[19] + V{0.296296296296296}*cell[1] + V{0.296296296296296}*cell[20] + V{0.296296296296296}*cell[21] + V{0.296296296296296}*cell[22] + V{0.296296296296296}*cell[23] + V{0.296296296296296}*cell[24] + V{0.296296296296296}*cell[25] + V{0.296296296296296}*cell[26] + V{0.296296296296296}*cell[2] + V{0.296296296296296}*cell[3] + V{0.296296296296296}*cell[4] + V{0.296296296296296}*cell[5] + V{0.296296296296296}*cell[6] + V{0.296296296296296}*cell[7] + V{0.296296296296296}*cell[8] + V{0.296296296296296}*cell[9] + V{0.296296296296296}) + cell[0] + V{0.296296296296296};
fNeq[1] = x70*(x87 + x89) + cell[1] + V{0.0740740740740741};
fNeq[2] = x70*(x94 + x96) + cell[2] + V{0.0740740740740741};
fNeq[3] = x70*(x98 + x99) + cell[3] + V{0.0740740740740741};
fNeq[4] = x100*(-x101*x118*x118 + x120) + cell[4] + V{0.0185185185185185};
fNeq[5] = x100*(x132 + x133) + cell[5] + V{0.0185185185185185};
fNeq[6] = x100*(-x101*x140*x140 + x119 + x98) + cell[6] + V{0.0185185185185185};
fNeq[7] = x100*(x119 + x146 + x97) + cell[7] + V{0.0185185185185185};
fNeq[8] = x100*(-x101*x151*x151 + x153) + cell[8] + V{0.0185185185185185};
fNeq[9] = x100*(x161 + x162) + cell[9] + V{0.0185185185185185};
fNeq[10] = x163*(-x101*x167*x167 + x120 + x98) + cell[10] + V{0.00462962962962963};
fNeq[11] = x163*(x120 + x172 + x97) + cell[11] + V{0.00462962962962963};
fNeq[12] = x163*(x133 + x174 + x98) + cell[12] + V{0.00462962962962963};
fNeq[13] = x163*(x133 + x176 + x97) + cell[13] + V{0.00462962962962963};
fNeq[14] = x70*(x86 + x89) + cell[14] + V{0.0740740740740741};
fNeq[15] = x70*(x93 + x96) + cell[15] + V{0.0740740740740741};
fNeq[16] = x70*(x97 + x99) + cell[16] + V{0.0740740740740741};
fNeq[17] = x100*(-x101*x180*x180 + x182) + cell[17] + V{0.0185185185185185};
fNeq[18] = x100*(x132 + x152 + x86) + cell[18] + V{0.0185185185185185};
fNeq[19] = x100*(-x101*x184*x184 + x181 + x97) + cell[19] + V{0.0185185185185185};
fNeq[20] = x100*(x146 + x181 + x98) + cell[20] + V{0.0185185185185185};
fNeq[21] = x100*(-x101*x186*x186 + x187 + x97) + cell[21] + V{0.0185185185185185};
fNeq[22] = x100*(x161 + x187 + x98) + cell[22] + V{0.0185185185185185};
fNeq[23] = x163*(-x101*x188*x188 + x182 + x97) + cell[23] + V{0.00462962962962963};
fNeq[24] = x163*(x172 + x182 + x98) + cell[24] + V{0.00462962962962963};
fNeq[25] = x163*(x162 + x174 + x86) + cell[25] + V{0.00462962962962963};
fNeq[26] = x163*(x153 + x176 + x86) + cell[26] + V{0.00462962962962963};

}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto bgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x27 = omega + V{-1};
auto x28 = u[0]*u[0];
auto x29 = V{1.5}*x28;
auto x30 = u[1]*u[1];
auto x31 = V{1.5}*x30;
auto x32 = u[2]*u[2];
auto x33 = V{1.5}*x32;
auto x34 = x31 + x33 + V{-1};
auto x35 = x29 + x34;
auto x36 = V{0.0740740740740741}*omega;
auto x37 = V{3}*u[0];
auto x38 = V{3}*x28;
auto x39 = V{3}*u[1];
auto x40 = V{3}*x30;
auto x41 = x29 + V{-1};
auto x42 = V{3}*u[2];
auto x43 = V{3}*x32;
auto x44 = V{0.0185185185185185}*omega;
auto x45 = u[0] + u[1];
auto x46 = V{4.5}*(x45*x45);
auto x47 = x35 + x37;
auto x48 = x39 + x47;
auto x49 = -u[0];
auto x50 = x49 + u[1];
auto x51 = -x39;
auto x52 = x47 + x51;
auto x53 = u[0] + u[2];
auto x54 = V{4.5}*(x53*x53);
auto x55 = -x42;
auto x56 = x49 + u[2];
auto x57 = u[1] + u[2];
auto x58 = V{4.5}*(x57*x57);
auto x59 = x35 + x39;
auto x60 = x42 + x59;
auto x61 = -u[1];
auto x62 = x61 + u[2];
auto x63 = V{0.00462962962962963}*omega;
auto x64 = x45 + u[2];
auto x65 = V{4.5}*(x64*x64);
auto x66 = x56 + x61;
auto x67 = -u[2];
auto x68 = x50 + x67;
auto x69 = -x37;
auto x70 = x50 + u[2];
auto x71 = -x31;
auto x72 = V{1} - x33;
auto x73 = x71 + x72;
auto x74 = -x29;
auto x75 = x39 + x74;
auto x76 = x42 + x73 + x75;
auto x77 = x37 + x73;
auto x78 = x42 + x74;
auto x79 = x75 + x77;
auto x80 = x61 + u[0];
auto x81 = x77 + x78;
auto x82 = x67 + u[0];
auto x83 = x35 + x42;
auto x84 = x67 + u[1];
auto x85 = x45 + x67;
auto x86 = x53 + x61;
auto x87 = x67 + x80;
cell[0] = -V{0.296296296296296}*omega*(rho*x35 + V{1}) - x27*cell[0];
cell[1] = -x27*cell[1] - x36*(rho*(x34 + x37 - x38) + V{1});
cell[2] = -x27*cell[2] - x36*(rho*(x33 + x39 - x40 + x41) + V{1});
cell[3] = -x27*cell[3] - x36*(rho*(x31 + x41 + x42 - x43) + V{1});
cell[4] = -x27*cell[4] - x44*(rho*(-x46 + x48) + V{1});
cell[5] = -x27*cell[5] - x44*(rho*(x52 - V{4.5}*x50*x50) + V{1});
cell[6] = -x27*cell[6] - x44*(rho*(x42 + x47 - x54) + V{1});
cell[7] = -x27*cell[7] - x44*(rho*(x47 + x55 - V{4.5}*x56*x56) + V{1});
cell[8] = -x27*cell[8] - x44*(rho*(-x58 + x60) + V{1});
cell[9] = -x27*cell[9] - x44*(rho*(x55 + x59 - V{4.5}*x62*x62) + V{1});
cell[10] = -x27*cell[10] - x63*(rho*(x42 + x48 - x65) + V{1});
cell[11] = -x27*cell[11] - x63*(rho*(x48 + x55 - V{4.5}*x66*x66) + V{1});
cell[12] = -x27*cell[12] - x63*(rho*(x42 + x52 - V{4.5}*x68*x68) + V{1});
cell[13] = -x27*cell[13] + x63*(rho*(x69 + x76 + V{4.5}*(x70*x70)) + V{-1});
cell[14] = -x27*cell[14] + x36*(rho*(x38 + x77) + V{-1});
cell[15] = -x27*cell[15] + x36*(rho*(x40 + x72 + x75) + V{-1});
cell[16] = -x27*cell[16] + x36*(rho*(x43 + x71 + x78 + V{1}) + V{-1});
cell[17] = -x27*cell[17] + x44*(rho*(x46 + x79) + V{-1});
cell[18] = -x27*cell[18] - x44*(rho*(x59 + x69 - V{4.5}*x80*x80) + V{1});
cell[19] = -x27*cell[19] + x44*(rho*(x54 + x81) + V{-1});
cell[20] = -x27*cell[20] - x44*(rho*(x69 + x83 - V{4.5}*x82*x82) + V{1});
cell[21] = -x27*cell[21] + x44*(rho*(x58 + x76) + V{-1});
cell[22] = -x27*cell[22] - x44*(rho*(x51 + x83 - V{4.5}*x84*x84) + V{1});
cell[23] = -x27*cell[23] + x63*(rho*(x42 + x65 + x79) + V{-1});
cell[24] = -x27*cell[24] + x63*(rho*(x55 + x79 + V{4.5}*(x85*x85)) + V{-1});
cell[25] = -x27*cell[25] + x63*(rho*(x51 + x81 + V{4.5}*(x86*x86)) + V{-1});
cell[26] = -x27*cell[26] - x63*(rho*(x60 + x69 - V{4.5}*x87*x87) + V{1});
return x28 + x30 + x32;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto adeBgkCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega)
{
auto x27 = omega + V{-1};
auto x28 = V{3}*u[0];
auto x29 = x28 + V{-1};
auto x30 = V{0.0740740740740741}*omega;
auto x31 = V{3}*u[1];
auto x32 = x31 + V{-1};
auto x33 = V{3}*u[2];
auto x34 = x29 + x31;
auto x35 = V{0.0185185185185185}*omega;
auto x36 = -x28;
auto x37 = x31 + V{1};
auto x38 = x36 + x37;
auto x39 = x29 + x33;
auto x40 = x33 + V{1};
auto x41 = -x31;
auto x42 = V{0.00462962962962963}*omega;
auto x43 = -x33;
auto x44 = x28 + V{1};
auto x45 = x31 + x44;
auto x46 = x41 + x44;
cell[0] = V{0.296296296296296}*omega*(rho + V{-1}) - x27*cell[0];
cell[1] = -x27*cell[1] - x30*(rho*x29 + V{1});
cell[2] = -x27*cell[2] - x30*(rho*x32 + V{1});
cell[3] = -x27*cell[3] - x30*(rho*(x33 + V{-1}) + V{1});
cell[4] = -x27*cell[4] - x35*(rho*x34 + V{1});
cell[5] = -x27*cell[5] + x35*(rho*x38 + V{-1});
cell[6] = -x27*cell[6] - x35*(rho*x39 + V{1});
cell[7] = -x27*cell[7] + x35*(rho*(x36 + x40) + V{-1});
cell[8] = -x27*cell[8] - x35*(rho*(x32 + x33) + V{1});
cell[9] = -x27*cell[9] + x35*(rho*(x40 + x41) + V{-1});
cell[10] = -x27*cell[10] - x42*(rho*(x33 + x34) + V{1});
cell[11] = -x27*cell[11] - x42*(rho*(x34 + x43) + V{1});
cell[12] = -x27*cell[12] - x42*(rho*(x39 + x41) + V{1});
cell[13] = -x27*cell[13] + x42*(rho*(x33 + x38) + V{-1});
cell[14] = -x27*cell[14] + x30*(rho*x44 + V{-1});
cell[15] = -x27*cell[15] + x30*(rho*x37 + V{-1});
cell[16] = -x27*cell[16] + x30*(rho*x40 + V{-1});
cell[17] = -x27*cell[17] + x35*(rho*x45 + V{-1});
cell[18] = -x27*cell[18] + x35*(rho*x46 + V{-1});
cell[19] = -x27*cell[19] + x35*(rho*(x33 + x44) + V{-1});
cell[20] = -x27*cell[20] + x35*(rho*(x43 + x44) + V{-1});
cell[21] = -x27*cell[21] + x35*(rho*(x33 + x37) + V{-1});
cell[22] = -x27*cell[22] + x35*(rho*(x37 + x43) + V{-1});
cell[23] = -x27*cell[23] + x42*(rho*(x33 + x45) + V{-1});
cell[24] = -x27*cell[24] + x42*(rho*(x43 + x45) + V{-1});
cell[25] = -x27*cell[25] + x42*(rho*(x33 + x46) + V{-1});
cell[26] = -x27*cell[26] + x42*(rho*(x43 + x46) + V{-1});
return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

template <typename CELL, typename PRESSURE, typename J, typename OMEGA, typename V=typename CELL::value_t>
static auto incBgkCollision(CELL& cell, PRESSURE& pressure, J& j, OMEGA& omega) any_platform
{
auto x27 = j[0]*j[0];
auto x28 = j[1]*j[1];
auto x29 = j[2]*j[2];
auto x30 = omega + V{-1};
auto x31 = V{0.111111111111111}*x28;
auto x32 = V{0.222222222222222}*j[0];
auto x33 = V{0.222222222222222}*x27;
auto x34 = V{0.111111111111111}*x29;
auto x35 = V{0.222222222222222}*pressure;
auto x36 = -x35;
auto x37 = x34 + x36 + V{0.0740740740740741};
auto x38 = V{0.111111111111111}*x27;
auto x39 = V{0.222222222222222}*j[1];
auto x40 = V{0.222222222222222}*x28;
auto x41 = V{0.222222222222222}*j[2];
auto x42 = V{0.222222222222222}*x29;
auto x43 = j[0] + j[1];
auto x44 = V{0.0833333333333333}*(x43*x43);
auto x45 = V{0.0555555555555556}*j[0];
auto x46 = V{0.0555555555555556}*j[1];
auto x47 = x45 + x46;
auto x48 = V{0.0277777777777778}*x27;
auto x49 = V{0.0277777777777778}*x28;
auto x50 = V{0.0277777777777778}*x29;
auto x51 = V{0.0555555555555556}*pressure;
auto x52 = x48 + x49 + x50 - x51 + V{0.0185185185185185};
auto x53 = -x46;
auto x54 = -j[0];
auto x55 = x54 + j[1];
auto x56 = x45 + x52;
auto x57 = V{0.0555555555555556}*j[2];
auto x58 = j[0] + j[2];
auto x59 = V{0.0833333333333333}*(x58*x58);
auto x60 = -x57;
auto x61 = x54 + j[2];
auto x62 = j[1] + j[2];
auto x63 = V{0.0833333333333333}*(x62*x62);
auto x64 = x46 + x52;
auto x65 = -j[1];
auto x66 = x65 + j[2];
auto x67 = x43 + j[2];
auto x68 = V{0.0208333333333333}*(x67*x67);
auto x69 = V{0.0138888888888889}*j[0];
auto x70 = V{0.0138888888888889}*j[1];
auto x71 = V{0.0138888888888889}*j[2];
auto x72 = x69 + x70 + x71;
auto x73 = V{0.00694444444444444}*x27;
auto x74 = V{0.00694444444444444}*x28;
auto x75 = V{0.00694444444444444}*x29;
auto x76 = V{0.0138888888888889}*pressure;
auto x77 = x73 + x74 + x75 - x76 + V{0.00462962962962963};
auto x78 = x61 + x65;
auto x79 = -x71;
auto x80 = x69 + x77;
auto x81 = x79 + x80;
auto x82 = -j[2];
auto x83 = x55 + x82;
auto x84 = -x70;
auto x85 = x71 + x84;
auto x86 = x55 + j[2];
auto x87 = -x31;
auto x88 = -x34 + x35 + V{-0.0740740740740741};
auto x89 = -x38;
auto x90 = -x48 - x49 - x50 + x51 + V{-0.0185185185185185};
auto x91 = -x45;
auto x92 = x65 + j[0];
auto x93 = x57 + x90;
auto x94 = x82 + j[0];
auto x95 = x52 + x57;
auto x96 = x82 + j[1];
auto x97 = x43 + x82;
auto x98 = -x69 + x77;
auto x99 = x58 + x65;
auto x100 = x70 + x98;
auto x101 = x82 + x92;
cell[0] = -omega*(-V{0.888888888888889}*pressure + V{0.444444444444444}*x27 + V{0.444444444444444}*x28 + V{0.444444444444444}*x29 + V{0.296296296296296}) - x30*cell[0];
cell[1] = -omega*(x31 + x32 - x33 + x37) - x30*cell[1];
cell[2] = -omega*(x37 + x38 + x39 - x40) - x30*cell[2];
cell[3] = -omega*(x31 + x36 + x38 + x41 - x42 + V{0.0740740740740741}) - x30*cell[3];
cell[4] = -omega*(-x44 + x47 + x52) - x30*cell[4];
cell[5] = -omega*(x53 + x56 - V{0.0833333333333333}*x55*x55) - x30*cell[5];
cell[6] = -omega*(x56 + x57 - x59) - x30*cell[6];
cell[7] = -omega*(x56 + x60 - V{0.0833333333333333}*x61*x61) - x30*cell[7];
cell[8] = -omega*(x57 - x63 + x64) - x30*cell[8];
cell[9] = -omega*(x60 + x64 - V{0.0833333333333333}*x66*x66) - x30*cell[9];
cell[10] = -omega*(-x68 + x72 + x77) - x30*cell[10];
cell[11] = -omega*(x70 + x81 - V{0.0208333333333333}*x78*x78) - x30*cell[11];
cell[12] = -omega*(x80 + x85 - V{0.0208333333333333}*x83*x83) - x30*cell[12];
cell[13] = -omega*(x81 + x84 - V{0.0208333333333333}*x86*x86) - x30*cell[13];
cell[14] = omega*(x32 + x33 + x87 + x88) - x30*cell[14];
cell[15] = omega*(x39 + x40 + x88 + x89) - x30*cell[15];
cell[16] = omega*(x35 + x41 + x42 + x87 + x89 + V{-0.0740740740740741}) - x30*cell[16];
cell[17] = omega*(x44 + x47 + x90) - x30*cell[17];
cell[18] = -omega*(x64 + x91 - V{0.0833333333333333}*x92*x92) - x30*cell[18];
cell[19] = omega*(x45 + x59 + x93) - x30*cell[19];
cell[20] = -omega*(x91 + x95 - V{0.0833333333333333}*x94*x94) - x30*cell[20];
cell[21] = omega*(x46 + x63 + x93) - x30*cell[21];
cell[22] = -omega*(x53 + x95 - V{0.0833333333333333}*x96*x96) - x30*cell[22];
cell[23] = omega*(x68 + x72 - x73 - x74 - x75 + x76 + V{-0.00462962962962963}) - x30*cell[23];
cell[24] = -omega*(x85 + x98 - V{0.0208333333333333}*x97*x97) - x30*cell[24];
cell[25] = -omega*(x100 + x79 - V{0.0208333333333333}*x99*x99) - x30*cell[25];
cell[26] = -omega*(x100 + x71 - V{0.0208333333333333}*x101*x101) - x30*cell[26];
return x27 + x28 + x29;
}

template <typename CELL, typename RHO, typename U, typename RATIORHO, typename OMEGA, typename V=typename CELL::value_t>
static auto constRhoBgkCollision(CELL& cell, RHO& rho, U& u, RATIORHO& ratioRho, OMEGA& omega) any_platform
{
auto x27 = omega + V{-1};
auto x28 = V{0.296296296296296}*rho;
auto x29 = u[0]*u[0];
auto x30 = V{1.5}*x29;
auto x31 = u[1]*u[1];
auto x32 = V{1.5}*x31;
auto x33 = u[2]*u[2];
auto x34 = V{1.5}*x33;
auto x35 = x32 + x34 + V{-1};
auto x36 = x30 + x35;
auto x37 = V{0.0740740740740741}*rho;
auto x38 = V{3}*u[0];
auto x39 = V{3}*x29;
auto x40 = x35 + x38 - x39;
auto x41 = ratioRho*x37;
auto x42 = V{3}*u[1];
auto x43 = V{3}*x31;
auto x44 = x30 + V{-1};
auto x45 = x34 + x42 - x43 + x44;
auto x46 = V{3}*u[2];
auto x47 = V{3}*x33;
auto x48 = x32 + x44 + x46 - x47;
auto x49 = V{0.0185185185185185}*rho;
auto x50 = u[0] + u[1];
auto x51 = V{4.5}*(x50*x50);
auto x52 = x36 + x38;
auto x53 = x42 + x52;
auto x54 = -x51 + x53;
auto x55 = ratioRho*x49;
auto x56 = -u[0];
auto x57 = x56 + u[1];
auto x58 = -x42;
auto x59 = x52 + x58;
auto x60 = x59 - V{4.5}*x57*x57;
auto x61 = u[0] + u[2];
auto x62 = V{4.5}*(x61*x61);
auto x63 = x46 + x52 - x62;
auto x64 = -x46;
auto x65 = x56 + u[2];
auto x66 = x52 + x64 - V{4.5}*x65*x65;
auto x67 = u[1] + u[2];
auto x68 = V{4.5}*(x67*x67);
auto x69 = x36 + x42;
auto x70 = x46 + x69;
auto x71 = -x68 + x70;
auto x72 = -u[1];
auto x73 = x72 + u[2];
auto x74 = x64 + x69 - V{4.5}*x73*x73;
auto x75 = V{0.00462962962962963}*rho;
auto x76 = x50 + u[2];
auto x77 = V{4.5}*(x76*x76);
auto x78 = x46 + x53 - x77;
auto x79 = ratioRho*x75;
auto x80 = x65 + x72;
auto x81 = x53 + x64 - V{4.5}*x80*x80;
auto x82 = -u[2];
auto x83 = x57 + x82;
auto x84 = x46 + x59 - V{4.5}*x83*x83;
auto x85 = -x38;
auto x86 = x57 + u[2];
auto x87 = -x32;
auto x88 = V{1} - x34;
auto x89 = x87 + x88;
auto x90 = -x30;
auto x91 = x42 + x90;
auto x92 = x46 + x89 + x91;
auto x93 = x85 + x92 + V{4.5}*(x86*x86);
auto x94 = x38 + x89;
auto x95 = x39 + x94;
auto x96 = x43 + x88 + x91;
auto x97 = x46 + x90;
auto x98 = x47 + x87 + x97 + V{1};
auto x99 = x91 + x94;
auto x100 = x51 + x99;
auto x101 = x72 + u[0];
auto x102 = x69 + x85 - V{4.5}*x101*x101;
auto x103 = x94 + x97;
auto x104 = x103 + x62;
auto x105 = x82 + u[0];
auto x106 = x36 + x46;
auto x107 = x106 + x85 - V{4.5}*x105*x105;
auto x108 = x68 + x92;
auto x109 = x82 + u[1];
auto x110 = x106 + x58 - V{4.5}*x109*x109;
auto x111 = x46 + x77 + x99;
auto x112 = x50 + x82;
auto x113 = x64 + x99 + V{4.5}*(x112*x112);
auto x114 = x61 + x72;
auto x115 = x103 + x58 + V{4.5}*(x114*x114);
auto x116 = x101 + x82;
auto x117 = x70 + x85 - V{4.5}*x116*x116;
cell[0] = -ratioRho*x28*x36 - x27*(x28*x36 + cell[0] + V{0.296296296296296}) + V{-0.296296296296296};
cell[1] = -x27*(x37*x40 + cell[1] + V{0.0740740740740741}) - x40*x41 + V{-0.0740740740740741};
cell[2] = -x27*(x37*x45 + cell[2] + V{0.0740740740740741}) - x41*x45 + V{-0.0740740740740741};
cell[3] = -x27*(x37*x48 + cell[3] + V{0.0740740740740741}) - x41*x48 + V{-0.0740740740740741};
cell[4] = -x27*(x49*x54 + cell[4] + V{0.0185185185185185}) - x54*x55 + V{-0.0185185185185185};
cell[5] = -x27*(x49*x60 + cell[5] + V{0.0185185185185185}) - x55*x60 + V{-0.0185185185185185};
cell[6] = -x27*(x49*x63 + cell[6] + V{0.0185185185185185}) - x55*x63 + V{-0.0185185185185185};
cell[7] = -x27*(x49*x66 + cell[7] + V{0.0185185185185185}) - x55*x66 + V{-0.0185185185185185};
cell[8] = -x27*(x49*x71 + cell[8] + V{0.0185185185185185}) - x55*x71 + V{-0.0185185185185185};
cell[9] = -x27*(x49*x74 + cell[9] + V{0.0185185185185185}) - x55*x74 + V{-0.0185185185185185};
cell[10] = -x27*(x75*x78 + cell[10] + V{0.00462962962962963}) - x78*x79 + V{-0.00462962962962963};
cell[11] = -x27*(x75*x81 + cell[11] + V{0.00462962962962963}) - x79*x81 + V{-0.00462962962962963};
cell[12] = -x27*(x75*x84 + cell[12] + V{0.00462962962962963}) - x79*x84 + V{-0.00462962962962963};
cell[13] = -x27*(-x75*x93 + cell[13] + V{0.00462962962962963}) + x79*x93 + V{-0.00462962962962963};
cell[14] = -x27*(-x37*x95 + cell[14] + V{0.0740740740740741}) + x41*x95 + V{-0.0740740740740741};
cell[15] = -x27*(-x37*x96 + cell[15] + V{0.0740740740740741}) + x41*x96 + V{-0.0740740740740741};
cell[16] = -x27*(-x37*x98 + cell[16] + V{0.0740740740740741}) + x41*x98 + V{-0.0740740740740741};
cell[17] = x100*x55 - x27*(-x100*x49 + cell[17] + V{0.0185185185185185}) + V{-0.0185185185185185};
cell[18] = -x102*x55 - x27*(x102*x49 + cell[18] + V{0.0185185185185185}) + V{-0.0185185185185185};
cell[19] = x104*x55 - x27*(-x104*x49 + cell[19] + V{0.0185185185185185}) + V{-0.0185185185185185};
cell[20] = -x107*x55 - x27*(x107*x49 + cell[20] + V{0.0185185185185185}) + V{-0.0185185185185185};
cell[21] = x108*x55 - x27*(-x108*x49 + cell[21] + V{0.0185185185185185}) + V{-0.0185185185185185};
cell[22] = -x110*x55 - x27*(x110*x49 + cell[22] + V{0.0185185185185185}) + V{-0.0185185185185185};
cell[23] = x111*x79 - x27*(-x111*x75 + cell[23] + V{0.00462962962962963}) + V{-0.00462962962962963};
cell[24] = x113*x79 - x27*(-x113*x75 + cell[24] + V{0.00462962962962963}) + V{-0.00462962962962963};
cell[25] = x115*x79 - x27*(-x115*x75 + cell[25] + V{0.00462962962962963}) + V{-0.00462962962962963};
cell[26] = -x117*x79 - x27*(x117*x75 + cell[26] + V{0.00462962962962963}) + V{-0.00462962962962963};
return x29 + x31 + x33;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, OMEGA& omega) any_platform
{
auto x27 = omega + V{-1};
auto x28 = V{3}*u[2];
auto x29 = -x28;
auto x30 = V{3}*u[1];
auto x31 = -x30;
auto x32 = V{3}*u[0];
auto x33 = V{1} - x32;
auto x34 = x31 + x33;
auto x35 = x29 + x34;
auto x36 = V{0.00102880658436214}*rho;
auto x37 = x35*x36;
auto x38 = V{0.0164609053497942}*rho;
auto x39 = x28 + x33;
auto x40 = x31 + x39;
auto x41 = x36*x40;
auto x42 = V{0.00411522633744856}*rho;
auto x43 = V{0.222222222222222}*cell[23];
auto x44 = V{0.222222222222222}*cell[24];
auto x45 = -V{0.222222222222222}*cell[10];
auto x46 = V{0.222222222222222}*cell[11];
auto x47 = x32 + V{1};
auto x48 = x30 + x47;
auto x49 = x28 + x48;
auto x50 = -x36*x49;
auto x51 = x29 + x48;
auto x52 = x36*x51;
auto x53 = -x42*x48 + x43 + x44 + x45 - x46 + x50 - x52 + V{0.222222222222222}*cell[17] - V{0.222222222222222}*cell[4];
auto x54 = x34*x42 + x41 + x53;
auto x55 = x30 + x33;
auto x56 = x29 + x55;
auto x57 = x36*x56;
auto x58 = x29 + x33;
auto x59 = V{0.222222222222222}*cell[25];
auto x60 = V{0.222222222222222}*cell[12];
auto x61 = x28 + x47;
auto x62 = x31 + x61;
auto x63 = x36*x62;
auto x64 = -x42*x61 + x59 - x60 - x63 + V{0.222222222222222}*cell[19] - V{0.222222222222222}*cell[6];
auto x65 = x42*x58 + x57 + x64;
auto x66 = V{0.222222222222222}*cell[18];
auto x67 = V{0.222222222222222}*cell[20];
auto x68 = V{0.222222222222222}*cell[26];
auto x69 = V{0.222222222222222}*cell[5];
auto x70 = V{0.222222222222222}*cell[7];
auto x71 = V{0.222222222222222}*cell[13];
auto x72 = x28 + x55;
auto x73 = x36*x72;
auto x74 = x42*x55;
auto x75 = x39*x42;
auto x76 = x31 + x47;
auto x77 = x29 + x76;
auto x78 = x36*x77;
auto x79 = x42*x76;
auto x80 = x29 + x47;
auto x81 = x42*x80;
auto x82 = -x38*x47 + x66 + x67 + x68 - x69 - x70 - x71 + x73 + x74 + x75 - x78 - x79 - x81 + V{0.222222222222222}*cell[14] - V{0.222222222222222}*cell[1] + V{6.16790569236198e-18};
auto x83 = x32 + V{-1};
auto x84 = V{0.0740740740740741}*rho;
auto x85 = x31 + V{1};
auto x86 = x29 + x85;
auto x87 = x30 + V{1};
auto x88 = x28 + x87;
auto x89 = -x42*x88 - x68 + x71 - x73 + x78 + V{0.222222222222222}*cell[21] - V{0.222222222222222}*cell[8];
auto x90 = x37 + x42*x86 + x89;
auto x91 = V{0.222222222222222}*cell[22];
auto x92 = V{0.222222222222222}*cell[9];
auto x93 = x28 + x85;
auto x94 = x42*x93;
auto x95 = x29 + x87;
auto x96 = x42*x95;
auto x97 = -x38*x87 - x59 + x60 + x63 - x66 + x69 - x74 + x79 + x91 - x92 + x94 - x96 + V{0.222222222222222}*cell[15] - V{0.222222222222222}*cell[2];
auto x98 = x30 + V{-1};
auto x99 = x29 + V{1};
auto x100 = x28 + V{1};
auto x101 = -x100*x38 + x43 - x44 + x45 + x46 + x50 + x52 - x67 + x70 - x75 + x81 - x91 + x92 - x94 + x96 + V{0.222222222222222}*cell[16] - V{0.222222222222222}*cell[3];
auto x102 = x28 + V{-1};
auto x103 = V{0.00051440329218107}*rho;
auto x104 = V{0.00205761316872428}*rho;
auto x105 = x103*x35;
auto x106 = V{0.0555555555555556}*cell[15];
auto x107 = V{0.111111111111111}*cell[23];
auto x108 = V{0.0555555555555556}*cell[2];
auto x109 = -x108;
auto x110 = -V{0.111111111111111}*cell[10];
auto x111 = -x103*x49;
auto x112 = x42*x87;
auto x113 = -x112;
auto x114 = V{0.0555555555555556}*cell[19];
auto x115 = V{0.0555555555555556}*cell[6];
auto x116 = x36*x61;
auto x117 = x114 - x115 - x116;
auto x118 = x106 + x107 + x109 + x110 + x111 + x113 + x117;
auto x119 = x105 + x118 + x36*x58 + x42*x85;
auto x120 = V{0.111111111111111}*cell[24];
auto x121 = -V{0.111111111111111}*cell[11];
auto x122 = -x103*x51;
auto x123 = V{0.0555555555555556}*cell[20];
auto x124 = V{0.0555555555555556}*cell[7];
auto x125 = x36*x39;
auto x126 = x36*x80;
auto x127 = x123 - x124 + x125 - x126;
auto x128 = x120 + x121 + x122 + x127;
auto x129 = V{0.0555555555555556}*cell[21];
auto x130 = V{0.0555555555555556}*cell[8];
auto x131 = x36*x88;
auto x132 = -x42*x47 + V{0.0555555555555556}*cell[14] - V{0.0555555555555556}*cell[1] + V{1.5419764230905e-18};
auto x133 = x129 - x130 - x131 + x132;
auto x134 = x133 + x33*x42 + x36*x86;
auto x135 = V{0.0555555555555556}*cell[22];
auto x136 = V{0.0555555555555556}*cell[9];
auto x137 = x36*x93;
auto x138 = x36*x95;
auto x139 = x135 - x136 + x137 - x138;
auto x140 = -x104*x48 + x139 + V{0.111111111111111}*cell[17] - V{0.111111111111111}*cell[4];
auto x141 = x30 + x83;
auto x142 = V{0.0185185185185185}*rho;
auto x143 = x28 + x83;
auto x144 = x143*x36;
auto x145 = -x144;
auto x146 = x42*x98;
auto x147 = x143 + x31;
auto x148 = x103*x147;
auto x149 = -x42*x83;
auto x150 = V{0.111111111111111}*cell[25];
auto x151 = V{0.111111111111111}*cell[12];
auto x152 = x103*x62;
auto x153 = -x135 + x136 - x137 + x138 + x150 - x151 - x152;
auto x154 = -x148 + x149 + x153;
auto x155 = V{0.111111111111111}*cell[26];
auto x156 = V{0.111111111111111}*cell[13];
auto x157 = x103*x72;
auto x158 = x28 + x98;
auto x159 = x158*x36;
auto x160 = x103*x77;
auto x161 = -x129 + x130 + x131 + x132 + x155 - x156 + x157 + x159 - x160;
auto x162 = x27*(x104*x55 - x104*x76 - x106 + x108 + x112 + x117 + x127 + x145 + x146 + x154 + x161 + V{0.111111111111111}*cell[18] - V{0.111111111111111}*cell[5]);
auto x163 = V{0.0555555555555556}*cell[16];
auto x164 = V{0.0555555555555556}*cell[3];
auto x165 = x100*x42;
auto x166 = -x36*x48 + V{0.0555555555555556}*cell[17] - V{0.0555555555555556}*cell[4];
auto x167 = x163 - x164 - x165 + x166;
auto x168 = x167 + x34*x36 + x42*x99;
auto x169 = V{0.0555555555555556}*cell[18];
auto x170 = V{0.0555555555555556}*cell[5];
auto x171 = x36*x55;
auto x172 = x36*x76;
auto x173 = x169 - x170 + x171 - x172;
auto x174 = -x104*x61 + x107 + x110 + x111 + x173 + V{0.111111111111111}*cell[19] - V{0.111111111111111}*cell[6];
auto x175 = x141 + x29;
auto x176 = -x103*x175;
auto x177 = -x141*x36;
auto x178 = x102*x42;
auto x179 = -x163 + x164 + x165 + x166 + x178;
auto x180 = x27*(x104*x39 - x104*x80 + x120 + x121 + x122 + x139 + x149 + x161 + x173 + x176 + x177 + x179 + V{0.111111111111111}*cell[20] - V{0.111111111111111}*cell[7]);
auto x181 = -x169 + x170 - x171 + x172;
auto x182 = -x104*x88 - x123 + x124 - x125 + x126 - x155 + x156 - x157 + x160 + x181 + V{0.111111111111111}*cell[21] - V{0.111111111111111}*cell[8];
auto x183 = -x146;
auto x184 = x177 + x183;
auto x185 = x128 + x176;
auto x186 = x27*(x104*x93 - x104*x95 + x106 + x109 + x113 - x114 + x115 + x116 + x144 + x148 - x150 + x151 + x152 + x179 + x181 + x184 + x185 + V{0.111111111111111}*cell[22] - V{0.111111111111111}*cell[9]);
auto x187 = V{0.000192901234567901}*rho;
auto x188 = V{6.43004115226337e-05}*rho;
auto x189 = x188*x56;
auto x190 = x188*x40;
auto x191 = x33*x36;
auto x192 = x36*x85;
auto x193 = V{0.0138888888888889}*cell[15];
auto x194 = V{0.0138888888888889}*cell[2];
auto x195 = x36*x87;
auto x196 = -x103*x48 + x193 - x194 - x195 + V{0.0277777777777778}*cell[17] - V{0.0277777777777778}*cell[4];
auto x197 = x103*x34 + x191 + x192 + x196;
auto x198 = x36*x99;
auto x199 = V{0.0138888888888889}*cell[16];
auto x200 = V{0.0138888888888889}*cell[3];
auto x201 = x100*x36;
auto x202 = -x103*x61 + x199 - x200 - x201 + V{0.0277777777777778}*cell[19] - V{0.0277777777777778}*cell[6];
auto x203 = x103*x58 + x198 + x202;
auto x204 = V{0.0138888888888889}*cell[13];
auto x205 = V{0.0277777777777778}*cell[21];
auto x206 = V{0.0138888888888889}*cell[26];
auto x207 = V{0.0277777777777778}*cell[8];
auto x208 = x188*x77;
auto x209 = x188*x72;
auto x210 = x103*x88;
auto x211 = V{0.0138888888888889}*cell[14];
auto x212 = V{0.0138888888888889}*cell[24];
auto x213 = V{0.0138888888888889}*cell[25];
auto x214 = -V{0.0138888888888889}*cell[1];
auto x215 = V{0.0138888888888889}*cell[11];
auto x216 = V{0.0138888888888889}*cell[12];
auto x217 = x188*x51;
auto x218 = x188*x62;
auto x219 = -x36*x47;
auto x220 = x211 + x212 + x213 + x214 - x215 - x216 - x217 - x218 + x219 + V{3.85494105772624e-19};
auto x221 = -x187*x49 + x204 + x205 - x206 - x207 + x208 - x209 - x210 + x220 - V{0.0416666666666667}*cell[10] + V{0.0416666666666667}*cell[23];
auto x222 = x141 + x28;
auto x223 = V{0.00462962962962963}*rho;
auto x224 = V{0.0138888888888889}*cell[23];
auto x225 = V{0.0138888888888889}*cell[10];
auto x226 = x188*x49;
auto x227 = -x204 + x206 - x208 + x209 + x211 + x214 + x219 + x224 - x225 - x226 + V{3.85494105772624e-19};
auto x228 = x188*x35 + x227;
auto x229 = V{0.0277777777777778}*cell[22];
auto x230 = V{0.0277777777777778}*cell[9];
auto x231 = x103*x93;
auto x232 = x103*x95;
auto x233 = x103*x39 - x103*x80 - x199 + x200 + x201 + V{0.0277777777777778}*cell[20] - V{0.0277777777777778}*cell[7];
auto x234 = -x187*x51 - x213 + x216 + x218 + x229 - x230 + x231 - x232 + x233 - V{0.0416666666666667}*cell[11] + V{0.0416666666666667}*cell[24];
auto x235 = x103*x55 - x103*x76 - x193 + x194 + x195 + V{0.0277777777777778}*cell[18] - V{0.0277777777777778}*cell[5];
auto x236 = -x187*x62 - x212 + x215 + x217 - x229 + x230 - x231 + x232 + x235 - V{0.0416666666666667}*cell[12] + V{0.0416666666666667}*cell[25];
auto x237 = x188*x222;
auto x238 = x103*x158;
auto x239 = x36*x98;
auto x240 = x102*x36;
auto x241 = -x36*x83;
auto x242 = x240 + x241;
auto x243 = x175*x188;
auto x244 = x147*x188;
auto x245 = -x243 - x244;
auto x246 = x27*(x187*x72 - x187*x77 - x205 + x207 + x210 + x220 - x224 + x225 + x226 + x233 + x235 + x237 + x238 + x239 + x242 + x245 - V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[26]);
auto x247 = -x222*x36;
auto x248 = x175*x36;
auto x249 = -x141*x42 - x248 + x53;
auto x250 = x147*x36;
auto x251 = -x143*x42 - x250 + x64;
auto x252 = -x158*x42 + x247 + x89;
auto x253 = -x103*x222;
auto x254 = x118 + x145 + x253;
auto x255 = x133 - x159;
auto x256 = x167 - x178;
auto x257 = -x103*x141 + x196 - x239;
auto x258 = -x103*x143 + x202 - x240 + x241;
auto x259 = x227 - x237;
cell[0] = V{0.296296296296296}*rho + V{-0.296296296296296};
cell[1] = x27*(x33*x38 + x37 + x54 + x65 + x82) - x83*x84 + V{-0.0740740740740741};
cell[2] = x27*(x38*x85 + x54 - x57 + x90 + x97) - x84*x98 + V{-0.0740740740740741};
cell[3] = -x102*x84 + x27*(x101 + x38*x99 - x41 + x65 + x90) + V{-0.0740740740740741};
cell[4] = -x141*x142 + x27*(x103*x40 + x104*x34 + x119 + x128 + x134 + x140) + V{-0.0185185185185185};
cell[5] = x142*x55 + x162 + V{-0.0185185185185185};
cell[6] = -x142*x143 + x27*(x103*x56 + x104*x58 + x105 + x134 + x153 + x168 + x174) + V{-0.0185185185185185};
cell[7] = x142*x39 + x180 + V{-0.0185185185185185};
cell[8] = -x142*x158 + x27*(x104*x86 + x119 + x168 + x182) + V{-0.0185185185185185};
cell[9] = x142*x93 + x186 + V{-0.0185185185185185};
cell[10] = -x222*x223 + x27*(x103*x86 + x187*x35 + x189 + x190 + x197 + x203 + x221) + V{-0.00462962962962963};
cell[11] = -x175*x223 + x27*(x187*x40 - x189 + x197 - x198 + x228 + x234) + V{-0.00462962962962963};
cell[12] = -x147*x223 + x27*(x187*x56 - x190 + x191 - x192 + x203 + x228 + x236) + V{-0.00462962962962963};
cell[13] = x223*x72 + x246 + V{-0.00462962962962963};
cell[14] = -x27*(x247 + x249 + x251 - x38*x83 + x82) + x47*x84 + V{-0.0740740740740741};
cell[15] = -x27*(x249 + x250 + x252 - x38*x98 + x97) + x84*x87 + V{-0.0740740740740741};
cell[16] = x100*x84 - x27*(x101 - x102*x38 + x248 + x251 + x252) + V{-0.0740740740740741};
cell[17] = x142*x48 - x27*(-x104*x141 + x140 + x149 + x183 + x185 + x254 + x255) + V{-0.0185185185185185};
cell[18] = x142*x76 - x162 + V{-0.0185185185185185};
cell[19] = x142*x61 - x27*(-x104*x143 + x154 + x174 + x177 + x253 + x255 + x256) + V{-0.0185185185185185};
cell[20] = x142*x80 - x180 + V{-0.0185185185185185};
cell[21] = x142*x88 - x27*(-x104*x158 + x182 + x184 + x254 + x256) + V{-0.0185185185185185};
cell[22] = x142*x95 - x186 + V{-0.0185185185185185};
cell[23] = x223*x49 - x27*(-x187*x222 + x221 - x238 + x245 + x257 + x258) + V{-0.00462962962962963};
cell[24] = x223*x51 - x27*(-x175*x187 + x234 + x242 + x244 + x257 + x259) + V{-0.00462962962962963};
cell[25] = x223*x62 - x27*(-x147*x187 + x236 + x239 + x243 + x258 + x259) + V{-0.00462962962962963};
cell[26] = x223*x77 - x246 + V{-0.00462962962962963};
return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
}

template <typename CELL, typename RHO, typename U, typename PI, typename OMEGA, typename V=typename CELL::value_t>
static auto rlbCollision(CELL& cell, RHO& rho, U& u, PI& pi, OMEGA& omega) any_platform
{
auto x27 = omega + V{-1};
auto x28 = u[0]*u[0];
auto x29 = V{1.5}*x28;
auto x30 = u[1]*u[1];
auto x31 = V{1.5}*x30;
auto x32 = u[2]*u[2];
auto x33 = V{1.5}*x32;
auto x34 = x31 + x33 + V{-1};
auto x35 = x29 + x34;
auto x36 = V{0.0740740740740741}*rho;
auto x37 = V{3}*u[0];
auto x38 = V{3}*x28;
auto x39 = V{0.111111111111111}*pi[3];
auto x40 = V{0.111111111111111}*pi[5];
auto x41 = x27*(x39 + x40 - V{0.222222222222222}*pi[0]) + V{-0.0740740740740741};
auto x42 = V{3}*u[1];
auto x43 = V{3}*x30;
auto x44 = x29 + V{-1};
auto x45 = V{0.111111111111111}*pi[0];
auto x46 = x27*(x40 + x45 - V{0.222222222222222}*pi[3]) + V{-0.0740740740740741};
auto x47 = V{3}*u[2];
auto x48 = V{3}*x32;
auto x49 = x27*(x39 + x45 - V{0.222222222222222}*pi[5]) + V{-0.0740740740740741};
auto x50 = V{0.166666666666667}*pi[1];
auto x51 = V{0.0555555555555556}*pi[0];
auto x52 = V{0.0555555555555556}*pi[3];
auto x53 = x51 + x52 - V{0.0277777777777778}*pi[5];
auto x54 = x27*(x50 + x53);
auto x55 = V{0.0185185185185185}*rho;
auto x56 = u[0] + u[1];
auto x57 = V{4.5}*(x56*x56);
auto x58 = x35 + x37;
auto x59 = x42 + x58;
auto x60 = -u[0];
auto x61 = x60 + u[1];
auto x62 = -x42;
auto x63 = x58 + x62;
auto x64 = x27*(-x50 + x53) + V{0.0185185185185185};
auto x65 = V{0.166666666666667}*pi[2];
auto x66 = V{0.0555555555555556}*pi[5];
auto x67 = x51 + x66 - V{0.0277777777777778}*pi[3];
auto x68 = x27*(x65 + x67);
auto x69 = u[0] + u[2];
auto x70 = V{4.5}*(x69*x69);
auto x71 = -x47;
auto x72 = x60 + u[2];
auto x73 = x27*(-x65 + x67) + V{0.0185185185185185};
auto x74 = V{0.166666666666667}*pi[4];
auto x75 = V{0.0277777777777778}*pi[0];
auto x76 = x27*(x52 + x66 + x74 - x75);
auto x77 = u[1] + u[2];
auto x78 = V{4.5}*(x77*x77);
auto x79 = x35 + x42;
auto x80 = x47 + x79;
auto x81 = -u[1];
auto x82 = x81 + u[2];
auto x83 = x27*(-x52 - x66 + x74 + x75) + V{-0.0185185185185185};
auto x84 = V{0.0416666666666667}*pi[2];
auto x85 = V{0.0416666666666667}*pi[4];
auto x86 = V{0.0416666666666667}*pi[1];
auto x87 = V{0.0138888888888889}*pi[0];
auto x88 = V{0.0138888888888889}*pi[3];
auto x89 = V{0.0138888888888889}*pi[5];
auto x90 = x86 + x87 + x88 + x89;
auto x91 = x27*(x84 + x85 + x90);
auto x92 = V{0.00462962962962963}*rho;
auto x93 = x56 + u[2];
auto x94 = V{4.5}*(x93*x93);
auto x95 = -x84;
auto x96 = -x85;
auto x97 = x27*(x90 + x95 + x96);
auto x98 = x72 + x81;
auto x99 = -x86 + x87 + x88 + x89;
auto x100 = x27*(x84 + x96 + x99);
auto x101 = -u[2];
auto x102 = x101 + x61;
auto x103 = x27*(x85 + x95 + x99);
auto x104 = -x37;
auto x105 = x61 + u[2];
auto x106 = -x31;
auto x107 = V{1} - x33;
auto x108 = x106 + x107;
auto x109 = -x29;
auto x110 = x109 + x42;
auto x111 = x108 + x110 + x47;
auto x112 = x108 + x37;
auto x113 = x109 + x47;
auto x114 = x110 + x112;
auto x115 = x81 + u[0];
auto x116 = x112 + x113;
auto x117 = x101 + u[0];
auto x118 = x35 + x47;
auto x119 = x101 + u[1];
auto x120 = x101 + x56;
auto x121 = x69 + x81;
auto x122 = x101 + x115;
cell[0] = -V{0.296296296296296}*rho*x35 + V{0.444444444444444}*x27*(pi[0] + pi[3] + pi[5]) + V{-0.296296296296296};
cell[1] = -x36*(x34 + x37 - x38) + x41;
cell[2] = -x36*(x33 + x42 - x43 + x44) + x46;
cell[3] = -x36*(x31 + x44 + x47 - x48) + x49;
cell[4] = -x54 - x55*(-x57 + x59) + V{-0.0185185185185185};
cell[5] = -x55*(x63 - V{4.5}*x61*x61) - x64;
cell[6] = -x55*(x47 + x58 - x70) - x68 + V{-0.0185185185185185};
cell[7] = -x55*(x58 + x71 - V{4.5}*x72*x72) - x73;
cell[8] = -x55*(-x78 + x80) - x76 + V{-0.0185185185185185};
cell[9] = -x55*(x71 + x79 - V{4.5}*x82*x82) + x83;
cell[10] = -x91 - x92*(x47 + x59 - x94) + V{-0.00462962962962963};
cell[11] = -x92*(x59 + x71 - V{4.5}*x98*x98) - x97 + V{-0.00462962962962963};
cell[12] = -x100 - x92*(x47 + x63 - V{4.5}*x102*x102) + V{-0.00462962962962963};
cell[13] = -x103 + x92*(x104 + x111 + V{4.5}*(x105*x105)) + V{-0.00462962962962963};
cell[14] = x36*(x112 + x38) + x41;
cell[15] = x36*(x107 + x110 + x43) + x46;
cell[16] = x36*(x106 + x113 + x48 + V{1}) + x49;
cell[17] = -x54 + x55*(x114 + x57) + V{-0.0185185185185185};
cell[18] = -x55*(x104 + x79 - V{4.5}*x115*x115) - x64;
cell[19] = x55*(x116 + x70) - x68 + V{-0.0185185185185185};
cell[20] = -x55*(x104 + x118 - V{4.5}*x117*x117) - x73;
cell[21] = x55*(x111 + x78) - x76 + V{-0.0185185185185185};
cell[22] = -x55*(x118 + x62 - V{4.5}*x119*x119) + x83;
cell[23] = -x91 + x92*(x114 + x47 + x94) + V{-0.00462962962962963};
cell[24] = x92*(x114 + x71 + V{4.5}*(x120*x120)) - x97 + V{-0.00462962962962963};
cell[25] = -x100 + x92*(x116 + x62 + V{4.5}*(x121*x121)) + V{-0.00462962962962963};
cell[26] = -x103 - x92*(x104 + x80 - V{4.5}*x122*x122) + V{-0.00462962962962963};
return x28 + x30 + x32;
}

template <typename CELL, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineEqFirstOrder(CELL& cell, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x27 = V{3}*newU[0];
auto x28 = x27 + V{-1};
auto x29 = V{3}*newU[1];
auto x30 = x29 + V{-1};
auto x31 = V{3}*newU[2];
auto x32 = x28 + x29;
auto x33 = -x27;
auto x34 = x29 + V{1};
auto x35 = x33 + x34;
auto x36 = x28 + x31;
auto x37 = x31 + V{1};
auto x38 = -x29;
auto x39 = -x31;
auto x40 = x27 + V{1};
auto x41 = x29 + x40;
auto x42 = x38 + x40;
cell[0] = V{0.296296296296296}*newRho + V{-0.296296296296296};
cell[1] = -V{0.0740740740740741}*newRho*x28 + V{-0.0740740740740741};
cell[2] = -V{0.0740740740740741}*newRho*x30 + V{-0.0740740740740741};
cell[3] = -V{0.0740740740740741}*newRho*(x31 + V{-1}) + V{-0.0740740740740741};
cell[4] = -V{0.0185185185185185}*newRho*x32 + V{-0.0185185185185185};
cell[5] = V{0.0185185185185185}*newRho*x35 + V{-0.0185185185185185};
cell[6] = -V{0.0185185185185185}*newRho*x36 + V{-0.0185185185185185};
cell[7] = V{0.0185185185185185}*newRho*(x33 + x37) + V{-0.0185185185185185};
cell[8] = -V{0.0185185185185185}*newRho*(x30 + x31) + V{-0.0185185185185185};
cell[9] = V{0.0185185185185185}*newRho*(x37 + x38) + V{-0.0185185185185185};
cell[10] = -V{0.00462962962962963}*newRho*(x31 + x32) + V{-0.00462962962962963};
cell[11] = -V{0.00462962962962963}*newRho*(x32 + x39) + V{-0.00462962962962963};
cell[12] = -V{0.00462962962962963}*newRho*(x36 + x38) + V{-0.00462962962962963};
cell[13] = V{0.00462962962962963}*newRho*(x31 + x35) + V{-0.00462962962962963};
cell[14] = V{0.0740740740740741}*newRho*x40 + V{-0.0740740740740741};
cell[15] = V{0.0740740740740741}*newRho*x34 + V{-0.0740740740740741};
cell[16] = V{0.0740740740740741}*newRho*x37 + V{-0.0740740740740741};
cell[17] = V{0.0185185185185185}*newRho*x41 + V{-0.0185185185185185};
cell[18] = V{0.0185185185185185}*newRho*x42 + V{-0.0185185185185185};
cell[19] = V{0.0185185185185185}*newRho*(x31 + x40) + V{-0.0185185185185185};
cell[20] = V{0.0185185185185185}*newRho*(x39 + x40) + V{-0.0185185185185185};
cell[21] = V{0.0185185185185185}*newRho*(x31 + x34) + V{-0.0185185185185185};
cell[22] = V{0.0185185185185185}*newRho*(x34 + x39) + V{-0.0185185185185185};
cell[23] = V{0.00462962962962963}*newRho*(x31 + x41) + V{-0.00462962962962963};
cell[24] = V{0.00462962962962963}*newRho*(x39 + x41) + V{-0.00462962962962963};
cell[25] = V{0.00462962962962963}*newRho*(x31 + x42) + V{-0.00462962962962963};
cell[26] = V{0.00462962962962963}*newRho*(x39 + x42) + V{-0.00462962962962963};

}

template <typename CELL, typename OLDRHO, typename OLDU, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
static void defineNEq(CELL& cell, OLDRHO& oldRho, OLDU& oldU, NEWRHO& newRho, NEWU& newU) any_platform
{
auto x27 = oldU[0]*oldU[0];
auto x28 = V{1.5}*x27;
auto x29 = oldU[1]*oldU[1];
auto x30 = V{1.5}*x29;
auto x31 = oldU[2]*oldU[2];
auto x32 = V{1.5}*x31;
auto x33 = x30 + x32 + V{-1};
auto x34 = x28 + x33;
auto x35 = newU[0]*newU[0];
auto x36 = V{1.5}*x35;
auto x37 = newU[1]*newU[1];
auto x38 = V{1.5}*x37;
auto x39 = newU[2]*newU[2];
auto x40 = V{1.5}*x39;
auto x41 = x38 + x40 + V{-1};
auto x42 = x36 + x41;
auto x43 = V{0.0740740740740741}*oldRho;
auto x44 = V{3}*oldU[0];
auto x45 = V{3}*x27;
auto x46 = V{0.0740740740740741}*newRho;
auto x47 = V{3}*newU[0];
auto x48 = V{3}*x35;
auto x49 = V{3}*oldU[1];
auto x50 = V{3}*x29;
auto x51 = x28 + V{-1};
auto x52 = V{3}*newU[1];
auto x53 = V{3}*x37;
auto x54 = x36 + V{-1};
auto x55 = V{3}*oldU[2];
auto x56 = V{3}*x31;
auto x57 = V{3}*newU[2];
auto x58 = V{3}*x39;
auto x59 = V{0.0185185185185185}*oldRho;
auto x60 = oldU[0] + oldU[1];
auto x61 = V{4.5}*(x60*x60);
auto x62 = x34 + x44;
auto x63 = x49 + x62;
auto x64 = V{0.0185185185185185}*newRho;
auto x65 = newU[0] + newU[1];
auto x66 = V{4.5}*(x65*x65);
auto x67 = x42 + x47;
auto x68 = x52 + x67;
auto x69 = oldU[0] - oldU[1];
auto x70 = -V{4.5}*x69*x69;
auto x71 = -x49;
auto x72 = x62 + x71;
auto x73 = newU[0] - newU[1];
auto x74 = -V{4.5}*x73*x73;
auto x75 = -x52;
auto x76 = x67 + x75;
auto x77 = oldU[0] + oldU[2];
auto x78 = V{4.5}*(x77*x77);
auto x79 = newU[0] + newU[2];
auto x80 = V{4.5}*(x79*x79);
auto x81 = -x55;
auto x82 = -oldU[2];
auto x83 = x82 + oldU[0];
auto x84 = -V{4.5}*x83*x83;
auto x85 = -x57;
auto x86 = -newU[2];
auto x87 = x86 + newU[0];
auto x88 = -V{4.5}*x87*x87;
auto x89 = oldU[1] + oldU[2];
auto x90 = V{4.5}*(x89*x89);
auto x91 = x34 + x49;
auto x92 = x55 + x91;
auto x93 = newU[1] + newU[2];
auto x94 = V{4.5}*(x93*x93);
auto x95 = x42 + x52;
auto x96 = x57 + x95;
auto x97 = x82 + oldU[1];
auto x98 = -V{4.5}*x97*x97;
auto x99 = x86 + newU[1];
auto x100 = -V{4.5}*x99*x99;
auto x101 = V{0.00462962962962963}*oldRho;
auto x102 = x60 + oldU[2];
auto x103 = V{4.5}*(x102*x102);
auto x104 = V{0.00462962962962963}*newRho;
auto x105 = x65 + newU[2];
auto x106 = V{4.5}*(x105*x105);
auto x107 = x60 + x82;
auto x108 = V{4.5}*(x107*x107);
auto x109 = x65 + x86;
auto x110 = V{4.5}*(x109*x109);
auto x111 = x69 + oldU[2];
auto x112 = V{4.5}*(x111*x111);
auto x113 = x73 + newU[2];
auto x114 = V{4.5}*(x113*x113);
auto x115 = -x47;
auto x116 = x93 - newU[0];
auto x117 = V{4.5}*(x116*x116);
auto x118 = -x38;
auto x119 = V{1} - x40;
auto x120 = x118 + x119;
auto x121 = -x36;
auto x122 = x121 + x52;
auto x123 = x120 + x122 + x57;
auto x124 = -x44;
auto x125 = x89 - oldU[0];
auto x126 = V{4.5}*(x125*x125);
auto x127 = -x30;
auto x128 = V{1} - x32;
auto x129 = x127 + x128;
auto x130 = -x28;
auto x131 = x130 + x49;
auto x132 = x129 + x131 + x55;
auto x133 = x120 + x47;
auto x134 = x129 + x44;
auto x135 = x121 + x57;
auto x136 = x130 + x55;
auto x137 = x122 + x133;
auto x138 = x131 + x134;
auto x139 = x133 + x135;
auto x140 = x134 + x136;
auto x141 = x34 + x55;
auto x142 = x42 + x57;
cell[0] = -V{0.296296296296296}*newRho*x42 + V{0.296296296296296}*oldRho*x34 + cell[0];
cell[1] = x43*(x33 + x44 - x45) - x46*(x41 + x47 - x48) + cell[1];
cell[2] = x43*(x32 + x49 - x50 + x51) - x46*(x40 + x52 - x53 + x54) + cell[2];
cell[3] = x43*(x30 + x51 + x55 - x56) - x46*(x38 + x54 + x57 - x58) + cell[3];
cell[4] = x59*(-x61 + x63) - x64*(-x66 + x68) + cell[4];
cell[5] = x59*(x70 + x72) - x64*(x74 + x76) + cell[5];
cell[6] = x59*(x55 + x62 - x78) - x64*(x57 + x67 - x80) + cell[6];
cell[7] = x59*(x62 + x81 + x84) - x64*(x67 + x85 + x88) + cell[7];
cell[8] = x59*(-x90 + x92) - x64*(-x94 + x96) + cell[8];
cell[9] = x59*(x81 + x91 + x98) - x64*(x100 + x85 + x95) + cell[9];
cell[10] = x101*(-x103 + x55 + x63) - x104*(-x106 + x57 + x68) + cell[10];
cell[11] = x101*(-x108 + x63 + x81) - x104*(-x110 + x68 + x85) + cell[11];
cell[12] = x101*(-x112 + x55 + x72) - x104*(-x114 + x57 + x76) + cell[12];
cell[13] = -x101*(x124 + x126 + x132) + x104*(x115 + x117 + x123) + cell[13];
cell[14] = -x43*(x134 + x45) + x46*(x133 + x48) + cell[14];
cell[15] = -x43*(x128 + x131 + x50) + x46*(x119 + x122 + x53) + cell[15];
cell[16] = -x43*(x127 + x136 + x56 + V{1}) + x46*(x118 + x135 + x58 + V{1}) + cell[16];
cell[17] = -x59*(x138 + x61) + x64*(x137 + x66) + cell[17];
cell[18] = x59*(x124 + x70 + x91) - x64*(x115 + x74 + x95) + cell[18];
cell[19] = -x59*(x140 + x78) + x64*(x139 + x80) + cell[19];
cell[20] = x59*(x124 + x141 + x84) - x64*(x115 + x142 + x88) + cell[20];
cell[21] = -x59*(x132 + x90) + x64*(x123 + x94) + cell[21];
cell[22] = x59*(x141 + x71 + x98) - x64*(x100 + x142 + x75) + cell[22];
cell[23] = -x101*(x103 + x138 + x55) + x104*(x106 + x137 + x57) + cell[23];
cell[24] = -x101*(x108 + x138 + x81) + x104*(x110 + x137 + x85) + cell[24];
cell[25] = -x101*(x112 + x140 + x71) + x104*(x114 + x139 + x75) + cell[25];
cell[26] = x101*(x124 - x126 + x92) - x104*(x115 - x117 + x96) + cell[26];

}

template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
static void defineNEqFromPi(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
{
auto x27 = u[0]*u[0];
auto x28 = V{1.5}*x27;
auto x29 = u[1]*u[1];
auto x30 = V{1.5}*x29;
auto x31 = u[2]*u[2];
auto x32 = V{1.5}*x31;
auto x33 = x30 + x32 + V{-1};
auto x34 = x28 + x33;
auto x35 = V{0.0740740740740741}*rho;
auto x36 = V{3}*u[0];
auto x37 = V{3}*x27;
auto x38 = -V{0.111111111111111}*pi[3];
auto x39 = -V{0.111111111111111}*pi[5] + V{-0.0740740740740741};
auto x40 = x38 + x39 + V{0.222222222222222}*pi[0];
auto x41 = V{3}*u[1];
auto x42 = V{3}*x29;
auto x43 = x28 + V{-1};
auto x44 = -V{0.111111111111111}*pi[0];
auto x45 = x39 + x44 + V{0.222222222222222}*pi[3];
auto x46 = V{3}*u[2];
auto x47 = V{3}*x31;
auto x48 = x38 + x44 + V{0.222222222222222}*pi[5] + V{-0.0740740740740741};
auto x49 = V{0.0185185185185185}*rho;
auto x50 = u[0] + u[1];
auto x51 = V{4.5}*(x50*x50);
auto x52 = x34 + x36;
auto x53 = x41 + x52;
auto x54 = V{0.166666666666667}*pi[1];
auto x55 = V{0.0555555555555556}*pi[3];
auto x56 = V{0.0555555555555556}*pi[0] + V{-0.0185185185185185};
auto x57 = x55 + x56 - V{0.0277777777777778}*pi[5];
auto x58 = x54 + x57;
auto x59 = -u[0];
auto x60 = x59 + u[1];
auto x61 = -x41;
auto x62 = x52 + x61;
auto x63 = -x54 + x57;
auto x64 = u[0] + u[2];
auto x65 = V{4.5}*(x64*x64);
auto x66 = V{0.166666666666667}*pi[2];
auto x67 = V{0.0555555555555556}*pi[5];
auto x68 = x56 + x67 - V{0.0277777777777778}*pi[3];
auto x69 = x66 + x68;
auto x70 = -x46;
auto x71 = x59 + u[2];
auto x72 = -x66 + x68;
auto x73 = u[1] + u[2];
auto x74 = V{4.5}*(x73*x73);
auto x75 = x34 + x41;
auto x76 = x46 + x75;
auto x77 = V{0.166666666666667}*pi[4];
auto x78 = x55 + x67 - V{0.0277777777777778}*pi[0] + V{-0.0185185185185185};
auto x79 = x77 + x78;
auto x80 = -u[1];
auto x81 = x80 + u[2];
auto x82 = -x77 + x78;
auto x83 = V{0.00462962962962963}*rho;
auto x84 = x50 + u[2];
auto x85 = V{4.5}*(x84*x84);
auto x86 = V{0.0416666666666667}*pi[2];
auto x87 = V{0.0416666666666667}*pi[4];
auto x88 = V{0.0416666666666667}*pi[1];
auto x89 = V{0.0138888888888889}*pi[0];
auto x90 = V{0.0138888888888889}*pi[3];
auto x91 = V{0.0138888888888889}*pi[5];
auto x92 = x88 + x89 + x90 + x91 + V{-0.00462962962962963};
auto x93 = x86 + x87 + x92;
auto x94 = x71 + x80;
auto x95 = -x86;
auto x96 = -x87;
auto x97 = x92 + x95 + x96;
auto x98 = -u[2];
auto x99 = x60 + x98;
auto x100 = -x88 + x89 + x90 + x91 + V{-0.00462962962962963};
auto x101 = x100 + x86 + x96;
auto x102 = -x36;
auto x103 = x60 + u[2];
auto x104 = -x30;
auto x105 = V{1} - x32;
auto x106 = x104 + x105;
auto x107 = -x28;
auto x108 = x107 + x41;
auto x109 = x106 + x108 + x46;
auto x110 = x100 + x87 + x95;
auto x111 = x106 + x36;
auto x112 = x107 + x46;
auto x113 = x108 + x111;
auto x114 = x80 + u[0];
auto x115 = x111 + x112;
auto x116 = x98 + u[0];
auto x117 = x34 + x46;
auto x118 = x98 + u[1];
auto x119 = x50 + x98;
auto x120 = x64 + x80;
auto x121 = x114 + x98;
cell[0] = -V{0.296296296296296}*rho*x34 - V{0.444444444444444}*pi[0] - V{0.444444444444444}*pi[3] - V{0.444444444444444}*pi[5] + V{-0.296296296296296};
cell[1] = -x35*(x33 + x36 - x37) + x40;
cell[2] = -x35*(x32 + x41 - x42 + x43) + x45;
cell[3] = -x35*(x30 + x43 + x46 - x47) + x48;
cell[4] = -x49*(-x51 + x53) + x58;
cell[5] = -x49*(x62 - V{4.5}*x60*x60) + x63;
cell[6] = -x49*(x46 + x52 - x65) + x69;
cell[7] = -x49*(x52 + x70 - V{4.5}*x71*x71) + x72;
cell[8] = -x49*(-x74 + x76) + x79;
cell[9] = -x49*(x70 + x75 - V{4.5}*x81*x81) + x82;
cell[10] = -x83*(x46 + x53 - x85) + x93;
cell[11] = -x83*(x53 + x70 - V{4.5}*x94*x94) + x97;
cell[12] = x101 - x83*(x46 + x62 - V{4.5}*x99*x99);
cell[13] = x110 + x83*(x102 + x109 + V{4.5}*(x103*x103));
cell[14] = x35*(x111 + x37) + x40;
cell[15] = x35*(x105 + x108 + x42) + x45;
cell[16] = x35*(x104 + x112 + x47 + V{1}) + x48;
cell[17] = x49*(x113 + x51) + x58;
cell[18] = -x49*(x102 + x75 - V{4.5}*x114*x114) + x63;
cell[19] = x49*(x115 + x65) + x69;
cell[20] = -x49*(x102 + x117 - V{4.5}*x116*x116) + x72;
cell[21] = x49*(x109 + x74) + x79;
cell[22] = -x49*(x117 + x61 - V{4.5}*x118*x118) + x82;
cell[23] = x83*(x113 + x46 + x85) + x93;
cell[24] = x83*(x113 + x70 + V{4.5}*(x119*x119)) + x97;
cell[25] = x101 + x83*(x115 + x61 + V{4.5}*(x120*x120));
cell[26] = x110 - x83*(x102 + x76 - V{4.5}*x121*x121);

}

template <typename CELL, typename FORCE, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell, FORCE& force) any_platform
{
auto x0 = V{1}*cell[12];
auto x1 = V{1}*cell[25];
auto x2 = -cell[23];
auto x3 = x2 + cell[10] + cell[12] - cell[19] - cell[25] + cell[6];
auto x4 = -cell[13] - cell[21] + cell[26] + cell[8];
auto x5 = cell[20] + cell[22] + cell[24] + cell[3];
auto x6 = x3 + x4 + x5 - cell[11] - cell[16] - cell[7] - cell[9];
auto x7 = cell[11] - cell[17] - cell[24] + cell[4];
auto x8 = cell[13] + cell[1] + cell[5] + cell[7];
auto x9 = x3 + x7 + x8 - cell[14] - cell[18] - cell[20] - cell[26];
auto x10 = cell[10] + cell[18] + cell[25] + cell[2] + cell[9];
auto x11 = x10 + x5 + x8 + cell[0] + cell[11] + cell[12] + cell[14] + cell[15] + cell[16] + cell[17] + cell[19] + cell[21] + cell[23] + cell[26] + cell[4] + cell[6] + cell[8];
auto x12 = V{1} / (x11 + V{1});
auto x13 = x12*(x11 + V{1});
auto x14 = V{0.5}*x13;
auto x15 = x12*x9;
auto x16 = V{1}*x15;
auto x17 = V{1}*cell[13];
auto x18 = V{1}*cell[26];
auto x19 = -V{1}*cell[10];
auto x20 = -V{1}*cell[23];
auto x21 = x17 + x18 + x19 + x20;
auto x22 = V{1}*cell[11];
auto x23 = V{1}*cell[24];
auto x24 = x22 + x23;
auto x25 = -x0 - x1 + x14*(x6*force[0] + x9*force[2]) + x16*x6 + x21 + x24 - V{1}*cell[19] + V{1}*cell[20] - V{1}*cell[6] + V{1}*cell[7];
auto x26 = x10 + x2 + x4 + x7 - cell[12] - cell[15] - cell[22] - cell[5];
auto x27 = x26*force[0] + x9*force[1];
auto x28 = x0 + x1;
auto x29 = V{2}*cell[13];
auto x30 = V{2}*cell[26];
auto x31 = V{2}*cell[11];
auto x32 = V{2}*cell[24];
auto x33 = V{1}*x13;
auto x34 = V{2}*x26;
auto x35 = -V{2}*cell[10] + V{2}*cell[12] - V{2}*cell[23] + V{2}*cell[25];
auto x36 = x26*force[2] + x6*force[1];
auto x37 = x12*x6;
auto x38 = V{0.666666666666667}*cell[10];
auto x39 = V{0.666666666666667}*cell[11];
auto x40 = V{0.666666666666667}*cell[12];
auto x41 = V{0.666666666666667}*cell[13];
auto x42 = V{0.666666666666667}*cell[23];
auto x43 = V{0.666666666666667}*cell[24];
auto x44 = V{0.666666666666667}*cell[25];
auto x45 = V{0.666666666666667}*cell[26];
auto x46 = -V{0.333333333333333}*cell[0];
auto x47 = x38 + x39 + x40 + x41 + x42 + x43 + x44 + x45 + x46 - V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5];
auto x48 = -V{0.333333333333333}*cell[15] + V{0.666666666666667}*cell[19] + V{0.666666666666667}*cell[20] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x49 = -x12*x9*x9 - x13*x9*force[0] + x47 + x48 + V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[21] - V{0.333333333333333}*cell[22] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9];
auto x50 = -V{0.333333333333333}*cell[14] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[21] + V{0.666666666666667}*cell[22] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
auto x51 = -x12*x26*x26 - x13*x26*force[1] + x47 + x50 + V{0.666666666666667}*cell[15] - V{0.333333333333333}*cell[19] - V{0.333333333333333}*cell[20] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7];
auto x52 = -x12*x6*x6 - x13*x6*force[2] + x38 + x39 + x40 + x41 + x42 + x43 + x44 + x45 + x46 + x48 + x50 + V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5];
return (x14*x27 + x16*x26 + x21 - x22 - x23 + x28 - V{1}*cell[17] + V{1}*cell[18] - V{1}*cell[4] + V{1}*cell[5])*(x15*x34 + x27*x33 + x29 + x30 - x31 - x32 + x35 - V{2}*cell[17] + V{2}*cell[18] - V{2}*cell[4] + V{2}*cell[5]) + (-x29 - x30 + x31 + x32 + x33*x36 + x34*x37 + x35 - V{2}*cell[21] + V{2}*cell[22] - V{2}*cell[8] + V{2}*cell[9])*(x14*x36 - x17 - x18 + x19 + x20 + x24 + V{1}*x26*x37 + x28 - V{1}*cell[21] + V{1}*cell[22] - V{1}*cell[8] + V{1}*cell[9]) + 2*(x25*x25) + x49*x49 + x51*x51 + x52*x52;
}

template <typename CELL, typename V=typename CELL::value_t>
static auto computePiNeqNormSqr(CELL& cell) any_platform
{
auto x0 = -cell[11];
auto x1 = -cell[12];
auto x2 = -cell[17];
auto x3 = -cell[24];
auto x4 = x2 + x3 + cell[18];
auto x5 = cell[10] + cell[11] + cell[4];
auto x6 = -cell[23];
auto x7 = -cell[13] - cell[21];
auto x8 = x6 + x7 + cell[26] + cell[8];
auto x9 = x1 + x4 + x5 + x8 - cell[15] - cell[22] + cell[25] + cell[2] - cell[5] + cell[9];
auto x10 = -cell[19];
auto x11 = -cell[25];
auto x12 = x10 + x11 + cell[7];
auto x13 = x6 - cell[26];
auto x14 = cell[12] + cell[6];
auto x15 = cell[13] + cell[1];
auto x16 = x12 + x13 + x14 + x15 + x2 + x3 + x5 - cell[14] - cell[18] - cell[20] + cell[5];
auto x17 = cell[12] + cell[25];
auto x18 = x17 + cell[5];
auto x19 = cell[11] + cell[24];
auto x20 = x19 + cell[20];
auto x21 = cell[22] + cell[9];
auto x22 = cell[10] + cell[3];
auto x23 = V{1} / (x15 + x18 + x20 + x21 + x22 + cell[0] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[19] + cell[21] + cell[23] + cell[26] + cell[2] + cell[4] + cell[6] + cell[7] + cell[8] + V{1});
auto x24 = x16*x23;
auto x25 = -cell[10];
auto x26 = x25 + x6 + cell[13] + cell[26];
auto x27 = x0 + x18 + x24*x9 + x26 + x4 - cell[4];
auto x28 = x0 + x10 + x11 + x14 + x22 + x8 - cell[16] + cell[20] + cell[22] + cell[24] - cell[7] - cell[9];
auto x29 = x1 + x12 + x20 + x24*x28 + x26 - cell[6];
auto x30 = x13 + x17 + x19 + x21 + x23*x28*x9 + x25 + x7 - cell[8];
auto x31 = V{1}*x23;
auto x32 = V{0.666666666666667}*cell[10];
auto x33 = V{0.666666666666667}*cell[11];
auto x34 = V{0.666666666666667}*cell[12];
auto x35 = V{0.666666666666667}*cell[13];
auto x36 = V{0.666666666666667}*cell[23];
auto x37 = V{0.666666666666667}*cell[24];
auto x38 = V{0.666666666666667}*cell[25];
auto x39 = V{0.666666666666667}*cell[26];
auto x40 = -V{0.333333333333333}*cell[0];
auto x41 = x32 + x33 + x34 + x35 + x36 + x37 + x38 + x39 + x40 - V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5];
auto x42 = -V{0.333333333333333}*cell[15] + V{0.666666666666667}*cell[19] + V{0.666666666666667}*cell[20] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x43 = -x31*x16*x16 + x41 + x42 + V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[21] - V{0.333333333333333}*cell[22] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9];
auto x44 = -V{0.333333333333333}*cell[14] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[21] + V{0.666666666666667}*cell[22] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
auto x45 = -x31*x9*x9 + x41 + x44 + V{0.666666666666667}*cell[15] - V{0.333333333333333}*cell[19] - V{0.333333333333333}*cell[20] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7];
auto x46 = -x31*x28*x28 + x32 + x33 + x34 + x35 + x36 + x37 + x38 + x39 + x40 + x42 + x44 + V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5];
return V{2}*(x27*x27) + V{2}*(x29*x29) + V{2}*(x30*x30) + x43*x43 + x45*x45 + x46*x46;
}

template <typename CELL, typename RHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
static void addExternalForce(CELL& cell, RHO& rho, U& u, OMEGA& omega, FORCE& force) any_platform
{
auto x27 = force[0]*u[0];
auto x28 = force[1]*u[1];
auto x29 = force[2]*u[2];
auto x30 = rho*(V{0.5}*omega + V{-1});
auto x31 = V{6}*u[0];
auto x32 = x31 + V{-3};
auto x33 = V{0.0740740740740741}*force[0];
auto x34 = V{0.222222222222222}*x28;
auto x35 = V{0.222222222222222}*x29;
auto x36 = x34 + x35;
auto x37 = V{6}*u[1];
auto x38 = x37 + V{-3};
auto x39 = V{0.0740740740740741}*force[1];
auto x40 = V{0.222222222222222}*x27;
auto x41 = x35 + x40;
auto x42 = V{6}*u[2];
auto x43 = x42 + V{-3};
auto x44 = V{0.0740740740740741}*force[2];
auto x45 = x34 + x40;
auto x46 = V{9}*u[1];
auto x47 = x32 + x46;
auto x48 = V{0.0185185185185185}*force[0];
auto x49 = V{9}*u[0];
auto x50 = x38 + x49;
auto x51 = V{0.0185185185185185}*force[1];
auto x52 = -V{0.0555555555555556}*x29;
auto x53 = -x49;
auto x54 = x37 + V{3};
auto x55 = x53 + x54;
auto x56 = V{3} - x31;
auto x57 = x46 + x56;
auto x58 = V{9}*u[2];
auto x59 = x32 + x58;
auto x60 = x43 + x49;
auto x61 = V{0.0185185185185185}*force[2];
auto x62 = -V{0.0555555555555556}*x28;
auto x63 = x42 + V{3};
auto x64 = x53 + x63;
auto x65 = -V{0.0555555555555556}*x27;
auto x66 = -x46;
auto x67 = x63 + x66;
auto x68 = -x37;
auto x69 = x68 + V{3};
auto x70 = x58 + x69;
auto x71 = V{0.00462962962962963}*x30;
auto x72 = -x58;
auto x73 = -x42;
auto x74 = x49 + V{-3};
auto x75 = x31 + V{3};
auto x76 = x46 + x75;
auto x77 = x49 + x54;
auto x78 = x66 + x75;
auto x79 = x49 + x69;
auto x80 = x49 + x63;
auto x81 = x73 + V{3};
auto x82 = x49 + x81;
cell[0] = V{0.888888888888889}*x30*(x27 + x28 + x29) + cell[0];
cell[1] = x30*(-x32*x33 + x36) + cell[1];
cell[2] = x30*(-x38*x39 + x41) + cell[2];
cell[3] = x30*(-x43*x44 + x45) + cell[3];
cell[4] = -x30*(x47*x48 + x50*x51 + x52) + cell[4];
cell[5] = -x30*(-x48*x57 + x51*x55 + x52) + cell[5];
cell[6] = -x30*(x48*x59 + x60*x61 + x62) + cell[6];
cell[7] = -x30*(-x48*(x56 + x58) + x61*x64 + x62) + cell[7];
cell[8] = -x30*(x51*(x38 + x58) + x61*(x43 + x46) + x65) + cell[8];
cell[9] = -x30*(-x51*x70 + x61*x67 + x65) + cell[9];
cell[10] = -x71*((x46 + x60)*force[2] + (x47 + x58)*force[0] + (x50 + x58)*force[1]) + cell[10];
cell[11] = -x71*((x47 + x72)*force[0] + (x50 + x72)*force[1] - (x46 + x73 + x74)*force[2]) + cell[11];
cell[12] = -x71*((x59 + x66)*force[0] + (x60 + x66)*force[2] - (x58 + x68 + x74)*force[1]) + cell[12];
cell[13] = -x71*((x46 + x64)*force[2] + (x55 + x58)*force[1] - (x57 + x58)*force[0]) + cell[13];
cell[14] = x30*(-x33*x75 + x36) + cell[14];
cell[15] = x30*(-x39*x54 + x41) + cell[15];
cell[16] = x30*(-x44*x63 + x45) + cell[16];
cell[17] = -x30*(x48*x76 + x51*x77 + x52) + cell[17];
cell[18] = -x30*(x48*x78 - x51*x79 + x52) + cell[18];
cell[19] = -x30*(x48*(x58 + x75) + x61*x80 + x62) + cell[19];
cell[20] = -x30*(x48*(x72 + x75) - x61*x82 + x62) + cell[20];
cell[21] = -x30*(x51*(x54 + x58) + x61*(x46 + x63) + x65) + cell[21];
cell[22] = -x30*(x51*(x54 + x72) - x61*(x46 + x81) + x65) + cell[22];
cell[23] = -x71*((x46 + x80)*force[2] + (x58 + x76)*force[0] + (x58 + x77)*force[1]) + cell[23];
cell[24] = -x71*(-(x46 + x82)*force[2] + (x72 + x76)*force[0] + (x72 + x77)*force[1]) + cell[24];
cell[25] = -x71*((x49 + x67)*force[2] - (x49 + x70)*force[1] + (x58 + x78)*force[0]) + cell[25];
cell[26] = -x71*(-(x66 + x82)*force[2] + (x72 + x78)*force[0] - (x72 + x79)*force[1]) + cell[26];

}

};

}

#endif

#endif
