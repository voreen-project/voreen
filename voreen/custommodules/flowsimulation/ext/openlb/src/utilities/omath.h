/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Jan E. Marquardt
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

#ifndef OLB_OMATH_H
#define OLB_OMATH_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <type_traits>

namespace olb {

namespace util {

// Pow
template <typename T, typename S>
any_platform inline auto pow(T base, S exp) -> std::enable_if_t<std::is_arithmetic_v<T>,decltype(std::pow(base,exp))>
{
  return std::pow(base, exp);
}

inline float powf(float base, float exp)
{
  return pow(base, exp);
}

inline long double powl(long double base, long double exp)
{
  return pow(base, exp);
}

//fmod
inline float fmod(float x, float y)
{
  return std::fmod(x, y);
}

inline float fmodf(float x, float y)
{
  return fmod(x, y);
}

inline double fmod(double x, double y)
{
  return std::fmod(x, y);
}

inline long double fmod(long double x, long double y)
{
  return std::fmod(x, y);
}

inline long double fmodl(long double x, long double y)
{
  return fmod(x, y);
}

inline float fmod(float x, int y)
{
  return std::fmod(x, y);
}

inline double fmod(double x, int y)
{
  return std::fmod(x, y);
}

inline long double fmod(long double x, int y)
{
  return std::fmod(x, y);
}

template <typename T, typename S>
inline auto fmod(T x, S y)
{
  return std::fmod(x, y);
}


// Exp
inline float exp(float arg)
{
  return std::exp(arg);
}

inline float expf(float arg)
{
  return exp(arg);
}

inline double exp(double arg)
{
  return std::exp(arg);
}

inline long double exp(long double arg)
{
  return std::exp(arg);
}

inline long double expl(long double arg)
{
  return exp(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> exp(T arg)
{
  return std::exp(arg);
}

// Log
inline float log(float arg)
{
  return std::log(arg);
}

inline float logf(float arg)
{
  return log(arg);
}

inline double log(double arg)
{
  return std::log(arg);
}

inline long double log(long double arg)
{
  return std::log(arg);
}

inline long double logl(long double arg)
{
  return log(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> log(T arg)
{
  return std::log(arg);
}

// Log10
inline float log10(float arg)
{
  return std::log10(arg);
}

inline float log10f(float arg)
{
  return log10(arg);
}

inline double log10(double arg)
{
  return std::log10(arg);
}

inline long double log10(long double arg)
{
  return std::log10(arg);
}

inline long double log10l(long double arg)
{
  return log10(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> log10(T arg)
{
  return std::log10(arg);
}

// Log2
inline float log2(float arg)
{
  return std::log2(arg);
}

inline float log2f(float arg)
{
  return log2(arg);
}

inline double log2(double arg)
{
  return std::log2(arg);
}

inline long double log2(long double arg)
{
  return std::log2(arg);
}

inline long double log2l(long double arg)
{
  return log2(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> log2(T arg)
{
  return std::log2(arg);
}

// Log1p
inline float log1p(float arg)
{
  return std::log1p(arg);
}

inline float log1pf(float arg)
{
  return log1p(arg);
}

inline double log1p(double arg)
{
  return std::log1p(arg);
}

inline long double log1p(long double arg)
{
  return std::log1p(arg);
}

inline long double log1pl(long double arg)
{
  return log1p(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> log1p(T arg)
{
  return std::log1p(arg);
}

// Sqrt
inline float sqrtf(float arg)
{
  return sqrt(arg);
}

inline long double sqrtl(long double arg)
{
  return sqrt(arg);
}

template <typename T>
inline std::enable_if_t<std::is_floating_point_v<T>, T> sqrt(T arg) any_platform
{
  return std::sqrt(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> sqrt(T arg) any_platform
{
  return std::sqrt(arg);
}

// Sin
inline float sin(float arg)
{
  return std::sin(arg);
}

inline float sinf(float arg)
{
  return sin(arg);
}

inline double sin(double arg)
{
  return std::sin(arg);
}

inline long double sin(long double arg)
{
  return std::sin(arg);
}

inline long double sinl(long double arg)
{
  return sin(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> sin(T arg)
{
  return std::sin(arg);
}

// Sinh
inline float sinh(float arg)
{
  return std::sinh(arg);
}

inline float sinhf(float arg)
{
  return sinh(arg);
}

inline double sinh(double arg)
{
  return std::sinh(arg);
}

inline long double sinh(long double arg)
{
  return std::sinh(arg);
}

inline long double sinhl(long double arg)
{
  return sinh(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> sinh(T arg)
{
  return std::sinh(arg);
}

// Cos
inline float cos(float arg)
{
  return std::cos(arg);
}

inline float cosf(float arg)
{
  return cos(arg);
}

inline double cos(double arg)
{
  return std::cos(arg);
}

inline long double cos(long double arg)
{
  return std::cos(arg);
}

inline long double cosl(long double arg)
{
  return cos(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> cos(T arg)
{
  return std::cos(arg);
}

// Cosh
inline float cosh(float arg)
{
  return std::cosh(arg);
}

inline float coshf(float arg)
{
  return cosh(arg);
}

inline double cosh(double arg)
{
  return std::cosh(arg);
}

inline long double cosh(long double arg)
{
  return std::cosh(arg);
}

inline long double coshl(long double arg)
{
  return cosh(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> cosh(T arg)
{
  return std::cosh(arg);
}

// Tan
inline float tan(float arg)
{
  return std::tan(arg);
}

inline float tanf(float arg)
{
  return tan(arg);
}

inline double tan(double arg)
{
  return std::tan(arg);
}

inline long double tan(long double arg)
{
  return std::tan(arg);
}

inline long double tanl(long double arg)
{
  return tan(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> tan(T arg)
{
  return std::tan(arg);
}

// Tanh
inline float tanh(float arg)
{
  return std::tanh(arg);
}

inline float tanhf(float arg)
{
  return tanh(arg);
}

inline double tanh(double arg)
{
  return std::tanh(arg);
}

inline long double tanh(long double arg)
{
  return std::tanh(arg);
}

inline long double tanhl(long double arg)
{
  return tanh(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> tanh(T arg)
{
  return std::tanh(arg);
}

// Asin
inline float asin(float arg)
{
  return std::asin(arg);
}

inline float asinf(float arg)
{
  return asin(arg);
}

inline double asin(double arg)
{
  return std::asin(arg);
}

inline long double asin(long double arg)
{
  return std::asin(arg);
}

inline long double asinl(long double arg)
{
  return asin(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> asin(T arg)
{
  return std::asin(arg);
}

// Asinh
inline float asinh(float arg)
{
  return std::asinh(arg);
}

inline float asinhf(float arg)
{
  return asinh(arg);
}

inline double asinh(double arg)
{
  return std::asinh(arg);
}

inline long double asinh(long double arg)
{
  return std::asinh(arg);
}

inline long double asinhl(long double arg)
{
  return asinh(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> asinh(T arg)
{
  return std::asinh(arg);
}

// Acos
inline float acos(float arg)
{
  return std::acos(arg);
}

inline float acosf(float arg)
{
  return acos(arg);
}

inline double acos(double arg)
{
  return std::acos(arg);
}

inline long double acos(long double arg)
{
  return std::acos(arg);
}

inline long double acosl(long double arg)
{
  return acos(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> acos(T arg)
{
  return std::acos(arg);
}

// Acosh
inline float acosh(float arg)
{
  return std::acosh(arg);
}

inline float acoshf(float arg)
{
  return acosh(arg);
}

inline double acosh(double arg)
{
  return std::acosh(arg);
}

inline long double acosh(long double arg)
{
  return std::acosh(arg);
}

inline long double acoshl(long double arg)
{
  return acosh(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> acosh(T arg)
{
  return std::acosh(arg);
}

// Atan
inline float atan(float arg)
{
  return std::atan(arg);
}

inline float atanf(float arg)
{
  return atan(arg);
}

inline double atan(double arg)
{
  return std::atan(arg);
}

inline long double atan(long double arg)
{
  return std::atan(arg);
}

inline long double atanl(long double arg)
{
  return atan(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> atan(T arg)
{
  return std::atan(arg);
}

// Atan2
inline float atan2(float y, float x)
{
  return std::atan2(y, x);
}

inline float atan2f(float y, float x)
{
  return atan2(y, x);
}

inline double atan2(double y, double x)
{
  return std::atan2(y, x);
}

inline long double atan2(long double y, long double x)
{
  return std::atan2(y, x);
}

inline long double atan2l(long double y, long double x)
{
  return atan2(y, x);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> atan2(T y, T x)
{
  return std::atan2(y, x);
}

// Atanh
inline float atanh(float arg)
{
  return std::atanh(arg);
}

inline float atanhf(float arg)
{
  return atanh(arg);
}

inline double atanh(double arg)
{
  return std::atanh(arg);
}

inline long double atanh(long double arg)
{
  return std::atanh(arg);
}

inline long double atanhl(long double arg)
{
  return atanh(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> atanh(T arg)
{
  return std::atanh(arg);
}


// Rounding
// floor
inline float floor(float arg)
{
  return std::floor(arg);
}

inline double floor(double arg)
{
  return std::floor(arg);
}

inline long double floor(long double arg)
{
  return std::floor(arg);
}

inline float floorf(float arg)
{
  return floor(arg);
}

inline long double floorl(long double arg)
{
  return floor(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> floor(T arg)
{
  return std::floor(arg);
}

// ceil
inline float ceil(float arg)
{
  return std::ceil(arg);
}

inline double ceil(double arg)
{
  return std::ceil(arg);
}

inline long double ceil(long double arg)
{
  return std::ceil(arg);
}

inline float ceilf(float arg)
{
  return ceil(arg);
}

inline long double ceill(long double arg)
{
  return ceil(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> ceil(T arg)
{
  return std::ceil(arg);
}

// trunc
inline float trunc(float arg)
{
  return std::trunc(arg);
}

inline float truncf(float arg)
{
  return trunc(arg);
}

inline double trunc(double arg)
{
  return std::trunc(arg);
}

inline long double trunc(long double arg)
{
  return std::trunc(arg);
}

inline long double truncl(long double arg)
{
  return trunc(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> trunc(T arg)
{
  return std::trunc(arg);
}

// round
template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> round(T arg)
{
  return std::round(arg);
}

inline float round(float arg)
{
  return std::round(arg);
}

inline double round(double arg)
{
  return std::round(arg);
}

inline long double round(long double arg)
{
  return std::round(arg);
}

inline float roundf(float arg)
{
  return round(arg);
}

inline double roundf(double arg)
{
  return round(arg);
}

inline long double roundl(long double arg)
{
  return round(arg);
}

template <typename T>
inline std::enable_if_t<std::is_arithmetic_v<T>, long> lround(T arg)
{
  return std::lround(arg);
}

inline long lroundf(float arg)
{
  return lround(arg);
}

inline long lroundl(long double arg)
{
  return lround(arg);
}

template <typename T>
inline std::enable_if_t<std::is_arithmetic_v<T>, long long> llround(T arg)
{
  return std::llround(arg);
}

inline long long llroundf(float arg)
{
  return llround(arg);
}

inline long long llroundl(long double arg)
{
  return llround(arg);
}


// Absolute
template <typename T>
inline std::enable_if_t<std::is_arithmetic_v<T>, T> abs(T arg) any_platform
{
  return std::abs(arg);
}

inline long labs(long arg)
{
  return abs(arg);
}

inline long long llabs(long long arg)
{
  return abs(arg);
}

inline float fabsf(float arg)
{
  return fabs(arg);
}

inline long double fabs(long double arg)
{
  return std::fabs(arg);
}

inline long double fabsl(long double arg)
{
  return fabs(arg);
}

template <typename T>
inline std::enable_if_t<std::is_floating_point_v<T>, T> fabs(T arg) any_platform
{
  return std::fabs(arg);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, double> fabs(T arg) any_platform
{
  return std::fabs(arg);
}

} // namespace util

} // namespace olb

#endif
