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

#ifndef SIMD_PACK_512_H_
#define SIMD_PACK_512_H_

#include <immintrin.h>

#include <cstdint>
#include <type_traits>

namespace olb {

namespace cpu {

namespace simd {

template <typename T> class Mask;
template <typename T> class Pack;

template <>
class Mask<double> {
private:
  __mmask8 _reg;

public:
  using storage_t = std::uint8_t;
  static constexpr unsigned storage_size = 8;

  static storage_t encode(bool* value)
  {
    storage_t mask = value[0];
    for (unsigned j=1; j < storage_size; ++j) {
      mask |= value[j] << j;
    }
    return mask;
  }

  Mask(bool b0, bool b1, bool b2, bool b3, bool b4, bool b5, bool b6, bool b7):
    _reg(std::uint16_t(b0 | b1<<1 | b2<<2 | b3<<3 | b4<<4 | b5<<5 | b6<<6 | b7<<7)) { }

  Mask(std::uint8_t* ptr):
    _reg(_load_mask16(reinterpret_cast<std::uint16_t*>(ptr))) { }

  Mask(storage_t* ptr, std::size_t iCell):
    Mask(ptr + iCell / storage_size) { }

  Mask(__mmask8 reg):
    _reg(reg) { }

  operator __mmask8()
  {
    return _reg;
  }

  __mmask8 neg() const
  {
    return _knot_mask8(_reg);
  }

  operator bool() const
  {
    const std::uint8_t* value = reinterpret_cast<const std::uint8_t*>(&_reg);
    return value[0] != 0;
  }
};

template <>
class Mask<float> {
private:
  __mmask16 _reg;

public:
  using storage_t = std::uint16_t;
  static constexpr unsigned storage_size = 16;

  static storage_t encode(bool* value)
  {
    storage_t mask = value[0];
    for (unsigned j=1; j < storage_size; ++j) {
      mask |= value[j] << j;
    }
    return mask;
  }

  Mask(std::uint16_t* ptr):
    _reg(_load_mask16(ptr)) { }

  Mask(storage_t* ptr, std::size_t iCell):
    Mask(ptr + iCell / storage_size) { }

  Mask(__mmask16 reg):
    _reg(reg) { }

  operator __mmask16()
  {
    return _reg;
  }

  __mmask16 neg() const
  {
    return _knot_mask16(_reg);
  }

  operator bool() const
  {
    const std::uint16_t* value = reinterpret_cast<const std::uint16_t*>(&_reg);
    return value[0] != 0;
  }
};


template <>
class Pack<double> : public SimdBase {
private:
  __m512d _reg;

public:
  using mask_t = Mask<double>;
  using index_t = std::uint32_t;

  static constexpr std::size_t size = 8;

  Pack() = default;

  Pack(__m512d reg):
    _reg(reg) { }

  Pack(double val):
    Pack(_mm512_set1_pd(val)) { }

  Pack(int val):
    Pack(static_cast<double>(val)) { }

  Pack(double a, double b, double c, double d, double e, double f, double g, double h):
    Pack(_mm512_set_pd(h,g,f,e,d,c,b,a)) { }

  Pack(const double* ptr):
    Pack(_mm512_loadu_pd(ptr)) { }

  Pack(const double* ptr, const index_t* idx):
    Pack(_mm512_i32gather_pd(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(idx)), ptr, sizeof(double))) { }

  operator __m512d()
  {
    return _reg;
  }

  Pack& operator=(Pack rhs)
  {
    _reg = rhs._reg;
    return *this;
  }

  double operator[](unsigned i) const
  {
    return reinterpret_cast<const double*>(&_reg)[i];
  }

  Pack operator+(Pack rhs) const
  {
    return Pack(_mm512_add_pd(_reg, rhs));
  }

  Pack& operator+=(Pack rhs)
  {
    _reg = _mm512_add_pd(_reg, rhs);
    return *this;
  }

  Pack operator-(Pack rhs) const
  {
    return Pack(_mm512_sub_pd(_reg, rhs));
  }

  Pack operator-=(Pack rhs)
  {
    _reg = _mm512_sub_pd(_reg, rhs);
    return *this;
  }

  Pack operator*(Pack rhs) const
  {
    return Pack(_mm512_mul_pd(_reg, rhs));
  }

  Pack& operator*=(Pack rhs)
  {
    _reg = _mm512_mul_pd(_reg, rhs);
    return *this;
  }

  Pack operator/(Pack rhs) const
  {
    return Pack(_mm512_div_pd(_reg, rhs));
  }

  Pack operator/=(Pack rhs)
  {
    _reg = _mm512_div_pd(_reg, rhs);
    return *this;
  }

  Pack operator-() const
  {
    return *this * Pack(-1);
  }

  Pack sqrt() const
  {
    return _mm512_sqrt_pd(_reg);
  }
};

template <>
class Pack<float> : public SimdBase {
private:
  __m512 _reg;

public:
  using mask_t = Mask<float>;
  using index_t = std::uint32_t;

  static constexpr std::size_t size = 16;

  Pack() = default;

  Pack(__m512 reg):
    _reg(reg) { }

  Pack(float val):
    Pack(_mm512_set1_ps(val)) { }

  Pack(double val):
    Pack(static_cast<float>(val)) { }

  Pack(int val):
    Pack(static_cast<float>(val)) { }

  Pack(float a, float b, float c, float d, float e, float f, float g, float h, float i, float j, float k, float l, float m, float n, float o, float p):
    Pack(_mm512_set_ps(p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)) { }

  Pack(const float* ptr):
    Pack(_mm512_loadu_ps(ptr)) { }

  Pack(const float* ptr, const index_t* idx):
    Pack(_mm512_i32gather_ps(_mm512_loadu_si512(reinterpret_cast<const __m512i*>(idx)), ptr, sizeof(float))) { }

  operator __m512()
  {
    return _reg;
  }

  Pack& operator=(Pack rhs)
  {
    _reg = rhs._reg;
    return *this;
  }

  float operator[](unsigned i) const
  {
    return reinterpret_cast<const float*>(&_reg)[i];
  }

  Pack operator+(Pack rhs) const
  {
    return Pack(_mm512_add_ps(_reg, rhs));
  }

  Pack& operator+=(Pack rhs)
  {
    _reg = _mm512_add_ps(_reg, rhs);
    return *this;
  }

  Pack operator-(Pack rhs) const
  {
    return Pack(_mm512_sub_ps(_reg, rhs));
  }

  Pack& operator-=(Pack rhs)
  {
    _reg = _mm512_sub_ps(_reg, rhs);
    return *this;
  }

  Pack operator*(Pack rhs) const
  {
    return Pack(_mm512_mul_ps(_reg, rhs));
  }

  Pack& operator*=(Pack rhs)
  {
    _reg = _mm512_mul_ps(_reg, rhs);
    return *this;
  }

  Pack operator/(Pack rhs) const
  {
    return Pack(_mm512_div_ps(_reg, rhs));
  }

  Pack& operator/=(Pack rhs)
  {
    _reg = _mm512_div_ps(_reg, rhs);
    return *this;
  }

  Pack operator-() const
  {
    return *this * Pack(-1);
  }

  __m512 sqrt() const
  {
    return _mm512_sqrt_ps(_reg);
  }
};


template <typename T>
Pack<T> pow(Pack<T> base, Pack<T> exp)
{
  if constexpr (std::is_same_v<T,double>) {
    return _mm512_pow_pd(base, exp);
  } else {
    return _mm512_pow_ps(base, exp);
  }
}

template <typename T>
Pack<T> min(Pack<T> rhs, Pack<T> lhs)
{
  if constexpr (std::is_same_v<T,double>) {
    return _mm512_min_pd(rhs, lhs);
  } else {
    return _mm512_min_ps(rhs, lhs);
  }
}

template <typename T>
Pack<T> max(Pack<T> rhs, Pack<T> lhs)
{
  if constexpr (std::is_same_v<T,double>) {
    return _mm512_max_pd(rhs, lhs);
  } else {
    return _mm512_max_ps(rhs, lhs);
  }
}

template <typename T>
Pack<T> fabs(Pack<T> x)
{
  if constexpr (std::is_same_v<T,double>) {
    return _mm512_abs_pd(x);
  } else {
    return _mm512_abs_ps(x);
  }
}

template <typename T>
void maskstore(T* target, Mask<T> mask, Pack<T> value);

template <>
void maskstore<double>(double* target, Mask<double> mask, Pack<double> value)
{
  _mm512_mask_storeu_pd(target, mask, value);
}

template <>
void maskstore<float>(float* target, Mask<float> mask, Pack<float> value)
{
  _mm512_mask_storeu_ps(target, mask, value);
}


template <typename T>
void store(T* target, Pack<T> value);

template <>
void store<double>(double* target, Pack<double> value)
{
  _mm512_storeu_pd(target, value);
}

template <>
void store<float>(float* target, Pack<float> value)
{
  _mm512_storeu_ps(target, value);
}


template <typename T>
void store(T* target, Pack<T> value, const typename Pack<T>::index_t* indices);

template <>
void store<double>(double* target, Pack<double> value, const Pack<double>::index_t* indices)
{
  _mm512_i32scatter_pd(target, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(indices)), value, sizeof(double));
}


template <>
void store<float>(float* target, Pack<float> value, const Pack<float>::index_t* indices)
{
  _mm512_i32scatter_ps(target, _mm512_loadu_si512(reinterpret_cast<const __m512i*>(indices)), value, sizeof(float));
}

}

}

}

#endif
