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

#ifndef SIMD_PACK_256_H_
#define SIMD_PACK_256_H_

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
  __m256i _reg;

public:
  using storage_t = std::uint64_t;
  static constexpr unsigned storage_size = 1;

  static constexpr storage_t true_v  = 1l << 63;
  static constexpr storage_t false_v = 0l;

  static storage_t encode(bool value)
  {
    return value ? true_v : false_v;
  }

  static storage_t encode(bool* value)
  {
    return encode(*value);
  }

  Mask(bool a, bool b, bool c, bool d):
    _reg(_mm256_set_epi64x(encode(d),encode(c),encode(b),encode(a))) { }

  Mask(std::uint64_t a, std::uint64_t b, std::uint64_t c, std::uint64_t d):
    _reg(_mm256_set_epi64x(d,c,b,a)) { }

  Mask(std::uint64_t* ptr):
    _reg(_mm256_loadu_si256(reinterpret_cast<__m256i*>(ptr))) { }

  Mask(storage_t* ptr, std::size_t iCell):
    Mask(ptr + iCell) { }

  Mask(bool* ptr):
    Mask(ptr[0],ptr[1],ptr[2],ptr[3]) { }

  Mask(bool* ptr, std::size_t iCell):
    Mask(ptr + iCell) { }

  Mask(__m256i reg):
    _reg(reg) { }

  operator __m256i()
  {
    return _reg;
  }

  __m256i neg() const
  {
    return _mm256_sub_epi64(_mm256_set1_epi64x(true_v), _reg);
  }

  operator bool() const
  {
    const std::uint64_t* values = reinterpret_cast<const std::uint64_t*>(&_reg);
    return values[0] == true_v
        || values[1] == true_v
        || values[2] == true_v
        || values[3] == true_v;
  }
};

template <>
class Mask<float> {
private:
  __m256i _reg;

public:
  using storage_t = std::uint32_t;
  static constexpr unsigned storage_size = 1;

  static constexpr storage_t true_v  = 1 << 31;
  static constexpr storage_t false_v = 0;

  static storage_t encode(bool value)
  {
    return value ? true_v : false_v;
  }

  static storage_t encode(bool* value)
  {
    return encode(*value);
  }

  Mask(bool a, bool b, bool c, bool d, bool e, bool f, bool g, bool h):
    _reg(_mm256_set_epi32(encode(h),encode(g),encode(f),encode(e),encode(d),encode(c),encode(b),encode(a))) { }

  Mask(storage_t* ptr):
    _reg(_mm256_loadu_si256(reinterpret_cast<__m256i*>(ptr))) { }

  Mask(storage_t* ptr, std::size_t iCell):
    Mask(ptr + iCell) { }

  Mask(bool* ptr):
    Mask(ptr[0],ptr[1],ptr[2],ptr[3],ptr[4],ptr[5],ptr[6],ptr[7]) { }

  Mask(bool* ptr, std::size_t iCell):
    Mask(ptr + iCell) { }

  Mask(__m256i reg):
    _reg(reg) { }

  operator __m256i()
  {
    return _reg;
  }

  __m256i neg() const
  {
    return _mm256_sub_epi32(_mm256_set1_epi32(true_v), _reg);
  }

  operator bool() const
  {
    const std::uint32_t* values = reinterpret_cast<const std::uint32_t*>(&_reg);
    return values[0] == true_v
        || values[1] == true_v
        || values[2] == true_v
        || values[3] == true_v
        || values[4] == true_v
        || values[5] == true_v
        || values[6] == true_v
        || values[7] == true_v;
  }
};

template <typename T> class Pack;

template <>
class Pack<double> : public SimdBase {
private:
  __m256d _reg;

public:
  using mask_t = Mask<double>;
  using index_t = std::uint32_t;

  static constexpr std::size_t size = 4;

  Pack() = default;

  Pack(__m256d reg):
    _reg(reg) { }

  Pack(double val):
    Pack(_mm256_set1_pd(val)) { }

  Pack(int val):
    Pack(static_cast<double>(val)) { }

  Pack(std::size_t val):
    Pack(static_cast<double>(val)) { }

  Pack(double a, double b, double c, double d):
    Pack(_mm256_set_pd(d,c,b,a)) { }

  Pack(const double* ptr):
    Pack(_mm256_loadu_pd(ptr)) { }

  Pack(const double* ptr, const index_t* idx):
    Pack(_mm256_i32gather_pd(ptr, _mm_loadu_si128(reinterpret_cast<const __m128i*>(idx)), sizeof(double))) { }

  operator __m256d()
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

  double& operator[](unsigned i)
  {
    return reinterpret_cast<double*>(&_reg)[i];
  }

  Pack operator+(Pack rhs) const
  {
    return Pack(_mm256_add_pd(_reg, rhs));
  }

  Pack& operator+=(Pack rhs)
  {
    _reg = _mm256_add_pd(_reg, rhs);
    return *this;
  }

  Pack operator-(Pack rhs) const
  {
    return Pack(_mm256_sub_pd(_reg, rhs));
  }

  Pack operator-=(Pack rhs)
  {
    _reg = _mm256_sub_pd(_reg, rhs);
    return *this;
  }

  Pack operator*(Pack rhs) const
  {
    return Pack(_mm256_mul_pd(_reg, rhs));
  }

  Pack& operator*=(Pack rhs)
  {
    _reg = _mm256_mul_pd(_reg, rhs);
    return *this;
  }

  Pack operator/(Pack rhs) const
  {
    return Pack(_mm256_div_pd(_reg, rhs));
  }

  Pack& operator/=(Pack rhs)
  {
    _reg = _mm256_div_pd(_reg, rhs);
    return *this;
  }

  Pack operator-() const
  {
    return *this * Pack(-1);
  }

  Pack sqrt() const
  {
    return _mm256_sqrt_pd(_reg);
  }
};

template <>
class Pack<float> : public SimdBase {
private:
  __m256 _reg;

public:
  using mask_t = Mask<float>;
  using index_t = std::uint32_t;

  static constexpr std::size_t size = 8;

  Pack() = default;

  Pack(__m256 reg):
    _reg(reg) { }

  Pack(float val):
    Pack(_mm256_set1_ps(val)) { }

  Pack(double val):
    Pack(static_cast<float>(val)) { }

  Pack(int val):
    Pack(static_cast<float>(val)) { }

  Pack(std::size_t val):
    Pack(static_cast<float>(val)) { }

  Pack(float a, float b, float c, float d, float e, float f, float g, float h):
    Pack(_mm256_set_ps(h,g,f,e,d,c,b,a)) { }

  Pack(const float* ptr):
    Pack(_mm256_loadu_ps(ptr)) { }

  Pack(const float* ptr, const index_t* idx):
    Pack(_mm256_i32gather_ps(ptr, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(idx)), sizeof(float))) { }

  operator __m256()
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

  float& operator[](unsigned i)
  {
    return reinterpret_cast<float*>(&_reg)[i];
  }

  Pack operator+(Pack rhs) const
  {
    return Pack(_mm256_add_ps(_reg, rhs));
  }

  Pack& operator+=(Pack rhs)
  {
    _reg = _mm256_add_ps(_reg, rhs);
    return *this;
  }

  Pack operator-(Pack rhs) const
  {
    return Pack(_mm256_sub_ps(_reg, rhs));
  }

  Pack& operator-=(Pack rhs)
  {
    _reg = _mm256_sub_ps(_reg, rhs);
    return *this;
  }

  Pack operator*(Pack rhs) const
  {
    return Pack(_mm256_mul_ps(_reg, rhs));
  }

  Pack& operator*=(Pack rhs)
  {
    _reg = _mm256_mul_ps(_reg, rhs);
    return *this;
  }

  Pack operator/(Pack rhs) const
  {
    return Pack(_mm256_div_ps(_reg, rhs));
  }

  Pack& operator/=(Pack rhs)
  {
    _reg = _mm256_div_ps(_reg, rhs);
    return *this;
  }

  Pack operator-() const
  {
    return *this * Pack(-1);
  }

  __m256 sqrt()
  {
    return _mm256_sqrt_ps(_reg);
  }
};


template <typename T>
Pack<T> pow(Pack<T> base, Pack<T> exp)
{
  // TODO: Replace by more efficient implementation
  Pack<T> result;
  for (unsigned i=0; i < Pack<T>::size; ++i) {
    result[i] = util::pow(base[i], exp[i]);
  }
  return result;
}

template <typename T>
Pack<T> min(Pack<T> rhs, Pack<T> lhs)
{
  if constexpr (std::is_same_v<T,double>) {
    return _mm256_min_pd(rhs, lhs);
  } else {
    return _mm256_min_ps(rhs, lhs);
  }
}

template <typename T>
Pack<T> max(Pack<T> rhs, Pack<T> lhs)
{
  if constexpr (std::is_same_v<T,double>) {
    return _mm256_max_pd(rhs, lhs);
  } else {
    return _mm256_max_ps(rhs, lhs);
  }
}

template <typename T>
Pack<T> fabs(Pack<T> x)
{
  return max(x, -x);
}

template <typename T>
void maskstore(T* target, Mask<T> mask, Pack<T> value);

template <>
void maskstore<double>(double* target, Mask<double> mask, Pack<double> value)
{
  _mm256_maskstore_pd(target, mask, value);
}

template <>
void maskstore<float>(float* target, Mask<float> mask, Pack<float> value)
{
  _mm256_maskstore_ps(target, mask, value);
}


template <typename T>
void store(T* target, Pack<T> value);

template <>
void store<double>(double* target, Pack<double> value)
{
  _mm256_storeu_pd(target, value);
}

template <>
void store<float>(float* target, Pack<float> value)
{
  _mm256_storeu_ps(target, value);
}


template <typename T>
void store(T* target, Pack<T> value, const typename Pack<T>::index_t* indices);

template <>
void store<double>(double* target, Pack<double> value, const Pack<double>::index_t* indices)
{
#ifdef __AVX512F__
  _mm256_i32scatter_pd(target, _mm_loadu_si128(reinterpret_cast<const __m128i*>(indices)), value, sizeof(double));
#else
  __m256d reg = value;
  for (unsigned i=0; i < simd::Pack<double>::size; ++i) {
    target[indices[i]] = reg[i];
  }
#endif
}

template <>
void store<float>(float* target, Pack<float> value, const Pack<float>::index_t* indices)
{
#ifdef __AVX512F__
  _mm256_i32scatter_ps(target, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(indices)), value, sizeof(float));
#else
  __m256 reg = value;
  for (unsigned i=0; i < simd::Pack<float>::size; ++i) {
    target[indices[i]] = reg[i];
  }
#endif
}

}

}

}

#endif
