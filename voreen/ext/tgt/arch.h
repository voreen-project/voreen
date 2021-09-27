/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2021 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#ifndef TGT_ARCH_H
#define TGT_ARCH_H

#include "tgt/types.h"

/// Use the following macros to surround structs/functions in order to compile
/// them with the desired cpu features (e.g., "avx2"). Then, inside, you can use
/// intrinsics specific to these architectures (see #include <immintrin.h>) and
/// use the functions below for runtime dispatch. Something like:
///
/// BEGIN_STRUCT_WITH_TARGET_FEATURE("avx2")
/// struct Foo {
///     __m256 bar;
/// };
/// END_STRUCT_WITH_TARGET_FEATURE()
///
/// BEGIN_FUNC_WITH_TARGET_FEATURE("avx2")
/// void myfunc_avx2(float foo) {
///     Foo { _mm_set1_ps(foo) };
///     // ...
/// }
/// END_FUNC_WITH_TARGET_FEATURE()
///
/// void myfunc_generic(float foo) {
///     if(tgt::cpuFeatureCheck(tgt::CPU_FEATURE_X86_AVX2)) {
///         myfunc_avx2(foo);
///     } else {
///         // Fallback...
///     }
/// }

/// ----------------------------------------------------------------------------
/// BEGIN_FUNC_WITH_TARGET macro definition ------------------------------------
/// ----------------------------------------------------------------------------
#define DO_PRAGMA_(x) _Pragma (#x)
#define DO_PRAGMA(x) DO_PRAGMA_(x)
#if defined(__clang__)
#define BEGIN_FUNC_WITH_TARGET_FEATURE(_target)\
    DO_PRAGMA(clang attribute push (__attribute__((target(_target))), apply_to=function))
#define END_FUNC_WITH_TARGET_FEATURE()\
    DO_PRAGMA(clang attribute pop)
#elif defined(__GNUC__)
#define BEGIN_FUNC_WITH_TARGET_FEATURE(_target)\
    DO_PRAGMA(GCC push_options)\
    DO_PRAGMA(GCC target(_target))
#define END_FUNC_WITH_TARGET_FEATURE()\
    DO_PRAGMA(GCC pop_options)
#elif defined(_MSC_VER)
// This macro is intended for use with avx intrinsics, but MSVC compiles them
// fine without activating the target (and apparently there is no way to
// activate a target for a single function).
#define BEGIN_FUNC_WITH_TARGET_FEATURE(_target)
#define END_FUNC_WITH_TARGET_FEATURE()
#else
    #warning Compiler could not be detected. Function cannot be compiled for the desired target.
#endif

/// ----------------------------------------------------------------------------
/// BEGIN_STRUCT_WITH_TARGET macro definition ----------------------------------
/// ----------------------------------------------------------------------------
#if defined(__clang__)
// Nothing required for clang
#define BEGIN_STRUCT_WITH_TARGET_FEATURE(_target)
#define END_STRUCT_WITH_TARGET_FEATURE()
#elif defined(__GNUC__)
#define BEGIN_STRUCT_WITH_TARGET_FEATURE(_target)\
    DO_PRAGMA(GCC push_options)\
    DO_PRAGMA(GCC target(_target))
#define END_STRUCT_WITH_TARGET_FEATURE()\
    DO_PRAGMA(GCC pop_options)
#elif defined(_MSC_VER)
// This macro is intended for use with avx intrinsics, but MSVC compiles them
// fine without activating the target (and apparently there is no way to
// activate a target for a single function).
#define BEGIN_STRUCT_WITH_TARGET_FEATURE(_target)
#define END_STRUCT_WITH_TARGET_FEATURE()
#else
    #warning Compiler could not be detected. Struct cannot be compiled for the desired target.
#endif

#if defined(__x86_64__) || defined(_M_AMD64) || defined(_M_X64)
#  define CPU_ARCH_X86_64
#elif defined(__i686__) || defined(__i586__) || defined(__i486__) || defined(__i386__) || defined(__i386) || defined(_M_IX86) || defined(_X86_) || defined(__THW_INTEL__)
#  define CPU_ARCH_X86
#elif defined(__arm__) || defined(_M_ARM)
#  define CPU_ARCH_ARM
#elif defined(__aarch64__)
#  define CPU_ARCH_ARM64
#endif

// Taken and adapted from "portable snippets" by Evan Nemerson (licensed CC0).
namespace tgt {

enum PSnipCPUFeature {
  CPU_FEATURE_NONE                = 0,

  CPU_FEATURE_CPU_MASK            = 0x1f000000,
  CPU_FEATURE_X86                 = 0x01000000,
  CPU_FEATURE_ARM                 = 0x04000000,

  /* x86 CPU features are constructed as:
   *
   *   (CPU_FEATURE_X86 | (eax << 16) | (ret_reg << 8) | (bit_position)
   *
   * For example, SSE3 is determined by the fist bit in the ECX
   * register for a CPUID call with EAX=1, so we get:
   *
   *   CPU_FEATURE_X86 | (1 << 16) | (2 << 8) | (0) = 0x01010200
   *
   * We should have information for inputs of EAX=0-7 w/ ECX=0.
   */
  CPU_FEATURE_X86_FPU             = 0x01010300,
  CPU_FEATURE_X86_VME             = 0x01010301,
  CPU_FEATURE_X86_DE              = 0x01010302,
  CPU_FEATURE_X86_PSE             = 0x01010303,
  CPU_FEATURE_X86_TSC             = 0x01010304,
  CPU_FEATURE_X86_MSR             = 0x01010305,
  CPU_FEATURE_X86_PAE             = 0x01010306,
  CPU_FEATURE_X86_MCE             = 0x01010307,
  CPU_FEATURE_X86_CX8             = 0x01010308,
  CPU_FEATURE_X86_APIC            = 0x01010309,
  CPU_FEATURE_X86_SEP             = 0x0101030b,
  CPU_FEATURE_X86_MTRR            = 0x0101030c,
  CPU_FEATURE_X86_PGE             = 0x0101030d,
  CPU_FEATURE_X86_MCA             = 0x0101030e,
  CPU_FEATURE_X86_CMOV            = 0x0101030f,
  CPU_FEATURE_X86_PAT             = 0x01010310,
  CPU_FEATURE_X86_PSE_36          = 0x01010311,
  CPU_FEATURE_X86_PSN             = 0x01010312,
  CPU_FEATURE_X86_CLFSH           = 0x01010313,
  CPU_FEATURE_X86_DS              = 0x01010314,
  CPU_FEATURE_X86_ACPI            = 0x01010316,
  CPU_FEATURE_X86_MMX             = 0x01010317,
  CPU_FEATURE_X86_FXSR            = 0x01010318,
  CPU_FEATURE_X86_SSE             = 0x01010319,
  CPU_FEATURE_X86_SSE2            = 0x0101031a,
  CPU_FEATURE_X86_SS              = 0x0101031b,
  CPU_FEATURE_X86_HTT             = 0x0101031c,
  CPU_FEATURE_X86_TM              = 0x0101031d,
  CPU_FEATURE_X86_IA64            = 0x0101031e,
  CPU_FEATURE_X86_PBE             = 0x0101031f,

  CPU_FEATURE_X86_SSE3            = 0x01010200,
  CPU_FEATURE_X86_PCLMULQDQ       = 0x01010201,
  CPU_FEATURE_X86_DTES64          = 0x01010202,
  CPU_FEATURE_X86_MONITOR         = 0x01010203,
  CPU_FEATURE_X86_DS_CPL          = 0x01010204,
  CPU_FEATURE_X86_VMX             = 0x01010205,
  CPU_FEATURE_X86_SMX             = 0x01010206,
  CPU_FEATURE_X86_EST             = 0x01010207,
  CPU_FEATURE_X86_TM2             = 0x01010208,
  CPU_FEATURE_X86_SSSE3           = 0x01010209,
  CPU_FEATURE_X86_CNXT_ID         = 0x0101020a,
  CPU_FEATURE_X86_SDBG            = 0x0101020b,
  CPU_FEATURE_X86_FMA             = 0x0101020c,
  CPU_FEATURE_X86_CX16            = 0x0101020d,
  CPU_FEATURE_X86_XTPR            = 0x0101020e,
  CPU_FEATURE_X86_PDCM            = 0x0101020f,
  CPU_FEATURE_X86_PCID            = 0x01010211,
  CPU_FEATURE_X86_DCA             = 0x01010212,
  CPU_FEATURE_X86_SSE4_1          = 0x01010213,
  CPU_FEATURE_X86_SSE4_2          = 0x01010214,
  CPU_FEATURE_X86_X2APIC          = 0x01010215,
  CPU_FEATURE_X86_MOVBE           = 0x01010216,
  CPU_FEATURE_X86_POPCNT          = 0x01010217,
  CPU_FEATURE_X86_TSC_DEADLINE    = 0x01010218,
  CPU_FEATURE_X86_AES             = 0x01010219,
  CPU_FEATURE_X86_XSAVE           = 0x0101021a,
  CPU_FEATURE_X86_OSXSAVE         = 0x0101021b,
  CPU_FEATURE_X86_AVX             = 0x0101021c,
  CPU_FEATURE_X86_F16C            = 0x0101021d,
  CPU_FEATURE_X86_RDRND           = 0x0101021e,
  CPU_FEATURE_X86_HYPERVISOR      = 0x0101021f,

  CPU_FEATURE_X86_FSGSBASE        = 0x01070100,
  CPU_FEATURE_X86_TSC_ADJ         = 0x01070101,
  CPU_FEATURE_X86_SGX             = 0x01070102,
  CPU_FEATURE_X86_BMI1            = 0x01070103,
  CPU_FEATURE_X86_HLE             = 0x01070104,
  CPU_FEATURE_X86_AVX2            = 0x01070105,
  CPU_FEATURE_X86_SMEP            = 0x01070107,
  CPU_FEATURE_X86_BMI2            = 0x01070108,
  CPU_FEATURE_X86_ERMS            = 0x01070109,
  CPU_FEATURE_X86_INVPCID         = 0x0107010a,
  CPU_FEATURE_X86_RTM             = 0x0107010b,
  CPU_FEATURE_X86_PQM             = 0x0107010c,
  CPU_FEATURE_X86_MPX             = 0x0107010e,
  CPU_FEATURE_X86_PQE             = 0x0107010f,
  CPU_FEATURE_X86_AVX512F         = 0x01070110,
  CPU_FEATURE_X86_AVX512DQ        = 0x01070111,
  CPU_FEATURE_X86_RDSEED          = 0x01070112,
  CPU_FEATURE_X86_ADX             = 0x01070113,
  CPU_FEATURE_X86_SMAP            = 0x01070114,
  CPU_FEATURE_X86_AVX512IFMA      = 0x01070115,
  CPU_FEATURE_X86_PCOMMIT         = 0x01070116,
  CPU_FEATURE_X86_CLFLUSHOPT      = 0x01070117,
  CPU_FEATURE_X86_CLWB            = 0x01070118,
  CPU_FEATURE_X86_INTEL_PT        = 0x01070119,
  CPU_FEATURE_X86_AVX512PF        = 0x0107011a,
  CPU_FEATURE_X86_AVX512ER        = 0x0107011b,
  CPU_FEATURE_X86_AVX512CD        = 0x0107011c,
  CPU_FEATURE_X86_SHA             = 0x0107011d,
  CPU_FEATURE_X86_AVX512BW        = 0x0107011e,
  CPU_FEATURE_X86_AVX512VL        = 0x0107011f,

  CPU_FEATURE_X86_PREFETCHWT1     = 0x01070200,
  CPU_FEATURE_X86_AVX512VBMI      = 0x01070201,
  CPU_FEATURE_X86_UMIP            = 0x01070202,
  CPU_FEATURE_X86_PKU             = 0x01070203,
  CPU_FEATURE_X86_OSPKE           = 0x01070204,
  CPU_FEATURE_X86_AVX512VPOPCNTDQ = 0x0107020e,
  CPU_FEATURE_X86_RDPID           = 0x01070215,
  CPU_FEATURE_X86_SGX_LC          = 0x0107021e,

  CPU_FEATURE_X86_AVX512_4VNNIW   = 0x01070302,
  CPU_FEATURE_X86_AVX512_4FMAPS   = 0x01070303,

  CPU_FEATURE_ARM_SWP             = CPU_FEATURE_ARM | 1,
  CPU_FEATURE_ARM_HALF            = CPU_FEATURE_ARM | 2,
  CPU_FEATURE_ARM_THUMB           = CPU_FEATURE_ARM | 3,
  CPU_FEATURE_ARM_26BIT           = CPU_FEATURE_ARM | 4,
  CPU_FEATURE_ARM_FAST_MULT       = CPU_FEATURE_ARM | 5,
  CPU_FEATURE_ARM_FPA             = CPU_FEATURE_ARM | 6,
  CPU_FEATURE_ARM_VFP             = CPU_FEATURE_ARM | 7,
  CPU_FEATURE_ARM_EDSP            = CPU_FEATURE_ARM | 8,
  CPU_FEATURE_ARM_JAVA            = CPU_FEATURE_ARM | 9,
  CPU_FEATURE_ARM_IWMMXT          = CPU_FEATURE_ARM | 10,
  CPU_FEATURE_ARM_CRUNCH          = CPU_FEATURE_ARM | 11,
  CPU_FEATURE_ARM_THUMBEE         = CPU_FEATURE_ARM | 12,
  CPU_FEATURE_ARM_NEON            = CPU_FEATURE_ARM | 13,
  CPU_FEATURE_ARM_VFPV3           = CPU_FEATURE_ARM | 14,
  CPU_FEATURE_ARM_VFPV3D16        = CPU_FEATURE_ARM | 15,
  CPU_FEATURE_ARM_TLS             = CPU_FEATURE_ARM | 16,
  CPU_FEATURE_ARM_VFPV4           = CPU_FEATURE_ARM | 17,
  CPU_FEATURE_ARM_IDIVA           = CPU_FEATURE_ARM | 18,
  CPU_FEATURE_ARM_IDIVT           = CPU_FEATURE_ARM | 19,
  CPU_FEATURE_ARM_VFPD32          = CPU_FEATURE_ARM | 20,
  CPU_FEATURE_ARM_LPAE            = CPU_FEATURE_ARM | 21,
  CPU_FEATURE_ARM_EVTSTRM         = CPU_FEATURE_ARM | 22,

  CPU_FEATURE_ARM_AES             = CPU_FEATURE_ARM | 0x0100 | 1,
  CPU_FEATURE_ARM_PMULL           = CPU_FEATURE_ARM | 0x0100 | 2,
  CPU_FEATURE_ARM_SHA1            = CPU_FEATURE_ARM | 0x0100 | 3,
  CPU_FEATURE_ARM_SHA2            = CPU_FEATURE_ARM | 0x0100 | 4,
  CPU_FEATURE_ARM_CRC32           = CPU_FEATURE_ARM | 0x0100 | 5
};

TGT_API int cpuCount();
TGT_API int cpuFeatureCheck(enum PSnipCPUFeature  feature);
TGT_API int cpuFeatureCheckMany(enum PSnipCPUFeature* feature);

}

#endif // TGT_ARCH_H
