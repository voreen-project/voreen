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

#include "tgt/arch.h"

// Taken and adapted from "portable snippets" by Evan Nemerson (licensed CC0).

namespace tgt {

#define ONCE__BACKEND_ATOMIC  1
#define ONCE__BACKEND_PTHREAD 2
#define ONCE__BACKEND_NONE    3
#define ONCE__BACKEND_C11     11
#define ONCE__BACKEND_WIN32   32

#if defined(__has_include)
#  if __has_include(<pthread.h>)
#    include<pthread.h>
#  endif
#endif

#include <limits.h>

#if !defined(ONCE_BACKEND)
#  if defined(__STDC_NO_THREADS__) && __STDC_NO_THREADS__
#  elif defined(__EMSCRIPTEN__)
#  elif defined(__has_include)
#    if __has_include(<threads.h>)
#      if !defined(__GLIBC__) || (defined(__GLIBC__) && defined(ENABLE_PTHREADS))
#        include <threads.h>
#        define ONCE_BACKEND ONCE__BACKEND_C11
#      endif
#    endif
#  elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201102L) && !defined(__STDC_NO_THREADS__)
#    if (defined(__GLIBC__) && (__GLIBC__ < 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ < 16)))
/* stdc-predef.h didn't include __STDC_NO_THREADS__ until 2.16. */
#    else
#      include <threads.h>
#      define ONCE_BACKEND ONCE__BACKEND_C11
#    endif
#  endif
#endif

#if !defined(ONCE_BACKEND) && defined(_WIN32) && (!defined(WINVER) || (defined(WINVER) && (WINVER >= 0x0600)))
#  include <Windows.h>
#  define ONCE_BACKEND ONCE__BACKEND_WIN32
#endif

#if !defined(ONCE_BACKEND) && defined(PTHREAD_ONCE_INIT)
#  define ONCE_BACKEND ONCE__BACKEND_PTHREAD
#endif

#if !defined(ONCE_BACKEND)
#  include "../atomic/atomic.h"
#  if !defined(ATOMIC_NOT_FOUND)
#    define ONCE_BACKEND ONCE__BACKEND_ATOMIC
#  endif
#endif

#if !defined(ONCE_BACKEND)
#  error No once backend found.
#endif

#if defined(__GNUC__) && (__GNUC__ >= 3)
#  define ONCE__UNLIKELY(expr) __builtin_expect(!!(expr), !!0)
#else
#  define ONCE__UNLIKELY(expr) (!!(expr))
#endif

#if ONCE_BACKEND == ONCE__BACKEND_C11
#  define ONCE_INIT ONCE_FLAG_INIT
typedef once_flag once;
#  define once_call(flag, func) call_once(flag, func)
#elif ONCE_BACKEND == ONCE__BACKEND_PTHREAD
#  define ONCE_INIT PTHREAD_ONCE_INIT
typedef pthread_once_t once;
#  define once_call(flag, func) pthread_once(flag, func)
#elif ONCE_BACKEND == ONCE__BACKEND_WIN32
#  define ONCE_INIT INIT_ONCE_STATIC_INIT
typedef INIT_ONCE once;
static BOOL CALLBACK once__callback_wrap(INIT_ONCE* InitOnce, void* Parameter, void** Context) {
  (void) Context;
  (void) InitOnce;
#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable:4055)
#endif
  ((void (*)(void)) Parameter)();
#if defined(_MSC_VER)
#  pragma warning(pop)
#endif
  return !0;
}
#  if defined(_MSC_VER) && (_MSC_VER >= 1500)
#    define once_call(flag, func) \
  __pragma(warning(push)) \
  __pragma(warning(disable:4152)) \
  InitOnceExecuteOnce(flag, &once__callback_wrap, func, NULL) \
  __pragma(warning(pop))
#  else
#    define once_call(flag, func) InitOnceExecuteOnce(flag, &once__callback_wrap, func, NULL)
#  endif
#elif ONCE_BACKEND == ONCE__BACKEND_ATOMIC
#  define ONCE_INIT ATOMIC_VAR_INIT(0)
typedef atomic_int32 once;
static void once_call(once* flag, void (*func)(void)) {
  int32_t state = atomic_int32_load(flag);
  if (ONCE__UNLIKELY(state == 0)) {
    if (atomic_int32_compare_exchange(flag, &state, 1)) {
      func();
      atomic_int32_store(flag, 2);
    } else {
      do {
	/* Spin; another thread is calling the initialization
	   function. */
      } while (atomic_int32_load(flag) == 1);
    }
  }
}
#elif ONCE_BACKEND == ONCE__BACKEND_NONE
#  define ONCE_INIT 0
typedef int once;
static void once_call(once* flag, void (*func)(void)) {
  if (*flag == 0) {
    func();
    *flag = 1;
  }
}
#endif

#include <assert.h>

#if defined(_WIN32)
#  include <Windows.h>
#  define CPU__IMPL_WIN32
#elif defined(unix) || defined(__unix__) || defined(__unix)
#  include <unistd.h>
#  if defined(_SC_NPROCESSORS_ONLN) || defined(_SC_NPROC_ONLN)
#    define CPU__IMPL_SYSCONF
#  else
#    include <sys/sysctl.h>
#  endif
#endif

#if defined(CPU_ARCH_X86) || defined(CPU_ARCH_X86_64)
#  if defined(_MSC_VER)
static void cpu_getid(int func, int* data) {
  __cpuid(data, func);
}
#  else
static void cpu_getid(int func, int* data) {
  __asm__ ("cpuid"
	   : "=a" (data[0]), "=b" (data[1]), "=c" (data[2]), "=d" (data[3])
	   : "0" (func), "2" (0));
}
#  endif
#elif defined(CPU_ARCH_ARM) || defined(CPU_ARCH_ARM64)
#  if (defined(__GNUC__) && ((__GNUC__ > 2) || (__GNUC__ == 2 && __GNUC_MINOR__ >= 16)))
#    define CPU__IMPL_GETAUXVAL
#    include <sys/auxv.h>
#  endif
#endif

static once cpu_once = ONCE_INIT;

#if defined(CPU_ARCH_X86) || defined(CPU_ARCH_X86_64)
static unsigned int cpuinfo[8 * 4] = { 0, };
#elif defined(CPU_ARCH_ARM) || defined(CPU_ARCH_ARM_64)
static unsigned long cpuinfo[2] = { 0, };
#endif

static void cpu_init(void) {
#if defined(CPU_ARCH_X86) || defined(CPU_ARCH_X86_64)
  int i;
  for (i = 0 ; i < 8 ; i++) {
    cpu_getid(i, (int*) &(cpuinfo[i * 4]));
  }
#elif defined(CPU_ARCH_ARM) || defined(CPU_ARCH_ARM_64)
  cpuinfo[0] = getauxval (AT_HWCAP);
  cpuinfo[1] = getauxval (AT_HWCAP2);
#endif
}

int cpuFeatureCheck (enum PSnipCPUFeature feature) {
#if defined(CPU_ARCH_X86) || defined(CPU_ARCH_X86_64)
  unsigned int i, r, b;
#elif defined(CPU_ARCH_ARM) || defined(CPU_ARCH_ARM_64)
  unsigned long b, i;
#endif

#if defined(CPU_ARCH_X86) || defined(CPU_ARCH_X86_64)
  if ((feature & CPU_FEATURE_CPU_MASK) != CPU_FEATURE_X86)
    return 0;
#elif defined(CPU_ARCH_ARM) || defined(CPU_ARCH_ARM_64)
  if ((feature & CPU_FEATURE_CPU_MASK) != CPU_FEATURE_ARM)
    return 0;
#else
  return 0;
#endif

  feature = (enum PSnipCPUFeature) (feature & ~CPU_FEATURE_CPU_MASK);
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4152)
#endif
  once_call (&cpu_once, cpu_init);
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

#if defined(CPU_ARCH_X86) || defined(CPU_ARCH_X86_64)
  i = (feature >> 16) & 0xff;
  r = (feature >>  8) & 0xff;
  b = (feature      ) & 0xff;

  if (i > 7 || r > 3 || b > 31)
    return 0;

  return (cpuinfo[(i * 4) + r] >> b) & 1;
#elif defined(CPU_ARCH_ARM) || defined(CPU_ARCH_ARM_64)
  b = 1 << ((feature & 0xff) - 1);
  i = cpuinfo[(feature >> 0x08) & 0xff];
  return (cpuinfo[(feature >> 0x08) & 0xff] & b) == b;
#endif
}

int cpuFeatureCheckMany (enum PSnipCPUFeature* feature) {
  int n;

  for (n = 0 ; feature[n] != CPU_FEATURE_NONE ; n++)
    if (!cpuFeatureCheck(feature[n]))
      return 0;

  return 1;
}

int cpuCount (void) {
  static int count = 0;
  int c;

#if defined(_WIN32)
  DWORD_PTR lpProcessAffinityMask;
  DWORD_PTR lpSystemAffinityMask;
  int i;
#elif defined(CPU__IMPL_SYSCONF) && defined(HW_NCPU)
  int mib[2];
  size_t len;
#endif

  if (count != 0)
    return count;

#if defined(_WIN32)
  if (!GetProcessAffinityMask(GetCurrentProcess(), &lpProcessAffinityMask, &lpSystemAffinityMask)) {
    c = -1;
  } else {
    c = 0;
    for (i = 0 ; lpProcessAffinityMask != 0 ; lpProcessAffinityMask >>= 1)
      c += lpProcessAffinityMask & 1;
  }
#elif defined(_SC_NPROCESSORS_ONLN)
  c = sysconf (_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROC_ONLN)
  c = sysconf (_SC_NPROC_ONLN);
#elif defined(_hpux)
  c = mpctl(MPC_GETNUMSPUS, NULL, NULL);
#elif defined(HW_NCPU)
  c = 0;
  mib[0] = CTL_HW;
  mib[1] = HW_NCPU;
  len = sizeof(c);
  sysctl (mib, 2, &c, &len, NULL, 0);
#endif

  count = (c > 0) ? c : -1;

  return count;
}

} // namespace tgt
