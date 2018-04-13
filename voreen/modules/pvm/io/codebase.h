/*

This file is taken from V^3, the Versatile Volume Viewer.
Copyright (c) 2003-2007 by Stefan Roettger.

The volume viewer is licensed under the terms of the GPL (see
http://www.gnu.org/copyleft/ for more information on the license).
Any commercial use of the code or parts of it requires the explicit
permission of the author! The source code is distributed with
ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY etc.!

The author's contact address is:

   mailto:stefan@stereofx.org
   http://stereofx.org

*/

// (c) by Stefan Roettger

#ifndef CODEBASE_H
#define CODEBASE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#if defined(IRIX) || defined(LINUX) || defined(MACOSX)
#define UNIX
#endif

#ifdef _WIN32
#define WINOS
#endif

#include <time.h>
#ifdef UNIX
#include <sys/time.h>
#endif
#ifdef WINOS
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <winbase.h>
#endif

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "tgt/exception.h"
#include "voreen/core/utils/stringutils.h"

#ifndef NULL
#define NULL (0)
#endif

#define BOOLINT char

#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

//#define ERRORMSG() errormsg(__FILE__,__LINE__)

#define ERRORMSG() throw tgt::FileException("Fatal error in PVM reader (file " + std::string(__FILE__) \
                                            + " at line " + voreen::itos(__LINE__) + ")")

inline void errormsg(const char *file, int line) {
    fprintf(stderr, "fatal error in <%s> at line %d!\n", file, line);
    exit(EXIT_FAILURE);
}

#define PI (3.141593f)
#define RAD (PI/180.0f)

#ifndef MAXFLOAT
#define MAXFLOAT (FLT_MAX)
#endif

#undef ffloor
#define ffloor(x) floor((double)(x))
#undef fceil
#define fceil(x) ceil((double)(x))
#define ftrc(x) (int)ffloor(x)

inline double FABS(const double x) {
    return ((x < 0.0) ? -x : x);
}
#define fabs(x) FABS(x)

inline int min(const int a, const int b) {
    return ((a < b) ? a : b);
}
inline double FMIN(const double a, const double b) {
    return ((a < b) ? a : b);
}
#define fmin(a,b) FMIN(a,b)

inline int max(const int a, const int b) {
    return ((a > b) ? a : b);
}
inline double FMAX(const double a, const double b) {
    return ((a > b) ? a : b);
}
#define fmax(a,b) FMAX(a,b)

inline int sqr(const int x) {
    return (x*x);
}
inline double fsqr(const double x) {
    return (x*x);
}

#undef fsqrt
#define fsqrt(x) sqrt((double)(x))

#undef fsin
#define fsin(x) sin((double)(x))
#undef fcos
#define fcos(x) cos((double)(x))
#undef ftan
#define ftan(x) tan((double)(x))

#undef fasin
#define fasin(x) asin((double)(x))
#undef facos
#define facos(x) acos((double)(x))
#undef fatan
#define fatan(x) atan((double)(x))

#undef fexp
#define fexp(x) exp((double)(x))
#undef flog
#define flog(x) log((double)(x))
#undef fpow
#define fpow(x,y) pow((double)(x),(double)(y))

#ifdef UNIX
#define GETRANDOM() drand48()
#endif
#ifdef WINOS
#define GETRANDOM() ((double)rand()/RAND_MAX)
#endif

inline double GETTIME() {
#ifdef UNIX
    struct timeval t;
    gettimeofday(&t, NULL);
    return (t.tv_sec + t.tv_usec / 1.0E6);
#elif defined(WINOS)

    static int cpus = 0;
    if (cpus == 0) {
        SYSTEM_INFO SystemInfo;
        GetSystemInfo(&SystemInfo);
        cpus = SystemInfo.dwNumberOfProcessors;
    }
    if (cpus == 1) {
        LARGE_INTEGER freq, count;
        if (QueryPerformanceFrequency(&freq) == 0)
            ERRORMSG();
        QueryPerformanceCounter(&count);
        return ((double)count.QuadPart / freq.QuadPart);
    }
    return ((double)clock() / CLOCKS_PER_SEC);
#else
    #error Neither UNIX nor WINOS seems to be defined.
    return 0;
#endif
}

inline double gettime() {
    static double time;
    static BOOLINT settime = FALSE;

    if (!settime) {
        time = GETTIME();
        settime = TRUE;
    }

    return (GETTIME() - time);
}

inline void waitfor(double secs) {
#ifdef UNIX
    struct timespec dt, rt;
    dt.tv_sec = ftrc(secs);
    dt.tv_nsec = ftrc(1.0E9 * (secs - ftrc(secs)));
    while (nanosleep(&dt, &rt) != 0)
        dt = rt;
#else

    double time = gettime() + secs;
    while (gettime() < time)
        ;
#endif

}

inline double getclockticks() {
    static double clockticks;
    static BOOLINT setclockticks = FALSE;

    if (!setclockticks) {
        double time = gettime();
        while (time == gettime())
            ;
        clockticks = 1.0 / (gettime() - time);
        setclockticks = TRUE;
    }

    return (clockticks);
}

#ifdef WINOS

inline char *strdup(char *str) {
    size_t len;
    char *dup;
    if ((dup = (char *)malloc(len = strlen(str) + 1)) == NULL)
        ERRORMSG();
    memcpy(dup, str, len);
    return (dup);
}

inline int strcasecmp(char *str1, char *str2) {
    char *ptr1, *ptr2;
    for (ptr1 = str1, ptr2 = str2; tolower(*ptr1) == tolower(*ptr2) && *ptr1 != '\0' && *ptr2 != '\0'; ptr1++, ptr2++)
        ;
    return (*ptr2 - *ptr1);
}

#define snprintf _snprintf

#endif

#endif
