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

#include "cpucapabilities.h"
namespace tgt{
    // from http://stackoverflow.com/questions/6121792/how-to-check-if-a-cpu-supports-the-sse3-instruction-set
#ifdef _WIN32

    //  Windows
#include <intrin.h>
#define cpuid(info,x)    __cpuidex(info,x,0)

#else

    //  GCC Inline Assembly
    static void cpuid(int CPUInfo[4],int InfoType){
        __asm__ __volatile__ (
            "cpuid":
        "=a" (CPUInfo[0]),
            "=b" (CPUInfo[1]),
            "=c" (CPUInfo[2]),
            "=d" (CPUInfo[3]) :
        "a" (InfoType), "c" (0)
            );
    }

#endif
    bool CPUCapabilities::hasSSE(){
        return HW_SSE;
    }
    bool CPUCapabilities::hasSSE2(){
        return HW_SSE2;
    }
    bool CPUCapabilities::hasSSE3(){
        return HW_SSE3;
    }
    bool CPUCapabilities::hasSSSE3(){
        return HW_SSSE3;
    }
    bool CPUCapabilities::hasSSE41(){
        return HW_SSE41;
    }
    bool CPUCapabilities::hasSSE42(){
        return HW_SSE42;
    }
    bool CPUCapabilities::hasSSE4a(){
        return HW_SSE4a;
    }
    bool CPUCapabilities::hasAES(){
        return HW_AES;
    }
    bool CPUCapabilities::hasSHA(){
        return HW_SHA;
    }

    bool CPUCapabilities::hasAVX(){
        return HW_AVX;
    }
    bool CPUCapabilities::hasXOP(){
        return HW_XOP;
    }
    bool CPUCapabilities::hasFMA3(){
        return HW_FMA3;
    }
    bool CPUCapabilities::hasFMA4(){
        return HW_FMA4;
    }
    bool CPUCapabilities::hasAVX2(){
        return HW_AVX2;
    }

    bool CPUCapabilities::hasAVX512F(){
        return HW_AVX512F;
    }
    bool CPUCapabilities::hasAVX512CD(){
        return HW_AVX512CD;
    }
    bool CPUCapabilities::hasAVX512PF(){
        return HW_AVX512PF;
    }
    bool CPUCapabilities::hasAVX512ER(){
        return HW_AVX512ER;
    }
    bool CPUCapabilities::hasAVX512VL(){
        return HW_AVX512VL;
    }
    bool CPUCapabilities::hasAVX512BW(){
        return HW_AVX512BW;
    }
    bool CPUCapabilities::hasAVX512DQ(){
        return HW_AVX512DQ;
    }
    bool CPUCapabilities::hasAVX512IFMA(){
        return HW_AVX512IFMA;
    }
    bool CPUCapabilities::hasAVX512VBMI(){
        return HW_AVX512VBMI;
    }
    CPUCapabilities::CPUCapabilities(){
        static bool inizialized = false;
        if (inizialized){
            return;
        }else{
            inizialized = true;
        }
        int info[4];
        cpuid(info, 0);
        int nIds = info[0];

        cpuid(info, 0x80000000);
        unsigned nExIds = info[0];

        //  Detect Features
        if (nIds >= 0x00000001){
            cpuid(info,0x00000001);
            HW_MMX    = (info[3] & ((int)1 << 23)) != 0;
            HW_SSE    = (info[3] & ((int)1 << 25)) != 0;
            HW_SSE2   = (info[3] & ((int)1 << 26)) != 0;
            HW_SSE3   = (info[2] & ((int)1 <<  0)) != 0;

            HW_SSSE3  = (info[2] & ((int)1 <<  9)) != 0;
            HW_SSE41  = (info[2] & ((int)1 << 19)) != 0;
            HW_SSE42  = (info[2] & ((int)1 << 20)) != 0;
            HW_AES    = (info[2] & ((int)1 << 25)) != 0;

            HW_AVX    = (info[2] & ((int)1 << 28)) != 0;
            HW_FMA3   = (info[2] & ((int)1 << 12)) != 0;

            HW_RDRAND = (info[2] & ((int)1 << 30)) != 0;
        }
        if (nIds >= 0x00000007){
            cpuid(info,0x00000007);
            HW_AVX2   = (info[1] & ((int)1 <<  5)) != 0;

            HW_BMI1        = (info[1] & ((int)1 <<  3)) != 0;
            HW_BMI2        = (info[1] & ((int)1 <<  8)) != 0;
            HW_ADX         = (info[1] & ((int)1 << 19)) != 0;
            HW_SHA         = (info[1] & ((int)1 << 29)) != 0;
            HW_PREFETCHWT1 = (info[2] & ((int)1 <<  0)) != 0;

            HW_AVX512F     = (info[1] & ((int)1 << 16)) != 0;
            HW_AVX512CD    = (info[1] & ((int)1 << 28)) != 0;
            HW_AVX512PF    = (info[1] & ((int)1 << 26)) != 0;
            HW_AVX512ER    = (info[1] & ((int)1 << 27)) != 0;
            HW_AVX512VL    = (info[1] & ((int)1 << 31)) != 0;
            HW_AVX512BW    = (info[1] & ((int)1 << 30)) != 0;
            HW_AVX512DQ    = (info[1] & ((int)1 << 17)) != 0;
            HW_AVX512IFMA  = (info[1] & ((int)1 << 21)) != 0;
            HW_AVX512VBMI  = (info[2] & ((int)1 <<  1)) != 0;
        }
        if (nExIds >= 0x80000001){
            cpuid(info,0x80000001);
            HW_x64   = (info[3] & ((int)1 << 29)) != 0;
            HW_ABM   = (info[2] & ((int)1 <<  5)) != 0;
            HW_SSE4a = (info[2] & ((int)1 <<  6)) != 0;
            HW_FMA4  = (info[2] & ((int)1 << 16)) != 0;
            HW_XOP   = (info[2] & ((int)1 << 11)) != 0;
        }
    }


}
