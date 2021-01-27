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

#ifndef TGT_CPUID_H
#define TGT_CPUID_H
#include "tgt/singleton.h"



namespace tgt{
    #ifdef DLL_TEMPLATE_INST
    class CPUCapabilities;
    template class TGT_API Singleton<CPUCapabilities>;
    #endif

    class TGT_API CPUCapabilities: public tgt::Singleton<CPUCapabilities>{
    public:
        CPUCapabilities();
        bool hasSSE();
        bool hasSSE2();
        bool hasSSE3();
        bool hasSSSE3();
        bool hasSSE41();
        bool hasSSE42();
        bool hasSSE4a();
        bool hasAES();
        bool hasSHA();

        bool hasAVX();
        bool hasXOP();
        bool hasFMA3();
        bool hasFMA4();
        bool hasAVX2();

        bool hasAVX512F();
        bool hasAVX512CD();
        bool hasAVX512PF();
        bool hasAVX512ER();
        bool hasAVX512VL();
        bool hasAVX512BW();
        bool hasAVX512DQ();
        bool hasAVX512IFMA();
        bool hasAVX512VBMI();
    private:

        //  Misc.
        bool HW_MMX;
        bool HW_x64;
        bool HW_ABM;      // Advanced Bit Manipulation
        bool HW_RDRAND;
        bool HW_BMI1;
        bool HW_BMI2;
        bool HW_ADX;
        bool HW_PREFETCHWT1;

        //  SIMD: 128-bit
        bool HW_SSE;
        bool HW_SSE2;
        bool HW_SSE3;
        bool HW_SSSE3;
        bool HW_SSE41;
        bool HW_SSE42;
        bool HW_SSE4a;
        bool HW_AES;
        bool HW_SHA;

        //  SIMD: 256-bit
        bool HW_AVX;
        bool HW_XOP;
        bool HW_FMA3;
        bool HW_FMA4;
        bool HW_AVX2;

        //  SIMD: 512-bit
        bool HW_AVX512F;    //  AVX512 Foundation
        bool HW_AVX512CD;   //  AVX512 Conflict Detection
        bool HW_AVX512PF;   //  AVX512 Prefetch
        bool HW_AVX512ER;   //  AVX512 Exponential + Reciprocal
        bool HW_AVX512VL;   //  AVX512 Vector Length Extensions
        bool HW_AVX512BW;   //  AVX512 Byte + Word
        bool HW_AVX512DQ;   //  AVX512 Doubleword + Quadword
        bool HW_AVX512IFMA; //  AVX512 Integer 52-bit Fused Multiply-Add
        bool HW_AVX512VBMI; //  AVX512 Vector Byte Manipulation Instructions

    };
}
#endif
