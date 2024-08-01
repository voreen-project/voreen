/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#ifndef VRN_MEMORYINFO_H
#define VRN_MEMORYINFO_H

#include "voreen/core/voreencoreapi.h"

#include "tgt/types.h"

#include <string>

namespace voreen {

/**
 * Provides information about the total and available CPU RAM as well as
 * the CPU memory used by the current process.
 *
 * @note Currently only available on Win32.
 */
class VRN_CORE_API MemoryInfo {

public:

    /**
     * Returns total amount in bytes of the virtual memory the system provides,
     * i.e., the size of the physically installed RAM plus the swap memory size.
     */
    static uint64_t getTotalVirtualMemory();

    /// Returns the byte size of the system's virtual memory that is currently unused.
    static uint64_t getAvailableVirtualMemory();

    /// Returns the byte size of the system's virtual memory that is currently used the process.
    static uint64_t getVirtualMemoryUsedByCurrentProcess();

    /// Returns the size of the physically installed CPU RAM in bytes.
    static uint64_t getTotalPhysicalMemory();

    /// Returns the byte size of the physically installed CPU RAM that currently unused.
    static uint64_t getAvailablePhysicalMemory();

    /// Returns the byte size of the physically installed CPU RAM that is currently used by the process.
    static uint64_t getPhysicalMemoryUsedByCurrentProcess();

    /// Returns the machine's total virtual and physical CPU RAM as info string.
    static std::string getTotalMemoryAsString();

    /// Returns the currently available virtual and physical CPU RAM as info string.
    static std::string getAvailableMemoryAsString();

    /// Returns the virtual and physical CPU RAM usage of the current process as info string.
    static std::string getProcessMemoryUsageAsString();

protected:

    static const std::string loggerCat_;

};

} // namespace

#endif // VRN_MEMORYINFO_H
