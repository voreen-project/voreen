/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/core/version.h"

#include "tgt/logmanager.h"
#include <sstream>

namespace voreen {

const std::string VoreenVersion::getCompilerVersion() {
    std::stringstream ver;
#if defined(WIN32) && defined(_MSC_VER)
    #if _MSC_VER == 1900
        ver << "Microsoft Visual C++ 2015";
    #elif _MSC_VER >= 1910
        ver << "Microsoft Visual C++ 2017";
    #else
        ver << "Unknown Microsoft Visual C++ (_MSC_VER=" << _MSC_VER << ")";
    #endif
#elif defined(__GNUC__)
    ver << "gcc " << __VERSION__;
#else
    ver << "Unknown compiler";
#endif

    // detect CPU architecture the application was compiled with
#if defined(__amd64__) || defined(__x86_64__) || defined(_M_X64)
    ver << " [x86_64]";
#elif defined(__i386__) || defined(_M_IX86) || defined(_X86_) || defined(__INTEL__)
    ver << " [x86]";
#elif defined(__powerpc__)
    ver << " [PowerPC]";
#else
    ver << " [unknown CPU architecture]";
#endif

    return ver.str();
}

const std::string VoreenVersion::getVersion() {
    return "5.1.1";
}

const std::string VoreenVersion::getRevision() {
#include "src/core/gen_revision.inc"
}

const std::string VoreenVersion::getCopyright() {
    return "Copyright (C) 2005-2019 University of Münster, Germany, \nDepartment of Computer Science.";
}

void VoreenVersion::logAll(const std::string& loggerCat) {
    LINFOC(loggerCat, "Version: " << getVersion());
    if (getRevision() != "")
        LINFOC(loggerCat, "Revision: " << getRevision());
    LINFOC(loggerCat, "Compiler: " << getCompilerVersion());
}

} // namespace voreen
