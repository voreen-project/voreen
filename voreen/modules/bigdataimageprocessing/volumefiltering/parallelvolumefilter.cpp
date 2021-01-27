/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "parallelvolumefilter.h"
namespace voreen {

// FilterValue -----------------------------------------------------------------------

// Some specializations for FilterValue1D


// VisualStudio does not conform with the c++ standard and wants the definition of static
// member specializations as part of the declaration somehow, even though the standard says:
//
// Every program shall contain exactly one definition of every non-inline function or object that is used in that program; no diagnostic required.
//
// (See http://stackoverflow.com/questions/2342550/static-member-initialization-for-specialized-template-class)
//
// The following definition is therefore only done in the .cpp file for non-windows platforms.
#ifndef WIN32
template<>
const size_t ParallelFilterValue<float>::dim = 1;
#endif


} // namespace voreen
