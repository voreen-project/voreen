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

#include "randomwalkermodule.h"

#include "processors/octreewalker.h"
#include "processors/randomwalker.h"
#include "processors/randomwalkeranalyzer.h"
#include "processors/supervoxelwalker.h"

#include "processors/rwmultilabelloopinitializer.h"
#include "processors/rwmultilabelloopfinalizer.h"

#ifdef VRN_RW_USE_MAGMA
#include "magma_v2.h"
#endif


namespace voreen {

RandomWalkerModule::RandomWalkerModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Random Walker");
    setGuiName("Random Walker");

    registerSerializableType(new OctreeWalker());
    registerSerializableType(new RandomWalker());
    registerSerializableType(new RandomWalkerAnalyzer());
    registerSerializableType(new SuperVoxelWalker());

    registerSerializableType(new RWMultiLabelLoopInitializer());
    registerSerializableType(new RWMultiLabelLoopFinalizer());

#ifdef VRN_RW_USE_MAGMA
    magma_init();
#endif
}

} // namespace
