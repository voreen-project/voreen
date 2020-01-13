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

/**
 * @file
 *
 * Precompiled Header File for voreenlib.
 *
 * The precompiled header must contain code which is *stable* and *static*
 * throughout the project.
 *
 */

/* Add C includes here */

#ifdef WIN32
#include <windows.h>
#endif

#if defined(VRN_OPENGL_COMPATIBILITY_PROFILE) && !defined(__APPLE__)
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif


#if defined __cplusplus
/* Add C++ includes here */

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "tgt/matrix.h"
#include "tgt/quaternion.h"
#include "tgt/shadermanager.h"
#include "tgt/texturemanager.h"
#include "tgt/types.h"
#include "tgt/vector.h"
#include "tgt/assert.h"

#include "voreen/core/voreenobject.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreencoreapi.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/utils/stringutils.h"

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/voxeltypeproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/eventproperty.h"

#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/genericport.h"
#include "voreen/core/ports/textport.h"

#endif
