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

#ifndef VRN_NWESTYLE_MATERIAL_DARK_H
#define VRN_NWESTYLE_MATERIAL_DARK_H

#include "nwestyle_material.h"

namespace voreen {

/**
 * Material-like design for Voreen.
 */
class NWEStyle_Material_Dark : public NWEStyle_Material {
public:
    NWEStyle_Material_Dark(NetworkEditor* networkeditor);
    ~NWEStyle_Material_Dark();

    /*********************************************************************
     *                       General Color Defines
     ********************************************************************/
    virtual QBrush getBackgroundBrush() const;
    virtual QColor getSelectionColor() const;
    virtual QColor getHoverColor() const;
    virtual QColor getProcessorColor1() const;
    virtual QColor getPortArrowColor() const;
    virtual QColor getPropertyLinkArrowColor() const;
    virtual QColor getPortSizeLinkArrowColor() const;

    virtual void ProcessorGI_initializePaintSettings(ProcessorGraphicsItem* item);
};

} //namespace voreen

#endif // VRN_NWESTYLE_CLASSIC_PRINT_H

