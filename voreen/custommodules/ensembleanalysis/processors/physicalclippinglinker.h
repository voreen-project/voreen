/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_PHYSICALCLIPPINGLINKER_H
#define VRN_PHYSICALCLIPPINGLINKER_H

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/boundingboxproperty.h"

namespace voreen {

class VRN_CORE_API PhysicalClippingLinker : public Processor {
public:
    PhysicalClippingLinker();

    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Utility";                }
    virtual std::string getClassName() const  { return "PhysicalClippingLinker"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;  }

protected:
    virtual void setDescriptions() {
        setDescription("Can be used to link a clipping region in voxel space to physical space.");
    }

    virtual void process();

    /// sets the max values of the clipping regions
    virtual void inputVolumeChanged();

    /// carries over changes from voxel space region to physical space region
    virtual void clipRegionVoxelChanged();

    /// carries over changes from physical space region to voxel space region
    virtual void clipRegionPhysicalChanged();

    static const std::string loggerCat_; ///< category used in logging

private:
    VolumePort inport_;

    IntBoundingBoxProperty clipRegionVoxel_;
    FloatBoundingBoxProperty clipRegionPhysical_;
};

}   //namespace

#endif
