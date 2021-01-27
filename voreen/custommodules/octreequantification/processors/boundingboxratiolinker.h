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

#ifndef VRN_BOUNDINGBOXRATIOLINKER_H
#define VRN_BOUNDINGBOXRATIOLINKER_H

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/boundingboxproperty.h"

namespace voreen {

class VRN_CORE_API BoundingBoxRatioLinker : public Processor {
public:
    BoundingBoxRatioLinker();

    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Utility";              }
    virtual std::string getClassName() const  { return "BoundingBoxRatioLinker";   }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

protected:
    virtual void setDescriptions() {
        setDescription("Can be used to link the clipping regions for two data sets with different dimensions (e.g. to adapt the clipping region for a downsampled version of the data set).\
                        Will synchronize the ratio of the clipping axes using the volume dimensions.");
    }

    virtual void process();

    /// sets the max values of the clipping regions 
    virtual void inputVolumesChanged();

    /// carries over changes from region A to region B
    virtual void clipRegionAChanged();

    /// carries over changes from region B to region A
    virtual void clipRegionBChanged();

    static const std::string loggerCat_; ///< category used in logging

private:
    VolumePort inportA_;
    VolumePort inportB_;
    
    IntBoundingBoxProperty clipRegionA_;
    IntBoundingBoxProperty clipRegionB_;
};

}   //namespace

#endif 
