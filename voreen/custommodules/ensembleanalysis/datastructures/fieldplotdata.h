/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_FIELDPLOTDATA_H
#define VRN_FIELDPLOTDATA_H

#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

class Volume;

class FieldPlotData : public Serializable
{
public:

    FieldPlotData(size_t width, size_t height, size_t numSlices);
    FieldPlotData(Volume* volume); // Creates a local copy! TODO: temporarily, remove
    ~FieldPlotData();

    //add connection between two points x1 and x2 with corresponding values v1 and v2
    void drawConnection(size_t x1, size_t x2, float v1, float v2, size_t sliceNumber);
    //add single value v at certain x
    void putSingleMass(size_t x, float v, size_t sliceNumber);

    /// Retrieve plot data prepared for the use of a transfer function
    Volume* getVolume() const;

    size_t getWidth()  const;
    size_t getHeight() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    //the content of the plot (one float value per pixel)
    Volume* plotData_;
    VolumeRAM_Float* representation_;
};

}

#endif
