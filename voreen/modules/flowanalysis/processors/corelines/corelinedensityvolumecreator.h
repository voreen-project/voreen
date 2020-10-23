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

#ifndef VRN_CORELINEDENSITYVOLUMECREATOR_H
#define VRN_CORELINEDENSITYVOLUMECREATOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

class CorelineDensityVolumeCreator : public Processor {
public:
    CorelineDensityVolumeCreator();
    virtual Processor *create() const { return new CorelineDensityVolumeCreator(); }
    virtual std::string getClassName() const { return "CorelineDensityVolumeCreator"; }
    virtual std::string getCategory() const { return "Volume Processing"; }

    static void Process(const std::vector<std::vector<tgt::vec3>>& corelines, VolumeRAM_Float& outBinaryVolume);

protected:
    virtual void process();

private:
    GeometryPort _inCorelines;
    VolumePort _inVolume;
    VolumePort _out;
};

} // namespace voreen

#endif