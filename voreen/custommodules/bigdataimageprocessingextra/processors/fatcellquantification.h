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

#ifndef VRN_FATCELLQUANTIFICATION_H
#define VRN_FATCELLQUANTIFICATION_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/volumeport.h"


#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"

#include "voreen/core/properties/filedialogproperty.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "modules/plotting/ports/plotport.h"

namespace voreen {

/**
 * Allows the quantification of the volume of a single (binary) segmentation or the quantification of two segmentation volumes (including the density of one within the other). 
 * The density is computed per slice and the output is given as plot data.
 */
class VRN_CORE_API FatCellQuantification : public Processor {

public:
    FatCellQuantification();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "FatCellQuantification";     }
    virtual std::string getCategory() const  { return "Quantification";         }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

    // progress is shown during quantification
    virtual bool usesExpensiveComputation() const { return true;    }

    virtual bool isEndProcessor() const {return true;   }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Specialized quantification of haegerling fat cells");
    }

    virtual void process();


    //BoolProperty useClipRegion_;
    //IntBoundingBoxProperty clipRegion_;

    FileDialogProperty csvSaveFile_;
    //ButtonProperty saveToCsv_;

    // quantification results
   // size_t numVoxelsTotal_;
    //size_t numVoxelsInOne_;
   // size_t numVoxelsInTwo_;
   // size_t numVoxelsInBoth_;

    /// plotting port for quantification results
    //PlotPort quantificationPlot_;

    VolumePort firstSegmentationVolume_;
    //VolumePort secondSegmentationVolume_;

    static const std::string loggerCat_;

private:

    //void useClipRegionChanged();

    void adjustToInputVolumes();
};


} // namespace

#endif
