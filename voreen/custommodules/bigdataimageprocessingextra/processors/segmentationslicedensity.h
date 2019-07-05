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

#ifndef VRN_SEGMENTATIONSLICEDENSITY_H
#define VRN_SEGMENTATIONSLICEDENSITY_H

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
class VRN_CORE_API SegmentationSliceDensity : public Processor {

public:
    SegmentationSliceDensity();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "SegmentationSliceDensity";     }
    virtual std::string getCategory() const  { return "Quantification";         }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

    // progress is shown during quantification
    virtual bool usesExpensiveComputation() const { return true;    }

    virtual bool isEndProcessor() const {return true;   }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Quantifies the density of a segmentation in each slice (i.e., relative number of foreground voxels per slice) along the z-axis and outputs plot data containing the number for each of the slices. <br> \
                <b> Caution: <\\b> The processor only supports float and unsigned integer (i.e., unsigned 8-bit, 16-bit, or 32-bit) volumes. For an unsigned integer volume, 0 is interpreted as background, and every value not equal to 0 as foreground. A float volume is interpreted as a probability volume, i.e., values < 0.5 are interpreted as background, and values >= 0.5 are interpreted as foreground.<br><br> \
                The processor only allows single-channel volumes.");
    }

    virtual void process();
    //virtual void initialize();

    virtual void computeQuantification();

    //virtual void exportToCSV();

    /// starts the quantification for all time steps
    ButtonProperty startComputation_;

    /// progress in application mode
    ProgressProperty progressProperty_;

    //BoolProperty useClipRegion_;
    //IntBoundingBoxProperty clipRegion_;

    //FileDialogProperty csvSaveFile_;
    //ButtonProperty saveToCsv_;

    // quantification results
   // size_t numVoxelsTotal_;
    //size_t numVoxelsInOne_;
   // size_t numVoxelsInTwo_;
   // size_t numVoxelsInBoth_;

    /// plotting port for quantification results
    PlotPort quantificationPlot_;

    VolumePort firstSegmentationVolume_;
    VolumePort secondSegmentationVolume_;

    static const std::string loggerCat_;

private:

    //void useClipRegionChanged();

    void adjustToInputVolumes();

    template<typename T>
    double genericQuantification(const VolumeRAM* slice, const VolumeRAM* referenceSlice) {

        const VolumeAtomic<T>* genericSlice = 0;
        const VolumeAtomic<T>* genericReferenceSlice = 0;

        // try to cast to generic type
        if (slice) {
            genericSlice = dynamic_cast<const VolumeAtomic<T>*>(slice);
            if (!genericSlice)
                return -1.0;
        }
        else 
            return -1.0;


        if (referenceSlice) {
            genericReferenceSlice = dynamic_cast<const VolumeAtomic<T>*>(referenceSlice);
            if (!genericReferenceSlice)
                return -1.0;
        }

        size_t numForegroundVoxels = 0;
        size_t numVoxels = 0;

        for (size_t y = 0; y < slice->getDimensions().y; ++y) {
            for (size_t x = 0; x < slice->getDimensions().x; ++x) {
            
                T voxel = genericSlice->voxel(x,y,0);

                if (voxel != 0) 
                    numForegroundVoxels++;


                if (genericReferenceSlice) {
                    if (genericReferenceSlice->voxel(x,y,0) != 0)
                        numVoxels++;
                }
                else
                    numVoxels++;
            }
        }

        double density = static_cast<double>(numForegroundVoxels) / static_cast<double>(numVoxels);

        return density;
    }

};


template<>
double SegmentationSliceDensity::genericQuantification<float>(const VolumeRAM* slice, const VolumeRAM* referenceSlice);

} // namespace

#endif
