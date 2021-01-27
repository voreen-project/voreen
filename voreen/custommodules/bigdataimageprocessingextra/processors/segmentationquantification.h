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

#ifndef VRN_SEGMENTATIONQUANTIFICATION_H
#define VRN_SEGMENTATIONQUANTIFICATION_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/volumeport.h"


#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"

#include "voreen/core/properties/filedialogproperty.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

/**
 * Allows the quantification of the volume of a single (binary) segmentation or the quantification of two segmentationi volumes (including the density of one within the other). 
 */
class VRN_CORE_API SegmentationQuantification : public Processor {

public:
    SegmentationQuantification();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "SegmentationQuantification";     }
    virtual std::string getCategory() const  { return "Quantification";         }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

    // progress is shown during quantification
    virtual bool usesExpensiveComputation() const { return true;    }

    virtual bool isEndProcessor() const {return true;   }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        // TODO: document supported volume types and how to deal with non-binary volumes
        setDescription("Quantifies the volume of one or two segmentation (i.e., binary) volumes. If two volumes are given as input, the amount of overlap is also computed. \
                The result can be exported to a CSV file. <br> \
                <b> Caution: <\\b> The processor only supports float and unsigned integer (i.e., unsigned 8-bit, 16-bit, or 32-bit) volumes. For an unsigned integer volume, 0 is interpreted as background, and every value not equal to 0 as foreground. A float volume is interpreted as a probability volume, i.e., values < 0.5 are interpreted as background, and values >= 0.5 are interpreted as foreground.<br><br> \
                For two input volumes, the dimensions and data type have to be the same. The processor only allows single-channel volumes.");
    }

    virtual void process();
    virtual void initialize();

    virtual void computeQuantification();

    virtual void exportToCSV();

    /// starts the quantification for all time steps
    ButtonProperty startComputation_;

    /// progress in application mode
    ProgressProperty progressProperty_;

    BoolProperty useClipRegion_;
    IntBoundingBoxProperty clipRegion_;

    FileDialogProperty csvSaveFile_;
    ButtonProperty saveToCsv_;

    // quantification results
    size_t numVoxelsTotal_;
    size_t numVoxelsInOne_;
    size_t numVoxelsInTwo_;
    size_t numVoxelsInBoth_;

    /// plotting port for quantification results
    //PlotPort quantificationPlot_;

    VolumePort firstSegmentationVolume_;
    VolumePort secondSegmentationVolume_;

    static const std::string loggerCat_;

private:

    void useClipRegionChanged();

    void adjustToInputVolumes();

    template<typename T>
    bool genericQuantification(const VolumeRAM* slice1, const VolumeRAM* slice2, size_t& numVoxelsTotal, size_t& numVoxelsInOne, size_t& numVoxelsInTwo, size_t& numVoxelsInBoth, tgt::svec3 llf, tgt::svec3 urb) {

        const VolumeAtomic<T>* genericSlice1 = 0;
        const VolumeAtomic<T>* genericSlice2 = 0;

        // try to cast to generic type
        if (slice1) {
            genericSlice1 = dynamic_cast<const VolumeAtomic<T>*>(slice1);
            if (!genericSlice1)
                return false;
        }

        if (slice2) {
            genericSlice2 = dynamic_cast<const VolumeAtomic<T>*>(slice2);
            if (!genericSlice2)
                return false;
        }

        for (size_t y = llf.y; y <= urb.y; ++y) {
            for (size_t x = llf.x; x <= urb.x; ++x) {
                
                T voxel1 = 0;
                T voxel2 = 0;
                if (genericSlice1)
                    voxel1 = genericSlice1->voxel(x,y,0);
                if (genericSlice2)
                    voxel2 = genericSlice2->voxel(x,y,0);

                numVoxelsTotal++;

                if (voxel1 != 0) 
                    numVoxelsInOne++;
                
                if (voxel2 != 0)
                    numVoxelsInTwo++;

                if (voxel1 != 0 && voxel2 != 0)
                    numVoxelsInBoth++;
            }
        }

        return true;
    }

};


template<>
bool SegmentationQuantification::genericQuantification<float>(const VolumeRAM* slice1, const VolumeRAM* slice2, size_t& numVoxelsTotal, size_t& numVoxelsInOne, size_t& numVoxelsInTwo, size_t& numVoxelsInBoth, tgt::svec3 llf, tgt::svec3 urb);

} // namespace

#endif
