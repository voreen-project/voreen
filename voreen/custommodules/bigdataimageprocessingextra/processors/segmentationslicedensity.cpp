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

#include "segmentationslicedensity.h"

#include "voreen/core/datastructures/volume/volumedisk.h"

#include "tgt/stopwatch.h"

#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"
#include "modules/plotting/datastructures/plotcell.h"

namespace voreen {

const std::string SegmentationSliceDensity::loggerCat_("voreen.SegmentationSliceDensity");

SegmentationSliceDensity::SegmentationSliceDensity()
    : Processor()
    , firstSegmentationVolume_(Port::INPORT, "firstsegmentation", "First Segmentation Volume", false)
    , secondSegmentationVolume_(Port::INPORT, "secondsegmentation", "Second (Reference) Segmentation Volume", false)
    //, useClipRegion_("useClipRegion", "Use Clip Region", false)
    //, clipRegion_("clipRegion", "Clip Region", tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(1)), tgt::ivec3(0), tgt::ivec3(1))
    , startComputation_("startComputation", "Start Computation")
    , progressProperty_("progressProperty", "Quantification Progress")
    //, csvSaveFile_("csvFileProp", "CSV Export Path", "CSV Export Path", ".", "Comma seperated values (*.csv)", FileDialogProperty::SAVE_FILE)
    //, saveToCsv_("savetocsv", "Save to CSV")
    , quantificationPlot_(Port::OUTPORT, "QuantificationPlot", "Quantification Plot", true,  Processor::VALID)
{
    addPort(firstSegmentationVolume_);
    addPort(secondSegmentationVolume_);
    
    addPort(quantificationPlot_);

    addProperty(startComputation_);

    addProperty(progressProperty_);
    addProgressBar(&progressProperty_);

    ON_CHANGE(firstSegmentationVolume_, SegmentationSliceDensity, adjustToInputVolumes);

    ON_CHANGE(startComputation_, SegmentationSliceDensity, computeQuantification);
}

Processor* SegmentationSliceDensity::create() const {
    return new SegmentationSliceDensity();
}

bool SegmentationSliceDensity::isReady() const {
    return isInitialized() && (firstSegmentationVolume_.isReady());
}

void SegmentationSliceDensity::process() {
    // nothing

}

void SegmentationSliceDensity::computeQuantification() {

    setProgress(0.f);
    quantificationPlot_.clear();
    
    // get volume
    const VolumeBase* volume1 = firstSegmentationVolume_.getData();
    if (!volume1) {
        LERROR("No input volume!");
        return;
    }

    std::string format = "unknown";

    // check number of channels, dimensions, and data type of the volume
    if (volume1->getNumChannels() != 1) {
        LERROR("volume has more than one channel!");
        return;
    }
    else {
        // one channel -> check if data type is supported
        if (!(volume1->getFormat() == "uint8" 
                || volume1->getFormat() == "uint16" 
                || volume1->getFormat() == "uint32" 
                || volume1->getFormat() == "float")) {
            LERROR("Unsupported format: only uint8, uint16, uint32, or float are supported");
            return;
        }
        else
            format = volume1->getFormat();
    }

    const VolumeBase* volume2 = 0;
    if (secondSegmentationVolume_.isReady() && secondSegmentationVolume_.getData()) {
        volume2 = secondSegmentationVolume_.getData();
        if (!volume2) {
            LERROR("No input volume for volume 2!");
            return;
        }
        if (volume2->getNumChannels() != 1) {
            LERROR("volume 2 has multiple channels!");
            return;
        }
        if (volume2->getFormat() != volume1->getFormat()) {
            LERROR("Volume 1 and Volume 2 have different formats!");
            return;
        }
        if (volume1->getSpacing() != volume2->getSpacing()) {
            LERROR("Spacing differs between data sets!");
            return;
        }
        if (volume1->getDimensions() != volume2->getDimensions()) {
            LERROR("Volume dimensions are not the same in volume 1 and 2");
            return;
        }
    }

    tgt::Stopwatch timer;
    timer.start();
    
    // result vector
    std::vector<double> results(volume1->getDimensions().z);

    // we do not know how large a single slice is, so we only load one slice at a time for each volume
    for (size_t z = 0; z < volume1->getDimensions().z; ++z) {

        VolumeRAM* slice1 = volume1->getSlice(z);
        if (!slice1) {
            LERROR("Could not create slice " << z << " for first volume!");
            timer.stop();
            return;
        }

        // if needed: get slice of second volume
        VolumeRAM* slice2 = 0;
        if (volume2) {
            slice2 = volume2->getSlice(z);
            if (!slice2) {
                LERROR("Could not create slice " << z << " for second volume!");
                timer.stop();
                return;
            }
        }

        double result = 0.0;

        // quantify slice according to the data type
        if (format == "uint8") {
            result = genericQuantification<uint8_t>(slice1, slice2);
        }
        else if (format == "uint16") {
            result = genericQuantification<uint16_t>(slice1, slice2);
        }
        else if (format == "uint32") {
            result = genericQuantification<uint32_t>(slice1, slice2);
        }
        else if (format == "float") {
            result = genericQuantification<float>(slice1, slice2);
        }

        delete slice1;
        delete slice2;

        if (result < 0.0) {
            LERROR("Error: could not quantify slice " << z);
            setProgress(0.f);
            timer.stop();
            return;
        }
        else
            results.at(z) = result;

        setProgress(std::min(0.99f, static_cast<float>(z) / static_cast<float>(volume1->getDimensions().z)));
    }

    timer.stop();
    LINFO("Quantification Time: " << (float) timer.getRuntime() / 1000.f << " seconds");

    PlotData* plotData = new PlotData(0, 2);
    plotData->setColumnLabel(0, "Slice Number");
    plotData->setColumnLabel(1, "Foreground Density");

    for (size_t i = 0; i < results.size(); ++i) {
        std::vector<PlotCellValue> v(2);
        v.at(0) = PlotCellValue(static_cast<plot_t>(i));
        v.at(1) = PlotCellValue(static_cast<plot_t>(results.at(i)));

        bool inserted1 = plotData->insert(v);
        if (!inserted1) {
            LERROR("Could not insert data into plot");
            delete plotData;
            plotData = 0;
            break;
        }
    }

    quantificationPlot_.setData(plotData, true);

    setProgress(1.f);
}

void SegmentationSliceDensity::adjustToInputVolumes() {
    if (!isInitialized())
        return;

    // clear the output if a new volume is present
    quantificationPlot_.clear();
}

/*void SegmentationQuantification::exportToCSV() {
    std::ofstream file(csvSaveFile_.get());
    file << "Number_of_voxels,voxels_in_volume1,voxels_in_volume2,voxels_in_both" << std::endl;
    file << numVoxelsTotal_ << "," << numVoxelsInOne_ << "," << numVoxelsInTwo_ << "," << numVoxelsInBoth_ << std::endl;

    file.close();
}*/

template<>
double SegmentationSliceDensity::genericQuantification<float>(const VolumeRAM* slice, const VolumeRAM* referenceSlice) {

    const VolumeAtomic<float>* genericSlice = 0;
    const VolumeAtomic<float>* genericReferenceSlice = 0;

    // try to cast to generic type
    if (slice) {
        genericSlice = dynamic_cast<const VolumeAtomic<float>*>(slice);
        if (!genericSlice)
            return -1.0;
    }
    else
        return -1.0;

    if (referenceSlice) {
        genericReferenceSlice = dynamic_cast<const VolumeAtomic<float>*>(referenceSlice);
        if (!genericReferenceSlice)
            return -1.0;
    }

    size_t numForegroundVoxels = 0;
    size_t numVoxels = 0;

    for (size_t y = 0; y < slice->getDimensions().y; ++y) {
        for (size_t x = 0; x < slice->getDimensions().x; ++x) {
            
            float voxel = genericSlice->voxel(x,y,0);
 
            if (voxel >= 0.5f) 
                numForegroundVoxels++;           

            if (genericReferenceSlice) {
                if (genericReferenceSlice->voxel(x,y,0) >= 0.5f)
                    numVoxels++;
            }
            else
                numVoxels++;
        }
    }

    double density = static_cast<double>(numForegroundVoxels) / static_cast<double>(numVoxels);

    return density;
}

} // namespace
