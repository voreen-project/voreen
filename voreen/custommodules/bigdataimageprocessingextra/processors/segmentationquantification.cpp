/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "segmentationquantification.h"

#include "voreen/core/datastructures/volume/volumedisk.h"

#include "tgt/stopwatch.h"

namespace voreen {

const std::string SegmentationQuantification::loggerCat_("voreen.SegmentationQuantification");

SegmentationQuantification::SegmentationQuantification()
    : Processor()
    , firstSegmentationVolume_(Port::INPORT, "firstsegmentation", "First Segmentation Volume", false)
    , secondSegmentationVolume_(Port::INPORT, "secondsegmentation", "Second Segmentation Volume", false)
    , useClipRegion_("useClipRegion", "Use Clip Region", false)
    , clipRegion_("clipRegion", "Clip Region", tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(1)), tgt::ivec3(0), tgt::ivec3(1))
    , startComputation_("startComputation", "Start Computation")
    , progressProperty_("progressProperty", "Quantification Progress")
    , csvSaveFile_("csvFileProp", "CSV Export Path", "CSV Export Path", ".", "Comma seperated values (*.csv)", FileDialogProperty::SAVE_FILE)
    , saveToCsv_("savetocsv", "Save to CSV")
    , numVoxelsTotal_(0)
    , numVoxelsInOne_(0)
    , numVoxelsInTwo_(0)
    , numVoxelsInBoth_(0)
    //, channel_("channel", "Quantification Channel", 0, 0, 0, Processor::INVALID_RESULT, NumericProperty<int>::DYNAMIC, Property::LOD_ADVANCED)
    //, quantificationPlot_(Port::OUTPORT, "QuantificationPlot", "Quantification Score Plot")
{
    addPort(firstSegmentationVolume_);
    addPort(secondSegmentationVolume_);
    //addPort(quantificationPlot_);

    addProperty(useClipRegion_);
    clipRegion_.setVisibleFlag(false);
    addProperty(clipRegion_);
    addProperty(startComputation_);

    addProperty(progressProperty_);
    addProgressBar(&progressProperty_);

    addProperty(csvSaveFile_);
    addProperty(saveToCsv_);
    ON_CHANGE(saveToCsv_, SegmentationQuantification, exportToCSV);
    saveToCsv_.setReadOnlyFlag(true);

    ON_CHANGE(firstSegmentationVolume_, SegmentationQuantification, adjustToInputVolumes);
    ON_CHANGE(secondSegmentationVolume_, SegmentationQuantification, adjustToInputVolumes);

    ON_CHANGE(useClipRegion_, SegmentationQuantification, useClipRegionChanged);

    ON_CHANGE(startComputation_, SegmentationQuantification, computeQuantification);
}

Processor* SegmentationQuantification::create() const {
    return new SegmentationQuantification();
}

bool SegmentationQuantification::isReady() const {
    return isInitialized() && (firstSegmentationVolume_.isReady() || secondSegmentationVolume_.isReady());
}

void SegmentationQuantification::process() {
    // nothing

}

void SegmentationQuantification::initialize() {
    Processor::initialize();

    adjustToInputVolumes();
}

void SegmentationQuantification::computeQuantification() {

    setProgress(0.f);
    //quantificationPlot_.clear();
    
    // get volumes
    const VolumeBase* volume1 = firstSegmentationVolume_.getData();
    const VolumeBase* volume2 = secondSegmentationVolume_.getData();
    std::string format = "unknown";

    // check number of channels, dimensions, and data type of the volumes
    if (volume1 && volume1->getNumChannels() != 1) {
        LERROR("First volume has more than one channel!");
        return;
    }
    else if (volume1) {
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


    if (volume2 && volume2->getNumChannels() != 1) {
        LERROR("Second volume has more than one channel!");
        return;
    }
    else if (volume2) {
        // one channel -> check if data type is supported
        if (!(volume2->getFormat() == "uint8" 
                || volume2->getFormat() == "uint16" 
                || volume2->getFormat() == "uint32" 
                || volume2->getFormat() == "float")) {
            LERROR("Unsupported format: only uint8, uint16, uint32, or float are supported");
            return;
        }
        else
            format = volume2->getFormat();
    }

    if (volume1 && volume2) {
        if (!(volume1->getDimensions() == volume2->getDimensions())) {
            LERROR("Volumes are not of the same size!");
            return;
        }

        if (!(volume1->getBaseType() == volume2->getBaseType() && volume1->getFormat() == volume2->getFormat())) {
            LERROR("Volumes do not have the same data type!");
            return;
        }

        if (!(volume1->getSpacing() == volume2->getSpacing())) {
            LERROR("Volumes do not have the same spacing!");
            return;
        }
    }

    // result values
    numVoxelsTotal_ = 0;
    numVoxelsInOne_ = 0;
    numVoxelsInTwo_ = 0;
    numVoxelsInBoth_ = 0;

    tgt::Stopwatch timer;
    timer.start();

    // get the quantification dimensions
    tgt::svec3 llf, urb;
    llf = tgt::svec3::zero;
    if (volume1) 
        urb = volume1->getDimensions() - tgt::svec3::one;
    
    if (volume2) 
        urb = volume2->getDimensions() - tgt::svec3::one;

    // if the clip region is used, crop our bounds
    if (useClipRegion_.get()) {
        llf = tgt::max(llf, tgt::svec3(clipRegion_.get().getLLF()));
        urb = tgt::min(urb, tgt::svec3(clipRegion_.get().getURB()));
    }
    
    // we do not know how large a single slice is, so we only load one slice at a time for each volume
    for (size_t z = llf.z; z <= urb.z; ++z) {

        VolumeRAM* slice1 = 0;
        VolumeRAM* slice2 = 0;

        if (volume1) {
            slice1 = volume1->getSlice(z);
            if (!slice1) {
                LERROR("Could not create slice " << z << " for first volume!");
                return;
            }
        }

        if (volume2) {
            slice2 = volume2->getSlice(z);
            if (!slice2) {
                LERROR("Could not create slice " << z << " for second volume!");
                delete slice1;
                return;
            }
        }

        bool sliceSuccess = true;

        // quantify slice according to the data type
        if (format == "uint8") {
            sliceSuccess = genericQuantification<uint8_t>(slice1, slice2, numVoxelsTotal_, numVoxelsInOne_, numVoxelsInTwo_, numVoxelsInBoth_, llf, urb);
        }
        else if (format == "uint16") {
            sliceSuccess = genericQuantification<uint16_t>(slice1, slice2, numVoxelsTotal_, numVoxelsInOne_, numVoxelsInTwo_, numVoxelsInBoth_, llf, urb);
        }
        else if (format == "uint32") {
            sliceSuccess = genericQuantification<uint32_t>(slice1, slice2, numVoxelsTotal_, numVoxelsInOne_, numVoxelsInTwo_, numVoxelsInBoth_, llf, urb);
        }
        else if (format == "float") {
            sliceSuccess = genericQuantification<float>(slice1, slice2, numVoxelsTotal_, numVoxelsInOne_, numVoxelsInTwo_, numVoxelsInBoth_, llf, urb);
        }

        delete slice1;
        delete slice2;

        if (!sliceSuccess) {
            LERROR("Error: could not quantify slice " << z);
            setProgress(0.f);
            timer.stop();
            return;
        }

        setProgress(std::min(0.99f, static_cast<float>(z - llf.z) / static_cast<float>(urb.z - llf.z)));
    }

    // TODO: compute volume with spacing, create plots     

    setProgress(1.f);

    timer.stop();
    LINFO("Quantification Time: " << (float) timer.getRuntime() / 1000.f << " seconds");

    saveToCsv_.setReadOnlyFlag(false);
    
    /*
    PlotData* plotData = new PlotData(0, 4);
    plotData->setColumnLabel(0, "Time Step");
    plotData->setColumnLabel(1, "Raw sum of voxels");
    plotData->setColumnLabel(2, "Weighted absolute sum");
    plotData->setColumnLabel(3, "Weighted normalized sum (quantification score)");

    PlotData* normalizedPlotData = new PlotData(0,2);
    normalizedPlotData->setColumnLabel(0, "Time Step");
    normalizedPlotData->setColumnLabel(1, "Weighted normalized sum (quantification score)");

    for (size_t i = static_cast<size_t>(quantificationFrame_.get().x); i <= static_cast<size_t>(quantificationFrame_.get().y); ++i) {
        std::vector<PlotCellValue> v(4);
        v.at(0) = PlotCellValue(static_cast<plot_t>(i));
        v.at(1) = PlotCellValue(static_cast<plot_t>(sumResults.at(i)));
        v.at(2) = PlotCellValue(static_cast<plot_t>(absoluteResults.at(i)));
        v.at(3) = PlotCellValue(static_cast<plot_t>(normalizedResults.at(i)));

        std::vector<PlotCellValue> v2(2);
        v2.at(0) = PlotCellValue(static_cast<plot_t>(i));
        v2.at(1) = PlotCellValue(static_cast<plot_t>(normalizedResults.at(i)));

        bool inserted1 = plotData->insert(v);
        bool inserted2 = normalizedPlotData->insert(v2);
        if (!inserted1 || !inserted2) {
            LERROR("Could not insert data into plot");
            delete plotData;
            delete normalizedPlotData;
            plotData = 0;
            normalizedPlotData = 0;
            break;
        }
    }

    quantificationPlot_.setData(plotData, true);
    normalizedQuantificationPlot_.setData(normalizedPlotData, true);*/
}

void SegmentationQuantification::adjustToInputVolumes() {
    if (!isInitialized())
        return;

    // new input volumes -> reset results
    numVoxelsTotal_ = 0;
    numVoxelsInOne_ = 0;
    numVoxelsInTwo_ = 0;
    numVoxelsInBoth_ = 0;   
    saveToCsv_.setReadOnlyFlag(true);

    // adjust clipping area to input
    const VolumeBase* firstVol = firstSegmentationVolume_.getData();
    const VolumeBase* secondVol = secondSegmentationVolume_.getData();
    if (firstVol && secondVol) {
        // check if dimensions are the same and adjust clipping area
        if (firstVol->getDimensions() == secondVol->getDimensions())
            clipRegion_.setMaxValue(tgt::ivec3(firstVol->getDimensions() - tgt::svec3::one));
    }
    else if (firstVol)
        clipRegion_.setMaxValue(tgt::ivec3(firstVol->getDimensions() - tgt::svec3::one));
    else if (secondVol)
        clipRegion_.setMaxValue(tgt::ivec3(secondVol->getDimensions() - tgt::svec3::one));

}

void SegmentationQuantification::useClipRegionChanged() {
    clipRegion_.setVisibleFlag(useClipRegion_.get());
}

void SegmentationQuantification::exportToCSV() {
    std::ofstream file(csvSaveFile_.get());
    file << "Number_of_voxels,voxels_in_volume1,voxels_in_volume2,voxels_in_both" << std::endl;
    file << numVoxelsTotal_ << "," << numVoxelsInOne_ << "," << numVoxelsInTwo_ << "," << numVoxelsInBoth_ << std::endl;

    file.close();
}

template<>
bool SegmentationQuantification::genericQuantification<float>(const VolumeRAM* slice1, const VolumeRAM* slice2, size_t& numVoxelsTotal, size_t& numVoxelsInOne, size_t& numVoxelsInTwo, size_t& numVoxelsInBoth, tgt::svec3 llf, tgt::svec3 urb) {

    const VolumeAtomic<float>* genericSlice1 = 0;
    const VolumeAtomic<float>* genericSlice2 = 0;

    // try to cast to generic type
    if (slice1) {
        genericSlice1 = dynamic_cast<const VolumeAtomic<float>*>(slice1);
        if (!genericSlice1)
            return false;
    }

    if (slice2) {
        genericSlice2 = dynamic_cast<const VolumeAtomic<float>*>(slice2);
        if (!genericSlice2)
            return false;
    }

    for (size_t y = llf.y; y <= urb.y; ++y) {
        for (size_t x = llf.x; x <= urb.x; ++x) {
            
            float voxel1 = 0;
            float voxel2 = 0;
            if (genericSlice1)
                voxel1 = genericSlice1->voxel(x,y,0);
            if (genericSlice2)
                voxel2 = genericSlice2->voxel(x,y,0);

            numVoxelsTotal++;

            if (voxel1 >= 0.5f) 
                numVoxelsInOne++;
            
            if (voxel2 >= 0.f)
                numVoxelsInTwo++;

            if (voxel1 >= 0.5f && voxel2 >= 0.5f)
                numVoxelsInBoth++;
        }
    }

    return true;
}

} // namespace
