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

#include "timeseriesextraction.h"

#include "../modules/plotting/datastructures/plotcell.h"
#include "../modules/ensembleanalysis/datastructures/ensembledataset.h"
#include <vector>
#include <map>
#include <string>


namespace voreen {
    
TimeseriesExtraction::TimeseriesExtraction()
    : Processor()
    , ensembleInport_(Port::INPORT, "ensemble_inport", "Ensemble-Data Inport")
    , samplePointInport_(Port::INPORT, "samplepoint_inport", "Samplepoint Inport")
    , outport_(Port::OUTPORT, "outport", "PlotData Outport")
    , fieldProp_("field_prop", "Selected Field")
    , useAgvSampling_("useAvgSampling", "Avg. sampling", false)
    , sphereDiameter_("sphereDiameter", "Sphere Diameter", 1.f, 0.1f, 5.f)
    , timeUnitProp_("timeUnitProp", "Time unit")
    , selTimeUnitProp_("selTimeUnitProp", "Selected time unit")
    , unitProp_("unitProp", "Unit")
    , selUnitProp_("selUnitProp", "Selected unit")
{
    // Ports
    addPort(ensembleInport_);
        ON_CHANGE(ensembleInport_, TimeseriesExtraction, ensembleChanged);
    addPort(samplePointInport_);
    addPort(outport_);

    // Properties
    addProperty(fieldProp_);
    addProperty(useAgvSampling_);
    addProperty(sphereDiameter_);
    sphereDiameter_.setStepping(0.1f);
    sphereDiameter_.setNumDecimals(1);
    addProperty(timeUnitProp_);
    addProperty(selTimeUnitProp_);
    selTimeUnitProp_.setEditable(false);
    timeUnitProp_.addOption("seconds", "Seconds", 1);
    timeUnitProp_.addOption("minutes", "Minutes", 60);
    timeUnitProp_.setDefaultValue("seconds");

    addProperty(unitProp_);
    addProperty(selUnitProp_);
    selUnitProp_.setEditable(false);
    unitProp_.addOption("cpscc", "CpsCC", "CpsCC");
    unitProp_.addOption("idml", "%ID/mL", "%ID/mL");
}

Processor* TimeseriesExtraction::create() const {
    return new TimeseriesExtraction();
}

std::string TimeseriesExtraction::getClassName() const {
    return "TimeseriesExtraction";
}

std::string TimeseriesExtraction::getCategory() const {
    return "Ensemble Processing";
}

void TimeseriesExtraction::setDescriptions() {
    setDescription("Processor which extracts a timeseries from ensemble-data" 
                   " at a given point.");
}

void TimeseriesExtraction::process() {

    // Update selected time-unit
    selTimeUnitProp_.set(timeUnitProp_.getDescription());
    // Update selected unit
    selUnitProp_.set(unitProp_.getDescription());

    // No Sample-Points --> Clear plot and return
    const EnsembleSamplePointMapping* samplePointMapping = samplePointInport_.getData();
    if (samplePointMapping == nullptr) {
        outport_.setData(0, true);
        return;
    }

    // No Ensemble-Members or no field selected --> Clear plot and return
    const EnsembleDataset* dataset = ensembleInport_.getData();
    if (dataset == 0 || dataset->getMembers().size() == 0 
                     || fieldProp_.getSelectedIndex() == -1) { 
        outport_.setData(0, true);
        return;
    }

    const Option<std::string>& selField = fieldProp_.getOptions()[fieldProp_.getSelectedIndex()];

    // Vector containing the column-names 
    std::vector<std::string> col_names;
    col_names.push_back("Time");

    // Map with key = timestep, value: vector with samples for this timestep (order same as col_names!)
    std::map<float, std::vector<PlotCellValue>> rowMap;

    // Add a vector for every timestep which already contains the value of the key-cloumn (time)
    // TODO: What if not all members share the same timesteps?
    float max_time = -1;
    for (const TimeStep& ts : dataset->getMembers()[0].getTimeSteps()) {
        std::vector<PlotCellValue> row;
        if (ts.getTime() > max_time) {
            max_time = ts.getTime();
        }
        const float timeConverted = convertTime(ts.getTime());
        row.push_back(PlotCellValue(timeConverted));
        rowMap.insert(std::pair<float, std::vector<PlotCellValue>>(ts.getTime(), row));
    }

    // Iterate over all members
    for (EnsembleMember member : dataset->getMembers()) {
        const std::string memberName = member.getName();
        const POIList* points = samplePointMapping->getPOIList(memberName);

        // No Samplepoints --> skip
        if (points->getPoints().size() == 0) {
            continue;
        }

        // Total counts in the volume (used to calculate %ID/mL)
        // Use last timestep
        double total_counts = -1;
        if (unitProp_.get() == "idml" && max_time != -1) {
            const size_t idx = member.getTimeStep(max_time);
            const VolumeBase* volumeBase = member.getTimeSteps()[idx].getVolume(selField.key_);
            const VolumeRAMRepresentationLock volumeData(volumeBase);
            const RealWorldMapping& rwm = volumeBase->getRealWorldMapping();
            // Use the last frame to calculate total number of counts
            total_counts = getTotalCount(volumeData, rwm, volumeBase->getSpacing());
        }

        // Iterate over all timesteps
        for (int i = 0; i < member.getTimeSteps().size(); ++i) {
            const TimeStep& timestep = member.getTimeSteps()[i];

            const VolumeBase* volumeBase = timestep.getVolume(selField.key_);
            VolumeRAMRepresentationLock volumeData(volumeBase);
            RealWorldMapping rwm = volumeBase->getRealWorldMapping();

            // Iterate over all sample-points
            for (POIPoint sp : points->getPoints()) {
                // Convert from World- to Voxel-Coordinates
                tgt::vec3 posWorld = sp.position_;
                tgt::vec3 posVoxel = volumeBase->getWorldToVoxelMatrix() * posWorld;

                // Add column label (only once for every Member and samplepoint)
                if (i == 0) {
                    POIGroupID gID = sp.group_;
                    std::string gName = points->getGroupName(gID);
                    POIPointID pID = sp.id_;
                    col_names.push_back(gName + "::" + std::to_string(pID) + "::" + memberName);
                }
                
                float valueRw;
                if (useAgvSampling_.get()) {
                    // Use average sampling
                    const tgt::vec3 spacing = volumeBase->getSpacing();
                    const float avg = getAvgSample(posVoxel, volumeData, spacing);
                    valueRw = rwm.normalizedToRealWorld(avg);
                }
                else {
                    // Use normal sampling
                    const float valueNorm = volumeData->getVoxelNormalizedLinear(posVoxel);
                    valueRw = rwm.normalizedToRealWorld(valueNorm);
                }

                // Unit conversion ( CpsCC requires no conversion):
                if (unitProp_.get() == "idml") { 
                    const float cpsCC_roi = valueRw;
                    valueRw = (cpsCC_roi / total_counts) * 100.0;
                }

                rowMap[timestep.getTime()].push_back(valueRw);

            }
        }
    }

    // Setup plot (One column for each sample-point + Key-column)
    const int numCols = col_names.size();
    PlotData* newData = new PlotData(1, numCols - 1);

    // Set column-names
    for (int i = 0; i < numCols; ++i) {
        newData->setColumnLabel(i, col_names[i]);
    }

    // Add Rows
    for (auto const& entry : rowMap) {
        std::vector<PlotCellValue> row = entry.second;
        newData->insert(row);
    }
    
    outport_.setData(newData, true);
    
}

void TimeseriesExtraction::ensembleChanged() {
    // Clear Propertys
    fieldProp_.reset();

    if (!ensembleInport_.isReady())
        return;

    const EnsembleDataset* dataset = ensembleInport_.getData();

    // Populate field-property
    fieldProp_.blockCallbacks(true);
    for (const std::string& fieldName : dataset->getCommonFieldNames()){
        if(!fieldProp_.hasKey(fieldName))
            fieldProp_.addOption(fieldName, fieldName, fieldName);
    } 
    fieldProp_.blockCallbacks(false);    
}

float TimeseriesExtraction::getAvgSample(const tgt::vec3& sp, const VolumeRAMRepresentationLock& volumeData, const tgt::vec3& spacing) {
    float radius = sphereDiameter_.get() / 2.0f;

    const float delta = 0.05f;
    
    // Construct lower-front-left-corner and upper-back-right-corner of a Cube with 
    // side length 'diameter' which is centered around sp.
    const tgt::svec3 lfl_corner(static_cast<size_t>(floor(sp.x - radius)), 
                                static_cast<size_t>(floor(sp.y - radius)), 
                                static_cast<size_t>(floor(sp.z - radius)));
    const tgt::vec3 ubr_corner(static_cast<size_t>(ceil(sp.x + radius)),
                               static_cast<size_t>(ceil(sp.y + radius)),
                               static_cast<size_t>(ceil(sp.z + radius)));

    float sum = 0.0f;
    size_t voxel_count = 0;

    // Check for all voxels inside the imaginary cube, wether they are inside an imaginary sphere
    // centered around sp with a diameter of 'diameter'
    for (size_t x = lfl_corner.x; x <= ubr_corner.x; ++x) {
        for (size_t y = lfl_corner.y; y <= ubr_corner.y; ++y) {
            for (size_t z = lfl_corner.z; z <= ubr_corner.z; ++z) {
                // Position of current voxel
                tgt::vec3 vox(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
                // Adjust for voxel-spacing
                tgt::vec3 diff = (vox - sp) * spacing;
                // Distance from voxel to sp
                const float dist = tgt::length(diff);
                if (dist <= radius) {
                    // Voxel is inside sphere
                    sum += volumeData->getVoxelNormalized(vox);
                    ++voxel_count;
                }
            }
        }
    }

    if (voxel_count == 0) {
        // No voxels inside sphere
        LWARNING("Average sampling used 0 voxels!!! Returning value -1 instead. Try to increase the diameter.");
        return 0;
    }
    else {
        float avg = sum / voxel_count;
        return avg;
    }
}

float voreen::TimeseriesExtraction::convertTime(float timeInSeconds) {
    const float div = timeUnitProp_.getValue();
    return timeInSeconds / div;
}

const float TimeseriesExtraction::getTotalCount(const VolumeRAMRepresentationLock& volumeData, const RealWorldMapping& rwm,
    const tgt::vec3& spacing) const {
    // Get total number of counts from PET-Volume:
    float total = 0;
    const tgt::svec3 dims = volumeData->getDimensions();
    for (size_t x = 0; x < dims.x; ++x) {
        for (size_t y = 0; y < dims.y; ++y) {
            for (size_t z = 0; z < dims.z; ++z) {
                const float valNorm = volumeData->getVoxelNormalized(tgt::svec3(x, y, z));
                const float valRw = rwm.normalizedToRealWorld(valNorm);
                total += valRw;
            }
        }
    }

    // Multiply counts (CpsCC = Bq/mL) by voxel-volume (mL) to get Cps
    const float voxel_volume = (spacing.x * spacing.y * spacing.z) / 1000.0f; // Spacing: mm
    total = total * voxel_volume;

    return total;
}
} // namespace
