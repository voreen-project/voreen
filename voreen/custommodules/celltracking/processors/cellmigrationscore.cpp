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

#include "cellmigrationscore.h"

#include "voreen/core/datastructures/volume/volumedisk.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "tgt/stopwatch.h"

namespace voreen {

const std::string CellMigrationScore::loggerCat_("voreen.celltracking.CellMigrationScore");

CellMigrationScore::CellMigrationScore()
    : Processor()
    , inport_(Port::INPORT, "volumecollection", "VolumeList Input", false)
    , quantificationFrame_("quantificationFrame", "Quantification Time Step Interval")
    , quantificationCenter_("quantificationCenter", "Quantification Center (Voxel Coordinates)", tgt::ivec3::one, tgt::ivec3::zero, tgt::ivec3::one, Processor::INVALID_RESULT, NumericProperty<tgt::ivec3>::DYNAMIC) 
    , radius_("quantificationRadius", "Quantification Radius (in um)", 1, 1, 1, Processor::INVALID_RESULT, NumericProperty<int>::DYNAMIC)
    , useClipRegion_("useClipRegion", "Use Clip Region", false)
    , clipRegion_("clipRegion", "Clip Region", tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(1)), tgt::ivec3(0), tgt::ivec3(1))
    , startComputation_("startComputation", "Start Computation")
    , progressProperty_("progressProperty", "Quantification Progress")
    , channel_("channel", "Quantification Channel", 0, 0, 0, Processor::INVALID_RESULT, NumericProperty<int>::DYNAMIC, Property::LOD_ADVANCED)
    , cylinderSubdivisions_("cylindersubdivisions", "Cylinder Subdivisions in Rendering", 20, 3, 200, 
            Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_DEBUG)
    , centerOutport_(Port::OUTPORT, "geometry.pointlist", "Quantification Center PointList Output")
    , areaOutport_(Port::OUTPORT, "geometry.roi", "Quantification Region Outport")
    , quantificationPlot_(Port::OUTPORT, "QuantificationPlot", "Quantification Score Plot")
    , normalizedQuantificationPlot_(Port::OUTPORT, "NormalizedQuantificationPlot", "Normalized Weighted Score Plot")
{
    addPort(inport_);
    addPort(centerOutport_);
    //addPort(lineOutport_);
    addPort(areaOutport_);
    addPort(quantificationPlot_);
    addPort(normalizedQuantificationPlot_);

    addProperty(quantificationFrame_);
    addProperty(quantificationCenter_);
    addProperty(radius_);
    addProperty(useClipRegion_);
    clipRegion_.setVisibleFlag(false);
    addProperty(clipRegion_);
    addProperty(startComputation_);

    addProperty(progressProperty_);
    addProgressBar(&progressProperty_);

    addProperty(channel_);
    addProperty(cylinderSubdivisions_);

    ON_CHANGE(useClipRegion_, CellMigrationScore, useClipRegionChanged);

    ON_CHANGE(startComputation_, CellMigrationScore, computeQuantification);
}

Processor* CellMigrationScore::create() const {
    return new CellMigrationScore();
}

bool CellMigrationScore::isReady() const {
    return isInitialized() && inport_.isReady();
}

void CellMigrationScore::process() {
    // nothing

    centerOutport_.clear();
    areaOutport_.clear();

    if (inport_.getData() && !inport_.getData()->empty()) {

        if (!checkVolumeList(inport_.getData())) {
            LERROR("No valid time series");
            return;
        }
    
        // 1. output quantification center

        // since we checked that, we can assume that every time step has the same transformation matrix...
        tgt::mat4 m = inport_.getData()->first()->getVoxelToWorldMatrix();
        PointListGeometryVec3* pointList = new PointListGeometryVec3();
        pointList->addPoint( (m * tgt::vec4(tgt::vec3(quantificationCenter_.get()), 1.f)).xyz() );
        centerOutport_.setData(pointList, true);

        // 2. output quantification region
        GlMeshGeometryUInt16Normal* mesh = new GlMeshGeometryUInt16Normal();
        //mesh->setSphereGeometry(static_cast<float>(radius_.get()) / 1000.f, (m * tgt::vec4(tgt::vec3(quantificationCenter_.get()), 1.f)).xyz(), tgt::vec4(1.f, 0.f,0.f,1.f), 10);

        // get the start and the end point of the cylinder, which are the begin and end of the volume in z-direction with the xy-coordinate of the quantification center
        tgt::vec3 start = (m * tgt::vec4(quantificationCenter_.get().x, quantificationCenter_.get().y, -0.5f, 1.f)).xyz();
        tgt::vec3 end = (m * tgt::vec4(quantificationCenter_.get().x, quantificationCenter_.get().y, static_cast<float>(inport_.getData()->first()->getDimensions().z - 1) + 0.5f, 1.f)).xyz(); 

        //calculate correct rotation matrix:
        tgt::vec3 rotz = normalize(end - start);
        tgt::vec3 roty = normalize(tgt::vec3(rotz.y, -rotz.z, 0.0f));
        tgt::vec3 rotx = tgt::cross(roty, rotz);

        float mm[16];

        mm[0] = rotx.x;
        mm[1] = rotx.y;
        mm[2] = rotx.z;
        mm[3] = 0.0f;

        mm[4] = roty.x;
        mm[5] = roty.y;
        mm[6] = roty.z;
        mm[7] = 0.0f;

        mm[8] = rotz.x;
        mm[9] = rotz.y;
        mm[10] = rotz.z;
        mm[11] = 0.0f;

        mm[12] = 0.0f;
        mm[13] = 0.0f;
        mm[14] = 0.0f;
        mm[15] = 1.0f;

        tgt::mat4 matrix(mm);
        matrix = tgt::mat4::createTranslation(start) * tgt::transpose(matrix);

        float l = length(start - end);

        mesh->setCylinderGeometry(tgt::vec4(1.f, 0.f,0.f,1.f), static_cast<float>(radius_.get()) / 1000.f, static_cast<float>(radius_.get()) / 1000.f, l, cylinderSubdivisions_.get(), 1, false, false);
        mesh->setTransformationMatrix(matrix);

        areaOutport_.setData(mesh, true); 
    }
}

void CellMigrationScore::initialize() {
    Processor::initialize();

    adjustPropertiesToInput();
}

bool CellMigrationScore::checkVolumeList(const VolumeList* collection) const {
    if (!collection || collection->empty()) {
        LERROR("Empty Volume List");
        return false;
    }

    // check transformation, spacing, dimensions and data type
    tgt::mat4 trafo = collection->first()->getVoxelToWorldMatrix();
    tgt::vec3 spacing = collection->first()->getSpacing();
    tgt::vec3 offset = collection->first()->getOffset();
    size_t numChannels = collection->first()->getNumChannels();
    std::string voxelFormat = collection->first()->getFormat(); 
    tgt::svec3 volDim = collection->first()->getDimensions();

    /*if (!collection->first()->hasRepresentation<VolumeDisk>()) {
        LERROR("First volume does not have a VolumeDisk representation!");
        return false;
    }*/

    for (size_t i = 1; i < collection->size(); ++i) {
        const VolumeBase* vol = collection->at(i);

        /*if (!vol->hasRepresentation<VolumeDisk>()) {
            LERROR("Volume " << i << " has no disk representation");
            return false;
        }*/

        bool valid = true;

        if (!(trafo == vol->getVoxelToWorldMatrix())) {
            LERROR("volumes differ in transformation matrix!");
            valid = false;
        }

        if (!(spacing == vol->getSpacing())) {
            LERROR("volumes differ in spacing!");
            valid = false;
        }

        if (!(offset == vol->getOffset())) {
            LERROR("volumes differ in offset!");
            valid = false;
        }

        if (!(numChannels == vol->getNumChannels())) {
            LERROR("volumes differ in number of channels!");
            valid = false;
        }

        if (!(voxelFormat == vol->getFormat())) {
            LERROR("volumes differ in voxel format!");
            valid = false;
        }

        if (!(volDim == vol->getDimensions())) {
            LERROR("volumes differ in dimensions!");
            valid = false;
        }
        
        if (!valid)
            return false;
    }

    // output multi-channel warning
    if (numChannels != 1)
        LWARNING("Volumes contain more than one channel - quantification will only be performed in the selected channel!");  

    return true;
}

void CellMigrationScore::computeQuantification() {

    setProgress(0.f);
    quantificationPlot_.clear();

    LINFO("Checking input data...");

    // check if the data is not empty and all volumes have the same parameters
    const VolumeList* collection = inport_.getData();
    int max = ((collection != 0) ? static_cast<int>(collection->size()) : 0);

    if (!collection || collection->empty()) {
        LERROR("No time steps found - no quantification computed");
        return;
    }

    if (!checkVolumeList(collection)) {
        LERROR("Cannot start quantification, volume list does not represent a time series");
        return;
    }

    // check if every volume has a volume disk representation
    /*for (size_t i = 0; i < collection->size(); ++i) {
        if (!collection->at(i)->hasRepresentation<VolumeDisk>()) {
            LERROR("No volume disk representation available for volumes in list, cannot compute quantification");
            return;
        }
    }*/

    LINFO("Starting quantification...");

    tgt::Stopwatch timer;
    timer.start();

    // result vector
    std::vector<double> sumResults(collection->size());
    std::vector<double> absoluteResults(collection->size());
    std::vector<double> normalizedResults(collection->size());
    
    // set time steps before and after the interval to -1
    for (size_t i = 0; i < static_cast<size_t>(quantificationFrame_.get().x); ++i) {
        sumResults.at(i) = -1.0;
        absoluteResults.at(i) = -1.0;
        normalizedResults.at(i) = -1.0;
    }

    for (size_t i = static_cast<size_t>(quantificationFrame_.get().y + 1); i < collection->size(); ++i) {
        sumResults.at(i) = -1.0;
        absoluteResults.at(i) = -1.0;
        normalizedResults.at(i) = -1.0;
    } 

    // compute sigma for gaussian distance weight
    double sigma = static_cast<double>(radius_.get()) / 2.0;

    // iterate over time steps and quantify the intensities 
    for (size_t i = static_cast<size_t>(quantificationFrame_.get().x); i < std::min(collection->size(), static_cast<size_t>(quantificationFrame_.get().y + 1)); ++i) {
        
        tgtAssert(static_cast<size_t>(channel_.get()) < collection->at(i)->getNumChannels(), "Not enough channels");

        /*if (!collection->at(i)->hasRepresentation<VolumeDisk>()) {
            // should not happen due to test in the beginning
            LERROR("Timestep " << i << " does not have a VolumeDisk representation. Setting score to -1");
            continue;
        }*/

        // get volume dimensions and bounding box of quantification
        tgt::svec3 volDim = collection->at(i)->getDimensions();
        tgt::vec3 spacing = collection->at(i)->getSpacing();
        tgt::dvec3 dSpacing(spacing);
    
        // compute the quantification radius (which is given in um) in voxels by taking into account the spacing
        tgt::ivec3 radius = tgt::iround(tgt::vec3(radius_.get()) / (spacing * 1000.f));
        tgt::svec3 llf = tgt::svec3(tgt::max(quantificationCenter_.get() - radius, tgt::ivec3::zero));
        // we do not want to take the radius in z-direction
        llf.z = 0;
        tgt::svec3 urb = tgt::min(tgt::svec3(quantificationCenter_.get() + radius), volDim - tgt::svec3::one);
        //see above
        urb.z = volDim.z - 1;

        // if the clip region is used, crop our bounds
        if (useClipRegion_.get()) {
            llf = tgt::max(llf, tgt::svec3(clipRegion_.get().getLLF()));
            urb = tgt::min(urb, tgt::svec3(clipRegion_.get().getURB()));
        }
        
        double score = 0.0;
        double sumScore = 0.0;
             
        // get VolumeRAM -> memory manager takes care of memory usage and will remove old VolumeRAM representations if necessary
        const VolumeRAM* vRam = collection->at(i)->getRepresentation<VolumeRAM>();
        if (!vRam) {
            LERROR("Could not load VolumeRAM from volume " << i << " - aborting quantification!");
            setProgress(0.f);
            return;
        }
        
        for (size_t z = llf.z; z <= urb.z; ++z) {
            /*const VolumeRAM* slice = collection->at(i)->getSlice(z); //vd->loadSlices(z,z);
            if (!slice) {
                LERROR("Could not load slice " << z << " from volume " << i << " - aborting quantification!");
                setProgress(0.f);
                return;
            }*/
            
            for (size_t y = llf.y; y <= urb.y; ++y) {
                for (size_t x = llf.x; x <= urb.x; ++x) {
                   
                    // compute distance to cylinder center in um and check if the voxel is inside the cylinder
                    double distance = tgt::length( (tgt::dvec2(x,y) - tgt::dvec2(quantificationCenter_.get().x, quantificationCenter_.get().y)) * (tgt::dvec2(dSpacing.x, dSpacing.y) * 1000.0) );
                    if (distance <= static_cast<double>(radius_.get())) {
                        // compute distance weight
                        double distanceWeight = std::exp(-0.5 * std::pow(distance / sigma, 2.0));

                        // intensity value is weighted by the voxel spacing (i.e., "size" of each voxel)
                        // sum values are just adding the voxels 
                        double value = static_cast<double>(vRam->getVoxelNormalized(x,y,z,static_cast<size_t>(channel_.get()))) * tgt::hmul(dSpacing * tgt::dvec3(1000.0));
                        sumScore += value;
                        score += value * distanceWeight; // the actual score weights the voxels using a gaussian distance function
                    }
                }
            }
            //delete slice;
        }
        sumResults.at(i) = sumScore;
        absoluteResults.at(i) = score;
             
        setProgress(std::min(0.99f, static_cast<float>(i - quantificationFrame_.get().x + 1) / static_cast<float>(quantificationFrame_.get().y - quantificationFrame_.get().x + 1)));
    }

    // compute normalized results by normalizing with the first frame
    for (size_t i = static_cast<size_t>(quantificationFrame_.get().x); (i < collection->size()) && (i <= static_cast<size_t>(quantificationFrame_.get().y)); ++i) {
        normalizedResults.at(i) = absoluteResults.at(i) / absoluteResults.at(static_cast<size_t>(quantificationFrame_.get().x));
    }
    
    setProgress(1.f);

    timer.stop();
    LINFO("Quantification Time: " << (float) timer.getRuntime() / 1000.f << " seconds");

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
    normalizedQuantificationPlot_.setData(normalizedPlotData, true);
}

void CellMigrationScore::adjustPropertiesToInput() {
    //if (!outport_.isInitialized())
    //    return;

    const VolumeList* collection = inport_.getData();
    int max = ((collection != 0) ? static_cast<int>(collection->size()) : 0);

    if (collection && !collection->empty()) {
        quantificationFrame_.setMaxValue(max - 1);
        // only adjust the center and radius max values if a valid time series is present
        if (checkVolumeList(inport_.getData())) {
            quantificationCenter_.setMaxValue(tgt::ivec3(collection->first()->getDimensions() - tgt::svec3::one));
            radius_.setMaxValue(tgt::max(collection->first()->getDimensions() - tgt::svec3::one));
            // adjust number of channels and clip region
            channel_.setMaxValue(collection->first()->getNumChannels() - 1);
            clipRegion_.setMaxValue(tgt::ivec3(collection->first()->getDimensions() - tgt::svec3::one));
        }
    }

}

void CellMigrationScore::useClipRegionChanged() {
    clipRegion_.setVisibleFlag(useClipRegion_.get());
}

} // namespace
