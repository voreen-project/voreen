/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#include "octreesegmentationquantification.h"

#include <queue>
#include <iostream>
#include <string>
#include <sstream>

#include "tgt/filesystem.h"
#include "tgt/bounds.h"
#include "tgt/stopwatch.h"

#include "voreen/core/datastructures/octree/octreeutils.h"

#include "../util/octreequantificationthread.h"
#include "../util/octreequantificationnodequeue.h"

#include "voreen/core/voreenapplication.h"

namespace voreen {

const std::string OctreeSegmentationQuantification::loggerCat_("voreen.OctreeQuantification.OctreeSegmentationQuantification");

OctreeSegmentationQuantification::OctreeSegmentationQuantification()
        : Processor()
        // input ports
        , volumeInport_(Port::INPORT, "volumeinport", "Volume Inport (Octree required)")
        , roi1Port_(Port::INPORT, "roi1", "Region of Interest 1 Input")
        , roi2Port_(Port::INPORT, "roi2", "Region of Interest 2 Input")        
        , roi3Port_(Port::INPORT, "roi3", "Region of Interest 3 Input")        
        , roi4Port_(Port::INPORT, "roi4", "Region of Interest 4 Input")         
        // textual information
        //, quantificationTextPort_(Port::OUTPORT, "QuantificationTextPort", "Quantification Text Information")
        // plot of the textual information (for export)
        //, quantificationInfoPlot_(Port::OUTPORT, "QuantificationInfoPlot", "Quantification Info Plot")
        // TODO: rework
        // absolute histograms
/*        , histogramPlotPortChannelA_(Port::OUTPORT, "HistogramPlotPortChannelA", "Histogram Data Channel A")
        , histogramPlotPortChannelB_(Port::OUTPORT, "HistogramPlotPortChannelB", "Histogram Data Channel B")
        , histogramPlotPortChannelC_(Port::OUTPORT, "HistogramPlotPortChannelC", "Histogram Data Channel C")
        , thresholdedHistogramPlotPortChannelA_(Port::OUTPORT, "ThresholdedHistogramPlotPortChannelA", "Thresholded Histogram Data Channel A")
        , thresholdedHistogramPlotPortChannelB_(Port::OUTPORT, "ThresholdedHistogramPlotPortChannelB", "Thresholded Histogram Data Channel B")
        , thresholdedHistogramPlotPortChannelC_(Port::OUTPORT, "ThresholdedHistogramPlotPortChannelC", "Thresholded Histogram Data Channel C")
        // relative histograms
        , relativeHistogramPlotPortChannelA_(Port::OUTPORT, "RelativeHistogramPlotPortChannelA", "Relative Histogram Data Channel A")
        , relativeHistogramPlotPortChannelB_(Port::OUTPORT, "RelativeHistogramPlotPortChannelB", "Relative Histogram Data Channel B")
        , relativeHistogramPlotPortChannelC_(Port::OUTPORT, "RelativeHistogramPlotPortChannelC", "Relative Histogram Data Channel C")
        , thresholdedRelativeHistogramPlotPortChannelA_(Port::OUTPORT, "ThresholdedRelativeHistogramPlotPortChannelA", "Thresholded Relative Histogram Data Channel A")
        , thresholdedRelativeHistogramPlotPortChannelB_(Port::OUTPORT, "ThresholdedRelativeHistogramPlotPortChannelB", "Thresholded Relative Histogram Data Channel B")
        , thresholdedRelativeHistogramPlotPortChannelC_(Port::OUTPORT, "ThresholdedRelativeHistogramPlotPortChannelC", "Thresholded Relative Histogram Data Channel C")
        // surface plots  
        , scatterABHistogramPort_(Port::OUTPORT, "ScatterABHistogram", "Scatter AB Histogram")
        , scatterCBHistogramPort_(Port::OUTPORT, "ScatterCBHistogram", "Scatter CB Histogram")
        , thresholdedScatterABHistogramPort_(Port::OUTPORT, "ThresholdedScatterABHistogram", "Thresholded Scatter AB Histogram")
        , thresholdedScatterCBHistogramPort_(Port::OUTPORT, "ThresholdedScatterCBHistogram", "Thresholded Scatter CB Histogram")
        // surface plots relative to number of voxels in region
        , scatterABRelativeHistogramPort_(Port::OUTPORT, "ScatterABRelativeHistogram", "Scatter AB Relative Histogram")
        , scatterCBRelativeHistogramPort_(Port::OUTPORT, "ScatterCBRelativeHistogram", "Scatter CB Relative Histogram")
        , thresholdedScatterABRelativeHistogramPort_(Port::OUTPORT, "ThresholdedScatterABRelativeHistogram", "Thresholded Scatter AB Relative Histogram")
        , thresholdedScatterCBRelativeHistogramPort_(Port::OUTPORT, "ThresholdedScatterCBRelativeHistogram", "Thresholded Scatter CB Relative Histogram")
*/
        // properties
        , quantificationProgressProperty_("quantificationprogressproperty", "Quantification Progress")
        , simMode_("segmentationinterpolationmode", "Segmentation Voxel Interpolation")
        , computeQuantification_("computequantification", "Compute Quantification")
        , exportFile_("exportfile", "Data Export (CSV)", "Save quantification results...", VoreenApplication::app()->getUserDataPath(), "Commaseparated Files(*.csv)", FileDialogProperty::SAVE_FILE, Processor::VALID)
        , exportButton_("exportbutton", "Export Data")
        // other internal variables
        // TODO: rework
/*        , numVoxelsInVolume_(0)
        , numVoxelsInOuterStructure_(0)
        , numVoxelsInInnerStructure_(0)
*/
        , threadVoxelCounter_(0)
        , currentResults_(1,1)
{
    // add volume input
    volumeInport_.onChange(MemberFunctionCallback<OctreeSegmentationQuantification>(this, &OctreeSegmentationQuantification::inputVolumeChanged));
    addPort(volumeInport_);

    roi1Port_.onChange(MemberFunctionCallback<OctreeSegmentationQuantification>(this, &OctreeSegmentationQuantification::segmentationChanged));
    addPort(roi1Port_);

    roi2Port_.onChange(MemberFunctionCallback<OctreeSegmentationQuantification>(this, &OctreeSegmentationQuantification::segmentationChanged));
    addPort(roi2Port_);

    roi3Port_.onChange(MemberFunctionCallback<OctreeSegmentationQuantification>(this, &OctreeSegmentationQuantification::segmentationChanged));
    addPort(roi3Port_);

    roi4Port_.onChange(MemberFunctionCallback<OctreeSegmentationQuantification>(this, &OctreeSegmentationQuantification::segmentationChanged));
    addPort(roi4Port_);

    //addPort(quantificationTextPort_);
    //addPort(quantificationInfoPlot_);

/*    addPort(histogramPlotPortChannelA_);
    addPort(histogramPlotPortChannelB_);
    addPort(histogramPlotPortChannelC_);

    addPort(thresholdedHistogramPlotPortChannelA_);
    addPort(thresholdedHistogramPlotPortChannelB_);
    addPort(thresholdedHistogramPlotPortChannelC_);

    addPort(relativeHistogramPlotPortChannelA_);
    addPort(relativeHistogramPlotPortChannelB_);
    addPort(relativeHistogramPlotPortChannelC_);

    addPort(thresholdedRelativeHistogramPlotPortChannelA_);
    addPort(thresholdedRelativeHistogramPlotPortChannelB_);
    addPort(thresholdedRelativeHistogramPlotPortChannelC_);

    addPort(scatterABHistogramPort_);
    addPort(scatterCBHistogramPort_);
    addPort(thresholdedScatterABHistogramPort_);
    addPort(thresholdedScatterCBHistogramPort_);

    addPort(scatterABRelativeHistogramPort_);
    addPort(scatterCBRelativeHistogramPort_);
    addPort(thresholdedScatterABRelativeHistogramPort_);
    addPort(thresholdedScatterCBRelativeHistogramPort_);*/

    // add progress bar
    addProperty(quantificationProgressProperty_);
    addProgressBar(&quantificationProgressProperty_);

    // add interpolation quality property
    simMode_.addOption("nearest", "Nearest", SIM_NEAREST);
    simMode_.addOption("linear", "Linear", SIM_LINEAR);
    simMode_.addOption("cubic", "Cubic", SIM_CUBIC);
    addProperty(simMode_);
    
    computeQuantification_.onClick(MemberFunctionCallback<OctreeSegmentationQuantification>(this, &OctreeSegmentationQuantification::computeQuantification));
    addProperty(computeQuantification_);

    exportFile_.setGroupID("dataexport");
    exportButton_.setGroupID("dataexport");
    addProperty(exportFile_);
    addProperty(exportButton_);
    setPropertyGroupGuiName("dataexport", "Data Export");

    exportButton_.onClick(MemberFunctionCallback<OctreeSegmentationQuantification>(this, &OctreeSegmentationQuantification::exportData));
}

Processor* OctreeSegmentationQuantification::create() const {
    return new OctreeSegmentationQuantification();
}

bool OctreeSegmentationQuantification::isReady() const {
    if (!isInitialized())
        return false;

    bool someRoiReady = roi1Port_.isReady() || roi2Port_.isReady() || roi3Port_.isReady() || roi4Port_.isReady();
    
    if (!volumeInport_.isReady() || !someRoiReady)
        return false;

    return true;
}

void OctreeSegmentationQuantification::inputVolumeChanged() {

    setProgress(0.f);

    // invalidate / delete current data if present
    clearQuantificationData();
}

void OctreeSegmentationQuantification::segmentationChanged() {
    // invalidate / delete current data if present
    clearQuantificationData();
}

void OctreeSegmentationQuantification::process() {
    // do nothing
}

void OctreeSegmentationQuantification::clearQuantificationData() {
    currentResults_.clear();
}


void OctreeSegmentationQuantification::computeQuantification() {

    deactivateProperties();
    setProgress(0.f);

    // invalidate / delete current data if present
    clearQuantificationData();

    if (!volumeInport_.getData() || !volumeInport_.getData()->hasRepresentation<VolumeOctree>()) {
        LERROR("No input volume or no octree representation in input volume");
        reactivateProperties();
        return;
    }

    bool someRoiReady = roi1Port_.isReady() || roi2Port_.isReady() || roi3Port_.isReady() || roi4Port_.isReady();
    
    if (!someRoiReady) {
        LERROR("No region of interest input");
        reactivateProperties();
        return;
    }

    // get volume octree
    const VolumeOctree* volumeOctree = volumeInport_.getData()->getRepresentation<VolumeOctree>();
    size_t numVoxels = volumeInport_.getData()->getNumVoxels();
    tgt::svec3 volumeURB = volumeInport_.getData()->getDimensions() - tgt::svec3::one;

    // check for ROI volumes
    tgt::bvec4 roiPresent = tgt::bvec4( (roi1Port_.getData() != 0),
                                        (roi2Port_.getData() != 0),
                                        (roi3Port_.getData() != 0),  
                                        (roi4Port_.getData() != 0)
                                      );
    // get ROI volumes
    std::vector<const VolumeRAM*> roiVolumes;
    size_t maxRoi = 0; 
    for (size_t i = 0; i < 4; ++i) {
        if (roiPresent[i])
            maxRoi = i;
    }
    roiVolumes.push_back(roi1Port_.getData() ? roi1Port_.getData()->getRepresentation<VolumeRAM>() : 0);
    if (maxRoi > 0)
        roiVolumes.push_back(roi2Port_.getData() ? roi2Port_.getData()->getRepresentation<VolumeRAM>() : 0);
    if (maxRoi > 1)
        roiVolumes.push_back(roi3Port_.getData() ? roi3Port_.getData()->getRepresentation<VolumeRAM>() : 0);
    if (maxRoi > 2)
        roiVolumes.push_back(roi4Port_.getData() ? roi4Port_.getData()->getRepresentation<VolumeRAM>() : 0);

    // if a ROI is missing in between, we still use it, but set all of its values to 0 
    size_t numROIs = roiVolumes.size();

    // compute ratios for determining voxel coordinates within the downsampled segmentation ROIs
    std::vector<tgt::vec3> voxelCoordinateRatios;
    for (size_t i = 0; i < numROIs; ++i) {
        if (roiVolumes.at(i) != 0) 
            voxelCoordinateRatios.push_back(tgt::vec3(roiVolumes.at(i)->getDimensions() - tgt::svec3::one) / tgt::vec3(volumeURB));
        else
            voxelCoordinateRatios.push_back(tgt::vec3::zero);
    }

    // check number of channels
    size_t numChannels = volumeInport_.getData()->getNumChannels();

    // clear quantification results / create new results with the current parameters
    currentResults_ = OctreeQuantificationResults(numChannels, numROIs);

    SegmentationInterpolationMode simMode = simMode_.getValue();

    threadVoxelCounter_ = 0;

    /*
     * Actual Quantification
     */
    tgt::Stopwatch stopwatch;
    LINFO("Starting Quantification");
    stopwatch.start();

    // get number of CPU cores to start the right amount of threads
    size_t numCores = static_cast<size_t>(boost::thread::hardware_concurrency());
    if (numCores == 0) {
        numCores = 1;
        LWARNING("Could not detect number of CPU cores. Using (slow) single threaded quantification.");
    }
    else 
        LINFO("Starting " << numCores << " threads for quantification.");

    // queue for nodes to be scheduled for quantification (limit to pause traversal if necessary to leave the quantification threads more processing time)
    OctreeQuantificationNodeQueue quantificationQueue(numCores * 8);

    // mutex for preventing race conditions when threads write the results of one block into the global results
    boost::mutex globalValueMutex;

    // now create the threads and start their execution
    std::vector<OctreeQuantificationThread*> quantificationThreads;
    for (size_t threadCount = 0; threadCount < numCores; ++threadCount) {
        // TODO: change parameters for threads
        OctreeQuantificationThread* thread = new OctreeQuantificationThread(this, quantificationQueue, volumeOctree, roiVolumes, voxelCoordinateRatios, 
                                                                            simMode, volumeURB, numChannels, globalValueMutex);
      
        quantificationThreads.push_back(thread);
        thread->run();
    }

    // start traversal at the octree root with the full octree dimensions
    const VolumeOctreeNode* root = volumeOctree->getRootNode();
    std::stack<std::pair< tgt::IntBounds, const VolumeOctreeNode* > > nodeStack;
    nodeStack.push(std::make_pair(tgt::IntBounds(tgt::ivec3(static_cast<int>(0)), tgt::ivec3(volumeOctree->getOctreeDim() - tgt::svec3::one)), root)); 

    // actual traversal
    while(!nodeStack.empty()) {
        // get the current node and its bounds and pop it from the stack
        const VolumeOctreeNode* currentNode = nodeStack.top().second;
        tgt::IntBounds currentBounds = nodeStack.top().first;
        nodeStack.pop();
    
        if (!currentNode->inVolume())     // dismiss nodes that are outside the volume 
            continue;
        else if (!currentNode->isLeaf() && !currentNode->isHomogeneous()) { // -> current node has to be traversed further
            // perform a depth-first-search by pushing the children on the stack
            tgt::ivec3 llf[] = { currentBounds.getLLF(), (tgt::ivec3(currentBounds.getURB()) - tgt::ivec3(currentBounds.getLLF())) / 2 + tgt::ivec3(currentBounds.getLLF()) + tgt::ivec3::one };
            tgt::ivec3 urb[] = { (tgt::ivec3(currentBounds.getURB()) - tgt::ivec3(currentBounds.getLLF())) / 2 + tgt::ivec3(currentBounds.getLLF()), currentBounds.getURB() };

            for (size_t z = 0; z <= 1; ++z) {
                for (size_t y = 0; y <= 1; ++y) {
                    for (size_t x = 0; x <= 1; ++x) {
                        tgt::ivec3 childLLF = tgt::ivec3(llf[x].x, llf[y].y, llf[z].z);
                        tgt::ivec3 childURB = tgt::ivec3(urb[x].x, urb[y].y, urb[z].z);
                        nodeStack.push(std::make_pair(tgt::IntBounds(childLLF, childURB), currentNode->children_[z*4 + y*2 + x]));
                    }
                }
            }
            continue;
        }
        else {
            // at this point the current node is either a homogeneous node or a leaf -> schedule it for quantification
            quantificationQueue.push(std::make_pair(currentBounds, currentNode));
        }

        // update progress bar after each block
        float p = std::min(static_cast<float>(threadVoxelCounter_) / static_cast<float>(numVoxels), 0.99f);
        setProgress(p);
    }

    // wait for all threads to finish and then delete them
    for (size_t t = 0; t < quantificationThreads.size(); ++t) {
        OctreeQuantificationThread* thread = quantificationThreads.at(t);
        thread->interruptAndJoin();
        delete thread;
        // update progress bar after each thread is finished
        float p = std::min(static_cast<float>(threadVoxelCounter_) / static_cast<float>(numVoxels), 0.99f);
        setProgress(p);
    }

    if (currentResults_.getNumVoxelsInVolume() != threadVoxelCounter_) {
        LERROR("Number of voxels visited by threads does not equal number of voxels visited in complete quantification");
        clearQuantificationData();
        setProgress(0.f);
        return;
    }
  
    if (currentResults_.getNumVoxelsInVolume() != numVoxels) {
        LERROR("Number of visited voxels does not equal number of voxels volume: " << currentResults_.getNumVoxelsInVolume() << " vs. " << numVoxels << " - something must have gone terribly wrong :(");
        clearQuantificationData();
        setProgress(0.f);
        return;
    }

    // TODO: output information using plotting etc.
/*

    std::stringstream quantificationTextPortStream;
    quantificationTextPortStream << "Quantification Info:" << std::endl;
    quantificationTextPortStream << "-------------------" << std::endl;
    quantificationTextPortStream << "Number of voxels in volume: " << numVoxelsInVolume_ << std::endl;
    quantificationTextPortStream << "Number of voxels inside the outer segmented structure: " << numVoxelsInOuterStructure_ << std::endl;
    quantificationTextPortStream << "Number of voxels inside the inner segmented structure: " << numVoxelsInInnerStructure_ << std::endl;
    quantificationTextPortStream << "Quantification Channel A: Data Channel " << quantificationChannelA_.get() << std::endl;
    quantificationTextPortStream << "Quantification Channel B: Data Channel " << quantificationChannelB_.get() << std::endl;
    quantificationTextPortStream << "Quantification Channel C: Data Channel " << quantificationChannelC_.get() << std::endl;
    quantificationTextPort_.setData(quantificationTextPortStream.str());


    PlotData* plotHistogramChannelA = new PlotData(0, 4);
    plotHistogramChannelA->setColumnLabel(0, "Data Value (Channel A)");
    plotHistogramChannelA->setColumnLabel(1, "In Volume");
    plotHistogramChannelA->setColumnLabel(2, "In Outer Structure");
    plotHistogramChannelA->setColumnLabel(3, "In Inner Structure");

    for (size_t i = 0; i < channelAHistogramVolume_.size(); ++i) {
        std::vector<PlotCellValue> v(4);
        v.at(0) = PlotCellValue(static_cast<plot_t>(i));
        v.at(1) = PlotCellValue(static_cast<plot_t>(channelAHistogramVolume_.at(i)));
        v.at(2) = PlotCellValue(static_cast<plot_t>(channelAHistogramOuterStructure_.at(i)));
        v.at(3) = PlotCellValue(static_cast<plot_t>(channelAHistogramInnerStructure_.at(i)));
        bool inserted = plotHistogramChannelA->insert(v);
        if (!inserted) {
            LERROR("Could not insert data in channel A");
            delete plotHistogramChannelA;
            plotHistogramChannelA = 0;
            break;
        }
    }
    histogramPlotPortChannelA_.setData(plotHistogramChannelA, true);

    PlotData* plotHistogramChannelB = new PlotData(0, 4);
    plotHistogramChannelB->setColumnLabel(0, "Data Value (Channel B)");
    plotHistogramChannelB->setColumnLabel(1, "In Volume");
    plotHistogramChannelB->setColumnLabel(2, "In Outer Structure");
    plotHistogramChannelB->setColumnLabel(3, "In Inner Structure");

    for (size_t i = 0; i < channelBHistogramVolume_.size(); ++i) {
        std::vector<PlotCellValue> v(4);
        v.at(0) = PlotCellValue(static_cast<plot_t>(i));
        v.at(1) = PlotCellValue(static_cast<plot_t>(channelBHistogramVolume_.at(i)));
        v.at(2) = PlotCellValue(static_cast<plot_t>(channelBHistogramOuterStructure_.at(i)));
        v.at(3) = PlotCellValue(static_cast<plot_t>(channelBHistogramInnerStructure_.at(i)));
        bool inserted = plotHistogramChannelB->insert(v);
        if (!inserted) {
            LERROR("Could not insert data in channel B");
            delete plotHistogramChannelB;
            plotHistogramChannelB = 0;
            break;
        }
    }
    histogramPlotPortChannelB_.setData(plotHistogramChannelB, true);

    PlotData* plotHistogramChannelC = new PlotData(0, 4);
    plotHistogramChannelC->setColumnLabel(0, "Data Value (Channel C)");
    plotHistogramChannelC->setColumnLabel(1, "In Volume");
    plotHistogramChannelC->setColumnLabel(2, "In Outer Structure");
    plotHistogramChannelC->setColumnLabel(3, "In Inner Structure");

    for (size_t i = 0; i < channelCHistogramVolume_.size(); ++i) {
        std::vector<PlotCellValue> v(4);
        v.at(0) = PlotCellValue(static_cast<plot_t>(i));
        v.at(1) = PlotCellValue(static_cast<plot_t>(channelCHistogramVolume_.at(i)));
        v.at(2) = PlotCellValue(static_cast<plot_t>(channelCHistogramOuterStructure_.at(i)));
        v.at(3) = PlotCellValue(static_cast<plot_t>(channelCHistogramInnerStructure_.at(i)));
        bool inserted = plotHistogramChannelC->insert(v);
        if (!inserted) {
            LERROR("Could not insert data in channel C");
            delete plotHistogramChannelC;
            plotHistogramChannelC = 0;
            break;
        }
    }
    histogramPlotPortChannelC_.setData(plotHistogramChannelC, true);

    // relative histograms

    PlotData* relativePlotHistogramChannelA = new PlotData(0, 4);
    relativePlotHistogramChannelA->setColumnLabel(0, "Data Value (Channel A)");
    relativePlotHistogramChannelA->setColumnLabel(1, "In Volume (relative in to voxels in volume in %)");
    relativePlotHistogramChannelA->setColumnLabel(2, "In Outer Structure (relative in to voxels in outer structure in %)");
    relativePlotHistogramChannelA->setColumnLabel(3, "In Inner Structure (relative in to voxels in inner structure in %)");

    for (size_t i = 0; i < channelAHistogramVolume_.size(); ++i) {
        std::vector<PlotCellValue> v(4);
        v.at(0) = PlotCellValue(static_cast<plot_t>(i));
        v.at(1) = PlotCellValue(static_cast<plot_t>(static_cast<float>(channelAHistogramVolume_.at(i)) / (static_cast<float>(numVoxelsInVolume_) / 100.f)));
        v.at(2) = PlotCellValue(static_cast<plot_t>(static_cast<float>(channelAHistogramOuterStructure_.at(i)) / (static_cast<float>(numVoxelsInOuterStructure_) / 100.f)));
        v.at(3) = PlotCellValue(static_cast<plot_t>(static_cast<float>(channelAHistogramInnerStructure_.at(i)) / (static_cast<float>(numVoxelsInInnerStructure_) / 100.f)));
        bool inserted = relativePlotHistogramChannelA->insert(v);
        if (!inserted) {
            LERROR("Could not insert relative histogram data in channel A");
            delete relativePlotHistogramChannelA;
            relativePlotHistogramChannelA = 0;
            break;
        }
    }
    relativeHistogramPlotPortChannelA_.setData(relativePlotHistogramChannelA, true);

    PlotData* relativePlotHistogramChannelB = new PlotData(0, 4);
    relativePlotHistogramChannelB->setColumnLabel(0, "Data Value (Channel B)");
    relativePlotHistogramChannelB->setColumnLabel(1, "In Volume (relative in to voxels in volume in %)");
    relativePlotHistogramChannelB->setColumnLabel(2, "In Outer Structure (relative in to voxels in outer structure in %)");
    relativePlotHistogramChannelB->setColumnLabel(3, "In Inner Structure (relative in to voxels in inner structure in %)");

    for (size_t i = 0; i < channelBHistogramVolume_.size(); ++i) {
        std::vector<PlotCellValue> v(4);
        v.at(0) = PlotCellValue(static_cast<plot_t>(i));
        v.at(1) = PlotCellValue(static_cast<plot_t>(static_cast<float>(channelBHistogramVolume_.at(i)) / (static_cast<float>(numVoxelsInVolume_) / 100.f)));
        v.at(2) = PlotCellValue(static_cast<plot_t>(static_cast<float>(channelBHistogramOuterStructure_.at(i)) / (static_cast<float>(numVoxelsInOuterStructure_) / 100.f)));
        v.at(3) = PlotCellValue(static_cast<plot_t>(static_cast<float>(channelBHistogramInnerStructure_.at(i)) / (static_cast<float>(numVoxelsInInnerStructure_) / 100.f)));
        bool inserted = relativePlotHistogramChannelB->insert(v);
        if (!inserted) {
            LERROR("Could not insert relative histogram data in channel B");
            delete relativePlotHistogramChannelB;
            relativePlotHistogramChannelB = 0;
            break;
        }
    }
    relativeHistogramPlotPortChannelB_.setData(relativePlotHistogramChannelB, true);

    PlotData* relativePlotHistogramChannelC = new PlotData(0, 4);
    relativePlotHistogramChannelC->setColumnLabel(0, "Data Value (Channel C)");
    relativePlotHistogramChannelC->setColumnLabel(1, "In Volume (relative in to voxels in volume in %)");
    relativePlotHistogramChannelC->setColumnLabel(2, "In Outer Structure (relative in to voxels in outer structure in %)");
    relativePlotHistogramChannelC->setColumnLabel(3, "In Inner Structure (relative in to voxels in inner structure in %)");

    for (size_t i = 0; i < channelCHistogramVolume_.size(); ++i) {
        std::vector<PlotCellValue> v(4);
        v.at(0) = PlotCellValue(static_cast<plot_t>(i));
        v.at(1) = PlotCellValue(static_cast<plot_t>(static_cast<float>(channelCHistogramVolume_.at(i)) / (static_cast<float>(numVoxelsInVolume_) / 100.f)));
        v.at(2) = PlotCellValue(static_cast<plot_t>(static_cast<float>(channelCHistogramOuterStructure_.at(i)) / (static_cast<float>(numVoxelsInOuterStructure_) / 100.f)));
        v.at(3) = PlotCellValue(static_cast<plot_t>(static_cast<float>(channelCHistogramInnerStructure_.at(i)) / (static_cast<float>(numVoxelsInInnerStructure_) / 100.f)));
        bool inserted = relativePlotHistogramChannelC->insert(v);
        if (!inserted) {
            LERROR("Could not insert relative histogram data in channel C");
            delete relativePlotHistogramChannelC;
            relativePlotHistogramChannelC = 0;
            break;
        }
    }
    relativeHistogramPlotPortChannelC_.setData(relativePlotHistogramChannelC, true);
    
    // 2D/3D plot data
    PlotData* scatterABPlotData = new PlotData(0, 5);
    scatterABPlotData->setColumnLabel(0, "Channel A (percent of max)");
    scatterABPlotData->setColumnLabel(1, "Channel B (percent of max)");
    scatterABPlotData->setColumnLabel(2, "Num Voxels (Volume)");   
    scatterABPlotData->setColumnLabel(3, "Num Voxels (Outer Structure)");       
    scatterABPlotData->setColumnLabel(4, "Num Voxels (Inner Structure)");
    for (size_t i = 0; i < 1001; ++i) {
        for (size_t j = 0; j < 1001; ++j) {
            std::vector<PlotCellValue> v(5);
            v.at(0) = PlotCellValue(static_cast<plot_t>(i) / 10.0);
            v.at(1) = PlotCellValue(static_cast<plot_t>(j) / 10.0);
            v.at(2) = PlotCellValue(static_cast<plot_t>(scatterABVolume_.at(i*1001 + j)));
            v.at(3) = PlotCellValue(static_cast<plot_t>(scatterABOuterStructure_.at(i*1001 + j)));
            v.at(4) = PlotCellValue(static_cast<plot_t>(scatterABInnerStructure_.at(i*1001 + j)));
            bool inserted = scatterABPlotData->insert(v);
            if (!inserted) {
                LERROR("Could not insert scatter data for channels A and B");
                delete scatterABPlotData;
                scatterABPlotData = 0;
                break;
            }
        }
    }
    scatterABHistogramPort_.setData(scatterABPlotData);

    PlotData* scatterCBPlotData = new PlotData(0, 5);
    scatterCBPlotData->setColumnLabel(0, "Channel C (percent of max)");
    scatterCBPlotData->setColumnLabel(1, "Channel B (percent of max)");
    scatterCBPlotData->setColumnLabel(2, "Num Voxels (Volume)");   
    scatterCBPlotData->setColumnLabel(3, "Num Voxels (Outer Structure)");       
    scatterCBPlotData->setColumnLabel(4, "Num Voxels (Inner Structure)");
    for (size_t i = 0; i < 1001; ++i) {
        for (size_t j = 0; j < 1001; ++j) {
            std::vector<PlotCellValue> v(5);
            v.at(0) = PlotCellValue(static_cast<plot_t>(i) / 10.0);
            v.at(1) = PlotCellValue(static_cast<plot_t>(j) / 10.0);
            v.at(2) = PlotCellValue(static_cast<plot_t>(scatterCBVolume_.at(i*1001 + j)));
            v.at(3) = PlotCellValue(static_cast<plot_t>(scatterCBOuterStructure_.at(i*1001 + j)));
            v.at(4) = PlotCellValue(static_cast<plot_t>(scatterCBInnerStructure_.at(i*1001 + j)));
            bool inserted = scatterCBPlotData->insert(v);
            if (!inserted) {
                LERROR("Could not insert scatter data for channels C and B");
                delete scatterCBPlotData;
                scatterCBPlotData = 0;
                break;
            }
        }
    }
    scatterCBHistogramPort_.setData(scatterCBPlotData);

    // relative 2D/3D plot data
    PlotData* scatterABRelativePlotData = new PlotData(0, 5);
    scatterABRelativePlotData->setColumnLabel(0, "Channel A (percent of max)");
    scatterABRelativePlotData->setColumnLabel(1, "Channel B (percent of max)");
    scatterABRelativePlotData->setColumnLabel(2, "Amount of Volume (%)");   
    scatterABRelativePlotData->setColumnLabel(3, "Amount of Outer Structure (%)");       
    scatterABRelativePlotData->setColumnLabel(4, "Amount of Inner Structure (%)");
    for (size_t i = 0; i < 1001; ++i) {
        for (size_t j = 0; j < 1001; ++j) {
            std::vector<PlotCellValue> v(5);
            v.at(0) = PlotCellValue(static_cast<plot_t>(i) / 10.0);
            v.at(1) = PlotCellValue(static_cast<plot_t>(j) / 10.0);
            v.at(2) = PlotCellValue(static_cast<plot_t>(static_cast<float>(scatterABVolume_.at(i*1001 + j)) / static_cast<float>(numVoxelsInVolume_) * 100.f ));
            v.at(3) = PlotCellValue(static_cast<plot_t>(static_cast<float>(scatterABOuterStructure_.at(i*1001 + j)) / static_cast<float>(numVoxelsInOuterStructure_) * 100.f ));
            v.at(4) = PlotCellValue(static_cast<plot_t>(static_cast<float>(scatterABInnerStructure_.at(i*1001 + j)) / static_cast<float>(numVoxelsInInnerStructure_) * 100.f ));
            bool inserted = scatterABRelativePlotData->insert(v);
            if (!inserted) {
                LERROR("Could not insert relative scatter data for channels A and B");
                delete scatterABRelativePlotData;
                scatterABRelativePlotData = 0;
                break;
            }
        }
    }
    scatterABRelativeHistogramPort_.setData(scatterABRelativePlotData);

    PlotData* scatterCBRelativePlotData = new PlotData(0, 5);
    scatterCBRelativePlotData->setColumnLabel(0, "Channel C (percent of max)");
    scatterCBRelativePlotData->setColumnLabel(1, "Channel B (percent of max)");
    scatterCBRelativePlotData->setColumnLabel(2, "Amount of Volume (%)");   
    scatterCBRelativePlotData->setColumnLabel(3, "Amount of Outer Structure (%)");       
    scatterCBRelativePlotData->setColumnLabel(4, "Amount of Inner Structure (%)");
    for (size_t i = 0; i < 1001; ++i) {
        for (size_t j = 0; j < 1001; ++j) {
            std::vector<PlotCellValue> v(5);
            v.at(0) = PlotCellValue(static_cast<plot_t>(i) / 10.0);
            v.at(1) = PlotCellValue(static_cast<plot_t>(j) / 10.0);
            v.at(2) = PlotCellValue(static_cast<plot_t>(static_cast<float>(scatterCBVolume_.at(i*1001 + j)) / static_cast<float>(numVoxelsInVolume_) * 100.f));
            v.at(3) = PlotCellValue(static_cast<plot_t>(static_cast<float>(scatterCBOuterStructure_.at(i*1001 + j)) / static_cast<float>(numVoxelsInOuterStructure_) * 100.f));
            v.at(4) = PlotCellValue(static_cast<plot_t>(static_cast<float>(scatterCBInnerStructure_.at(i*1001 + j)) / static_cast<float>(numVoxelsInInnerStructure_) * 100.f));
            bool inserted = scatterCBRelativePlotData->insert(v);
            if (!inserted) {
                LERROR("Could not insert relative scatter data for channels C and B");
                delete scatterCBRelativePlotData;
                scatterCBRelativePlotData = 0;
                break;
            }
        }
    }
    scatterCBRelativeHistogramPort_.setData(scatterCBRelativePlotData);
*/
    stopwatch.stop();
    LINFO("Quantification finished. Runtime: " << static_cast<double>(stopwatch.getRuntime()) / 1000.0 << " seconds");

    setProgress(1.f);
    reactivateProperties();
}

void OctreeSegmentationQuantification::exportData() {

    deactivateProperties();
    setProgress(0.f);

    if (!volumeInport_.getData() || !volumeInport_.getData()->hasRepresentation<VolumeOctree>()) {
        LERROR("No input volume or no octree representation in input volume");
        reactivateProperties();
        return;
    }

    bool someRoiReady = roi1Port_.isReady() || roi2Port_.isReady() || roi3Port_.isReady() || roi4Port_.isReady();
    
    if (!someRoiReady) {
        LERROR("No region of interest input");
        reactivateProperties();
        return;
    }

    if (currentResults_.getNumVoxelsInVolume() == 0) {
        LERROR("No voxels quantified");
        reactivateProperties();
        return;
    }

    std::string exportPath = tgt::FileSystem::cleanupPath(exportFile_.get());
    std::ofstream exportFile(exportPath);
    if (!exportFile.good()) {
        LERROR("Cannot write to " << exportPath);
        exportFile.close();
        reactivateProperties();
        return;
    }


    // write number of voxels
    exportFile << "Number of Voxels in Volume" << std::endl;
    exportFile << currentResults_.getNumVoxelsInVolume() << std::endl;

    // write number of voxels in the different ROIs
    for (size_t i = 0; i < currentResults_.getNumVoxelsInROIs().size(); ++i) {
        exportFile << "Num Voxels in ROI " << i;
        if (i < currentResults_.getNumVoxelsInROIs().size() - 1)
            exportFile << ",";
        else
            exportFile << std::endl;
    }

    for (size_t i = 0; i < currentResults_.getNumVoxelsInROIs().size(); ++i) {
        exportFile << currentResults_.getNumVoxelsInROIs().at(i);
        if (i < currentResults_.getNumVoxelsInROIs().size() - 1)
            exportFile << ",";
        else
            exportFile << std::endl;
    }

    // write histograms for each pair (channel, ROI)
    std::vector<std::vector<std::vector<size_t> > >& histograms = currentResults_.getIntensityHistogramsInROIs();
    
    // write header
    exportFile << "Intensity (uint16)" << ",";
    for (size_t channel = 0; channel < currentResults_.getNumChannels(); ++channel) {
        for (size_t roi = 0; roi < currentResults_.getNumROIs(); ++roi) {
            exportFile << "ROI " << roi << " - Channel " << channel;
            if (!(channel == currentResults_.getNumChannels() - 1 && roi == currentResults_.getNumROIs() - 1))
                exportFile << ",";
            else
                exportFile << std::endl;
        }
    }

    // write histograms (TODO: very inefficient traversal...)
    size_t maxIntensity = (size_t) std::numeric_limits<uint16_t>::max();
    for (size_t intensity = 0; intensity <= maxIntensity; ++intensity) {
        exportFile << intensity << ",";
        for (size_t channel = 0; channel < currentResults_.getNumChannels(); ++channel) {
            for (size_t roi = 0; roi < currentResults_.getNumROIs(); ++roi) {
                exportFile << histograms.at(channel).at(roi).at(intensity);
                if (!(channel == currentResults_.getNumChannels() - 1 && roi == currentResults_.getNumROIs() - 1))
                    exportFile << ",";
                else
                    exportFile << std::endl;
            }
        }
        setProgress(std::min(0.99f, static_cast<float>(intensity) / static_cast<float>(maxIntensity)));
    }

    exportFile.close();

    LINFO("Finished data export");

    setProgress(1.f);        
    reactivateProperties();


/*    

     


    std::string dataBaseDirectory = exportParentFolder_.get() + "/" + tgt::FileSystem::baseName(volumeInport_.getData()->getOrigin().getFilename());
    std::string segmentationDataDirectory = dataBaseDirectory + "/" + "segmentation/";
    std::string quantificationDataDirectory = dataBaseDirectory + "/" + "quantification/";
    std::string textBaseDirectory = dataBaseDirectory + "/";
    tgt::FileSystem::createDirectory(dataBaseDirectory);
    tgt::FileSystem::createDirectory(segmentationDataDirectory);
    tgt::FileSystem::createDirectory(quantificationDataDirectory);

    // Export general quantification info
    std::string quantInfoPath = tgt::FileSystem::cleanupPath(textBaseDirectory + "QuantificationInformation.csv");
    std::ofstream quantInfoExportFile(quantInfoPath.c_str());
    quantInfoExportFile << "Number of Voxels in Volume" << ";";
    quantInfoExportFile << "Number of Voxels in Outer Structure"<< ";";
    quantInfoExportFile << "Number of Voxels in Inner Structure"<< ";";
    quantInfoExportFile << "Channel A Data Channel"<< ";";
    quantInfoExportFile << "Channel B Data Channel"<< ";";
    quantInfoExportFile << "Channel C Data Channel"<< ";";
    quantInfoExportFile << "Channel A Max Intensity Value"<< ";";
    quantInfoExportFile << "Channel B Max Intensity Value"<< ";";
    quantInfoExportFile << "Channel C Max Intensity Value"<< ";";
    quantInfoExportFile << "Channel A Intensity Threshold Min Value"<< ";";
    quantInfoExportFile << "Channel A Intensity Threshold Max Value"<< ";";
    quantInfoExportFile << "Channel B Intensity Threshold Min Value"<< ";";
    quantInfoExportFile << "Channel B Intensity Threshold Max Value"<< ";";
    quantInfoExportFile << "Channel C Intensity Threshold Min Value"<< ";";
    quantInfoExportFile << "Channel C Intensity Threshold Max Value"<< "\n";
    
    quantInfoExportFile << numVoxelsInVolume_ << ";";
    quantInfoExportFile << numVoxelsInOuterStructure_ << ";";
    quantInfoExportFile << numVoxelsInInnerStructure_ << ";";
    quantInfoExportFile << quantificationChannelA_.get() << ";";
    quantInfoExportFile << quantificationChannelB_.get() << ";";
    quantInfoExportFile << quantificationChannelC_.get() << ";";
    quantInfoExportFile << quantificationMaxValueA_.get() << ";";
    quantInfoExportFile << quantificationMaxValueB_.get() << ";";
    quantInfoExportFile << quantificationMaxValueC_.get() << ";";
    quantInfoExportFile << quantificationAbsoluteThresholdA_.get().x << ";";
    quantInfoExportFile << quantificationAbsoluteThresholdA_.get().y << ";";
    quantInfoExportFile << quantificationAbsoluteThresholdB_.get().x << ";";
    quantInfoExportFile << quantificationAbsoluteThresholdB_.get().y << ";";
    quantInfoExportFile << quantificationAbsoluteThresholdC_.get().x << ";";
    quantInfoExportFile << quantificationAbsoluteThresholdC_.get().y << "\n";

    quantInfoExportFile.close();
    
    // Export Histograms Channel A
    std::string histAPath = tgt::FileSystem::cleanupPath(quantificationDataDirectory + "HistogramChannelA.csv");
    std::ofstream histAExportFile(histAPath.c_str());
    histAExportFile << "Data Value (Channel A)" << ";";
    histAExportFile << "In Volume" << ";";
    histAExportFile << "In Outer Structure" << ";";
    histAExportFile << "In Inner Structure" << "\n";

    for (size_t i = 0; i < channelAHistogramVolume_.size(); ++i) {
        histAExportFile << i << ";";
        histAExportFile << channelAHistogramVolume_.at(i) << ";";
        histAExportFile << channelAHistogramOuterStructure_.at(i) << ";";
        histAExportFile << channelAHistogramInnerStructure_.at(i) << "\n";
    }
    histAExportFile.close();

    // Export Histograms Channel B
    std::string histBPath = tgt::FileSystem::cleanupPath(quantificationDataDirectory + "HistogramChannelB.csv");
    std::ofstream histBExportFile(histBPath.c_str());
    histBExportFile << "Data Value (Channel B)" << ";";
    histBExportFile << "In Volume" << ";";
    histBExportFile << "In Outer Structure" << ";";
    histBExportFile << "In Inner Structure" << "\n";

    for (size_t i = 0; i < channelBHistogramVolume_.size(); ++i) {
        histBExportFile << i << ";";
        histBExportFile << channelBHistogramVolume_.at(i) << ";";
        histBExportFile << channelBHistogramOuterStructure_.at(i) << ";";
        histBExportFile << channelBHistogramInnerStructure_.at(i) << "\n";
    }
    histBExportFile.close();

    // Export Histograms Channel C
    std::string histCPath = tgt::FileSystem::cleanupPath(quantificationDataDirectory + "HistogramChannelC.csv");
    std::ofstream histCExportFile(histCPath.c_str());
    histCExportFile << "Data Value (Channel C)" << ";";
    histCExportFile << "In Volume" << ";";
    histCExportFile << "In Outer Structure" << ";";
    histCExportFile << "In Inner Structure" << "\n";

    for (size_t i = 0; i < channelCHistogramVolume_.size(); ++i) {
        histCExportFile << i << ";";
        histCExportFile << channelCHistogramVolume_.at(i) << ";";
        histCExportFile << channelCHistogramOuterStructure_.at(i) << ";";
        histCExportFile << channelCHistogramInnerStructure_.at(i) << "\n";
    }
    histCExportFile.close();

    // Export Thresholded Histograms Channel A
    std::string histAThreshPath = tgt::FileSystem::cleanupPath(quantificationDataDirectory + "ThresholdedHistogramChannelA.csv");
    std::ofstream histAThreshExportFile(histAThreshPath.c_str());
    histAThreshExportFile << "Data Value (Channel A)" << ";";
    histAThreshExportFile << "In Volume" << ";";
    histAThreshExportFile << "In Outer Structure" << ";";
    histAThreshExportFile << "In Inner Structure" << "\n";

    size_t aMin = static_cast<size_t>(quantificationAbsoluteThresholdA_.get().x);
    size_t aMax = std::min(static_cast<size_t>(quantificationAbsoluteThresholdA_.get().y) + 1, channelAHistogramVolume_.size());
    for (size_t i = aMin; i < aMax; ++i) {
        histAThreshExportFile << i << ";";
        histAThreshExportFile << channelAHistogramVolume_.at(i) << ";";
        histAThreshExportFile << channelAHistogramOuterStructure_.at(i) << ";";
        histAThreshExportFile << channelAHistogramInnerStructure_.at(i) << "\n";
    }
    histAThreshExportFile.close();

    // Export Thresholded Histograms Channel B
    std::string histBThreshPath = tgt::FileSystem::cleanupPath(quantificationDataDirectory + "ThresholdedHistogramChannelB.csv");
    std::ofstream histBThreshExportFile(histBThreshPath.c_str());
    histBThreshExportFile << "Data Value (Channel B)" << ";";
    histBThreshExportFile << "In Volume" << ";";
    histBThreshExportFile << "In Outer Structure" << ";";
    histBThreshExportFile << "In Inner Structure" << "\n";

    size_t bMin = static_cast<size_t>(quantificationAbsoluteThresholdB_.get().x);
    size_t bMax = std::min(static_cast<size_t>(quantificationAbsoluteThresholdB_.get().y) + 1, channelBHistogramVolume_.size());
    for (size_t i = bMin; i < bMax; ++i) {
        histBThreshExportFile << i << ";";
        histBThreshExportFile << channelBHistogramVolume_.at(i) << ";";
        histBThreshExportFile << channelBHistogramOuterStructure_.at(i) << ";";
        histBThreshExportFile << channelBHistogramInnerStructure_.at(i) << "\n";
    }
    histBThreshExportFile.close();
        
    // Export Thresholded Histograms Channel C
    std::string histCThreshPath = tgt::FileSystem::cleanupPath(quantificationDataDirectory + "ThresholdedHistogramChannelC.csv");
    std::ofstream histCThreshExportFile(histCThreshPath.c_str());
    histCThreshExportFile << "Data Value (Channel C)" << ";";
    histCThreshExportFile << "In Volume" << ";";
    histCThreshExportFile << "In Outer Structure" << ";";
    histCThreshExportFile << "In Inner Structure" << "\n";

    size_t cMin = static_cast<size_t>(quantificationAbsoluteThresholdC_.get().x);
    size_t cMax = std::min(static_cast<size_t>(quantificationAbsoluteThresholdC_.get().y) + 1, channelCHistogramVolume_.size());
    for (size_t i = cMin; i < cMax; ++i) {
        histCThreshExportFile << i << ";";
        histCThreshExportFile << channelCHistogramVolume_.at(i) << ";";
        histCThreshExportFile << channelCHistogramOuterStructure_.at(i) << ";";
        histCThreshExportFile << channelCHistogramInnerStructure_.at(i) << "\n";
    }
    histCThreshExportFile.close();

     // Export Correlation Data A-B
    std::string scatterABPath = tgt::FileSystem::cleanupPath(quantificationDataDirectory + "CorrelationChannelAB.csv");
    std::ofstream scatterABExportFile(scatterABPath.c_str());
    scatterABExportFile << "Channel A (percent of max)" << ";";
    scatterABExportFile << "Channel B (percent of max)" << ";";
    scatterABExportFile << "Num Voxels (Volume)" << ";";   
    scatterABExportFile << "Num Voxels (Outer Structure)" << ";";
    scatterABExportFile << "Num Voxels (Inner Structure)" << "\n";
    for (size_t i = 0; i < 1001; ++i) {
        for (size_t j = 0; j < 1001; ++j) {
            scatterABExportFile << (static_cast<float>(i) / 10.f) << ";";
            scatterABExportFile << (static_cast<float>(j) / 10.f) << ";";
            scatterABExportFile << scatterABVolume_.at(i*1001 + j) << ";";
            scatterABExportFile << scatterABOuterStructure_.at(i*1001 + j) << ";";
            scatterABExportFile << scatterABInnerStructure_.at(i*1001 + j) << "\n";
        }
    }
    scatterABExportFile.close();

    // Export Correlation Data C-B
    std::string scatterCBPath = tgt::FileSystem::cleanupPath(quantificationDataDirectory + "CorrelationChannelCB.csv");
    std::ofstream scatterCBExportFile(scatterCBPath.c_str());
    scatterCBExportFile << "Channel C (percent of max)" << ";";
    scatterCBExportFile << "Channel B (percent of max)" << ";";
    scatterCBExportFile << "Num Voxels (Volume)" << ";";   
    scatterCBExportFile << "Num Voxels (Outer Structure)" << ";";
    scatterCBExportFile << "Num Voxels (Inner Structure)" << "\n";
    for (size_t i = 0; i < 1001; ++i) {
        for (size_t j = 0; j < 1001; ++j) {
            scatterCBExportFile << (static_cast<float>(i) / 10.f) << ";";
            scatterCBExportFile << (static_cast<float>(j) / 10.f) << ";";
            scatterCBExportFile << scatterCBVolume_.at(i*1001 + j) << ";";
            scatterCBExportFile << scatterCBOuterStructure_.at(i*1001 + j) << ";";
            scatterCBExportFile << scatterCBInnerStructure_.at(i*1001 + j) << "\n";
        }
    }
    scatterCBExportFile.close();

    //Compute Thresholds for percent-of-max scatter
    size_t scatterAThresholdMin = static_cast<size_t>(tgt::iround(quantificationRelativeThresholdA_.get().x * 10.f));
    size_t scatterAThresholdMax = std::min(static_cast<size_t>(tgt::iround(quantificationRelativeThresholdA_.get().y * 10.f)), static_cast<size_t>(1000));
    size_t scatterBThresholdMin = static_cast<size_t>(tgt::iround(quantificationRelativeThresholdB_.get().x * 10.f));
    size_t scatterBThresholdMax = std::min(static_cast<size_t>(tgt::iround(quantificationRelativeThresholdB_.get().y * 10.f)), static_cast<size_t>(1000));
    size_t scatterCThresholdMin = static_cast<size_t>(tgt::iround(quantificationRelativeThresholdC_.get().x * 10.f));
    size_t scatterCThresholdMax = std::min(static_cast<size_t>(tgt::iround(quantificationRelativeThresholdC_.get().y * 10.f)), static_cast<size_t>(1000));

    // Export Thresholded Correlation Data A-B
    std::string scatterABThreshPath = tgt::FileSystem::cleanupPath(quantificationDataDirectory + "ThresholdedCorrelationChannelAB.csv");
    std::ofstream scatterABThreshExportFile(scatterABThreshPath.c_str());
    scatterABThreshExportFile << "Channel A (percent of max)" << ";";
    scatterABThreshExportFile << "Channel B (percent of max)" << ";";
    scatterABThreshExportFile << "Num Voxels (Volume)" << ";";   
    scatterABThreshExportFile << "Num Voxels (Outer Structure)" << ";";
    scatterABThreshExportFile << "Num Voxels (Inner Structure)" << "\n";
    for (size_t i = scatterAThresholdMin; i <= scatterAThresholdMax; ++i) {
        for (size_t j = scatterBThresholdMin; j <= scatterBThresholdMax; ++j) {
            scatterABThreshExportFile << (static_cast<float>(i) / 10.f) << ";";
            scatterABThreshExportFile << (static_cast<float>(j) / 10.f) << ";";
            scatterABThreshExportFile << scatterABVolume_.at(i*1001 + j) << ";";
            scatterABThreshExportFile << scatterABOuterStructure_.at(i*1001 + j) << ";";
            scatterABThreshExportFile << scatterABInnerStructure_.at(i*1001 + j) << "\n";
        }
    }
    scatterABThreshExportFile.close();

    // Export Thresholded Correlation Data C-B
    std::string scatterCBThreshPath = tgt::FileSystem::cleanupPath(quantificationDataDirectory + "ThresholdedCorrelationChannelCB.csv");
    std::ofstream scatterCBThreshExportFile(scatterCBThreshPath.c_str());
    scatterCBThreshExportFile << "Channel C (percent of max)" << ";";
    scatterCBThreshExportFile << "Channel B (percent of max)" << ";";
    scatterCBThreshExportFile << "Num Voxels (Volume)" << ";";   
    scatterCBThreshExportFile << "Num Voxels (Outer Structure)" << ";";
    scatterCBThreshExportFile << "Num Voxels (Inner Structure)" << "\n";
    for (size_t i = scatterCThresholdMin; i <= scatterCThresholdMax; ++i) {
        for (size_t j = scatterBThresholdMin; j <= scatterBThresholdMax; ++j) {
            scatterCBThreshExportFile << (static_cast<float>(i) / 10.f) << ";";
            scatterCBThreshExportFile << (static_cast<float>(j) / 10.f) << ";";
            scatterCBThreshExportFile << scatterCBVolume_.at(i*1001 + j) << ";";
            scatterCBThreshExportFile << scatterCBOuterStructure_.at(i*1001 + j) << ";";
            scatterCBThreshExportFile << scatterCBInnerStructure_.at(i*1001 + j) << "\n";
        }
    }
    scatterCBThreshExportFile.close();

    */
}


void OctreeSegmentationQuantification::deactivateProperties() {
    simMode_.setReadOnlyFlag(true);
    computeQuantification_.setReadOnlyFlag(true); 

    exportFile_.setReadOnlyFlag(true);
    exportButton_.setReadOnlyFlag(true);
}

void OctreeSegmentationQuantification::reactivateProperties() {
    simMode_.setReadOnlyFlag(false);
    computeQuantification_.setReadOnlyFlag(false); 

    exportFile_.setReadOnlyFlag(false);
    exportButton_.setReadOnlyFlag(false);
}

}   // namespace
