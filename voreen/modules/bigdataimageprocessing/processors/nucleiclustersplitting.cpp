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

#include "nucleiclustersplitting.h"

#include "voreen/core/voreenapplication.h"

#include "../util/clustersplittingthread.h"

#include "tgt/stopwatch.h"

#include "voreen/core/datastructures/volume/volumedisk.h"

namespace voreen {

const std::string NucleiClusterSplitting::loggerCat_("voreen.nucleusdetection.NucleiClusterSplitting");

NucleiClusterSplitting::NucleiClusterSplitting()
    : VolumeProcessor()
    , inportImage_(Port::INPORT, "volumehandle.input.image", "Volume Input Image")
    , inportLabels_(Port::INPORT, "volumehandle.input.cca", "Volume Input Connected Component Labeling")
    //, outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", true, voreen::Processor::VALID)
    , centroidOutport_(Port::OUTPORT, "centroid.output", "Nuclei Centroids Output", true, voreen::Processor::VALID)
    , connectedComponentsDescription_("connectedComponentsDescription", "Connected Components List (CSV File)", "Open CSV input file", VoreenApplication::app()->getUserDataPath(), "Commaseparated Files(*.csv);;Textfile (*.txt);;All TextFiles (*.csv *.txt)", FileDialogProperty::OPEN_FILE,Processor::VALID)
    , kernelSize_("KernelSize", "Kernel Size")
    , minSeedRank_("MinSeedRank", "Min Seed Rank", 1, 1, 65535)
    , deriveSigmaFromKernelSize_("DeriveSigmaFromKernelSize", "Derive sigma from kernel size", false)
    , sigma_("sigma","sigma", 0.6f, 0.1f, 20.0f)
    , expandMarkers_("expandMarkers", "Expand Markers", true)
    //, enableProcessing_("enableProcessing", "Enable")
    , smoothForSeedDetection_("smoothforseeddetection", "Smooth before Seeddetection", true)
    , maskBeforeSmoothing_("maskbeforesmoothing", "Mask before Smoothing", true)
    , progressProperty_("progressProperty", "Progress")
    , compute_("compute", "Compute")
    , forceComputation_(false)
    , threadComponentCounter_(0)
{
    addPort(inportImage_);
    addPort(inportLabels_);
    //addPort(outport_);
    addPort(centroidOutport_);

    addProperty(connectedComponentsDescription_);

    // properties for pre-processing (smoothing)
    smoothForSeedDetection_.setGroupID("preprocessing");
    addProperty(smoothForSeedDetection_);

    kernelSize_.addOption("3",  "3x3x3",    3);
    kernelSize_.addOption("5",  "5x5x5",    5);
    kernelSize_.addOption("7",  "7x7x7",    7);
    kernelSize_.addOption("9",  "9x9x9",    9);
    //kernelSize_.addOption("11", "11x11x11", 11);
    //kernelSize_.addOption("13", "13x13x13", 13);
    //kernelSize_.addOption("15", "15x15x15", 15);
    kernelSize_.selectByKey("5");
    kernelSize_.setGroupID("preprocessing");
    addProperty(kernelSize_);

    deriveSigmaFromKernelSize_.setGroupID("preprocessing");
    addProperty(deriveSigmaFromKernelSize_);
    
    sigma_.setGroupID("preprocessing");
    addProperty(sigma_);

    maskBeforeSmoothing_.setGroupID("preprocessing");
    addProperty(maskBeforeSmoothing_);

    setPropertyGroupGuiName("preprocessing", "Pre-processing (Gaussian Smoothing)");

    ON_CHANGE(kernelSize_, NucleiClusterSplitting, adjustPropertyVisibility);
    ON_CHANGE(deriveSigmaFromKernelSize_, NucleiClusterSplitting, adjustPropertyVisibility);
    ON_CHANGE(smoothForSeedDetection_, NucleiClusterSplitting, adjustPropertyVisibility);

    // properties for marker extension
    expandMarkers_.setGroupID("markerexpansion");
    addProperty(expandMarkers_);
    setPropertyGroupGuiName("markerexpansion", "Marker Expansion");
    addProperty(minSeedRank_);

    // progress and compute button
    addProperty(progressProperty_);
    addProgressBar(&progressProperty_);

    ON_CHANGE(compute_, NucleiClusterSplitting, performClusterSplitting);
    addProperty(compute_);
}

Processor* NucleiClusterSplitting::create() const {
    return new NucleiClusterSplitting();
}

void NucleiClusterSplitting::process() {
    // reset status
    setProgress(0.f);

    if (forceComputation_) {
        forceComputation_ = false;
        segment();
    }
}

void NucleiClusterSplitting::segment() {
    tgtAssert(inportImage_.hasData() && inportLabels_.hasData(), "Inport has not data");

    centroidOutport_.clear();
    setProgress(0.f);
    threadComponentCounter_ = 0;
    std::vector< std::vector<tgt::vec3> > emptyList;
    nucleiSegmentList_.setData(emptyList);

    // check input data
    const VolumeBase* image = inportImage_.getData();
    const VolumeBase* labels = inportLabels_.getData();

    if (!image || !image->hasRepresentation<VolumeDisk>()) {
        LERROR("No input or no disk representation for image input");
        return;
    }

    if (image->getNumChannels() != 1 || image->getFormat() != "uint16") {
        LERROR("wrong image format: expected single-channel uint16");
        return;
    }

    if (!labels || !labels->hasRepresentation<VolumeDisk>()) {
        LERROR("No input or no disk representation for image input");
        return;
    }

    if (labels->getNumChannels() != 1 || labels->getFormat() != "uint32") {
        LERROR("wrong connected component labeling format: expected single-channel uint32");
        return;
    }

    if (image->getDimensions() != labels->getDimensions()) {
        LERROR("Dimensions of input image and label volume are different");
        return;
    }

    const VolumeDisk* imageDisk = image->getRepresentation<VolumeDisk>();
    const VolumeDisk* labelDisk = labels->getRepresentation<VolumeDisk>();

    // TODO: check for VolumeDiskHDF5 or for functionality to load bricks

    // check and load CSV file 
    std::vector<ConnectedComponentQueue::ConnectedComponent> connectedComponents = readDescriptionFile();
    if (connectedComponents.empty()) {
        LERROR("Empty list of connected components - not processing");
        return;
    }

    // check if llf <= urb for all connected components and if urb is within volume dimensions
    for (auto& i : connectedComponents) {
        if (tgt::hor(tgt::greaterThan(i.llf_, i.urb_))) {
            LERROR("LLF > URB in connected component list");
            return;
        }
        if (tgt::hor(tgt::greaterThanEqual(i.urb_, image->getDimensions()))) {
            LERROR("Connected components in the list are outside of volume dimensions!");
            return;
        }
    }

    // deactivate all settings properties
    deactivateSettings();

    tgt::Stopwatch stopwatch;
    LINFO("Starting nuclei cluster splitting...");
    stopwatch.start();

    // get number of CPU cores to start the right amount of threads
    size_t numCores = static_cast<size_t>(boost::thread::hardware_concurrency());
    if (numCores == 0) {
        numCores = 1;
        LWARNING("Could not detect number of CPU cores. Using (slow) single threaded quantification.");
    }
    else 
        LINFO("Starting " << numCores << " threads for cluster splitting...");

    // queue for connected components to be scheduled for splitting (limited to stay within a reasonable amount of main memory)
    ConnectedComponentQueue componentQueue(numCores);

    // mutex for preventing race conditions when threads write the results
    boost::mutex globalValueMutex;

    // now create the threads and start their execution
    std::vector<ClusterSplittingThread*> workerThreads;
    for (size_t threadCount = 0; threadCount < numCores; ++threadCount) {
        ClusterSplittingThread* thread = new ClusterSplittingThread(this, componentQueue, globalValueMutex, smoothForSeedDetection_.get(), kernelSize_.getValue(), sigma_.get(), maskBeforeSmoothing_.get(), expandMarkers_.get(), static_cast<size_t>(minSeedRank_.get()));
        workerThreads.push_back(thread);
        thread->run();
    }
    
    tgt::svec3 volDim = image->getDimensions();

    // iterate over connected components
    for (auto component : connectedComponents) {

        // load original image and label block from files (one additional border voxel if possible)
        component.llf_ = tgt::svec3(component.llf_.x == 0 ? 0 : component.llf_.x - 1, 
                                    component.llf_.y == 0 ? 0 : component.llf_.y - 1,
                                    component.llf_.z == 0 ? 0 : component.llf_.z - 1);
        component.urb_ = tgt::svec3(component.urb_.x == (volDim.x - 1) ? (volDim.x - 1) : component.urb_.x + 1, 
                                    component.urb_.y == (volDim.y - 1) ? (volDim.y - 1) : component.urb_.y + 1,
                                    component.urb_.z == (volDim.z - 1) ? (volDim.z - 1) : component.urb_.z + 1);

        tgt::svec3 dim = component.urb_ - component.llf_ + tgt::svec3::one;

        // set spacing and offset (leave offset at zero, since it is not required)
        component.spacing_ = image->getSpacing();
        //component.offset_ = tgt::vec3(0.f);

        VolumeRAM* imageRAM = 0; 
        VolumeRAM* labelRAM = 0;

        try {
            imageRAM = imageDisk->loadBrick(component.llf_, dim);
        }
        catch (tgt::Exception& e) {
            for (size_t t = 0; t < workerThreads.size(); ++t) {
                ClusterSplittingThread* thread = workerThreads.at(t);
                thread->interruptAndJoin(); // this will write back the thread's results
                delete thread;
            }
            LERROR("Could not load brick [" << component.llf_ << ", " << component.urb_ << "] from image file");
            setProgress(0.f);
            reactivateSettings();
            return;
        }
        component.image_ = dynamic_cast<VolumeRAM_UInt16*>(imageRAM); 
        if (!component.image_) {
            for (size_t t = 0; t < workerThreads.size(); ++t) {
                ClusterSplittingThread* thread = workerThreads.at(t);
                thread->interruptAndJoin(); // this will write back the thread's results
                delete thread;
            }
            delete imageRAM;
            LERROR("Could not convert brick [" << component.llf_ << ", " << component.urb_ << "] from image file to uint16");
            setProgress(0.f);
            reactivateSettings();
            return;
        }

        try {
            labelRAM = labelDisk->loadBrick(component.llf_, dim);
        }
        catch (tgt::Exception& e) {
            for (size_t t = 0; t < workerThreads.size(); ++t) {
                ClusterSplittingThread* thread = workerThreads.at(t);
                thread->interruptAndJoin(); // this will write back the thread's results
                delete thread;
            }
            delete imageRAM;
            LERROR("Could not load brick [" << component.llf_ << ", " << component.urb_ << "] from connected components labeling file");
            setProgress(0.f);
            reactivateSettings();
            return;
        }
        component.labels_ = dynamic_cast<VolumeRAM_UInt32*>(labelRAM); 
        if (!component.labels_) {
            for (size_t t = 0; t < workerThreads.size(); ++t) {
                ClusterSplittingThread* thread = workerThreads.at(t);
                thread->interruptAndJoin(); // this will write back the thread's results
                delete thread;
            }
            delete imageRAM;
            LERROR("Could not convert brick [" << component.llf_ << ", " << component.urb_ << "] from connected component labeling file to uint32");
            setProgress(0.f);
            reactivateSettings();
            return;
        }

        // push the connected component on the queue
        componentQueue.push(component);

        // update progress bar after each block
        float p = std::min(static_cast<float>(threadComponentCounter_) / static_cast<float>(connectedComponents.size()), 0.99f);
        setProgress(p);
    }

    // wait for all threads to finish and then delete them
    for (size_t t = 0; t < workerThreads.size(); ++t) {
        ClusterSplittingThread* thread = workerThreads.at(t);
        thread->interruptAndJoin(); // this will write back the thread's results
        delete thread;
        // update progress bar after each thread is finished
        float p = std::min(static_cast<float>(threadComponentCounter_) / static_cast<float>(connectedComponents.size()), 0.99f);
        setProgress(p);
    }

    stopwatch.stop();
    LINFO("Cluster splitting finished. Runtime: " << static_cast<double>(stopwatch.getRuntime()) / 1000.0 << " seconds");

    // re-activate the settings
    reactivateSettings();

    // output number of detected cell nuclei
    LINFO("Detected " << nucleiSegmentList_.getNumPoints() << " cell nuclei in total");

    // set outout list
    centroidOutport_.setData(&nucleiSegmentList_, false);

    setProgress(1.f);
}

void NucleiClusterSplitting::deactivateSettings() {
    connectedComponentsDescription_.setReadOnlyFlag(true);
    smoothForSeedDetection_.setReadOnlyFlag(true);
    kernelSize_.setReadOnlyFlag(true);
    deriveSigmaFromKernelSize_.setReadOnlyFlag(true);
    sigma_.setReadOnlyFlag(true);
    expandMarkers_.setReadOnlyFlag(true);
    compute_.setReadOnlyFlag(true);
    minSeedRank_.setReadOnlyFlag(true);
}

void NucleiClusterSplitting::reactivateSettings() {
    // re-activate all settings properties 
    connectedComponentsDescription_.setReadOnlyFlag(false);
    smoothForSeedDetection_.setReadOnlyFlag(false);
    kernelSize_.setReadOnlyFlag(false);
    deriveSigmaFromKernelSize_.setReadOnlyFlag(false);
    sigma_.setReadOnlyFlag(false);
    expandMarkers_.setReadOnlyFlag(false);
    compute_.setReadOnlyFlag(false);
    minSeedRank_.setReadOnlyFlag(false);
    adjustPropertyVisibility();
}

std::vector<ConnectedComponentQueue::ConnectedComponent> NucleiClusterSplitting::readDescriptionFile() {
    std::vector<ConnectedComponentQueue::ConnectedComponent> result;
    std::string filename = connectedComponentsDescription_.get();
    if (filename.empty()) {
        LERROR("Filename ist empty.");
        return result;
    }    

    std::ifstream inFile;
    inFile.open(filename.c_str(), std::ios::in);
    if (inFile.fail()) {
        LERROR(" Unable to open connected components description file: " << filename);
        return result;
    }

    LINFO("Reading connected components description...");
    // read every line in the file
    size_t currentLine = 0;
    std::string line;
    while(inFile.good() && std::getline(inFile, line)) {
        currentLine++;
        //ignore empty lines
        if (line.empty())
            continue;
        // split the line at the delimiter        
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<size_t> row;
        while(std::getline(lineStream, cell, ',')) {
            size_t value;
            std::stringstream cellStream(cell);
            cellStream >> value;            
            row.push_back(value);
        }
        // check for right amount of values
        if (row.size() != 8) {
            LERROR("Error in line " << currentLine << " of file " << filename << ": wrong number of values (8 expected)");
            inFile.close();
            std::vector<ConnectedComponentQueue::ConnectedComponent> errorResult;
            return errorResult;
        }

        ConnectedComponentQueue::ConnectedComponent c;
        c.id_ = row.at(0);
        c.volume_ = row.at(1);
        c.llf_ = tgt::svec3(row.at(2), row.at(3), row.at(4));
        c.urb_ = tgt::svec3(row.at(5), row.at(6), row.at(7));
        c.image_ = 0;
        c.labels_ = 0;
        c.spacing_ = tgt::vec3::one;
        c.offset_ = tgt::vec3::zero;

        result.push_back(c);
    }

    inFile.close();

    LINFO("Read " << result.size() << " connected components from description file");
    return result;
}

void NucleiClusterSplitting::adjustPropertyVisibility() {
    if (!isInitialized())
        return;
    sigma_.setVisibleFlag(smoothForSeedDetection_.get());
    maskBeforeSmoothing_.setVisibleFlag(smoothForSeedDetection_.get());
    deriveSigmaFromKernelSize_.setVisibleFlag(smoothForSeedDetection_.get());
    kernelSize_.setVisibleFlag(smoothForSeedDetection_.get());
    sigma_.setReadOnlyFlag(deriveSigmaFromKernelSize_.get());
    if (deriveSigmaFromKernelSize_.get()) {
        size_t kernelRadius = static_cast<size_t>(kernelSize_.getValue()) / 2;
        float sigma = static_cast<float>(kernelRadius) / 2.5f; 
        sigma_.set(sigma);
    }
}

void NucleiClusterSplitting::performClusterSplitting() {
    if (isInitialized() && isReady()) {
        forceComputation_ = true;
        invalidate();
    }
}


}   // namespace
