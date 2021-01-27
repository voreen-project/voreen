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

#ifndef VRN_OCTREESEGMENTATIONQUANTIFICATION_H
#define VRN_OCTREESEGMENTATIONQUANTIFICATION_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/textport.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/progressproperty.h"

#include "voreen/core/datastructures/octree/volumeoctree.h"

// plotting stuff
#include "../../../modules/plotting/datastructures/plotdata.h"
#include "../../../modules/plotting/ports/plotport.h"
#include "../../../modules/plotting/datastructures/plotcell.h"

// results
#include "../util/octreequantificationresults.h"

namespace voreen {

class VRN_CORE_API OctreeSegmentationQuantification : public Processor {
public:
    OctreeSegmentationQuantification();

    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Quantification";              }
    virtual std::string getClassName() const  { return "OctreeSegmentationQuantification";   }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

    virtual bool usesExpensiveComputation() const { return true; }

    virtual bool isReady() const;

protected:

    friend class OctreeQuantificationThread;

    enum SegmentationInterpolationMode {
        SIM_NEAREST,
        SIM_LINEAR,
        SIM_CUBIC
    };

    virtual void setDescriptions() {
        setDescription("This processor receives a data volume (for which an octree representation is required) and up to four regions of interest (ROIs) in form of (lower resolution) probability volumes (i.e., float volumes with the data range [0, 1]) for which a RAM representation is requested. The volume of each input ROI is computed in the octree data volume by locally upsampling the probability volumes (and thresholding them with 0.5). Additionally, the intensity histogram of each octree channel is computed for each of the ROIs separately.");
    }

    virtual void process();

    /**
     * Invalidates / removes the current quantification and sets the channel settings according to the new volume
     */
    void inputVolumeChanged();
    
    /**
     * Invalidates / removes the current quantification
     */
    void segmentationChanged();

    void computeQuantification();

    /// clears all data and ports
    void clearQuantificationData();

    /// exports quantification data to a CSV file  
    void exportData();

    void deactivateProperties();
    void reactivateProperties();


    /// inport for the (octree) volume for which the quantification should be computed
    VolumePort volumeInport_;

    // ports for probability volumes
    VolumePort roi1Port_;
    VolumePort roi2Port_;
    VolumePort roi3Port_;
    VolumePort roi4Port_;

    // analysis information
    //TextPort quantificationTextPort_;

    // plotting ports (TODO: rework)
/*
    // information as plot data
    PlotPort quantificationInfoPlot_;

    // histogram channel data 
    PlotPort histogramPlotPortChannelA_;
    PlotPort histogramPlotPortChannelB_;
    PlotPort histogramPlotPortChannelC_;

    // thresholded channel data
    PlotPort thresholdedHistogramPlotPortChannelA_;
    PlotPort thresholdedHistogramPlotPortChannelB_;
    PlotPort thresholdedHistogramPlotPortChannelC_;

    // histogram channel data 
    PlotPort relativeHistogramPlotPortChannelA_;
    PlotPort relativeHistogramPlotPortChannelB_;
    PlotPort relativeHistogramPlotPortChannelC_;

    // thresholded channel data
    PlotPort thresholdedRelativeHistogramPlotPortChannelA_;
    PlotPort thresholdedRelativeHistogramPlotPortChannelB_;
    PlotPort thresholdedRelativeHistogramPlotPortChannelC_;

    // surface plots
    PlotPort scatterABHistogramPort_;
    PlotPort scatterCBHistogramPort_;
    PlotPort thresholdedScatterABHistogramPort_;
    PlotPort thresholdedScatterCBHistogramPort_;

    // relative surface plots
    PlotPort scatterABRelativeHistogramPort_;
    PlotPort scatterCBRelativeHistogramPort_;
    PlotPort thresholdedScatterABRelativeHistogramPort_;
    PlotPort thresholdedScatterCBRelativeHistogramPort_;
*/

    /// determines how the voxel access into the segmentation is interpolated (linear or cubic)
    OptionProperty<SegmentationInterpolationMode> simMode_;
    
    ButtonProperty computeQuantification_; 

    ProgressProperty quantificationProgressProperty_;

    // export data properties
    FileDialogProperty exportFile_;
    ButtonProperty exportButton_;


    // quantification data which is stored within the processor
    OctreeQuantificationResults currentResults_;
        
    size_t threadVoxelCounter_;


    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif 
