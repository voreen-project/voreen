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

#ifndef VRN_NUCLEIPOSITIONQUANTIFICATION_H
#define VRN_NUCLEIPOSITIONQUANTIFICATION_H

#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/ports/geometryport.h"

//#include "voreen/core/properties/optionproperty.h"
//#include "voreen/core/properties/floatproperty.h"
//#include "voreen/core/properties/intproperty.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/buttonproperty.h"


//#include "../util/connectedcomponentqueue.h"

namespace voreen {

/**
 * Determines positions of points corresponding to cell nuclei centroids with respect to a foreground segmentation.
 * Outputs two groups: within the foreground and outside of the foreground.
 *
 * For the nuclei centroids outside the region of interest, the distance to the region of interest is also computed.
 * The results can also be written to a CSV file. 
 */
class VRN_CORE_API NucleiPositionQuantification : public Processor {

public:
    NucleiPositionQuantification();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "NucleiPositionQuantification"; }
    virtual std::string getCategory() const   { return "Quantification";      }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING;       }

    virtual bool isReady() const;

    virtual bool usesExpensiveComputation() const { return true; }
    virtual bool isEndProcessor() const { return true; }

protected:


    virtual void setDescriptions() {
        setDescription("Receives an original volume, a PointSegmentList containing cell nuclei centroids, and a (potentially downsampled) binary volume (= region of interest). The input points are separated into two output groups, one within the region of interest, and one outside of it. The resolution of the binary volume can be different (e.g., lower), but the volumes have to contain the same region, since a coordinate transformation is computed based on the volume dimensions of the two input images. \
                \n\n For points outside the region of interest, the distance to the roi is computed. This assumes that the spacing in the segmentation volume is uniform!!!");

        inportImage_.setDescription("Orginal input volume");
        inportSegmentation_.setDescription("Binary foreground volume (preferably uint8), spacing has to be uniform for the distance computation to be correct.");
        centroidInport_.setDescription("Points representing the nuclei centroid positions in voxel coordinates of the original volume.");
    }

    virtual void process();

    virtual void performQuantification();

    /**
     * Sets all settings properties as read-only.
     */
    virtual void deactivateSettings();

    /**
     * Re-activates properties after setting them to read-only
     */
    virtual void reactivateSettings();

    //virtual void adjustPropertyVisibility();

    VolumePort inportImage_;            ///< input of original image 
    VolumePort inportSegmentation_;     ///< input of binary image

    GeometryPort centroidInport_;       ///< centroids of detected cell nuclei as PointSegmentList

    // output of two point lists
    GeometryPort insideOutport_;
    GeometryPort outsideOutport_;

    //BoolProperty exportAsCSV_;
    FileDialogProperty csvFile_;
    ProgressProperty progressProperty_;
    ButtonProperty compute_;

    bool forceComputation_;

    static const std::string loggerCat_;
};

}   //namespace

#endif // VRN_NUCLEIPOSITIONQUANTIFICATION_H
