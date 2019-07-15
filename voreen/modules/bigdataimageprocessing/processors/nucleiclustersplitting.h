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

#ifndef VRN_NUCLEICLUSTERSPLITTING_H
#define VRN_NUCLEICLUSTERSPLITTING_H

#include "voreen/core/processors/volumeprocessor.h"

#include "voreen/core/ports/geometryport.h"

#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"

#include "voreen/core/properties/buttonproperty.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

#include "../util/connectedcomponentqueue.h"

namespace voreen {

/**
 * Performs segmentation of touching cell nuclei in fluorescence microscopy images.
 */
class VRN_CORE_API NucleiClusterSplitting : public VolumeProcessor {

public:
    NucleiClusterSplitting();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "NucleiClusterSplitting"; }
    virtual std::string getCategory() const   { return "Volume Processing";      }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING;       }

    virtual bool usesExpensiveComputation() const { return true; }
    virtual bool isEndProcessor() const { return true; }

protected:

    friend class ClusterSplittingThread;

    virtual void setDescriptions() {
        setDescription("Segmentation of touching cell nuclei in ultramicroscopy images. Realizes the method proposed in Scherzinger et al. - \"Automated Segmentation of Immunostained Cell Nuclei in 3D Ultramicroscopy Images\". The processor does only perform the cluster splitting (i.e., expects a foreground segmentation in form of a connected component labeling) and processes each cluster separately by only loading the current cluster into main memory in a blockwise fashion. The processor does not output a labeled segmentation image, but a PointSegmentList of the cell nuclei centroids, where each segment corresponds to a single cluster.");

        inportImage_.setDescription("Original input image. Is expected to contain only one channel, data type must be uint16. A VolumeDiskHDF5 representation is expected to allow loading arbitrary blocks into main memory.");
        inportLabels_.setDescription("Imnage containing the connected component labeling of the foreground segmentation of the input image. Is expected to contain only one channel, data type must be uint32. The dimensions must match the input image. This data set is used both as a foreground mask and to identify each individual connected component.");
        connectedComponentsDescription_.setDescription("CSV file which contains the description of the connected components, i.e., their ID, volume, LLF and URB voxel. Is written out by the ConnectedComponentLabeling processor for large image data.");
        centroidOutport_.setDescription("Output of cell nuclei centroids as a PointSegmentList. Each segment corresponds to a single cluster and each point in the segment corresponds to a cell nucleus centroid within that cluster.");
    }

    virtual void process();

    /// callback for compute-button press
    virtual void performClusterSplitting();

    /// actual segmentation function
    virtual void segment();

    /**
     * Sets all settings properties as read-only (called when starting the cluster splitting).
     */
    virtual void deactivateSettings();

    /**
     * Re-activates properties after setting them to read-only (e.g., when cluster splitting is finished)
     */
    virtual void reactivateSettings();

    virtual void adjustPropertyVisibility();

    /**
     * Helper function which reads in the CSV files containing the connected components descriptions.
     * Returns an empty vector on failure.
     */
    virtual std::vector<ConnectedComponentQueue::ConnectedComponent> readDescriptionFile();

    VolumePort inportImage_;        ///< inport of original ultramicroscopy image (one channel, uint16, VolumeDiskHDF5 representation)
    VolumePort inportLabels_;   ///< inport of a connected component labeling image (one channel, uint32, VolumeDiskHDF5 representation, same dimensions as original image)

    GeometryPort centroidOutport_;  ///< centroids of detected cell nuclei as PointSegmentList (one segment per cluster)

    // TODO: optionally output a complete segmentation image?!
    //VolumePort outport_;

    FileDialogProperty connectedComponentsDescription_; ///< description file of the connected components, i.e., their IDs, volume, LLF, URB

    // properties for Gaussian smoothing (pre-processing)
    BoolProperty smoothForSeedDetection_;
    BoolProperty maskBeforeSmoothing_;
    IntOptionProperty kernelSize_;
    BoolProperty deriveSigmaFromKernelSize_;
    FloatProperty sigma_;

    BoolProperty expandMarkers_;    ///< perform the marker extension?

    ProgressProperty progressProperty_;
    ButtonProperty compute_;

    bool forceComputation_;

    static const std::string loggerCat_;

    size_t threadComponentCounter_;     ///< counter for keeping track of progress despite using several worker threads
    PointSegmentListGeometryVec3 nucleiSegmentList_; ///< result which is computed by worker threads

    IntProperty minSeedRank_;
};


}   //namespace

#endif // VRN_NUCLEICLUSTERSPLITTING_H
