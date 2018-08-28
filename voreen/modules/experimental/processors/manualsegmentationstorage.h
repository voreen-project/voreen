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

#ifndef VRN_MANUALSEGMENTATIONSTORAGE_H
#define VRN_MANUALSEGMENTATIONSTORAGE_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/genericcoprocessorport.h"

#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"


namespace voreen {

/**
 * This class stores a manually made segmentation. It stores a uint16_t volume which assignes
 * a id to each voxel in the volume that should be segmented. The id 0 has the semantic, that
 * this voxel is not segmented yet.
 *
 * A already stared segmentation can be loaded in again via a volume port. It need to have
 * exactly the same size as the volume to segment and need to be in uint16 format.
 *
 * One or multiple ManuelSegmentation processors should be connected to this processor with the
 * coprocessor port.
 *
 * The storages also controls (not user editable) the colors of the rendering. These are used
 * by the coprocessors for painting the brush and so on. The is also a Transfer Function for
 * linking with the same colors. It works only for the lower id's and is only useful with
 * nearest filtering.
 */
class ManualSegmentationStorage : public Processor{
public:
    /************************************************************************/
    /*                                                                      */
    /*                          Coprocessor Interface                       */
    /*                                                                      */
    /************************************************************************/
    /**
     * Tell the processor the id of any painting operation with this. So it can
     * maintain the highest used id in the current segmentation effiently.
     * @param id An id used in painting the volume. Any id is acceptable and
     *           it does NOT need to be the highest id.
     */
    void updateHighestId(int id);
    /**
     * Gets the volume with the current state of segmentation.
     *
     * After changing it EVERY id painted into it should be told this processor
     * with updateHighestId(int id) and invalidateSegmentation() NEEDS to be called
     */
    VolumeRAM_UInt16* getSegmentationVolume();
    /**
     *  Invalidates the Segmentation and invalidates the outport. This need to be
     * called after changing the volume.
     */
    void invalidateSegmentation();
    /**
     * Gets a color for a segmented id
     */
    tgt::vec4 getColorForId(int id);



public:
    ManualSegmentationStorage();
    ~ManualSegmentationStorage();

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "ManualSegmentationStorage"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    virtual Processor* create() const { return new ManualSegmentationStorage(); }

    virtual void process();

    virtual bool isReady() const;

    virtual void initialize();

    virtual void deinitialize();

protected:
    virtual void setDescriptions();
    /**
     * Callback for a new volume
     */
    void newVolume();
    /**
     * Recreates the volume. Either when a new volume to segment is avalible or when a old segmentation is loaded
     */
    void recreateVolume();

    /**
     * Clears the segmentation and sets all voxels to 0
     */
    void clearSegmentation();
private:
    VolumePort volumeInport_; ///< volume to segment
    VolumePort segmentationInport_; ///< optional volume with already existing (partial) segmentation

    VolumePort volumeOutport_; ///< The segmentation
    GenericCoProcessorPort<ManualSegmentationStorage> coprocessorPort_; ///< Port for ManualSegmentation processors.

    IntProperty hightestUsedId_; ///< the highest used id in the current segmentation

    TransFunc1DKeysProperty transferFunc_; ///< A transfunction for the segmentation to be linked

    ButtonProperty clearSegmentationButton_; ///< Button for clearing the segmentation

    bool shouldRecreateVolume_; ///< set of the volume containing the segmentation needs to be recreated
    bool shouldSetOutport_; ///< the segmentation volume changed and the outport need to be changed in the next process call

    VolumeRAM_UInt16* volumeram_; ///< VolumeRam for the segmentation
    Volume* volume_; ///< Volume for the segmentation

    std::vector<tgt::vec4> colors_; ///< colors for id's used in a round robin fashion


    static const std::string loggerCat_;
};

} // namespace voreen

#endif
