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

#include "manualsegmentationstorage.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/utils/stringutils.h"

#include "tgt/textureunit.h"
#include "tgt/immediatemode/immediatemode.h"
#include "tgt/event/mouseevent.h"

namespace voreen {

const std::string ManualSegmentationStorage::loggerCat_("voreen.ManualSegmentationStorage");

ManualSegmentationStorage::ManualSegmentationStorage()
    : Processor()
    , volumeInport_      (Port::INPORT,  "volumeInport",       "Volume Inport")
    , segmentationInport_(Port::INPORT,  "segmentationInport", "Optional Presegmentation")
    , volumeOutport_     (Port::OUTPORT, "volumeOutport",      "Segmentation Outport")
    , coprocessorPort_            (Port::OUTPORT, "coprozessorPort",    "Coprozessor Port", true)
    , hightestUsedId_("highestUsedId", "Highest Used Id", 0, 0, std::numeric_limits<uint16_t>::max())
    , transferFunc_("transferFunction", "Transfer Function")
    , clearSegmentationButton_("clearSegmentationButton", "Clear Segmentation")
    , shouldRecreateVolume_(true)
    , shouldSetOutport_(true)
    , volumeram_(0)
    , volume_(0)
{
    addPort(volumeInport_);
    addPort(segmentationInport_);
    addPort(volumeOutport_);
    addPort(coprocessorPort_);

    addProperty(hightestUsedId_);
    hightestUsedId_.setReadOnlyFlag(true);

    addProperty(transferFunc_);
    transferFunc_.setReadOnlyFlag(true);

    addProperty(clearSegmentationButton_);

    ON_CHANGE(volumeInport_, ManualSegmentationStorage, newVolume);
    ON_CHANGE(segmentationInport_, ManualSegmentationStorage, newVolume);


    colors_.push_back(tgt::vec4(1.0f, 0.0f, 0.0f, 1.0f));
    colors_.push_back(tgt::vec4(0.0f, 1.0f, 0.0f, 1.0f));
    colors_.push_back(tgt::vec4(0.0f, 0.0f, 1.0f, 1.0f));
    colors_.push_back(tgt::vec4(1.0f, 1.0f, 0.0f, 1.0f));
    colors_.push_back(tgt::vec4(1.0f, 0.0f, 1.0f, 1.0f));
    colors_.push_back(tgt::vec4(0.0f, 1.0f, 1.0f, 1.0f));
    colors_.push_back(tgt::vec4(1.0f, 0.5f, 0.5f, 1.0f));
    colors_.push_back(tgt::vec4(0.5f, 1.0f, 0.5f, 1.0f));
    colors_.push_back(tgt::vec4(0.5f, 0.5f, 1.0f, 1.0f));

    ON_CHANGE(clearSegmentationButton_, ManualSegmentationStorage, clearSegmentation)
}

ManualSegmentationStorage::~ManualSegmentationStorage() {
}

bool ManualSegmentationStorage::isReady() const{
    return volumeInport_.isReady() && volumeOutport_.isReady();
}

void ManualSegmentationStorage::initialize() {
    Processor::initialize();
    float max = static_cast<float>(std::numeric_limits<uint16_t>::max());
    TransFunc1DKeys* keys = transferFunc_.get();
    std::vector<TransFuncMappingKey*> keyvec;
    int numkeys = 1000;

    keyvec.push_back(new TransFuncMappingKey(0, tgt::vec4::zero));
    for(int i = 0; i != numkeys; i++){
        float intensity = static_cast<float>(i+1)/(numkeys+1);
        tgt::vec4 color = 255.0f*colors_[i%colors_.size()];
        keyvec.push_back(new TransFuncMappingKey(intensity, color));
    }
    keyvec.push_back(new TransFuncMappingKey(1.0, tgt::vec4::zero));

    keys->setDomain(0, static_cast<float>(numkeys)/std::numeric_limits<uint16_t>::max());

    keys->setKeys(keyvec);

    transferFunc_.invalidate();
}

void ManualSegmentationStorage::deinitialize() {

    delete volume_;
    volume_ = 0;

    Processor::deinitialize();
}

void ManualSegmentationStorage::process(){
    if (shouldRecreateVolume_){
        recreateVolume();
        shouldRecreateVolume_ = false;
        shouldSetOutport_ = true;
    }
    if (shouldSetOutport_){
        volume_->invalidate();
        volumeOutport_.setData(volume_, false);
        shouldSetOutport_ = false;
    }
}

void ManualSegmentationStorage::setDescriptions(){
    setDescription("Processor that manages data storage while doing a manual segmentation.<br/>"
        "Multiple painters can be connected to work on the same volume.");
    volumeInport_.setDescription("The volume that should be segmented.");
    segmentationInport_.setDescription("An optional input with an segmentation that should be worked on further.");
    volumeOutport_.setDescription("The segmentation.");
    coprocessorPort_.setDescription("Port for the ManualSegmentation processors which do the painting.");
    hightestUsedId_.setDescription("The highest used id in the current segmentation.");
    transferFunc_.setDescription(
        "A transfer function for the segmented data.<br/>"
        "It matches the brush colors while painting and need to be used with nearest filtering.");
    clearSegmentationButton_.setDescription("Clears the current segmentation data.");
}



void ManualSegmentationStorage::newVolume(){
    shouldRecreateVolume_ = true;
}


void ManualSegmentationStorage::recreateVolume(){
    const VolumeBase* inputVolume = volumeInport_.getData();
    const VolumeBase* segmentationVolume = segmentationInport_.getData();
    tgt::svec3 dim = inputVolume->getDimensions();
    size_t size = tgt::hmul(dim);

    delete volume_;

    volumeram_ = new VolumeRAM_UInt16(dim);
    volume_ = new Volume(volumeram_, inputVolume->getSpacing(), inputVolume->getOffset());

    hightestUsedId_.set(0);

    if (segmentationVolume){
        const VolumeRAM * inram = segmentationVolume->getRepresentation<VolumeRAM>();
        const VolumeRAM_UInt16 * inram16 = dynamic_cast<const VolumeRAM_UInt16*>(inram);

        if (!inram16){
            LERROR("Input segmentation need to be Uint16 volume ram");
            std::fill_n(static_cast<uint16_t*>(volumeram_->getData()), size, 0);
            return;
        }
        if (segmentationInport_.getData() && volumeInport_.getData()->getDimensions() != segmentationInport_.getData()->getDimensions()){
            LERROR("Segmentation inport and volume inport need the same dimensions, ignores segmentatation volume");
            std::fill_n(static_cast<uint16_t*>(volumeram_->getData()), size, 0);
            return;
        }

        // Use the supplied volume
        const uint16_t *indata = reinterpret_cast<const uint16_t*>(inram16->getData());
        memcpy(volumeram_->getData(), indata, sizeof(uint16_t)*size);

        uint16_t max = *std::max_element(indata, indata+size);
        hightestUsedId_.set(max);
    }else{
        std::fill_n(static_cast<uint16_t*>(volumeram_->getData()), size, 0);
    }
}

void ManualSegmentationStorage::updateHighestId(int id)
{
    if (hightestUsedId_.get() < id)
        hightestUsedId_.set(id);
}

VolumeRAM_UInt16* ManualSegmentationStorage::getSegmentationVolume()
{
    return volumeram_;
}

void ManualSegmentationStorage::invalidateSegmentation()
{
    shouldSetOutport_ = true;
    invalidate();
}

tgt::vec4 ManualSegmentationStorage::getColorForId(int id)
{
    if (id == 0)
        return tgt::vec4(1.0);
    else
        return colors_[(id-1)%colors_.size()];
}

void ManualSegmentationStorage::clearSegmentation()
{
    if (volumeram_){
        size_t size = tgt::hmul(volumeram_->getDimensions());
        std::fill_n(static_cast<uint16_t*>(volumeram_->getData()), size, 0);
        shouldSetOutport_ = true;
        invalidate();
    }
}


} // namespace voreen
