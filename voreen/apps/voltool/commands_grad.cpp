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

#include "commands_grad.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/datastructures/volume/volumecollection.h"

#include "tgt/exception.h"

namespace voreen {

CommandGrad::CommandGrad() /*:
    Command("--grad", "", "Calculate gradients. Writes 32-bit dataset.\n\
\t\tsimple: neighborhood of 6\n\
\t\t26: neighborhood of 26\n\
\t\tsobel: neighborhood of 26 using sobel mask\n\
\t\tsobelic: sobel using intensity check",
"<[simple|26|sobel|sobelic] IN[8|12|16] OUT32>", 3) */
{
//    loggerCat_ += "." + name_;
}

bool CommandGrad::checkParameters(const std::vector<std::string>& parameters) {
    /*std::set<std::string> set;
    set.insert("simple");
    set.insert("26");
    set.insert("sobel");
    set.insert("sobelic");
    return (parameters.size() == 3) && isValueInSet(parameters[0], set); */
    return false;
}

bool CommandGrad::execute(const std::vector<std::string>& parameters) {
    /*VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    //load volume dataset
    const VolumeBase* sourceDataset_;
    Volume* targetDataset_ = 0;

    VolumeCollection* volumeCollection = serializer->read(parameters[1]);
    sourceDataset_ = volumeCollection->first();

    if (parameters[0] == "simple")
        targetDataset_ = calcGradientsCentralDifferences<uint8_t>(sourceDataset_);
    else if (parameters[0] == "26")
        targetDataset_ = calcGradients26(sourceDataset_);
    else if (parameters[0] == "sobel")
        targetDataset_ = calcGradientsSobel<uint8_t>(sourceDataset_);

    if (targetDataset_) {
        VolumeSerializerPopulator volLoadPop;
        const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
        serializer->write(parameters[2], targetDataset_);
        delete targetDataset_;
    }
    else {
        LERROR("Failed to calculate target dataset!");
        delete sourceDataset_;
    delete volumeCollection;
        return false;
    }

    delete sourceDataset_;
    delete volumeCollection; */
    return true;
}

//-----------------------------------------------------------------------------

CommandFilterGrad::CommandFilterGrad() /*:
    Command("--filtergrad", "", "Filter gradients", "<[simple|mid|weighted|weightedic] TIMES IN OUT>", 4) */
{
//    loggerCat_ += "." + name_;
}

bool CommandFilterGrad::checkParameters(const std::vector<std::string>& parameters) {
    /*std::set<std::string> set;
    set.insert("simple");
    set.insert("mid");
    set.insert("weighted");
    set.insert("weightedic");
    return ((parameters.size() == 4) && isValueInSet(parameters[0], set)); */
    return false;
}

bool CommandFilterGrad::execute(const std::vector<std::string>& parameters) {
    /*VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    //load volume dataset
    const VolumeBase* sourceDataset_;
    Volume* targetDataset_ = 0;

    VolumeCollection* volumeCollection = serializer->read(parameters[2]);
    sourceDataset_ = volumeCollection->first();

    if (parameters[0] == "simple")
        targetDataset_ = filterGradients(sourceDataset_);
    else if (parameters[0] == "mid")
        targetDataset_ = filterGradientsMid(sourceDataset_);
    else if (parameters[0] == "weighted") {
        int times;
        times = cast<int>(parameters[1]);
        for (int i=0; i<times; ++i) {
            targetDataset_ = filterGradientsWeighted(sourceDataset_, false);
            sourceDataset_ = targetDataset_;
        }
    }
    else if (parameters[0] == "weightedic") {
        int times;
        times = cast<int>(parameters[1]);
        for (int i=0; i<times; ++i) {
            targetDataset_ = filterGradientsWeighted(sourceDataset_, true);
            sourceDataset_ = targetDataset_;
        }
    }
    else {
        delete sourceDataset_;
        throw tgt::Exception("Unknown filter!");
    }

    if (targetDataset_) {
        VolumeSerializerPopulator volLoadPop;
        const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
        serializer->write(parameters[3], targetDataset_);
        delete targetDataset_;
    }
    else {
        LERROR("Failed to calculate target dataset!");
        delete sourceDataset_;
        return false;
    }

    delete sourceDataset_;
    delete volumeCollection; */
    return true;
}

}   //namespace voreen
