/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "volumeseriessource.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/processors/processorwidgetfactory.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/io/textfilereader.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include "tgt/filesystem.h"

namespace voreen {

const std::string VolumeSeriesSource::loggerCat_("voreen.deprecated.VolumeSeriesSource");

VolumeSeriesSource::VolumeSeriesSource()
    : Processor()
    , volume_(0)
    , filename_("seriesfile", "Volume Series File", "Select Volume Series File",
                VoreenApplication::app()->getUserDataPath("volumes"), "Volume Series Files (*.sdat)",
                FileDialogProperty::OPEN_FILE, Processor::INVALID_RESULT)
    , step_("step", "Time Step", 0, 0, 1000)
    , outport_(Port::OUTPORT, "volumehandle.volumehandle", "volumehandle.volumehandle", false)
    , needUpload_(false)
{
    addProperty(filename_);
    filename_.onChange(MemberFunctionCallback<VolumeSeriesSource>(this, &VolumeSeriesSource::openSeries));

    addProperty(step_);
    step_.setTracking(false);
    step_.onChange(MemberFunctionCallback<VolumeSeriesSource>(this, &VolumeSeriesSource::loadStep));

    addPort(outport_);
}

VolumeSeriesSource::~VolumeSeriesSource() {
    delete volume_;
}

Processor* VolumeSeriesSource::create() const {
    return new VolumeSeriesSource();
}

void VolumeSeriesSource::process() {
    if (needUpload_) {
        if(volume_)
            volume_->makeRepresentationExclusive<VolumeRAM>();
        needUpload_ = false;
    }
}

void VolumeSeriesSource::initialize() {
    Processor::initialize();
    outport_.setData(volume_, false);
    loadStep();
}

void VolumeSeriesSource::loadStep() {
    if(!volume_)
        return;

    VolumeRAM* v = volume_->getWritableRepresentation<VolumeRAM>();
    if (!v)
        return;

    int step = step_.get();
    std::string filename;

    if (files_.size() > 0 && step < static_cast<int>(files_.size()))
        filename = files_[step];
    else
        return;

    std::ifstream f(filename.c_str(), std::ios_base::binary);
    if (!f) {
        LERROR("Could not open file: " << filename);
        return;
    }

    // read the volume as a raw blob, should be fast
    LINFO("Loading raw file " << filename);
    f.read(reinterpret_cast<char*>(v->getData()), v->getNumBytes());
    if (!f.good()) {
        LERROR("Reading from file failed: " << filename);
        return;
    }

    // Special handling for float volumes: normalize values to [0.0; 1.0]
    VolumeRAM_Float* vf = dynamic_cast<VolumeRAM_Float*>(v);
    if (vf && spreadMin_ != spreadMax_) {
        const size_t n = vf->getNumVoxels();

        // use spread values if available
        if (spreadMin_ != spreadMax_) {
            const float d = spreadMax_ - spreadMin_;
            float* voxel = vf->voxel();
            for (size_t i = 0; i < n; ++i)
                voxel[i] = (voxel[i] - spreadMin_) / d;
        } else {
            LINFO("Normalizing float data to [0.0; 1.0]. "
                  << "This might not be what you want, better define 'Spread: <min> <max>' in the .sdat file.");
            const float d = vf->max() - vf->min();
            const float p = vf->min();
            float* voxel = vf->voxel();
            for (size_t i = 0; i < n; ++i)
                voxel[i] = (voxel[i] - p) / d;
        }
        vf->invalidate();
    }

    needUpload_ = true;
    invalidate();
}

void VolumeSeriesSource::openSeries() {
    if(volume_)
        delete volume_;
    volume_ = 0;

    try {
        TextFileReader reader(filename_.get());
        if (!reader)
            throw tgt::FileNotFoundException("Reading sdat file failed", filename_.get());

        std::string type;
        std::istringstream args;
        std::string format, model;
        tgt::ivec3 resolution(0);
        tgt::vec3 sliceThickness(1.f);
        spreadMin_ = 0.f;
        spreadMax_ = 0.f;

        files_.clear();
        std::string blockfile;

        while (reader.getNextLine(type, args, false)) {
            if (type == "ObjectFileName:") {
                std::string s;
                args >> s;
                LDEBUG(type << " " << s);
                files_.push_back(tgt::FileSystem::dirName(filename_.get()) + "/" + s);
            } else if (type == "Resolution:") {
                args >> resolution[0];
                args >> resolution[1];
                args >> resolution[2];
                LDEBUG(type << " " << resolution[0] << " x " <<
                       resolution[1] << " x " << resolution[2]);
            } else if (type == "SliceThickness:") {
                args >> sliceThickness[0] >> sliceThickness[1] >> sliceThickness[2];
                LDEBUG(type << " " << sliceThickness[0] << " "
                       << sliceThickness[1] << " " << sliceThickness[2]);
            } else if (type == "Format:") {
                args >> format;
                LDEBUG(type << " " << format);
            } else if (type == "ObjectModel:") {
                args >> model;
                LDEBUG(type << " " << model);
            } else if (type == "Spread:") {
                args >> spreadMin_ >> spreadMax_;
                LDEBUG(type << " " << spreadMin_ << " " << spreadMax_);
            } else {
                LWARNING("Unknown type: " << type);
            }

            if (args.fail()) {
                LERROR("Format error");
            }
        }

        VolumeRAM* vol = 0;

        if (format == "UCHAR") {
            if(model == "I")
                vol = new VolumeRAM_UInt8(resolution);
            if(model == "RGB")
                vol = new VolumeRAM_3xUInt8(resolution);
            if(model == "RGBA")
                vol = new VolumeRAM_4xUInt8(resolution);
        }
        else if (format == "USHORT") {
            if(model == "I")
                vol = new VolumeRAM_UInt16(resolution);
            else if(model == "RGB")
                vol = new VolumeRAM_3xUInt16(resolution);
            else if(model == "RGBA")
                vol = new VolumeRAM_4xUInt16(resolution);
        }
        else if (format == "FLOAT") {
            if(model == "I")
                vol = new VolumeRAM_Float(resolution);
            else if(model == "RGB")
                vol = new VolumeRAM_3xFloat(resolution);
            else if(model == "RGBA")
                vol = new VolumeRAM_4xFloat(resolution);
        }
        else {
            LERROR("Unsupported format: " << format);
            return;
        }

        volume_ = new Volume(vol, sliceThickness, tgt::vec3(0.0f));
        oldVolumePosition(volume_);

        if (step_.get() >= static_cast<int>(files_.size()))
            step_.set(0);
        step_.setMaxValue(static_cast<int>(files_.size()) - 1);
    }
    catch (tgt::FileException) {
        LERROR("Loading failed");
    }
    catch (std::bad_alloc) {
        LERROR("Out of memory");
    }

    loadStep();
    outport_.setData(volume_, false);
}

} // namespace
