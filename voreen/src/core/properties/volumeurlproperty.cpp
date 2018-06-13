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

#include "voreen/core/properties/volumeurlproperty.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/properties/volumeinfoproperty.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/progressbar.h"
#include "voreen/core/voreenapplication.h"

#include "tgt/filesystem.h"

namespace voreen {

const std::string VolumeURLProperty::loggerCat_("voreen.VolumeURLProperty");

VolumeURLProperty::VolumeURLProperty(const std::string& id, const std::string& guiText,
                    const std::string& url, int invalidationLevel, Property::LevelOfDetail lod,
                    VoreenFileWatchListener::WatchMode watchMode)
    : StringProperty(id, guiText, url, invalidationLevel, lod)
    , VoreenFileWatchListener(watchMode)
    , volume_(0)
    , volumeOwner_(false)
    , infoProp_(0)
    , recursiveSet_(false)
    , progressBar_(0)
{}

VolumeURLProperty::VolumeURLProperty()
    : StringProperty()
    , VoreenFileWatchListener(VoreenFileWatchListener::OPTIONAL_OFF)
    , volume_(0)
    , volumeOwner_(false)
    , infoProp_(0)
    , recursiveSet_(false)
    , progressBar_(0)
{}

Property* VolumeURLProperty::create() const {
    return new VolumeURLProperty();
}

void VolumeURLProperty::fileActionCallback() {
    recursiveSet_ = true;
    try {
        loadVolume();
    }
    catch (const tgt::FileException& e) {

        // Remove file watch.
        removeWatch(VolumeURL(get()).getPath());

        // Hide progressbar.
        if (getProgressBar())
            getProgressBar()->hide();

        // Reset volume.
        setVolume(nullptr);

        // Output error message.
        LWARNING(e.what() << " The file might has been moved or deleted.");
    }
    recursiveSet_ = false;
}

void VolumeURLProperty::deinitialize() {
    std::string curURL = get();
    setVolume(0); //< also deletes the volume, if property is owner
    delete progressBar_;
    progressBar_ = 0;

    StringProperty::deinitialize();

    // restore URL (has been cleared by setVolume(0)
    set(curURL);
}

void VolumeURLProperty::set(const std::string& url) {

    if (getVolume() && getVolume()->getOrigin().getURL() != url)
        setVolume(0);

    // Remove old watch.
    bool recursive = recursiveSet_;
    if (!recursive) {
        recursiveSet_ = true;
        removeWatch(VolumeURL(get()).getPath());
    }

    // Set actual value (might trigger set by callback).
    StringProperty::set(url);

    // Add new watch.
    if (!recursive) {
        // Add watch if possible.
        std::string path = VolumeURL(get()).getPath();
        bool success = addWatch(path);
        if (!success) {
            LWARNING("Parent directory of " << tgt::FileSystem::fileName(path) << " does not exist. Resetting path.");
            setFileWatchEnabled(false);
            StringProperty::set("");
        }
        recursiveSet_ = false;
    }

    // invalidate even if the url is the same to fix problems with loading the same url again
    invalidate();
}

void VolumeURLProperty::setURL(const std::string& url) {
    set(url);
}

std::string VolumeURLProperty::getURL() const {
    return get();
}

void VolumeURLProperty::setVolume(VolumeBase* handle, bool owner) {

    if ((handle != volume_) && volumeOwner_) {
        delete volume_;
    }

    volume_ = handle;
    volumeOwner_ = owner;

    if(infoProp_)
        infoProp_->setVolume(handle);

    set(handle ? handle->getOrigin().getURL() : "");

    updateWidgets();
}

VolumeBase* VolumeURLProperty::getVolume() const {
    return volume_;
}

void VolumeURLProperty::loadVolume() {

    std::string url = get();
    if (url.empty()) {
        LWARNING("loadVolume(): empty URL");
        return;
    }

    ProgressBar* progressBar = getProgressBar();
    if (progressBar) {
        progressBar->setTitle("Loading volume");
        progressBar->setProgressMessage("Loading volume ...");
    }
    VolumeSerializerPopulator serializerPopulator(progressBar);

    VolumeBase* handle = serializerPopulator.getVolumeSerializer()->read(VolumeURL(url));

    if (progressBar)
        progressBar->hide();

    if (handle) {
        tgtAssert(handle, "No handle");
        setVolume(static_cast<Volume*>(handle));

        // property does take ownership of loaded handles
        volumeOwner_ = true;
    }
}

void VolumeURLProperty::clear() {
    if (volumeOwner_)
        delete volume_;
    volume_ = 0;
    volumeOwner_ = false;
    if(infoProp_)
        infoProp_->setVolume(0);
    removeWatch(VolumeURL(get()).getPath());
    StringProperty::set("");
}

void VolumeURLProperty::addInfoProperty(VolumeInfoProperty* pointer) {
    infoProp_ = pointer;
}

void VolumeURLProperty::serialize(Serializer& s) const {
    Property::serialize(s);
    VoreenFileWatchListener::serialize(s);
    s.serialize("urls", VolumeURL(get()));
}

void VolumeURLProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);
    VoreenFileWatchListener::deserialize(s);
    
    std::string value;
    try {
        VolumeURL tmp;
        s.deserialize("urls", tmp);
        value = tmp.getURL();
    }
    catch(SerializationNoSuchDataException) {
        s.removeLastError();
        //LINFO("trying old deserialization");
        //old deserialization
        std::string relativeURL;
        s.deserialize("url", relativeURL);
        if (relativeURL.empty())
            value = "";
        else {
            std::string basePath = tgt::FileSystem::dirName(s.getDocumentPath());
            value = VolumeURL::convertURLToAbsolutePath(relativeURL, basePath);
        }
    }

    // Set value.
    set(value);
}

ProgressBar* VolumeURLProperty::getProgressBar() {
    if (!progressBar_)
        progressBar_ = VoreenApplication::app()->createProgressDialog();
    return progressBar_;
}

} // namespace voreen
