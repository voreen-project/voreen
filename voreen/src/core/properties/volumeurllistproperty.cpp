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

#include "voreen/core/properties/volumeurllistproperty.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/progressbar.h"
#include "voreen/core/voreenapplication.h"

#include "tgt/logmanager.h"
#include "tgt/filesystem.h"

namespace voreen {

const std::string VolumeURLListProperty::loggerCat_("voreen.VolumeURLListProperty");

VolumeURLListProperty::VolumeURLListProperty(const std::string& id, const std::string& guiText,
                    const std::vector<std::string>& value, int invalidationLevel, Property::LevelOfDetail lod,
                    VoreenFileWatchListener::WatchMode watchMode)
    : TemplateProperty<std::vector<std::string> >(id, guiText, value, invalidationLevel, lod)
    , VoreenFileWatchListener(watchMode)
    , previewsVisible_(true)
    , progressBar_(0)
{}

VolumeURLListProperty::VolumeURLListProperty()
    : TemplateProperty<std::vector<std::string> >("", "", std::vector<std::string>(), Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , VoreenFileWatchListener()
    , progressBar_(0)
{}

Property* VolumeURLListProperty::create() const {
    return new VolumeURLListProperty();
}

void VolumeURLListProperty::fileActionCallback() {
    loadVolumes(false, false, true); // force reload!
}

void VolumeURLListProperty::deinitialize() {
    clear();
    delete progressBar_;
    progressBar_ = 0;

    TemplateProperty<std::vector<std::string> >::deinitialize();
}

void VolumeURLListProperty::set(const std::vector<std::string>& URLs) {
    clear();
    TemplateProperty<std::vector<std::string> >::set(URLs);
    for (const std::string& url : URLs)
        addWatch(VolumeURL(url).getPath());
}

void VolumeURLListProperty::setURLs(const std::vector<std::string>& URLs, bool selected /*=true*/) {
    set(URLs);
    setAllSelected(selected);

    invalidate();
}

void VolumeURLListProperty::addURL(const std::string& url, bool selected /*=true*/) {
    if (containsURL(url)) {
        LWARNING("Volume url already present, not adding to VolumeList");
        return;
    }

    value_.push_back(url);
    selectionMap_[url] = selected;
    addWatch(VolumeURL(url).getPath());

    invalidate();
}

const std::vector<std::string>& VolumeURLListProperty::getURLs() const {
    return value_;
}

bool VolumeURLListProperty::containsURL(const std::string& url) const {
    return std::find(value_.begin(), value_.end(), url) != value_.end();
}

void VolumeURLListProperty::removeURL(const std::string& url) {
    if (!containsURL(url)) {
        LWARNING("removeURL(): passed URL not contained by this property: " << url);
        return;
    }

    // remove file watch
    removeWatch(VolumeURL(url).getPath());

    // free corresponding volume, if owned by the property
    VolumeBase* handle = getVolume(url);
    if (handle)
        removeVolume(handle);

    // remove url from url list
    std::vector<std::string>::iterator it = std::find(value_.begin(), value_.end(), url);
    tgtAssert(it != value_.end(), "url not found");
    value_.erase(it);

    // remove url from maps
    handleMap_.erase(url);
    selectionMap_.erase(url);
    ownerMap_.erase(url);

    invalidate();
}

VolumeList* VolumeURLListProperty::getVolumes(bool selectedOnly /*= true*/) const {
    VolumeList* collection = new VolumeList();
    std::vector<std::string> urls = get();
    for (size_t i=0; i<urls.size(); i++) {
        std::string url = urls.at(i);
        if (handleMap_.find(url) != handleMap_.end()) {
            if (!selectedOnly || (selectionMap_.find(url) != selectionMap_.end() && selectionMap_.find(url)->second == true) ) {
                VolumeBase* handle = handleMap_.find(url)->second;
                tgtAssert(handle, "handleMap_ contains null pointer");
                collection->add(handle);
            }
        }
    }

    return collection;
}

VolumeBase* VolumeURLListProperty::getVolume(const std::string& url) const {
    if (!containsURL(url)) {
        LWARNING("getVolume(): passed URL not contained by this property: " << url);
        return 0;
    }

    VolumeBase* result = 0;
    if (handleMap_.find(url) != handleMap_.end()) {
        result = handleMap_.find(url)->second;
        tgtAssert(result, "handleMap_ contains null pointer");
    }

    return result;
}

void VolumeURLListProperty::loadVolume(const std::string& url, bool invalidateProperty /*= true */) {

    if (!containsURL(url)) {
        LWARNING("loadVolume(): passed URL not contained by this property: " << url);
        return;
    }

    // delete volume, if already loaded and owned by the property
    if (getVolume(url) && isOwner(url))
        delete getVolume(url);
    handleMap_.erase(url);
    ownerMap_.erase(url);

    ProgressBar* progressBar = getProgressBar();
    if (progressBar) {
        progressBar->setTitle("Loading volume");
        progressBar->setProgressMessage("Loading volume ...");
    }

    VolumeBase* handle = VolumeSerializerPopulator(progressBar).getVolumeSerializer()->read(VolumeURL(url));

    if (progressBar)
        progressBar->hide();

    if (handle) {
        // url may have been altered by loading routine
        if (url != handle->getOrigin().getURL()) {
            bool selected = isSelected(url);
            selectionMap_.erase(url);
            selectionMap_[handle->getOrigin().getURL()] = selected;

            for (size_t i=0; i<value_.size(); i++) {
                if (value_[i] == url) {
                    value_[i] = handle->getOrigin().getURL();
                    break;

                }
            }
        }

        handleMap_[handle->getOrigin().getURL()] = handle;
        ownerMap_[handle->getOrigin().getURL()] = true;
    }

    if (invalidateProperty)
        invalidate();
}

void VolumeURLListProperty::loadVolumes(bool selectedOnly /*= false*/, bool removeOnFailure /*= false*/, bool forceReload /*=false*/) {

    std::vector<std::string> failedURLs;

    for (size_t i=0; i<value_.size(); i++) {
        std::string url = value_[i];
        if (selectedOnly && !isSelected(url))
            continue;
        if (!forceReload && getVolume(url))
            continue;

        try {
            // load volume, but do not (yet) invalidate, otherwise the network is evaluated for every loaded volume individually
            loadVolume(url, false);
        }
        catch (tgt::FileException& e) {
            LERROR(e.what());
            failedURLs.push_back(url);
        }
        catch (std::bad_alloc&) {
            LERROR("bad allocation while loading volume: " << url);
            failedURLs.push_back(url);
        }
        catch (tgt::Exception& e) {
            LERROR("unknown exception while loading volume '" << url << "':" << e.what());
            failedURLs.push_back(url);
        }
    }

    if (removeOnFailure) {
        for (size_t i=0; i<failedURLs.size(); i++)
            removeURL(failedURLs.at(i));
    }

    invalidate();
}

void VolumeURLListProperty::addVolume(VolumeBase* handle, bool owner /*= false*/, bool selected /*= false*/) {
    tgtAssert(handle, "null pointer passed");

    std::string url = handle->getOrigin().getURL();
    if (!containsURL(url))
        addURL(url, selected);
    else {
        LWARNING("Volume already in list... removing old instance");
        //removeVolume(handle);
        removeURL(url);
        addURL(url, selected);
    }

    handleMap_[url] = handle;
    ownerMap_[url] = owner;

    invalidate();
}

void VolumeURLListProperty::removeVolume(VolumeBase* handle) {
    tgtAssert(handle, "null pointer passed");
    std::string url = handle->getOrigin().getURL();
    if (!containsURL(url)) {
        LWARNING("removeVolume(): passed handle's URL not contained by this property: " << url);
        return;
    }

    if ((getVolume(url) == handle) && isOwner(url)) {
        delete handle;
    }

    handleMap_.erase(url);
    ownerMap_[url] = false;

    invalidate();
}

void VolumeURLListProperty::setSelected(const std::string& url, bool selected) {
    if (!containsURL(url)) {
        LWARNING("setSelected(): passed URL not contained by this property: " << url);
        return;
    }

    if (isSelected(url) != selected) {
        selectionMap_[url] = selected;
        invalidate();
    }
}

void VolumeURLListProperty::setAllSelected(bool selected) {
    for (size_t i=0; i<value_.size(); i++)
        selectionMap_[value_.at(i)] = selected;

    invalidate();
}

bool VolumeURLListProperty::isSelected(const std::string& url) const {
    if (!containsURL(url)) {
        LWARNING("isURLSelected(): passed URL not contained by this property: " << url);
        return false;
    }

    return (selectionMap_.find(url) != selectionMap_.end()) && selectionMap_.at(url);
}

void VolumeURLListProperty::serialize(Serializer& s) const {
    Property::serialize(s);
    VoreenFileWatchListener::serialize(s);

    int numURLs = static_cast<int>(value_.size());
    std::map<std::string, bool> selectMap = selectionMap_;
    s.serialize("numURLs", numURLs);
    for(int i = 0; i < numURLs; i++) {
        s.serialize("VolumeURL_" + std::to_string(i), VolumeURL(value_[i]));
        bool selected = selectionMap_.find(value_[i]) != selectionMap_.end() && selectionMap_.at(value_[i]);
        s.serialize("VolumeSelected_" + std::to_string(i), selected);
    }
    
    s.serialize("previewsVisible", previewsVisible_);
}

void VolumeURLListProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);
    VoreenFileWatchListener::deserialize(s);

    s.optionalDeserialize("previewsVisible", previewsVisible_, true);

    try {
        int numURLs = 0;
        s.deserialize("numURLs", numURLs);
        std::vector<std::string> urlList(numURLs);
        std::map<std::string, bool> selectMap;
        VolumeURL tmp;
        for(int i = 0; i < numURLs; i++) {
            s.deserialize("VolumeURL_" + std::to_string(i), tmp);
            urlList[i] = tmp.getURL();
            s.deserialize("VolumeSelected_" + std::to_string(i), selectMap[tmp.getURL()]);

            // Add watch for each URL.
            addWatch(tmp.getPath());
        }
        value_ = urlList;
        selectionMap_ = selectMap;
    }
    catch(SerializationNoSuchDataException&) {
        s.removeLastError();
    
        //LINFO("trying old deserialization");
        //old deserialization
        std::vector<std::string> urlList = value_;
        std::map<std::string, bool> selectMap = selectionMap_;
        s.deserialize("VolumeURLs", urlList, "url");

        try {
            s.deserialize("Selection", selectMap, "entry", "url");
        }
        catch (SerializationException& e) {
            s.removeLastError();
            LWARNING("Failed to deserialize selection map: " << e.what());
        }

        // convert URLs to absolute path
        std::string basePath = tgt::FileSystem::dirName(s.getDocumentPath());
        for (size_t i=0; i<urlList.size(); i++) {
            std::string url = urlList[i];
            std::string urlConv = VolumeURL::convertURLToAbsolutePath(url, basePath);
            urlList[i] = urlConv;
            if (selectMap.find(url) != selectMap.end()) {
                bool selected = selectMap[url];
                selectMap.erase(url);
                selectMap.insert(std::pair<std::string, bool>(urlConv, selected));
            }

            addWatch(VolumeURL(urlConv).getPath());
        }

        value_ = urlList;
        selectionMap_ = selectMap;
    }
    
    invalidate();
}

void VolumeURLListProperty::clear() {
    std::vector<std::string> urls = getURLs();
    for (size_t i=0; i<urls.size(); i++)
        removeURL(urls.at(i));
}

bool VolumeURLListProperty::isOwner(const std::string& url) const {
    return (ownerMap_.find(url) != ownerMap_.end()) && (ownerMap_.find(url)->second == true);
}

ProgressBar* VolumeURLListProperty::getProgressBar() {
    if (!progressBar_)
        progressBar_ = VoreenApplication::app()->createProgressDialog();
    return progressBar_;
}

bool VolumeURLListProperty::getPreviewsVisible() {
    return previewsVisible_;
}

void VolumeURLListProperty::setPreviewsVisible(bool previewsVisible) {
    previewsVisible_ = previewsVisible;
}

} // namespace voreen
