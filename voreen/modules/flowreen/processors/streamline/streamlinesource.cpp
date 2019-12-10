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

#include "streamlinesource.h"

#include "../../datastructures/streamlinelist.h"

#include "modules/core/io/vvdformat.h"

#include "voreen/core/voreenapplication.h"
#include "tgt/filesystem.h"

#include <fstream>

namespace voreen {

const std::string StreamlineSource::loggerCat_("voreen.flowreen.StreamlineSource");

StreamlineSource::StreamlineSource()
    : Processor()
    // ports
    , streamlineOutport_(Port::OUTPORT, "streamlineOutport", "Streamline Output")
    , volumeOutport_(Port::OUTPORT, "volumeOutport", "Reference Volume Output")
    //properties
    , loadOptionProp_("loadOptionProp","Load Option:")
        //single
    , singleFileProp_("singleFileProp", "Load Streamlines:", "Select File...",
                      VoreenApplication::app()->getUserDataPath(), "Voreen streamline data (*.vsd)",
                      FileDialogProperty::OPEN_FILE,Processor::INVALID_PATH)
        //directory
    , directoryFileProp_("directoryFileProp", "Load Directory:", "Select Directory...",
                      VoreenApplication::app()->getUserDataPath(), "Voreen streamline data (*.vsd)",
                      FileDialogProperty::DIRECTORY,Processor::INVALID_PATH)
    , tableProp_("tableProp","",4)
    , currentListProp_("currentListProp","Current Output:",-1,-1,-1,Processor::INVALID_RESULT,NumericProperty<int>::DYNAMIC)
    // members
    , forceStreamlineLoadFlag_(true), updateOutputFlag_(false), firstDirectoryLoad_(true)
{
    //ports
    addPort(streamlineOutport_);
    addPort(volumeOutport_);
    //properties
    addProperty(loadOptionProp_);
        loadOptionProp_.addOption("single","Load Single File",LOAD_SINGLE);
        loadOptionProp_.addOption("directory","Load Entire Directory",LOAD_DIRECTORY);
        loadOptionProp_.onChange(MemberFunctionCallback<StreamlineSource>(this,&StreamlineSource::loadOptionOnChange));
        //single
    addProperty(singleFileProp_);
        singleFileProp_.setGroupID("single");
        //directory
    addProperty(directoryFileProp_);
        directoryFileProp_.setGroupID("directory");
    addProperty(tableProp_);
        tableProp_.setColumnLabel(0,"#");
        tableProp_.setColumnLabel(1,"Streamlines");
        tableProp_.setColumnLabel(2,"Ref.Volume");
        tableProp_.setColumnLabel(3,"Filename");
        tableProp_.onChange(MemberFunctionCallback<StreamlineSource>(this,&StreamlineSource::tableSelectionOnChange));
        tableProp_.setGroupID("directory");
    addProperty(currentListProp_);
        currentListProp_.onChange(MemberFunctionCallback<StreamlineSource>(this,&StreamlineSource::currentListOnChange));
        currentListProp_.setGroupID("directory");

    setPropertyGroupGuiName("single","Load Single File");
    setPropertyGroupGuiName("directory","Load Entire Directory");

    //set visibility
    loadOptionOnChange();
}

bool StreamlineSource::isReady() const {
    if(!isInitialized())
        return false;

    if(!streamlineOutport_.isReady())
        return false;

    return true;
}

void StreamlineSource::invalidate(int inv) {
    Processor::invalidate(inv);

    if(inv == Processor::INVALID_PATH && isInitialized()) {
        forceStreamlineLoadFlag_ = true;
    }
}

void StreamlineSource::process() {
    //load files if necessary
    if (forceStreamlineLoadFlag_){
        switch(loadOptionProp_.getValue()) {
        case LOAD_SINGLE:
            //loads a single streamline and updates the output
            loadSingleStreamlines();
            break;
        case LOAD_DIRECTORY:
            //loads all VSD files from a directory. updateOutputFlag_ is set true.
            loadDirectoryStreamlines();
            break;
        default:
            LERROR("Unknown load option");
            tgtAssert(false,"Should not get here");
        }
        //set flag back to false
        forceStreamlineLoadFlag_ = false;
    }

    //update output if we are in directory mode and files have been loaded
    if(updateOutputFlag_) {
        // The cast is ok here, since -1 translates to std::numeric_limits<size_t>:max(), which is totally fine here.
        if((loadOptionProp_.getValue() == LOAD_DIRECTORY) && (currentlyLoadedListsAndVolumes_.size() > static_cast<size_t>(currentListProp_.get()))) {
            streamlineOutport_.setData(currentlyLoadedListsAndVolumes_[currentListProp_.get()].first,false);
            volumeOutport_.setData(currentlyLoadedListsAndVolumes_[currentListProp_.get()].second,false);
        }
        updateOutputFlag_ = false;
    }
}

void StreamlineSource::loadSingleStreamlines() {
    if (!isInitialized())
        return;
    if (singleFileProp_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }
    std::string extension = FileSys.fileExtension(singleFileProp_.get(),true);
    if(extension.compare("vsd")) {
        LWARNING("Selected file extension is not *.vsd");
        return;
    }

    // load streamlines
    std::pair<StreamlineList*,Volume*> output = loadVSDFile(singleFileProp_.get());

    //reset property on error
    if(!output.first)
        singleFileProp_.set("");

    //set new streamlines
    streamlineOutport_.setData(output.first); volumeOutport_.setData(output.second);
}

void StreamlineSource::loadDirectoryStreamlines() {
    if (!isInitialized())
        return;
    //clear list
    for(size_t i = 1; i < currentlyLoadedListsAndVolumes_.size(); i++) {
        delete currentlyLoadedListsAndVolumes_[i].first;
        delete currentlyLoadedListsAndVolumes_[i].second;
    }
    currentlyLoadedListsAndVolumes_.clear();

    int lastSelectedRow = 0;
    if(firstDirectoryLoad_) {
        lastSelectedRow = currentListProp_.get();
        firstDirectoryLoad_ = false;
    }

    currentListProp_.setMinValue(-1);
    currentListProp_.setMaxValue(-1);
    tableProp_.reset();

    if (directoryFileProp_.get().empty()) {
        LWARNING("no directory specified");
        return;
    }
    if (!FileSys.dirExists(directoryFileProp_.get())) {
        LWARNING("directory no longer exists");
        return;
    }

    //load vsd files
    std::vector<std::string> files = FileSys.listFiles(directoryFileProp_.get(),true);
    int length = static_cast<int>(files.size() + 10 - (files.size() % 10));
    for(size_t i = 0; i < files.size(); i++) {
        if(!FileSys.fileExtension(files[i],true).compare("vsd")) {
            // load listand volume
            std::pair<StreamlineList*,Volume*> output = loadVSDFile(directoryFileProp_.get() + "/" + files[i]);
            if(output.first) {
                currentlyLoadedListsAndVolumes_.push_back(output);
                std::vector<std::string> row(4);
                row[0] = itos(currentlyLoadedListsAndVolumes_.size() - 1,length);
                row[1] = itos(output.first->getStreamlines().size());
                row[2] = (output.second ? "yes" : "no");
                row[3] = files[i];
                tableProp_.addRow(row);
            }
        }
    }

    if(currentlyLoadedListsAndVolumes_.empty()) {
        LWARNING("No streamlines could be loaded");
        streamlineOutport_.setData(0); volumeOutport_.setData(0);
        return;
    }

    //update max range of int slider
    currentListProp_.setMaxValue(currentlyLoadedListsAndVolumes_.size()-1);
    tableProp_.setSelectedRowIndex(lastSelectedRow);
    currentListProp_.setMinValue(0);
    //output must be updated
    updateOutputFlag_ = true;
}

//-------------------------------------------------------------------------------
//      Helpers
//-------------------------------------------------------------------------------
std::pair<StreamlineList*,Volume*> StreamlineSource::loadVSDFile(const std::string& absPath) {
        StreamlineList* list = new StreamlineList();
        Volume* volume = 0;
        std::vector<VvdObject> vec; ///< used for deserialization
    try {
        std::ifstream inFile;
        inFile.open(absPath.c_str());
        LINFO("Reading streamlines from file " << absPath);

        XmlDeserializer d;
        d.read(inFile);
        Deserializer deserializer(d);
        deserializer.deserialize("VSD-File", *list);

        //optional deserialize not imlpemented for std::vector :-(
        try {
            deserializer.deserialize("Magnitude-Volumes",vec,"Volume");
            volume = vec[0].createVolume();
        } catch (VoreenException& e) {
            //do nothing if deserialize fails
        }

        inFile.close();
    }
    catch(VoreenException& e) {
        LERROR(e.what());
        delete list;
        list = 0;
    }
    return std::make_pair(list,volume);
}

//-------------------------------------------------------------------------------
//      Callbacks
//-------------------------------------------------------------------------------
void StreamlineSource::loadOptionOnChange() {
    //clear output
    streamlineOutport_.setData(0); volumeOutport_.setData(0);
    //are we in single mode?
    bool single = (loadOptionProp_.getValue() == LOAD_SINGLE);
    if(single) {
        //clear old directory list
        for(size_t i = 1; i < currentlyLoadedListsAndVolumes_.size(); i++) {
            delete currentlyLoadedListsAndVolumes_[i].first;
            delete currentlyLoadedListsAndVolumes_[i].second;
        }
        currentlyLoadedListsAndVolumes_.clear();
        currentListProp_.blockCallbacks(true);
        currentListProp_.setMinValue(-1);
        currentListProp_.setMaxValue(-1);
        currentListProp_.blockCallbacks(false);
        tableProp_.reset();
        //force new load
        if(FileSys.exists(singleFileProp_.get()))
            forceStreamlineLoadFlag_ = true;
    } else {
        if(FileSys.dirExists(directoryFileProp_.get()))
            forceStreamlineLoadFlag_ = true;
    }

    //update visibility
    setPropertyGroupVisible("single",single);
    setPropertyGroupVisible("directory",!single);
}

void StreamlineSource::currentListOnChange() {
    updateOutputFlag_ = true;
    tableProp_.setSelectedRowIndex(currentListProp_.get());
}

void StreamlineSource::tableSelectionOnChange() {
    currentListProp_.set(tableProp_.getSelectedRowIndex());
}

} // namespace
