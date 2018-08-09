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

#include "voreen/core/datastructures/roi/roicollection.h"
#include "voreen/core/datastructures/roi/roiaggregation.h"

#include "voreen/core/io/serialization/serialization.h"

#include "tgt/filesystem.h"

#ifdef __unix__
#include <unistd.h>
#endif

namespace voreen {

const std::string ROICollection::loggerCat_ = "voreen.ROICollection";

ROICollection::~ROICollection() {
    clear();
}

void ROICollection::clear() {
    while(!rois_.empty()) {
        delete rois_.back();
        rois_.pop_back();
    }
}

void ROICollection::serialize(Serializer& s) const {
    s.serialize("ROIs", rois_, "ROI");
}

void ROICollection::deserialize(Deserializer& s) {
    clear();
    s.deserialize("ROIs", rois_, "ROI");
}

bool ROICollection::moveControlPoint(const ControlPoint* cp, tgt::vec3 to) {
    for(size_t i=0; i<rois_.size(); i++) {
        if(rois_[i]->moveControlPoint(cp, to))
            return true;
    }
    return false;
}

void ROICollection::removeROI(ROIBase* roi) {
    for(size_t i=0; i<rois_.size(); i++) {
        if(rois_[i] == roi) {
            rois_.erase(rois_.begin()+i);
            delete roi;
            return;
        }
        else if(removeROIrecursive(rois_[i], roi))
            return;
    }
}

void ROICollection::deregisterROI(ROIBase* roi) {
    for(size_t i=0; i<rois_.size(); i++) {
        if(rois_[i] == roi) {
            rois_.erase(rois_.begin()+i);
            return;
        }
        else if(removeROIrecursive(rois_[i], roi)) //TODO
            return;
    }
}

bool ROICollection::removeROIrecursive(ROIBase* p, ROIBase* roi) {
    std::vector<ROIBase*> c = p->getChildren();
    for(size_t i=0; i<c.size(); i++) {
        if(c[i] == roi) {
            if(dynamic_cast<ROIAggregation*>(p)) {
                static_cast<ROIAggregation*>(p)->removeChild(roi);
            }
            else {
                LERROR("Tried to remove child from non-aggregation ROI");
            }
            return true;
        }
        else if(removeROIrecursive(c[i], roi))
            return true;
    }
    return false;
}

void ROICollection::save(const std::string& filename) const {
    // serialize
    XmlSerializer s(FileSys.dirName(filename));
    s.setUseAttributes(true);
    s.serialize("ROICollection", *this);

    // write serialization data to temporary string stream
    std::ostringstream textStream;

    try {
        s.write(textStream);
        if (textStream.fail())
            throw SerializationException("Failed to write serialization data to string stream.");
    }
    catch (std::exception& e) {
        throw SerializationException("Failed to write serialization data to string stream: " + std::string(e.what()));
    }
    catch (...) {
        throw SerializationException("Failed to write serialization data to string stream (unknown exception).");
    }

    // Now we have a valid StringStream containing the serialization data.
    // => Open output file and write it to the file.
    // For added data security we write to a temporary file and afterwards move it into place
    // (which should be an atomic operation).
    const std::string tmpfilename = filename + ".tmp";
    std::fstream fileStream(tmpfilename.c_str(), std::ios_base::out);
    if (fileStream.fail())
        throw SerializationException("Failed to open file '" + tmpfilename + "' for writing.");

    try {
        fileStream << textStream.str();
    }
    catch (std::exception& e) {
        throw SerializationException("Failed to write serialization data stream to file '"
                                     + tmpfilename + "': " + std::string(e.what()));
    }
    catch (...) {
        throw SerializationException("Failed to write serialization data stream to file '"
                                     + tmpfilename + "' (unknown exception).");
    }
    fileStream.close();

    // Finally move the temporary file into place. It is important that this happens in-place,
    // without deleting the old file first.
    bool success;
#ifdef WIN32
    // rename() does not replace existing files on Windows, so we have to use this
    success = (MoveFileEx(tmpfilename.c_str(), filename.c_str(),
                          MOVEFILE_REPLACE_EXISTING | MOVEFILE_COPY_ALLOWED) != 0);

#else
    // atomic replace
    success = (rename(tmpfilename.c_str(), filename.c_str()) == 0);
#endif
    if (!success) {
#ifdef WIN32
        _unlink(tmpfilename.c_str()); // ignore failure here
#else
        unlink(tmpfilename.c_str()); // ignore failure here
#endif
        throw SerializationException("Failed to rename temporary file '" + tmpfilename + "' to '"
                                     + filename + "'");
    }
}

void ROICollection::load(const std::string& filename) {
    // open file for reading
    std::fstream fileStream(filename.c_str(), std::ios_base::in);
    if (fileStream.fail()) {
        throw SerializationException("Failed to open ROICollection file '" + tgt::FileSystem::absolutePath(filename) + "' for reading.");
    }

    // read data stream into deserializer
    XmlDeserializer d(FileSys.dirName(filename));
    d.setUseAttributes(true);
    try {
        d.read(fileStream);
    }
    catch (SerializationException& e) {
        throw SerializationException("Failed to read serialization data stream from ROICollection file '"
                                     + filename + "': " + e.what());
    }
    catch (...) {
        throw SerializationException("Failed to read serialization data stream from ROICollection file '"
                                     + filename + "' (unknown exception).");
    }

    // deserialize ROICollection from data stream
    try {
        d.deserialize("ROICollection", *this);
    }
    catch (std::exception& e) {
        throw SerializationException("Deserialization from ROICollection file '" + filename + "' failed: " + e.what());
    }
    catch (...) {
        throw SerializationException("Deserialization from ROICollection file '" + filename + "' failed (unknown exception).");
    }
}

} // namespace
