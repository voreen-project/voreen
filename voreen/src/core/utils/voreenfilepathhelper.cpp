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

#include "voreen/core/utils/voreenfilepathhelper.h"

#include "tgt/filesystem.h"

#include "voreen/core/voreenapplication.h"

namespace voreen {

VoreenFilePathHelper::VoreenFilePathHelper(const std::string& filePath)
    : path_(filePath)
{}

VoreenFilePathHelper::~VoreenFilePathHelper() {}

const std::string& VoreenFilePathHelper::getPath() const {
    return path_;
}

void VoreenFilePathHelper::serialize(Serializer& s) const {
    // An empty serialized value does not necessarily mean that it wasn't set, but could also mean that it was the
    // same path as the document path passed during serialization, which makes the relative path empty.  We need an
    // extra bool to remember if this was the case.
    if (path_.empty()) {
        s.serialize("noPathSet", true);
    } else {
        s.serialize("noPathSet", false);

        std::vector<std::string> pathVector;

        //add relative path
        pathVector.push_back(tgt::FileSystem::cleanupPath(
                             tgt::FileSystem::relativePath(
                                              path_,tgt::FileSystem::dirName(s.getDocumentPath())),false));

        //add voreen path
        if(VoreenApplication* vapp = VoreenApplication::app())
            pathVector.push_back(tgt::FileSystem::cleanupPath(
                                 tgt::FileSystem::relativePath(
                                                  path_,tgt::FileSystem::dirName(vapp->getBasePath())),false));

        //add absolute path
        // call to tgt::FileSystem::absolutePath fixes problems with links under linux
        std::string absPath = tgt::FileSystem::absolutePath(tgt::FileSystem::cleanupPath(path_));
        pathVector.push_back(tgt::FileSystem::cleanupPath(absPath,false));

        s.serialize("paths",pathVector);
    }
}

void VoreenFilePathHelper::deserialize(Deserializer& s) {
    path_ = "";
    // An empty serialized value does not necessarily mean that it wasn't set, but could also mean that it was the
    // same path as the document path passed during serialization, which makes the relative path empty.  We need an
    // extra bool to remember if this was the case.
    try {
        bool noPathSet;
        s.deserialize("noPathSet", noPathSet);
        if(noPathSet) {
            return;
        }
    }
    catch (SerializationNoSuchDataException&) {
        s.removeLastError();
        return;
    }

    //deserialize path values
    std::vector<std::string> pathVector;
    try {
        s.deserialize("paths", pathVector);
    }
    catch (SerializationNoSuchDataException&) {
        s.removeLastError();
        return;
    }

    //find valid path
    std::string tmpString;
    for(size_t i = 0; i < pathVector.size(); i++) {
        // check absolute path
        if(tgt::FileSystem::isAbsolutePath(pathVector[i])) {
            tmpString = tgt::FileSystem::cleanupPath(pathVector[i],true);
            if(tgt::FileSystem::fileExists(tmpString) || tgt::FileSystem::dirExists(tmpString)) {
                path_ = tmpString;
                return;
            }
        } else {
            // check relative path
            if (!s.getDocumentPath().empty()) {
                tmpString = tgt::FileSystem::cleanupPath(
                            tgt::FileSystem::absolutePath(tgt::FileSystem::dirName(s.getDocumentPath()) + "/" + pathVector[i]),true);
                if(tgt::FileSystem::fileExists(tmpString) || tgt::FileSystem::dirExists(tmpString)) {
                    path_ = tmpString;
                    return;
                }
            }
            //check relative voreen path
            if (VoreenApplication* vapp = VoreenApplication::app()) {
                tmpString = tgt::FileSystem::cleanupPath(
                            tgt::FileSystem::absolutePath(tgt::FileSystem::dirName(vapp->getBasePath()) + "/" + pathVector[i]),true);
                if(tgt::FileSystem::fileExists(tmpString) || tgt::FileSystem::dirExists(tmpString)) {
                    path_ = tmpString;
                    return;
                }
            }
        }
    }

    //no valid path found. Use first stored path as default
    if(!pathVector.empty()) {
        // check absolute path
        if(tgt::FileSystem::isAbsolutePath(pathVector[0])) {
            path_ = tgt::FileSystem::cleanupPath(pathVector[0],true);
            return;
        }
        // check relative path
        if (!s.getDocumentPath().empty()) {
            path_ = tgt::FileSystem::cleanupPath(
                    tgt::FileSystem::absolutePath(tgt::FileSystem::dirName(s.getDocumentPath()) + "/" + pathVector[0]),true);
            return;
        }
        //check relative voreen path
        if (VoreenApplication* vapp = VoreenApplication::app()) {
            path_ = tgt::FileSystem::cleanupPath(
                    tgt::FileSystem::absolutePath(tgt::FileSystem::dirName(vapp->getBasePath()) + "/" + pathVector[1]),true);
            return;
        }
    }
}

} // namespace
