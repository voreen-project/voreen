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

#include "voreen/core/utils/stringutils.h"

#include "bodypartscatalogbase.h"
#include "tgt/filesystem.h"

namespace voreen{

void BodyParts3DCatalogBase::connectBodyParts(const std::string& filename){
    tgt::File* file = FileSys.open(path_ + "/" + filename);
    std::string tmp, id, current;
    if(file)
        tmp = file->getLine();
    else
        throw tgt::CorruptedFileException("Cant read file " ,path_ + "/" + filename);
    std::vector<std::string> sections;
    BodyPart* currentNode;
    while(!file->eof()){
        tmp = file->getLine();
        sections = strSplit(tmp, '\t');
        if(sections.size() > 3){
            if(current != sections[0]){
                current = sections[0];
                id = sections[0].substr(sections[0].find_first_of("0123456789"),std::string::npos);
                currentNode = getBodyPart(static_cast<size_t>(findBodyPart(genericFromString<int>(id))));
            }
            if(currentNode){
                id = sections[2].substr(sections[2].find_first_of("0123456789"),std::string::npos);
                BodyPart* child = getBodyPart(static_cast<size_t>(findBodyPart(genericFromString<int>(id))));
                if(child){
                    if(child != currentNode){
                        currentNode->parts.push_back(child);
                        child->partOf = currentNode;
                    }else{
                        if(id.find("nsn") != std::string::npos){
                            child = getBodyPart(static_cast<size_t>(findBodyPart(genericFromString<int>(id)+1000000)));
                            if(child){
                                currentNode->parts.push_back(child);
                                child->partOf = currentNode;
                            }
                        }
                    }
                }
            }else{
                break;
            }
        }else{
            throw tgt::CorruptedFileException("Not Enough data in line  " + path_ + "/" + filename);
        }
    }
    file->close();
    delete file;
}

void BodyParts3DCatalogBase::sort(){
    bool swapped = true;
    while(swapped){
        swapped = false;
        for(int i = 1; i < catalog_.size() ; ++i){
            if(catalog_[i-1]->id > catalog_[i]->id){
                BodyPart* tmp = catalog_[i-1];
                catalog_[i-1] = catalog_[i];
                catalog_[i] = tmp;
                swapped = true;
            }
        }
    }
}

void BodyParts3DCatalogBase::setRenderState(BodyPart* part){
    if(part->parts.size() > 0){
        for(std::vector<BodyPart*>::iterator partIter = part->parts.begin() ; partIter != part->parts.end() ; ++partIter){
            setRenderState(*partIter);
        }
    }else{
        for(std::vector<BodyPartfile*>::iterator fileIter = part->files.begin(); fileIter != part->files.end() ; ++fileIter){
            (*fileIter)->render = part->checkState;
        }
    }
}

void BodyParts3DCatalogBase::resetRenderState(){
    for(std::vector<BodyPartfile*>::iterator fileIter = files_.begin() ; fileIter != files_.end(); ++fileIter){
        (*fileIter)->render = false;
    }
}

size_t BodyParts3DCatalogBase::findBodyPart(int id){
    size_t lowerlimit = 0;
    size_t upperlimit = catalog_.size()-1;
    if(catalog_[lowerlimit]->id == id)
        return lowerlimit;
    if(catalog_[upperlimit]->id == id)
        return upperlimit;
    while(lowerlimit < upperlimit){
        size_t middle = (lowerlimit + upperlimit) / 2;
        if(id < catalog_[middle]->id){
            upperlimit = middle;
        }else if(id > catalog_[middle]->id){
            lowerlimit = middle;
        }else{
            return middle;
        }
    }
    return -1;
}

BodyPart* BodyParts3DCatalogBase::getBodyPart(size_t index){
    if(index >= 0 && index < catalog_.size())
        return catalog_[index];
    return 0;
}

BodyPart* BodyParts3DCatalogBase::getHead(){
    return getBodyPart(head_);
}

std::vector<std::string> BodyParts3DCatalogBase::fileNameOutput(){
    resetRenderState();
    for(std::vector<BodyPart*>::iterator partIter = catalog_.begin() ; partIter != catalog_.end() ; ++partIter){
        setRenderState(*partIter);
    }
    std::vector<std::string> fileNameList;
    for(std::vector<BodyPartfile*>::iterator fileIter = files_.begin() ; fileIter != files_.end() ; ++fileIter){
        if((*fileIter)->render){
            fileNameList.push_back((*fileIter)->fileName);
        }
    }
    return fileNameList;
}

std::vector<tgt::vec4> BodyParts3DCatalogBase::colorOutput(){
    std::vector<tgt::vec4> colorList;

    for(std::vector<BodyPartfile*>::iterator fileIter = files_.begin() ; fileIter != files_.end() ; ++fileIter){
        if((*fileIter)->render){
            colorList.push_back((*fileIter)->color);
        }
    }
    return colorList;
}

std::vector<std::string> BodyParts3DCatalogBase::descriptionOutput(){
    std::vector<std::string> descList;
    for(std::vector<BodyPartfile*>::iterator fileIter = files_.begin() ; fileIter != files_.end() ; ++fileIter){
        if((*fileIter)->render){
            descList.push_back((*fileIter)->description);
        }
    }
    return descList;
}

}//namespace voreen
