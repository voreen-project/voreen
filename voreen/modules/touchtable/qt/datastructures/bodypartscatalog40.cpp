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

#include "bodypartscatalog40.h"
#include "tgt/filesystem.h"
#include "voreen/core/utils/stringutils.h"

namespace voreen{

BodyParts3DCatalog40::BodyParts3DCatalog40(const std::string& dir, bool colorVersion){
    path_ = dir;
    if(colorVersion){
        readBodyPartsColored("");
        path_ = path_.substr(0,path_.find_last_of('/'));
    }else{
        readBodyPartsColored("isa_parts_list_e.txt");
    }
    sort();
    connectBodyParts("isa_inclusion_relation_list.txt");
    addPrimitives("isa_element_parts.txt");
    setFileAttributes();
}

BodyParts3DCatalog40::BodyParts3DCatalog40(){

}

BodyParts3DCatalog40::~BodyParts3DCatalog40(){
    for(size_t i = 0; i < catalog_.size() ; ++i){
        delete catalog_[i];
        catalog_[i] = 0;
    }
}

void BodyParts3DCatalog40::readBodyParts(const std::string& filename){
    tgt::File* file = FileSys.open(path_ + "/" + filename);
    std::string tmp, id;
    if(file)
        tmp = file->getLine();
    else
        throw tgt::CorruptedFileException("Cant read file",path_ + "/" + filename);
    std::vector<std::string> sections;
    while(!file->eof()){
        tmp = file->getLine();
        sections = strSplit(tmp, '\t');
        if(sections.size() > 1){
            id = sections[0].substr(sections[0].find_first_of("0123456789"),std::string::npos);
            catalog_.push_back(new BodyPart(genericFromString<int>(id),sections[0],sections[2],false,tgt::vec4(0.5f,0.5f,0.5f,1.0f)));
        }else{
            throw tgt::CorruptedFileException("Not enough data in line");
        }
    }
    std::string::iterator it = catalog_[catalog_.size()-1]->description.end();
    catalog_[catalog_.size()-1]->description.erase(it-1,it);
    file->close();
    delete file;
}

void BodyParts3DCatalog40::readBodyPartsColored(const std::string& filename){
    tgt::File* file;
    if(filename == ""){
        file = FileSys.open(path_);
    }else{
        file = FileSys.open(path_ + "/" + filename);
    }
    std::string tmp, id;
    if(!file)
        throw tgt::CorruptedFileException("Cant read file",path_ + "/" + filename );

    std::vector<std::string> sections;
    while(!file->eof()){
        tmp = file->getLine();
        sections = strSplit(tmp, '\t');
        if(sections.size() > 1){
            id = sections[0].substr(sections[0].find_first_of("0123456789"),std::string::npos);
            tgt::vec4 color;
            std::vector<std::string> colors = strSplit(sections[2],' ');
            if(colors.size() > 3){
                color.x = genericFromString<float>(colors[0]);
                color.y = genericFromString<float>(colors[1]);
                color.z = genericFromString<float>(colors[2]);
                color.w = genericFromString<float>(colors[3]);
            }else{
                throw tgt::CorruptedFileException("Insufficient Color Data");
            }
            catalog_.push_back(new BodyPart(genericFromString<int>(id),sections[0],sections[1],false,color));
        }else{
            throw tgt::CorruptedFileException("Not enough data in line");
        }
    }
    std::string::iterator it = catalog_[catalog_.size()-1]->description.end();
    catalog_[catalog_.size()-1]->description.erase(it-1,it);
    file->close();
    delete file;

}


void BodyParts3DCatalog40::addPrimitives(const std::string& filename){
    tgt::File* file = FileSys.open(path_ + "/" + filename);
    std::string tmp, id, current;
    if(file)
        tmp = file->getLine();
    else
        throw tgt::CorruptedFileException("Cant read file ", path_ + "/" + filename);

    std::vector<std::string> sections;
    BodyPart* currentNode;
    while(!file->eof()){
        tmp = file->getLine();
        sections = strSplit(tmp, '\t');
        if(current != sections[0]){
            current = sections[0];
            id = sections[0].substr(sections[0].find_first_of("0123456789"),std::string::npos);
            currentNode = getBodyPart(findBodyPart(genericFromString<int>(id)));
        }
        if(BodyPartfile* bpfile = getFilePointer(sections[2])){
            currentNode->files.push_back(bpfile);
        }else{
            BodyPartfile* bpFile = new BodyPartfile(sections[2]);
            currentNode->files.push_back(bpFile);
            files_.push_back(bpFile);
        }
    }
    file->close();
    delete file;
}

void BodyParts3DCatalog40::exportCatalog(const std::string& fileName){ //baseclass
    std::ofstream stream;
    stream.open(fileName.c_str(),std::ios_base::out);
    for(int i = 0 ; i < catalog_.size() ; i++){
        std::string line;
        std::stringstream ss(std::stringstream::out);
        ss << catalog_[i]->color.x << " " << catalog_[i]->color.y << " " << catalog_[i]->color.z << " " << catalog_[i]->color.w;
        line = catalog_[i]->filename + "\t" + catalog_[i]->description + "\t" + ss.str() +"\n" ;
        stream.write(line.c_str(),line.size());
    }
    stream.close();
}

void BodyParts3DCatalog40::setFileAttributes(){
    for(size_t i = 0; i < catalog_.size() ; ++i){
        if(catalog_[i]->parts.size() == 0){
            for(size_t j = 0; j < catalog_[i]->files.size(); j++){
                catalog_[i]->files[j]->color = catalog_[i]->color;
                catalog_[i]->files[j]->description = catalog_[i]->description;
            }
        }
    }
}

BodyPartfile* BodyParts3DCatalog40::getFilePointer(const std::string& filename){
    for(size_t i = 0; i < files_.size() ; i++){
        if(files_.at(i)->fileName == filename){
            return files_.at(i);
        }
    }
    return 0;
}



}//namespace voreen
