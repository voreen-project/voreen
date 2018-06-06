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

#include "bodyparts3dsource.h"

#include "tgt/filesystem.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/processors/processorwidget.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen{
const std::string BodyParts3DSource::loggerCat_("voreen.touchtable.BodyParts3D");

BodyParts3DSource::BodyParts3DSource()
                            : Processor()
                            , version_(Version3_0)
                            , versionoption_("BodyPartsVersion","BodyParts Version")
                            , folderPath_("folderPathBodyParts3D","Open Folder","Open", VoreenApplication::app()->getUserDataPath(),"",FileDialogProperty::DIRECTORY)
                            , loadBodyParts_("loadBodyParts", "Load Body Parts")
                            , geometryPort_(Port::OUTPORT,"geometry.output","Geometry output")
                            , pathValid_("path.valid","path valid", false)
{
    loadBodyParts_.onChange(MemberFunctionCallback<BodyParts3DSource>(this, &BodyParts3DSource::updateWidget));
    versionoption_.onChange(MemberFunctionCallback<BodyParts3DSource>(this, &BodyParts3DSource::versionChanged));
    versionoption_.addOption("BodyParts3D 3.0","BodyParts3D 3.0", "BodyParts3D 3.0");
    versionoption_.addOption("BodyParts3D 3.0 HQ","BodyParts3D 3.0 HQ", "BodyParts3D 3.0 HQ");
    versionoption_.addOption("BodyParts3D 4.0","BodyParts3D 4.0", "BodyParts3D 4.0");
    addProperty(versionoption_);
    addProperty(folderPath_);
    addProperty(loadBodyParts_);
    addProperty(pathValid_);
    //pathValid_.setVisible(false);

    addPort(geometryPort_);
}

BodyParts3DSource::~BodyParts3DSource(){
    //free resources
    geometryPort_.clear();
}

bool BodyParts3DSource::isReady()const{
    return geometryPort_.isReady();
}

void BodyParts3DSource::setVersion(BodyPartsVersion v){
    switch(v){
    case Version3_0:
        versionoption_.selectByValue("BodyParts3D 3.0");
        break;
    case Version3_0HQ:
        versionoption_.selectByValue("BodyParts3D 3.0 HQ");
        break;
    case Version4_0:
        versionoption_.selectByValue("BodyParts3D 4.0");
        break;
    }
    version_ = v;
}

void BodyParts3DSource::setPathValid(bool b){
    pathValid_.set(b);
}

void BodyParts3DSource::versionChanged(){
    if(versionoption_.get() == "BodyParts3D 3.0"){
        version_ = Version3_0;
    }else if(versionoption_.get() == "BodyParts3D 3.0 HQ"){
        version_ = Version3_0HQ;
    }else{
        version_ = Version4_0;
    }
}

void BodyParts3DSource::updateWidget(){
    processorWidget_;
    if(getProcessorWidget())
        getProcessorWidget()->updateFromProcessor();
}

void BodyParts3DSource::initialize() {
    Processor::initialize();
    if(pathValid_.get())
        updateWidget();
}

TriangleMeshGeometryBodyParts3D* BodyParts3DSource::readPart(const tgt::vec4& color,const std::string& fileName){
    std::string folder;
    if(version_ == Version3_0){
            folder = folderPath_.get() + "\\BodyParts3D_3.0_obj_99\\" + fileName;
    }else if(version_ == Version3_0HQ){
            folder = folderPath_.get() + "\\BodyParts3D_3.0_obj_95\\" + fileName;
    }else{
        folder = folderPath_.get() + "\\isa_BP3D_4.0_obj_99\\" + fileName;
    }
    tgt::FileSystem sys;
    if(sys.fileExists(folder + ".bpvrn")){
        std::string url = folder + ".bpvrn";
        return readBpvrn(color, url);
    }
    else{
        return readObj(color, folder);
    }
}

TriangleMeshGeometryBodyParts3D* BodyParts3DSource::readObj(tgt::vec4 color, std::string& path){
    std::vector<tgt::vec3> vertices;
    std::vector<tgt::vec3> normals;

    TriangleMeshGeometryBodyParts3D* mesh = new TriangleMeshGeometryBodyParts3D();
    mesh->setColor(color);
    size_t buffersize = sizeof(tgt::vec3);

    //correct face offset
    vertices.push_back(tgt::vec3(0,0,0));
    normals.push_back(tgt::vec3(0,0,0));


    tgt::File* file = FileSys.open(path + ".obj");
    if(!file)
        throw tgt::CorruptedFileException("Cant read file",path + ".obj");

    //create new file in binaryformat
    std::ofstream stream;
    stream.open((path.append(".bpvrn")).c_str(),std::ios_base::out | std::ios_base::binary);

    while(!file->eof()){
        std::string line = file->getLine();
        std::vector<std::string> prt = strSplit(line, ' ');
        if(prt.size() > 0){
            if( prt[0][0] == '#' ){
                if(line.size() >= 128){ //due to a bug in getLine() you cant read a line with exactly 128 characters so you
                                        //need to check if a comment swallowed another line by accident
                    prt = strSplit(line, '\n');
                    if(prt.size() > 1){
                        line = prt[1];
                        prt = strSplit(line, ' ');
                    }
                }
            }
            if(prt[0] == "vn"){
                //read vector normal
                if(prt.size() >= 3){
                    normals.push_back(tgt::vec3(stof(prt[1]),stof(prt[2]),stof(prt[3])));
                }else{
                    throw tgt::CorruptedFileException("Not Enough Data");
                }
            }else if(prt[0] == "v"){
                //read vector position
                if(prt.size() >= 3){
                    vertices.push_back(tgt::vec3(stof(prt[1]),stof(prt[2]),stof(prt[3])));
                }else{
                    throw tgt::CorruptedFileException("Not Enough Data");
                }
            }else if(prt[0] == "f"){
                //read face (Triangle), write its data into the new .bpvrn file in binary data
                //and add it to geometry
                if(prt.size() >= 3){
                    std::vector<std::string> faceVertex = strSplit(prt[1],'/');
                    int index = stoi(faceVertex[0]);
                    VertexNormal vertex1(vertices[index],normals[index]);
                    stream.write(reinterpret_cast<char*>(&vertices[index]),buffersize);
                    stream.write(reinterpret_cast<char*>(&normals[index]),buffersize);
                    faceVertex = strSplit(prt[2],'/');
                    index = stoi(faceVertex[0]);
                    VertexNormal vertex2(vertices[index],normals[index]);
                    stream.write(reinterpret_cast<char*>(&vertices[index]),buffersize);
                    stream.write(reinterpret_cast<char*>(&normals[index]),buffersize);
                    faceVertex = strSplit(prt[3],'/');
                    index = stoi(faceVertex[0]);
                    VertexNormal vertex3(vertices[index],normals[index]);
                    stream.write(reinterpret_cast<char*>(&vertices[index]),buffersize);
                    stream.write(reinterpret_cast<char*>(&normals[index]),buffersize);

                    mesh->addTriangle(Triangle<VertexNormal>(vertex1,vertex2,vertex3));
                }else{
                    throw tgt::CorruptedFileException("Not Enough Data");
                }
            }
        }else{
            throw tgt::CorruptedFileException("Not Enough Data");
        }
    }
    stream.close();
    file->close();
    delete file;
    std::remove((path + ".obj").c_str());
    return mesh;
}


TriangleMeshGeometryBodyParts3D* BodyParts3DSource::readBpvrn(tgt::vec4 color, std::string& path){
    tgt::vec3 vertex;
    tgt::vec3 normal;
    size_t bufferSize = sizeof(tgt::vec3);
    std::ifstream stream;

    TriangleMeshGeometryBodyParts3D* mesh = new TriangleMeshGeometryBodyParts3D();
    mesh->setColor(color);
    stream.open(path.c_str(),std::ios_base::in | std::ios_base::binary);
    //read 12 bytes at a time and cast them into a tgt::vec3
    while(true){
        if(!stream.read(reinterpret_cast<char*>(&vertex), bufferSize))
            break;
        if(!stream.read(reinterpret_cast<char*>(&normal), bufferSize))
            break;
        VertexNormal vertex1(vertex,normal);

        if(!stream.read(reinterpret_cast<char*>(&vertex), bufferSize))
            break;
        if(!stream.read(reinterpret_cast<char*>(&normal), bufferSize))
            break;
        VertexNormal vertex2(vertex,normal);

        if(!stream.read(reinterpret_cast<char*>(&vertex), bufferSize))
            break;
        if(!stream.read(reinterpret_cast<char*>(&normal), bufferSize))
            break;
        VertexNormal vertex3(vertex,normal);

        mesh->addTriangle(Triangle<VertexNormal>(vertex1,vertex2,vertex3));
    }
    stream.close();
    return mesh;
}

void BodyParts3DSource::createMesh(const std::vector<std::string>& filenames, const std::vector<tgt::vec4>& colors, const std::vector<std::string>& descriptions){

    TriangleMeshGeometryCollectionBodyParts3D* meshCatalog = new TriangleMeshGeometryCollectionBodyParts3D();
    for(int i = 0; i < filenames.size() ; ++i){
        meshCatalog->addMesh(readPart(colors[i],filenames[i]));
        meshCatalog->getMesh(static_cast<size_t>(i))->setTag(descriptions[i]);//todo: re-add description to geometries
        setProgress(static_cast<float>(i)/static_cast<float>(filenames.size()-1));
    }
    tgt::Bounds bb = meshCatalog->getBoundingBox(false);
    float maxSideLength = std::max(std::abs(bb.getURB().x - bb.getLLF().x),
                        std::max(std::abs(bb.getURB().y - bb.getLLF().y),
                                std::abs(bb.getURB().z - bb.getLLF().z)));

    meshCatalog->setTransformationMatrix(tgt::mat4::createScale(tgt::vec3(2/maxSideLength))*tgt::mat4::createTranslation(-bb.center()));

    /*
    TriangleMeshGeometrySingleColor* mesh = new TriangleMeshGeometrySingleColor();
    float y =  (bb.getURB().y - bb.center().y)/1000.0f;
    mesh->addQuad(VertexNormal(3000.0f*tgt::vec3(-1.0f, y, -1.0f), tgt::vec3(0.0f,1.0f,0.0f)),
        VertexNormal(3000.0f*tgt::vec3(1.0f, y, -1.0f), tgt::vec3(0.0f,1.0f,0.0f)),
        VertexNormal(3000.0f*tgt::vec3(1.0f, y, 1.0f), tgt::vec3(0.0f,1.0f,0.0f)),
        VertexNormal(3000.0f*tgt::vec3(-1.0f, y, 1.0f), tgt::vec3(0.0f,1.0f,0.0f)));

    mesh->setColor(tgt::vec4(1.0f));
    meshCatalog->addMesh(mesh);
    */

    geometryPort_.clear();
    if(!meshCatalog->isEmpty())
        geometryPort_.setData(meshCatalog);
}

void BodyParts3DSource::process(){

}

} //namespace voreen
