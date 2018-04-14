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

#include "bodypartsport.h"



namespace voreen {

BodyPartsPort::BodyPartsPort(PortDirection direction, const std::string& name, const std::string& guiName,
    bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
    : GenericPort<std::vector<std::pair<TriangleMeshGeometrySingleColor*,std::string> > >(direction, name, guiName, allowMultipleConnections, invalidationLevel)
{}

tgt::col3 BodyPartsPort::getColorHint() const{
    return tgt::col3(128, 145, 128);
}
}

/*
namespace voreen{

const std::string TextListPort::loggerCat_("voreen.TextListPort");

TextListPort::TextListPort(PortDirection direction, const std::string& id, const std::string& guiName,
                    bool allowMultipleConnections,
                    Processor::InvalidationLevel invalidationLevel)
                    : Port(direction, id, guiName, allowMultipleConnections, invalidationLevel){
    portData_ = 0;
    owner_ = false;
}

TextListPort::~TextListPort(){

}

bool TextListPort::hasData() const{
    return portData_;
}

void TextListPort::forwardData() const{
    for(std::vector<Port*>::const_iterator it = forwardPorts_.begin(); it != forwardPorts_.end(); ++it){
        dynamic_cast<TextListPort*>(*it)->setData(getData());
    }
}

void TextListPort::setData(std::vector<std::string>* data, bool owner){
    if(owner_){
        delete portData_;
    }
    owner_ = owner;
    portData_ = data;
    forwardData();
}

std::vector<std::string>* TextListPort::getData()const{
    return portData_;
}

bool TextListPort::isReady()const{
    if (isOutport())
        return isConnected();
    else
        return (!getConnected().empty() && hasData() && areConditionsMet());
}

void TextListPort::clear(){
    if(!isOutport())
        LERROR("called clear() on inport");
    else{
        if(owner_){
            delete portData_;
        }else{
            portData_ = 0;
        }
    }
}

tgt::col3 TextListPort::getColorHint()const{
    return tgt::col3(0,255,255);
}
}*/
