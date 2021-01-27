/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "metadataextractor.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

using std::string;
using std::vector;

namespace voreen {

const std::string MetaDataExtractor::loggerCat_("voreen.base.MetaDataExtractor");

MetaDataExtractor::MetaDataExtractor()
    : Processor(),
    inport_(Port::INPORT, "volume.inport", "Volume Input"),
    outport_(Port::OUTPORT, "text", "Text Output", false),
    expressionProperty_("expression", "Expression"),
    dumpButtonProperty_("dumpButton", "Dump to Console", VALID),
    updateButtonProperty_("updateButton", "Update Meta Data List", VALID)
{
    addPort(inport_);
    addPort(outport_);
    addProperty(expressionProperty_);

    dumpButtonProperty_.onChange(MemberFunctionCallback<MetaDataExtractor>(this,&MetaDataExtractor::dumpToConsole));
    addProperty(dumpButtonProperty_);

    updateButtonProperty_.onChange(MemberFunctionCallback<MetaDataExtractor>(this,&MetaDataExtractor::updateMetaDataList));
    addProperty(updateButtonProperty_);
}

Processor* MetaDataExtractor::create() const {
    return new MetaDataExtractor();
}

bool MetaDataExtractor::isReady() const {
    if (!isInitialized())
        return false;

    return true;
}

void MetaDataExtractor::process() {

    std::stringstream output;

    updateMetaDataList();

    if (inport_.isReady())
        output << replaceMetaDataAndGetString();

    //this is a workaround for graphics problems when rendering empty lines...
    vector<string> lines = strSplit(output.str(),'\n');
    vector<string>::iterator iter;
    for (iter = lines.begin(); iter != lines.end(); ++iter) {
        if (*iter == "")
            iter->append(" ");
    }
    string result = strJoin(lines, "\n");

    if (outport_.isReady())
        outport_.setData(result);
}

void MetaDataExtractor::dumpToConsole() {
    if (inport_.isReady()) {
        LINFO("Console Dump:\n" + replaceMetaDataAndGetString());
    }
    else
        LWARNING("Cannot dump to console: volume.inport not ready!");
}

void MetaDataExtractor::updateMetaDataList() {

    //erase the items of the StringExpressionProperty
    expressionProperty_.eraseItems();

    if (inport_.isReady()) {
        //add some data contained in the volume
        //bits per voxel
        expressionProperty_.addPlaceHolder("Bits per Voxel", expressionProperty_.makePlaceHolder("volBitsPerVoxel"));
        //dimensions
        expressionProperty_.addPlaceHolder("Volume Dimensions", expressionProperty_.makePlaceHolder("volDimensions"));
        //memory size
        expressionProperty_.addPlaceHolder("Volume Memory Size", expressionProperty_.makePlaceHolder("volMemSize"));
        //number of voxels
        expressionProperty_.addPlaceHolder("Number of Voxels", expressionProperty_.makePlaceHolder("volNumVoxels"));
        //path to the volume
        expressionProperty_.addPlaceHolder("Path to Volume", expressionProperty_.makePlaceHolder("volPath"));

        //get all Meta Data keys and add corresponding items to the StringExpressionProptery
        vector<string> keys = inport_.getData()->getMetaDataKeys();

        vector<string>::const_iterator iter;
        for (iter = keys.begin(); iter != keys.end(); ++iter) {
            expressionProperty_.addPlaceHolder(*iter, expressionProperty_.makePlaceHolder(*iter));
        }
    }
}

std::string MetaDataExtractor::replaceMetaDataAndGetString() const {
    std::stringstream output;

    //get the Meta Data keys to be replaced
    std::set<string> placeholders = expressionProperty_.getPlaceholdersInText();
    //construct replacement map for these keys
    std::map<string, string> replacements;
    //add the volume information not contained in meta data
    std::stringstream dim;
    dim << inport_.getData()->getDimensions();
    replacements.insert(std::make_pair("volDimensions", dim.str()));
    replacements.insert(std::make_pair("volBitsPerVoxel", itos(inport_.getData()->getBytesPerVoxel() * 8)));
    replacements.insert(std::make_pair("volNumVoxels", itos(inport_.getData()->getNumVoxels())));
    replacements.insert(std::make_pair("volPath", inport_.getData()->getOrigin().getPath()));
    replacements.insert(std::make_pair("volFilename", inport_.getData()->getOrigin().getFilename()));
    replacements.insert(std::make_pair("volBasename", tgt::FileSystem::baseName(inport_.getData()->getOrigin().getPath())));

    //memory size
    std::stringstream memSize;
    size_t bytes = inport_.getData()->getBytesPerVoxel()*inport_.getData()->getNumVoxels();
    float mb = tgt::round(bytes/104857.6f) / 10.f;    //calculate mb with 0.1f precision
    float kb = tgt::round(bytes/102.4f) / 10.f;
    if (mb >= 0.5f) {
        memSize << mb << " MB";
    }
    else if (kb >= 0.5f) {
        memSize << kb << " kB";
    }
    else {
        memSize << bytes << " bytes";
    }
    replacements.insert(std::make_pair("volMemSize", memSize.str()));

    //iterate over placeholders to add relevant meta data to the map
    std::set<string>::const_iterator i;

    for (i = placeholders.begin(); i != placeholders.end(); ++i) {
        //check if dot operator is used
        if (i->find(".") == string::npos) {
            const MetaDataBase* m = inport_.getData()->getMetaData(*i);
            if (m) {
                //found MetaData -> add std::string representation to map
                replacements.insert(std::make_pair(*i, trim(m->toString())));
            }
        }
        else {
            //use of dot operator: split the string, get the MetaData and call component-wise toString-method
            size_t position = i->find(".");
            const MetaDataBase* m = inport_.getData()->getMetaData(i->substr(0,position));
            if (m) {
                //found MetaData -> add std::string representation with component argument to map
                replacements.insert(std::make_pair(*i, trim(m->toString(i->substr(position+1)))));
            }

            //check for Volume Dimensions (not MetaData, but may be queried component-wise
            if (i->substr(0,position) == "volDimensions") {
                if (i->substr(position+1) == "x")
                    replacements.insert(std::make_pair(*i, itos(inport_.getData()->getDimensions().x)));
                else if (i->substr(position+1) == "y")
                    replacements.insert(std::make_pair(*i, itos(inport_.getData()->getDimensions().y)));
                else if (i->substr(position+1) == "z")
                    replacements.insert(std::make_pair(*i, itos(inport_.getData()->getDimensions().z)));
                else
                    replacements.insert(std::make_pair(*i, dim.str()));
            }

            //check for Memory Size
            if (i->substr(0,position) == "volMemSize") {
                if (i->substr(position+1) == "B")
                    replacements.insert(std::make_pair(*i, itos(bytes)));
                else if (i->substr(position+1) == "kB")
                    replacements.insert(std::make_pair(*i, ftos(kb)));
                else if (i->substr(position+1) == "MB")
                    replacements.insert(std::make_pair(*i, ftos(mb)));
                else
                    replacements.insert(std::make_pair(*i, memSize.str()));
            }
        }
    }

    output << expressionProperty_.replacePlaceHoldersInText(replacements);

    return output.str();
}

}   // namespace
