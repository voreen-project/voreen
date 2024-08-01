/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

//header file
#include "cmhalodescriptor.h"
#include "voreen/core/utils/stringutils.h"

using std::string;
using std::vector;
//we are in namespace voreen
namespace voreen {

CMHaloDescriptor::CMHaloDescriptor()
    : Processor()
    , inport_(Port::INPORT, "haloport.input", "Halo Data Input")
    , outport_(Port::OUTPORT, "haloport.output", "Description of selected Halo")
    , selectedHaloIDProp_("selectedHaloIDProp", "ID of selected halo", 257, -1, 1000000)
    , expressionProperty_("expression", "Expression")
    , dumpButtonProperty_("dumpButton", "Dump to Console", VALID)
    , updateButtonProperty_("updateButton", "Update Meta Data List", VALID)
{
    addPort(inport_);
    addPort(outport_);


    //register properties
    addProperty(selectedHaloIDProp_);
    addProperty(expressionProperty_);
    addProperty(dumpButtonProperty_);
    addProperty(updateButtonProperty_);

    ON_CHANGE(dumpButtonProperty_, CMHaloDescriptor, dumpToConsole);
    ON_CHANGE(updateButtonProperty_, CMHaloDescriptor, updateMetaDataList);
}

CMHaloDescriptor::~CMHaloDescriptor(){
}

void CMHaloDescriptor::initialize() {
}

void CMHaloDescriptor::deinitialize() {
}

void CMHaloDescriptor::process() {

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
void CMHaloDescriptor::dumpToConsole() {
    if (inport_.isReady()) {
        LINFO("Console Dump:\n" + replaceMetaDataAndGetString());
    }
    else
        LWARNING("Cannot dump to console: volume.inport not ready!");
}
//void CMHaloDescriptor::process() {
//    std::stringstream text;
//    text << "ID:" << selectedHalo->ID << std::endl;
//    text << "M:" << selectedHalo->mass << std::endl;
//    text << "P:" << selectedHalo->pos << std::endl;
//    text << "V:" << selectedHalo->velocity << std::endl;
//    text << "A:" << selectedHalo->angularMomenta << std::endl;
//    text << "R:" << selectedHalo->radius << std::endl;
//    outport_.setData(text.str());
//}
void CMHaloDescriptor::updateMetaDataList() {

    //erase the items of the StringExpressionProperty
    expressionProperty_.eraseItems();

    if (inport_.isReady()) {
        expressionProperty_.addPlaceHolder("ID", expressionProperty_.makePlaceHolder("haloID"));
        expressionProperty_.addPlaceHolder("Original ID", expressionProperty_.makePlaceHolder("haloOrigID"));
        expressionProperty_.addPlaceHolder("Mass", expressionProperty_.makePlaceHolder("haloMass"));
        expressionProperty_.addPlaceHolder("Spin", expressionProperty_.makePlaceHolder("haloSpin"));
        expressionProperty_.addPlaceHolder("Position", expressionProperty_.makePlaceHolder("haloPos"));
        expressionProperty_.addPlaceHolder("Velocity", expressionProperty_.makePlaceHolder("haloVelocity"));
        expressionProperty_.addPlaceHolder("Angular Momenta", expressionProperty_.makePlaceHolder("haloAngularMomenta"));
        expressionProperty_.addPlaceHolder("Radius", expressionProperty_.makePlaceHolder("haloRadius"));
        expressionProperty_.addPlaceHolder("Scale Radius", expressionProperty_.makePlaceHolder("haloScaleRadius"));
        expressionProperty_.addPlaceHolder("Black Hole Mass", expressionProperty_.makePlaceHolder("haloBlackHoleMass"));
        expressionProperty_.addPlaceHolder("Black Hole Spin", expressionProperty_.makePlaceHolder("haloBlackHoleSpin"));
        expressionProperty_.addPlaceHolder("Spheroid Radius", expressionProperty_.makePlaceHolder("haloSpheroidRadius"));
        expressionProperty_.addPlaceHolder("Spheroid Mass Gas", expressionProperty_.makePlaceHolder("haloSpheroidMassGas"));
        expressionProperty_.addPlaceHolder("Spheroid Velocity", expressionProperty_.makePlaceHolder("haloSpheroidVelocity"));
        expressionProperty_.addPlaceHolder("Disk Radius", expressionProperty_.makePlaceHolder("haloDiskRadius"));
        expressionProperty_.addPlaceHolder("Disk Mass Gas", expressionProperty_.makePlaceHolder("haloDiskMassGas"));
        expressionProperty_.addPlaceHolder("Disk Velocity", expressionProperty_.makePlaceHolder("haloDiskVelocity"));
    }
}

template <typename T> std::string toString(const T& v) {
    std::ostringstream os;
    os<<v;
    return os.str();
}
std::string CMHaloDescriptor::replaceMetaDataAndGetString() const {
    std::stringstream output;

    const CMHalo* selectedHalo;
    try {
        selectedHalo = getSelectedHalo();
    } catch(...) {
        return "";
    }
    const CMMergerTree* tree = inport_.getData();
    if(!tree) {
        throw "no tree";
        return "";
    }

    //get the Meta Data keys to be replaced
    std::set<string> placeholders = expressionProperty_.getPlaceholdersInText();
    //construct replacement map for these keys
    std::map<string, string> replacements;
    //add the volume information not contained in meta data
    std::stringstream dim;
    replacements.insert(std::make_pair("haloID"             , toString(selectedHalo->ID              )                       ));
    replacements.insert(std::make_pair("haloOrigID"         , toString(selectedHalo->origID          )                       ));
    replacements.insert(std::make_pair("haloMass"           , toString(selectedHalo->mass            ) + " Msun/h"           ));
    replacements.insert(std::make_pair("haloSpin"           , toString(selectedHalo->spinParameter   )                       ));
    replacements.insert(std::make_pair("haloPos"            , toString(selectedHalo->pos             ) + " Mpc/h"            ));
    replacements.insert(std::make_pair("haloVelocity"       , toString(selectedHalo->velocity        ) + " km/s"             ));
    replacements.insert(std::make_pair("haloAngularMomenta" , toString(selectedHalo->angularMomenta  ) + " Msun/h*Mpc/h*km/s"));
    replacements.insert(std::make_pair("haloRadius"         , toString(selectedHalo->radius/1000     ) + " Mpc/h"            ));
    replacements.insert(std::make_pair("haloScaleRadius"    , toString(selectedHalo->scaleRadius/1000) + " Mpc/h"            ));
    if(tree->containsGalacticusData()) {
        replacements.insert(std::make_pair("haloBlackHoleMass"    , toString(selectedHalo->blackHoleMass   ) + " Msun/h"));
        replacements.insert(std::make_pair("haloBlackHoleSpin"    , toString(selectedHalo->blackHoleSpin   )            ));
        replacements.insert(std::make_pair("haloSpheroidRadius"   , toString(selectedHalo->spheroidRadius  ) + " Mpc/h" ));
        replacements.insert(std::make_pair("haloSpheroidMassGas"  , toString(selectedHalo->spheroidMassGas ) + " Msun/h"));
        replacements.insert(std::make_pair("haloSpheroidVelocity" , toString(selectedHalo->spheroidVelocity) + " km/s"  ));
        replacements.insert(std::make_pair("haloDiskRadius"       , toString(selectedHalo->diskRadius      ) + " Mpc/h" ));
        replacements.insert(std::make_pair("haloDiskMassGas"      , toString(selectedHalo->diskMassGas     ) + " Msun/h"));
        replacements.insert(std::make_pair("haloDiskVelocity"     , toString(selectedHalo->diskVelocity    ) + " km/s"  ));
    } else {
        replacements.insert(std::make_pair("haloBlackHoleMass"    , "N/A"));
        replacements.insert(std::make_pair("haloBlackHoleSpin"    , "N/A"));
        replacements.insert(std::make_pair("haloSpheroidRadius"   , "N/A"));
        replacements.insert(std::make_pair("haloSpheroidMassGas"  , "N/A"));
        replacements.insert(std::make_pair("haloSpheroidVelocity" , "N/A"));
        replacements.insert(std::make_pair("haloDiskRadius"       , "N/A"));
        replacements.insert(std::make_pair("haloDiskMassGas"      , "N/A"));
        replacements.insert(std::make_pair("haloDiskVelocity"     , "N/A"));
    }

    //iterate over placeholders to add relevant meta data to the map
    std::set<string>::const_iterator i;

    for (i = placeholders.begin(); i != placeholders.end(); ++i) {
        //check if dot operator is used
        if (i->find(".") != string::npos) {
            //use of dot operator: split the string, get the MetaData and call component-wise toString-method
            size_t position = i->find(".");
            //check for Volume Dimensions (not MetaData, but may be queried component-wise
#define ADD_VECTOR3D_COMPONENTS(name, haloComponent) \
            if (i->substr(0,position) == name) { \
                if (i->substr(position+1) == "x") \
                    replacements.insert(std::make_pair(*i, toString(haloComponent.x)));\
                else if (i->substr(position+1) == "y")\
                    replacements.insert(std::make_pair(*i, toString(haloComponent.y)));\
                else if (i->substr(position+1) == "z")\
                    replacements.insert(std::make_pair(*i, toString(haloComponent.z)));\
                else\
                    replacements.insert(std::make_pair(*i, toString(haloComponent)));\
            }
            ADD_VECTOR3D_COMPONENTS("haloPos", selectedHalo->pos);
            ADD_VECTOR3D_COMPONENTS("haloVelocity", selectedHalo->velocity);
            ADD_VECTOR3D_COMPONENTS("haloAngularMomenta", selectedHalo->angularMomenta);
#undef ADD_VECTOR3D_COMPONENTS
        }
    }

    output << expressionProperty_.replacePlaceHoldersInText(replacements);

    return output.str();
}
const CMHalo* CMHaloDescriptor::getSelectedHalo() const {
    const CMMergerTree* tree = inport_.getData();
    if(!tree) {
        throw "no tree";
    }
    const CMHalo* centerHalo = tree->haloByID(selectedHaloIDProp_.get());
    if(!centerHalo) {
        throw "no halo";
    }
    return centerHalo;
}

} // namespace
