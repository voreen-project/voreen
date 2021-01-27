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

#include "poicsvimport.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/io/serialization/serializable.h"

#include <fstream>
namespace voreen{
static int read_int(char** b){
    int sign = 1;
    int num = 0;

    if (**b == '-'){
        sign = -1;
        (*b)++;
    }

    if (!isdigit(**b)){
        throw VoreenException("expected digit");
    }

    while(isdigit(**b)){
        int d = **b-'0';
        num = num*10+d;
        (*b)++;
    }

    return sign*num;
}

static float read_float(char** b){
    float sign = 1.0;
    if (**b == '-'){
        sign =  -1.0;
        (*b)++;
    }

    float num = (float)read_int(b);
    if (**b != '.')
        return num;
    else
        (*b)++;

    float mul = 0.1f;
    while(isdigit(**b)){
        int d = **b-'0';
        num += mul*d;
        mul *= 0.1f;
        (*b)++;
    }
    return num*sign;
}

static void skip_spaces(char** b){
    while(**b == ' ')
        (*b)++;
}

static void read_char(char** b, char c){
    if (**b != c)
        throw VoreenException("Unexpected char");
    (*b)++;
}

void line_end(char** b){
    if (**b == 0)
        return;
    if (**b == '\r')
        read_char(b, '\r');
    read_char(b, '\n');
}

static void skip_line(char** b){
    while(**b && **b != '\r' && **b != '\n')
        (*b)++;
    line_end(b);

}

std::string read_string_to_nl(char** b){
    std::string str;
    while (**b && **b != '\r' && **b != '\n'){
        str.push_back(**b);
        (*b)++;
    }
    return str;
}

static void read_poi(char** b, POIList* pois){
    skip_spaces(b);
    if (!isdigit(**b)){
        skip_line(b);
        return;
    }else{
        int id = read_int(b); skip_spaces(b); read_char(b, ','); skip_spaces(b);
        tgt::vec3 pos;
        pos.x = read_float(b); skip_spaces(b); read_char(b, ','); skip_spaces(b);
        pos.y = read_float(b); skip_spaces(b); read_char(b, ','); skip_spaces(b);
        pos.z = read_float(b); skip_spaces(b); read_char(b, ','); skip_spaces(b);
        tgt::vec3 pos_vx;
        pos_vx.x = read_float(b); skip_spaces(b); read_char(b, ','); skip_spaces(b);
        pos_vx.y = read_float(b); skip_spaces(b); read_char(b, ','); skip_spaces(b);
        pos_vx.z = read_float(b); skip_spaces(b); read_char(b, ','); skip_spaces(b);
        std::string group = read_string_to_nl(b);
        skip_line(b);
        POIGroupID gid = pois->addGroup(group);
        pois->addPoint(pos, gid);
    }
}


POICSVImport::POICSVImport()
    : outport_(Port::OUTPORT, "outport", "POI Outport", true)
    , fileProp_("fileProp", "Path to load", "Path for POI loading", ".", "Comma seperated values (*.csv)", FileDialogProperty::OPEN_FILE)
    , autoLoad_("autoLoad", "Auto load on change")
    , loadButton_("loadButton", "Load")
    , shouldLoadOnProcess_(false)
{
    addPort(outport_);

    addProperty(fileProp_);
    addProperty(autoLoad_);
    addProperty(loadButton_);

    ON_CHANGE_LAMBDA(loadButton_, [this]{
        shouldLoadOnProcess_ = true;
    });
}

Processor* POICSVImport::create() const
{
    return new POICSVImport();
}

std::string POICSVImport::getClassName() const
{
    return "POICSVImport";
}

std::string POICSVImport::getCategory() const
{
    return "Points of Interest";
}

void POICSVImport::setDescriptions()
{
    setDescription("The POICSVImport processor imports data from a CSV file similar to the POISource "
                  "processor. The coordinates in voxel space are ignored and only the world space positions "
                  "are used.");
    fileProp_.setDescription("File to import data from.");
    autoLoad_.setDescription("Automatical data loading");
    loadButton_.setDescription("Load now.");
}

void POICSVImport::process()
{
    if (autoLoad_.get() || shouldLoadOnProcess_){
        shouldLoadOnProcess_ = false;
        std::string fileName = fileProp_.get();
        POIList *pois= new POIList;
        std::ifstream file(fileName);
        file.seekg(0, std::ios_base::end);
        std::vector<char> data;
        data.resize(file.tellg());
        file.seekg(0, std::ios_base::beg);
        file.read(data.data(), data.size());
        data.push_back(0);

        char* b = data.data();
        while(*b){
            read_poi(&b, pois);
        }
        outport_.setData(pois, true);
    }
}
}