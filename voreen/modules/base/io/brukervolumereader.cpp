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

#include "brukervolumereader.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "tgt/exception.h"
#include "tgt/vector.h"
#include "tgt/filesystem.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatortranspose.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

using tgt::vec3;
using tgt::ivec3;
using std::vector;
using std::string;

namespace voreen {

std::string JCampToken::toString() const {
    string st;
    switch(type_) {
        case OPENING_BRACKET:
            return "(";
        case CLOSING_BRACKET:
            return ")";
        case COMMA:
            return ",";
        case NUMBER:
            return "NUM[" + content_ + "]";
        case STRING:
            return "STRING[" + content_ + "]";
        case ENUM:
            return "ENUM[" + content_ + "]";
    }
    return "";
}

bool JCampToken::isInteger() const {
    if(content_.find('.') != string::npos)
        return false;
    if(content_.find('e') != string::npos)
        return false;

    return true;
}

int JCampToken::getInteger() const {
    if(type_ != NUMBER)
        throw JCampParserException("Not a number!");

    if(!isInteger())
        throw JCampParserException("Not an integer!");

    return stoi(content_);
}

float JCampToken::getFloat() const {
    if(type_ != NUMBER)
        throw JCampParserException("Not a number!");

    return stof(content_);
}

const std::string JCampElement::loggerCat_ = "voreen.base.JCampElement";

std::vector<JCampToken> JCampElement::tokenize(std::string data) {
    vector<JCampToken> tokens;
    string whitespaceList = " \t\n\r\0\x0B";
    string numStartList = "0123456789-";
    string numList = "0123456789-.e";
    string stopList = ",()<>";

    size_t pos = 0;
    while(pos < data.size()) {
        //skip whitespaces
        while(whitespaceList.find(data[pos]) != std::string::npos)
            pos++;

        if(data[pos] == '<') {
            //string:
            pos++;

            string curToken;
            while((pos < data.size()) && (data[pos] != '>')) {
                curToken += data[pos];
                pos++;
            }
            pos++;

            curToken = strReplaceAll(curToken, "\n", "");
            tokens.push_back(JCampToken(JCampToken::STRING, curToken));
        }
        else if(numStartList.find(data[pos]) != std::string::npos) {
            string curToken;
            curToken += data[pos];
            pos++;

            while( (pos < data.size()) && (numList.find(data[pos]) != std::string::npos)) {
                curToken += data[pos];
                pos++;
            }

            tokens.push_back(JCampToken(JCampToken::NUMBER, curToken));
        }
        else if(data[pos] == '(') {
            tokens.push_back(JCampToken(JCampToken::OPENING_BRACKET, "("));
            pos++;
        }
        else if(data[pos] == ')') {
            tokens.push_back(JCampToken(JCampToken::CLOSING_BRACKET, ")"));
            pos++;
        }
        else if(data[pos] == ',') {
            tokens.push_back(JCampToken(JCampToken::COMMA, ","));
            pos++;
        }
        else {
            //enum:
            string curToken;
            while( (pos < data.size()) && (whitespaceList.find(data[pos]) == std::string::npos) && (stopList.find(data[pos]) == std::string::npos)){
                curToken += data[pos];
                pos++;
            }

            tokens.push_back(JCampToken(JCampToken::ENUM, curToken));
        }

    }
    return tokens;
}

void JCampElement::addLine(std::string line) {
    lines_.push_back(line);
}

bool JCampElement::parse() {
    if(lines_.empty())
        return false;

    size_t eq = lines_[0].find('=');

    if(eq != string::npos) {
        name_ = lines_[0].substr(0, eq);

        if(name_.find("##$") != std::string::npos) {
            private_ = true;
            name_ = strReplaceAll(name_, "##$", "");
        }
        else if(name_.find("##") != std::string::npos) {
            private_ = false;
            name_ = strReplaceAll(name_, "##", "");
        }
        else {
            LERROR("Failed to parse first line of element (no ## found)");
            return false;
        }

        name_ = trim(name_);

        string firstLineData = lines_[0].substr(eq+1);
        firstLineData = trim(firstLineData);
        lines_[0] = firstLineData;
        string data = strJoin(lines_, "\n");

        vector<JCampToken> tokens = tokenize(data);

        if(!tokens.empty()) {
            if(tokens.size() == 1) {
                LDEBUG("Name: [" << name_ << "] : " << tokens[0].toString());
                data_.push_back(tokens[0]);
                return true;
            }
            else {
                if(tokens.front().getType() == JCampToken::OPENING_BRACKET) {
                    bool isStruct = true;
                    if(tokens.back().getType() == JCampToken::CLOSING_BRACKET) {
                        for(size_t i=1; (i+1)<tokens.size(); i++) {
                            if((tokens[i].getType() == JCampToken::OPENING_BRACKET) || (tokens[i].getType() == JCampToken::CLOSING_BRACKET) )
                                isStruct = false;
                        }
                    }
                    else
                        isStruct = false;


                    if(isStruct) {
                        for(size_t i=1; (i+1)<tokens.size(); i++) {
                            if((tokens[i].getType() == JCampToken::NUMBER) || (tokens[i].getType() == JCampToken::STRING) ||(tokens[i].getType() == JCampToken::ENUM) )
                                data_.push_back(tokens[i]);
                        }

                        dims_.push_back(static_cast<int>(data_.size()));
                        LDEBUG("Name: [" << name_ << "] : STRUCT(" << dims_[0] << ")");
                        return true;
                    }
                    else {
                        //array (or string)
                        size_t i=1;

                        std::stringstream strstr;
                        while( (i<tokens.size()) && (tokens[i].getType() != JCampToken::CLOSING_BRACKET) ) {
                            if(tokens[i].getType() == JCampToken::NUMBER) {
                                if(tokens[i].isInteger()) {
                                    dims_.push_back(tokens[i].getInteger());

                                    if(dims_.size() == 1)
                                        strstr << dims_.back();
                                    else
                                        strstr << "x" << dims_.back();
                                }
                                else {
                                    LERROR("Expected integer!");
                                    return false;
                                }
                            }
                            else if( (tokens[i].getType() != JCampToken::COMMA) && (tokens[i].getType() != JCampToken::CLOSING_BRACKET) ) {
                                LERROR("Unexpected Token in array length indicator: " << tokens[i].toString());
                                return false;
                            }

                            i++;
                        }
                        i++; //move past closing bracket

                        size_t expectedElements = 1;
                        for(size_t j=0; j<dims_.size(); j++) {
                            expectedElements *= dims_[j];
                        }

                        LDEBUG("Name: [" << name_ << "] : ARRAY(" << strstr.str() << ") expected: " << expectedElements);

                        //check if we have an array of strucs:
                        bool structs = false;
                        for(size_t j=i; j<tokens.size(); j++) {
                            switch(tokens[j].getType()) {
                                case JCampToken::OPENING_BRACKET:
                                case JCampToken::CLOSING_BRACKET:
                                case JCampToken::COMMA:
                                    structs = true;
                                    break;
                                case JCampToken::NUMBER:
                                case JCampToken::ENUM:
                                case JCampToken::STRING:
                                    break;
                            }
                        }
                        if(structs) {
                            std::vector<int> structDims;

                            while(i<tokens.size()) {
                                switch(tokens[i].getType()) {
                                    case JCampToken::OPENING_BRACKET:
                                        structDims.push_back(0); //new struct
                                        break;
                                    case JCampToken::CLOSING_BRACKET:
                                        break;
                                    case JCampToken::COMMA:
                                        break;
                                    case JCampToken::NUMBER:
                                    case JCampToken::ENUM:
                                    case JCampToken::STRING:
                                        if(!structDims.empty()) {
                                            structDims.back()++;
                                            data_.push_back(tokens[i]);
                                        }
                                        break;
                                }

                                i++;
                            }
                            int l = structDims[0];
                            for(size_t k=0; k<structDims.size(); k++) {
                                if(structDims[k] != l) {
                                    LWARNING("Mismatch in array of structs!");
                                    return false;
                                }
                            }
                            LDEBUG("Array of structs: " << structDims.size() << "x" << l);
                            //We have one more dimension:
                            //dims_.insert(dims_.begin(), l);
                            dims_.push_back(l);

                            expectedElements = 1;
                            for(size_t j=0; j<dims_.size(); j++) {
                                expectedElements *= dims_[j];
                            }
                        }
                        else {
                            //Check if we have an array of strings:
                            bool strings = true;
                            while(i<tokens.size()) {
                                switch(tokens[i].getType()) {
                                    case JCampToken::OPENING_BRACKET:
                                    case JCampToken::CLOSING_BRACKET:
                                    case JCampToken::COMMA:
                                        return false;
                                    case JCampToken::NUMBER:
                                    case JCampToken::ENUM:
                                        strings = false;
                                    case JCampToken::STRING:
                                        data_.push_back(tokens[i]);
                                }

                                i++;
                            }
                            if(strings) {
                                //the last dimension was actually a string length indicator:
                                dims_.pop_back();

                                expectedElements = 1;
                                for(size_t j=0; j<dims_.size(); j++) {
                                    expectedElements *= dims_[j];
                                }
                            }
                        }

                        if(data_.size() != expectedElements) {
                            LWARNING("Expected " << expectedElements << " elements, got " << (int)data_.size());

                            if(data_.size() < expectedElements)
                                return false;
                        }
                        return true;
                    }
                }
                else {
                    //Handle enums like "2D" and non-private string tags
                    dims_.clear();

                    if(!isPrivate()) {
                        data_.push_back(JCampToken(JCampToken::STRING, data));
                    }
                    else {
                        data_.push_back(JCampToken(JCampToken::ENUM, data));
                    }

                    LDEBUG("Name: [" << name_ << "] : " << data_[0].toString());
                    return true;
                }
            }
        }
        else
            return false;
    }
    else {
        LERROR("Failed to parse element! [" << lines_[0] << "]");
        return false;
    }
}

//-----------------------------------------------------------------------------

const std::string JCampParser::loggerCat_ = "voreen.base.JCampParser";

bool JCampParser::isCommentLine(const std::string& line) {
    if((line.length() >= 2) && (line[0] == '$') && (line[0] == '$'))
        return true;
    else
        return false;
}

bool JCampParser::isParameterLine(const std::string& line) {
    if((line.length() >= 2) && (line[0] == '#') && (line[1] == '#')/* && (line[2] == '$')*/) // "$" denotes private tag
        return true;
    else
        return false;
}

bool JCampParser::parseFile(std::string filename) {
    std::ifstream methodReader;
    methodReader.open(filename.c_str(), std::ios::in);

    if (methodReader.fail()) {
        LERROR("Could not open file " << filename);
        return false;
    }
    LINFO("Opened file " << filename);

    std::string line;

    vector<string> lines;
    // read each line from input file, skip comments
    while (methodReader.good() && std::getline(methodReader, line)) {
        if (line.empty())
            continue;

        if(isCommentLine(line))
            continue;

        lines.push_back(trim(line)); //trim line and add it
    }
    vector<JCampElement> elements;

    //split into elements:
    size_t l = 0;
    bool done = false;
    while(!done) {
        string line = lines[l];

        if(isParameterLine(line)) {
            JCampElement elem;
            elem.addLine(line);

            elements.push_back(elem);
        }
        else {
            if(!elements.empty())
                elements.back().addLine(line);
            else {
                LWARNING("Skipping line: " << line);
            }
        }

        l++;

        if(l >= lines.size())
            done = true;
    }

    for(size_t i=0; i<elements.size(); i++) {
        if(elements[i].parse())
            elements_[elements[i].getName()] = elements[i];
    }

    return true;
}

//-----------------------------------------------------------------------------

const std::string BrukerVolumeReader::loggerCat_ = "voreen.base.BrukerVolumeReader";

BrukerVolumeReader::BrukerVolumeReader(ProgressBar* progress /*= 0*/)
    : VolumeReader(progress)
{
    filenames_.push_back("2dseq");
}

VolumeBase* BrukerVolumeReader::read(const VolumeURL& origin) {
    VolumeBase* result = 0;

    int volumeId = -1;
    std::string tmp = origin.getSearchParameter("volumeId");
    if (! tmp.empty())
        volumeId = stoi(tmp);

    VolumeList* collection = read(origin.getPath(), volumeId);

    if (collection && collection->size() == 1) {
        result = collection->first();
    }
    else if (collection && collection->size() > 1) {
        delete collection;
        throw tgt::FileException("Only one volume expected", origin.getPath());
    }

    delete collection;

    return result;
}

VolumeList* BrukerVolumeReader::read(const std::string &url) {
    VolumeURL origin(url);
    int volumeId = -1;
    std::string tmp = origin.getSearchParameter("volumeId");
    if (! tmp.empty())
        volumeId = stoi(tmp);

    return read(url, volumeId);
}

VolumeList* BrukerVolumeReader::read(const std::string &url, int volumeId) {
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    // The 2dseq file only contains the raw data.
    // The meta information is to be found in several other files around:
    std::string path = tgt::FileSystem::dirName(fileName);
    std::string runPath = tgt::FileSystem::parentDir(tgt::FileSystem::parentDir(path));
    std::string recoFile = path + "/reco";
    std::string methodFile = runPath + "/method";
    std::string acqpFile = runPath + "/acqp";
    std::string visuParsFile = path + "/visu_pars";

    std::string objectFilename = fileName;
    vec3 sliceThickness(1.f);
    std::string format;
    std::string model;
    bool error = false;
    bool bigEndian = false;
    int numFrames = 1;
    int recoTranspose = -1;

    LINFO("BrukerVolumeReader: " << fileName);
    std::string type;
    std::istringstream args;

    // The Bruker file format is a mess: We need to collect all necessary data from at least
    // three different files. We start with the acqp (ACQuisition Parameters) file:
    LINFO("  parsing acqp file: " << acqpFile);
    JCampParser acqpParser;
    if(acqpParser.parseFile(acqpFile)) {
        try {
            numFrames = acqpParser.getInteger("NR");
            LDEBUG("nRuns: " << numFrames);
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }

        try {
            int list_size = acqpParser.getInteger("ACQ_ns_list_size");
            if (list_size != 1)
                LERROR("Unexpected value in field ACQ_ns_list_size: " << list_size << ". Results may not be as expected!");
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }
    }

    // Now parse the reco file
    LINFO("  parsing reco file: " << recoFile);

    JCampParser recoParser;
    if(recoParser.parseFile(recoFile)) {
        try {
            JCampElement e = recoParser.getElement("RECO_transposition");
            e.verifyNumDimensions(1);

            // Bruker file format allows to define transposition (axis swap) for each slice.
            // Currently we only support homogeneous transposition, meaning same transposition for each slice:
            int temp;
            for (size_t i = 0; i < e.getSize(0); ++i) {
                temp = e.getData(i).getInteger();
                if (recoTranspose == -1)
                    recoTranspose = temp;
                else if (recoTranspose != temp) {
                    LERROR("Inhomogeneous transposition values detected in reco file - this feature is currently not supported.");
                    recoTranspose = 0;
                    continue;
                }
            }
            LDEBUG("Transposition: " << recoTranspose);
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }

        try {
            std::string s = recoParser.getEnum("RECO_byte_order");

            if (s == "big-endian" || s == "bigendian" || s == "bigEndian")
                bigEndian = true;
            else if (s == "little-endian" || s == "littleendian" || s == "littleEndian")
                bigEndian = false;
            else
                LERROR("Unknown byte order: " << s);

            LDEBUG("byte order: " << bigEndian);
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }

    }

    // now parse the method file
    // DWI meta data for Bruker images is to be parsed from the 'method' file.
    // An official documentation can be found at http://filer.case.edu/vxs33/pvman/A/Docs/A06_MethDescr.pdf.
    //
    // When parsing DWI meta data from Bruker method files at least two fields seem to be important:
    //  - PVM_DwGradVec:    Diffusion Gradient Vector - Parameter of dimension N × 3 specifying
    //                      the diffusion gradient amplitudes in the x,y,z coordinate system.
    //                      This parameter is used in the pulse-program part of the DTI module in
    //                      the case of direct scaled switching.
    //  - PVM_DwEffBval:    Eff. B Value (PVM_DwEffBval) - The trace of the b-matrix (is implemented as an
    //                      array of size N. In homogenous isotropic media, this parameter may be used
    //                      for a mono exponential fit of the signal as a function of the b-value.
    //                      There are slight differences between the b-values as specified in PVM_DwBvalEach
    //                      and the values of this parameter since the effect of the imaging gradients
    //                      is now considered.

    std::vector<float> PVM_DwEffBval;
    std::vector<tgt::dvec3> PVM_DwGradVec;
    int numRepetitions = 1;

    JCampParser methodParser;
    if(methodParser.parseFile(methodFile)) {
        try {
            PVM_DwEffBval = methodParser.getFloatArray("PVM_DwEffBval");
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }

        try {
            numRepetitions = methodParser.getInteger("PVM_NRepetitions");
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }

        try {
            JCampElement e = methodParser.getElement("PVM_DwGradVec");
            e.verifyNumDimensions(2);
            size_t numEntries = e.getSize(0);
            e.verifySize(1, 3);

            for (size_t i = 0; i < numEntries; ++i) {
                tgt::dvec3 grad(0.f);
                grad.x = e.getData(i, 0).getFloat();
                grad.y = e.getData(i, 1).getFloat();
                grad.z = e.getData(i, 2).getFloat();
                PVM_DwGradVec.push_back(grad);
            }
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }
    }
    else {
        LERROR("Could not open gradient file..");
        error = true;
    }

    if (PVM_DwEffBval.size() != PVM_DwGradVec.size()) {
        LERROR("Number of b-values does not match number of gradients.");
        error = true;
    }

    // now parse the visu_pars file
    // DWI meta data for Bruker images is to be parsed from the 'method' file.
    // An official documentation can be found at http://filer.case.edu/vxs33/pvman/A/Docs/A06_MethDescr.pdf.
    //
    // When parsing DWI meta data from Bruker method files at least two fields seem to be important:
    //  - PVM_DwGradVec:    Diffusion Gradient Vector - Parameter of dimension N × 3 specifying
    //                      the diffusion gradient amplitudes in the x,y,z coordinate system.
    //                      This parameter is used in the pulse-program part of the DTI module in
    //                      the case of direct scaled switching.
    //  - PVM_DwEffBval:    Eff. B Value (PVM_DwEffBval) - The trace of the b-matrix (is implemented as an
    //                      array of size N. In homogenous isotropic media, this parameter may be used
    //                      for a mono exponential fit of the signal as a function of the b-value.
    //                      There are slight differences between the b-values as specified in PVM_DwBvalEach
    //                      and the values of this parameter since the effect of the imaging gradients
    //                      is now considered.
    std::vector<float> mappingOffset;
    std::vector<float> mappingScale;

    tgt::ivec2 resolution = tgt::ivec2::zero;
    JCampParser visuParsParser;
    if(visuParsParser.parseFile(visuParsFile)) {
        try {
            resolution = visuParsParser.getIVec2("VisuCoreSize");
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }
        try {
            format = visuParsParser.getEnum("VisuCoreWordType");

            if (format == "_8BIT_SGN_INT") {
                format = "CHAR";
                model = "I";
            } else if (format == "_8BIT_USGN_INT") {
                format = "UCHAR";
                model = "I";
            } else if (format == "_16BIT_SGN_INT") {
                format = "SHORT";
                model = "I";
            } else if (format == "_16BIT_USGN_INT") {
                format = "USHORT";
                model = "I";
            } else if (format == "_32BIT_USGN_INT") {
                format = "UINT";
                model = "I";
            } else if (format == "_32BIT_SGN_INT") {
                format = "INT";
                model = "I";
            } else if (format == "_32BIT_FLOAT") {
                format = "FLOAT";
                model = "I";
            } else{
                LERROR("Unknown data format: " << format);
            }
            LDEBUG("format: " << format);
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }

        try {
            sliceThickness[2] = visuParsParser.getFloat("VisuCoreFrameThickness");
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }

        try {
            sliceThickness.xy() = visuParsParser.getVec2("VisuCoreExtent");
            sliceThickness[0] /= static_cast<float>(resolution[0]);
            sliceThickness[1] /= static_cast<float>(resolution[1]);
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }
        LDEBUG("spacing: " << sliceThickness);

        try {
            mappingOffset = visuParsParser.getFloatArray("VisuCoreDataOffs");
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }
        try {
            mappingScale = visuParsParser.getFloatArray("VisuCoreDataSlope");
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
        }
    }
    else{
        LERROR("Could not open visu_pars file.");
        error = true;
    }

    if (PVM_DwEffBval.size() != PVM_DwGradVec.size()) {
        LERROR("Number of b-values does not match number of gradients.");
        error = true;
    }
    else {
        size_t tmp = PVM_DwEffBval.size();
        for (int rep = 1; rep < numRepetitions; ++rep) {
            for (size_t img = 0; img < tmp; ++img) {
                PVM_DwEffBval.push_back(PVM_DwEffBval[img]);
                PVM_DwGradVec.push_back(PVM_DwGradVec[img]);
            }
        }
    }

    vector<SliceCollection> sliceCols = listSliceCollections(url);

    if (!error) {
        // recoTransposition values define axis swapping (LHS/RHS coordinate system transformation):
        //   0: no swapping
        //   1: swap x/y axis
        //   2: swap y/z axis
        //   3: swap x/z axis
        if (recoTranspose == 1) {
            for (std::vector<tgt::dvec3>::iterator it = PVM_DwGradVec.begin(); it != PVM_DwGradVec.end(); ++it)
                std::swap(it->x, it->y);
        }
        else if (recoTranspose == 2) {
            for (std::vector<tgt::dvec3>::iterator it = PVM_DwGradVec.begin(); it != PVM_DwGradVec.end(); ++it)
                std::swap(it->y, it->z);
        }
        else if (recoTranspose == 3) {
            for (std::vector<tgt::dvec3>::iterator it = PVM_DwGradVec.begin(); it != PVM_DwGradVec.end(); ++it)
                std::swap(it->x, it->z);
        }

        // do we have a relative path?
        if ((objectFilename.substr(0,1) != "/")  && (objectFilename.substr(0,1) != "\\") &&
                (objectFilename.substr(1,2) != ":/") && (objectFilename.substr(1,2) != ":\\"))
        {
            size_t p = fileName.find_last_of("\\");
            if (p == std::string::npos)
                p = fileName.find_last_of("/");

            // construct path relative to dat file
            objectFilename = fileName.substr(0, p + 1) + objectFilename;
        }

        int start = 0;
        int end = numFrames;
        if (volumeId != -1) {
            if (volumeId >= numFrames)
                throw tgt::FileException("Specified volume id not in volume", fileName);

            start = volumeId;
            end = volumeId+1;
        }


        VolumeList* toReturn = new VolumeList();
        for (int frame = start; frame < end; ++frame) {
            VolumeRAM* vol = 0;

            if (getProgressBar()) {
                getProgressBar()->setTitle("Loading Volume");
                getProgressBar()->setProgressMessage("Loading volume: " + fileName);
            }

            ivec3 dims = ivec3(resolution, static_cast<int>(sliceCols[frame].getNumSlices()));
            if (format == "CHAR") {
                vol = new VolumeAtomic<int8_t>(dims);
            } else if (format == "UCHAR") {
                vol = new VolumeAtomic<uint8_t>(dims);
            } else if (format == "SHORT") {
                vol = new VolumeAtomic<int16_t>(dims);
            } else if (format == "USHORT") {
                vol = new VolumeAtomic<uint16_t>(dims);
            } else if (format == "UINT") {
                vol = new VolumeAtomic<uint32_t>(dims);
            } else if (format == "INT") {
                vol = new VolumeAtomic<int32_t>(dims);
            } else if (format == "FLOAT") {
                vol = new VolumeAtomic<float>(dims);
            } else {
                LERROR("Unknown Format");
            }

            std::fstream fin(objectFilename.c_str(), std::ios::in | std::ios::binary);
            if (!fin.good())
                throw tgt::IOException();

            size_t sliceSize = vol->getBytesPerVoxel() * resolution.x * resolution.y;
            for(size_t i=0; i<sliceCols[frame].getNumSlices(); i++) {
                fin.seekg(sliceSize * sliceCols[frame].getSlice(i).id_);
                fin.read(reinterpret_cast<char*>((size_t)vol->getData() + (i*sliceSize)) , sliceSize);
                if (getProgressBar())
                    getProgressBar()->setProgress(static_cast<float>(i) / static_cast<float>(sliceCols[frame].getNumSlices()));
            }
            if (getProgressBar())
                getProgressBar()->hide();

            if ( fin.eof() )
                throw tgt::CorruptedFileException();

            fin.close();

            if (vol) {
                VolumeURL origin(fileName);
                origin.addSearchParameter("volumeId", itos(frame));
                Volume* vh = new Volume(vol, sliceThickness, vec3(0.0f));

                //if (recoTranspose == 1) {
                    //Volume* vhTransposed = VolumeOperatorTranspose::APPLY_OP(vh, 0, 1);
                    //delete vh;
                    //vh = vhTransposed;
                //}
                //else if (recoTranspose == 2) {
                    //Volume* vhTransposed = VolumeOperatorTranspose::APPLY_OP(vh, 1, 2);
                    //delete vh;
                    //vh = vhTransposed;
                //}
                //else if (recoTranspose == 3) {
                    //Volume* vhTransposed = VolumeOperatorTranspose::APPLY_OP(vh, 0, 2);
                    //delete vh;
                    //vh = vhTransposed;
                //}

                // add metadata
                vh->setOrigin(origin);
                if (static_cast<size_t>(frame) < PVM_DwGradVec.size())
                    vh->setMetaDataValue<DVec3MetaData>("DiffusionGradientOrientation", PVM_DwGradVec[frame]);
                if (static_cast<size_t>(frame) < PVM_DwEffBval.size())
                    vh->setMetaDataValue<DoubleMetaData>("DiffusionBValue", PVM_DwEffBval[frame]);

                if ( (static_cast<size_t>(frame) < mappingOffset.size()) &&
                        (static_cast<size_t>(frame) < mappingScale.size()) ) {
                    float offset = mappingOffset[frame];
                    float scale = mappingScale[frame];

                    RealWorldMapping denormalize;
                    if (format == "CHAR") {
                        denormalize = RealWorldMapping::createDenormalizingMapping<int8_t>();
                    } else if (format == "UCHAR") {
                        denormalize = RealWorldMapping::createDenormalizingMapping<uint8_t>();
                    } else if (format == "SHORT") {
                        denormalize = RealWorldMapping::createDenormalizingMapping<int16_t>();
                    } else if (format == "USHORT") {
                        denormalize = RealWorldMapping::createDenormalizingMapping<uint16_t>();
                    } else if (format == "UINT") {
                        denormalize = RealWorldMapping::createDenormalizingMapping<uint32_t>();
                    } else if (format == "INT") {
                        denormalize = RealWorldMapping::createDenormalizingMapping<int32_t>();
                    } else if (format == "FLOAT") {
                        denormalize = RealWorldMapping::createDenormalizingMapping<float>();
                    } else{
                        LERROR("Unknown Format");
                    }

                    RealWorldMapping rwm(scale, offset, "");

                    vh->setRealWorldMapping(RealWorldMapping::combine(denormalize, rwm));
                    mapMetaData("VisuFGElemComment", "Description", vh, visuParsParser, frame, numFrames);
                }

                toReturn->add(vh);
            }
        }
        return toReturn;
    }
    else {
        throw tgt::CorruptedFileException("Failed to load file", fileName);
    }
}

VolumeReader* BrukerVolumeReader::create(ProgressBar* progress) const {
    return new BrukerVolumeReader(progress);
}

bool BrukerVolumeReader::mapTokenToMetaData(const std::string& voreenKey, Volume* vh, const JCampToken& t) {
    switch(t.getType()) {
        case JCampToken::STRING:
        case JCampToken::ENUM:
            vh->setMetaDataValue<StringMetaData>(voreenKey, t.getString());
            return true;
            break;
        case JCampToken::NUMBER:
            if(t.isInteger())
                vh->setMetaDataValue<IntMetaData>(voreenKey, t.getInteger());
            else
                vh->setMetaDataValue<FloatMetaData>(voreenKey, t.getFloat());

            return true;
            break;
        case JCampToken::OPENING_BRACKET:
        case JCampToken::CLOSING_BRACKET:
        case JCampToken::COMMA:
            LWARNING("Cannot map element!");
            return false;
    }
    return false;
}

bool BrukerVolumeReader::mapMetaData(const std::string& key, const std::string& voreenKey, Volume* vh, const JCampParser& parser, size_t id, size_t numFrames) {
    if(parser.containsElement(key)) {
        JCampElement e = parser.getElement(key);

        if(e.getNumDimensions() == 0) {
            JCampToken t = e.getData(0);

            return mapTokenToMetaData(voreenKey, vh, t);
        }
        else if(e.getNumDimensions() == 1) {
            if(e.getSize(0) == numFrames) {
                JCampToken t = e.getData(id);

                return mapTokenToMetaData(voreenKey, vh, t);
            }
            else
                return false;
        }
        else if(e.getNumDimensions() == 2) {
            if(e.getSize(0) == numFrames) {
                //TODO: map vec3...
                return false;
            }
            else
                return false;
        }
        else
            return false;
    }
    else {
        //LWARNING("key " << key << " not found!");
        return false;
    }
}

std::vector<VolumeURL> BrukerVolumeReader::listVolumes(const std::string& url) const {
    VolumeURL urlOrigin(url);
    std::string fileName = urlOrigin.getPath();

    std::vector<VolumeURL> result;

    //Build slice collections:
    vector<SliceCollection> sliceCols = listSliceCollections(url);

    for(size_t i=0; i<sliceCols.size(); i++) {
        VolumeURL origin("bruker", fileName);
        origin.addSearchParameter("volumeId", itos(i));
        origin.getMetaDataContainer().addMetaData("volumeId", new IntMetaData(static_cast<int>(i)));
        result.push_back(origin);
    }

    return result;
}

std::vector<BrukerVolumeReader::SliceCollection> BrukerVolumeReader::listSliceCollections(const std::string& url) const {
    VolumeURL urlOrigin(url);
    std::string fileName = urlOrigin.getPath();

    std::string path = tgt::FileSystem::dirName(fileName);
    std::string visuParsFile = path + "/visu_pars";

    JCampParser visuParsParser;
    if(visuParsParser.parseFile(visuParsFile)) {
        int numFrames = 0;
        try {
            numFrames = visuParsParser.getInteger("VisuCoreFrameCount");
            LINFO("Got " << numFrames << " frames.");
        }
        catch(JCampParserException& e) {
            LERROR(e.what());
            return vector<BrukerVolumeReader::SliceCollection>();
        }

        vector<FrameGroup> frameGroups;
        try {
            JCampElement e = visuParsParser.getElement("VisuFGOrderDesc");
            e.verifyNumDimensions(2);
            e.verifySize(1, 5);

            LINFO("Got " << e.getSize(0) << " grouping elements");
            for(size_t i=0; i<e.getSize(0); i++) {
                FrameGroup fg;
                fg.len_ = e.getData(i, 0).getInteger();
                fg.groupId_ = e.getData(i, 1).getString();
                fg.groupComment_ = e.getData(i, 2).getString();
                fg.valsStart_ = e.getData(i, 3).getInteger();
                fg.valsCnt_ = e.getData(i, 4).getInteger();

                frameGroups.push_back(fg);

                LINFO("Framegroup: " << fg.len_ << " " << fg.groupId_ << "(" << fg.groupComment_ << ") " << fg.valsStart_ << " " << fg.valsCnt_);
            }

        }
        catch(JCampParserException& e) {
            LERROR(e.what());
            return vector<BrukerVolumeReader::SliceCollection>();
        }

        int numVolumes = 1;
        int numSlices = 0;
        size_t sliceGroup = -1;
        for(size_t i=0; i<frameGroups.size(); i++) {
            if(frameGroups[i].groupId_ == "FG_SLICE") {
                numSlices = frameGroups[i].len_;
                sliceGroup = i;
            }
            else {
                numVolumes *= frameGroups[i].len_;
            }
        }
        LINFO("Found " << numVolumes << " volumes with " << numSlices << " slices.");

        //calculate group ids for all slices:
        vector<Slice> slices;
        for(int i=0; i<numFrames; i++) {
            size_t temp = i;
            slices.push_back(Slice(i));
            for(size_t j=0; j<frameGroups.size(); j++) {
                size_t c = temp % frameGroups[j].len_;

                slices.back().groupIds_.push_back(static_cast<int>(c));

                temp = temp / frameGroups[j].len_;
            }
        }

        //Build slice collections:
        vector<SliceCollection> sliceCols;
        for(int i=0; i<numVolumes; i++) {
            sliceCols.push_back(SliceCollection(numSlices));
        }

        for(int i=0; i<numFrames; i++) {
            for(size_t j=0; j<sliceCols.size(); j++) {
                if(sliceCols[j].hasMatchingGroupIds(slices[i].groupIds_, sliceGroup)) {
                    sliceCols[j].setSlice(slices[i].groupIds_[sliceGroup], slices[i]);
                    break;
                }
            }
        }

        std::cout << "\n";
        for(size_t i=0; i<sliceCols.size(); i++) {
            for(size_t j=0; j<sliceCols[i].getNumSlices(); j++) {
                std::cout << sliceCols[i].getSlice(j).id_ << " ";
            }
            std::cout << "\n";
        }
        return sliceCols;
    }
    return vector<BrukerVolumeReader::SliceCollection>();
}

} // namespace voreen
