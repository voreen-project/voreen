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

#include "plotdatasource.h"
#include "voreen/core/voreenapplication.h"
#include "../datastructures/plotdata.h"
#include "../datastructures/plotcell.h"

#include <limits>
#include "tgt/tgt_math.h"
#include <iostream>
#include <fstream>
#include <string>

namespace voreen {

const std::string PlotDataSource::loggerCat_("voreen.PlotDataSource");

PlotDataSource::PlotDataSource():
    Processor()
    , outPort_(Port::OUTPORT,"PlotDataOutPort")
    , inputFile_("inputFile", "CSV File", "Open CSV input file", VoreenApplication::app()->getUserDataPath(),
        "Commaseparated Files(*.csv);;Textfile (*.txt);;All TextFiles (*.csv *.txt)",
        FileDialogProperty::OPEN_FILE,Processor::VALID)
    , separator_("Separator","CSV Separator",";",Processor::VALID)
    , pData_(0)
    , countLine_("countLine","Set Header Lines",1,0,200,Processor::VALID)
    , countKeyColumn_("countKeyColumn","Set Key Columns",1,1,200,Processor::VALID)
    , constantOrder_("constantOrder","Constant Data Order",Processor::VALID)
    , recalculate_("recalculate","Reload CSV File")
{
    inputFile_.onChange(MemberFunctionCallback<PlotDataSource>(this, &PlotDataSource::recalculate));
    recalculate_.onChange(MemberFunctionCallback<PlotDataSource>(this, &PlotDataSource::recalculate));
    addProperty(&separator_);
    addProperty(&countLine_);
    addProperty(&countKeyColumn_);
    addProperty(&constantOrder_);
    addProperty(&inputFile_);
    addProperty(&recalculate_);
    addPort(&outPort_);

    inputFile_.setGroupID("FileSelection");
    recalculate_.setGroupID("FileSelection");

    separator_.setGroupID("Fileparams");
    countLine_.setGroupID("Fileparams");
    countKeyColumn_.setGroupID("Fileparams");
    constantOrder_.setGroupID("Fileparams");

    setPropertyGroupGuiName("FileSelection", "File Selection");
    setPropertyGroupGuiName("FileParams", "File Parameter");
}

PlotDataSource::~PlotDataSource() {
    // nothing to do here
}

Processor* PlotDataSource::create() const {
    return new PlotDataSource();
}

void PlotDataSource::process() {
    tgtAssert(pData_, "no plotdata object");

    if (pData_->getColumnCount() == 0)
        recalculate();
    Processor::setProgress(1.f);
}

void PlotDataSource::initialize() {
    Processor::initialize();

    pData_ = new PlotData(0,0);
    Processor::setProgress(1.f);
}

void PlotDataSource::deinitialize() {
    if (pData_ != outPort_.getData())
        delete pData_;
    outPort_.setData(0, true);
    pData_ = 0;

    Processor::deinitialize();
}

void PlotDataSource::recalculate() {
    Processor::setProgress(0.f);
    PlotData* newData = readCSVData();
    PlotData* oldData = pData_;
    pData_ = newData;
    Processor::setProgress(1.f);
    outPort_.setData(pData_);
    delete oldData;
}

PlotData* PlotDataSource::readCSVData() {
    PlotData* newData = new PlotData(0,0);
    if (!isInitialized())
        return newData;

    std::string filename = inputFile_.get();
    if (filename.empty()) {
        LERROR("        Filename ist empty.");
        return newData;
    }

    std::ifstream inFile;
    LINFO("        Open file: " <<  filename);
    inFile.open(filename.c_str(), std::ios::in);
    if (inFile.fail()) {
        LERROR("        Unable to open data file: " << filename);
        return newData;
    }
    // ProgressBar Position
    inFile.seekg(0, std::ios::end);
    size_t size = static_cast<size_t>(inFile.tellg());
    inFile.seekg(0);
    float position = 0;
    //char* buffer = new char [size];
    //inFile.read(buffer,size);
    //std::stringstream selfFile(buffer);
    std::string line;
    std::vector<std::string> myvector (0);
    size_t mass = 0;
    inFile.seekg(0);
    int k = 0;
    for(int j=0; j< countLine_.get(); j++){
        std::getline(inFile, line);
        position += (line.size()*1.f)/(size*1.f);
        if (j == countLine_.get() - 1){
            csvline_populate(myvector, line);
            if (!constantOrder_.get()){
                newData->reset(countKeyColumn_.get(),
                    static_cast<int>(myvector.size())-countKeyColumn_.get());
                k = 0;
            }
            else {
                newData->reset(1, static_cast<int>(myvector.size()));
                newData->setColumnLabel(0,"Index");
                k = 1;
            }
            for (size_t i = 0; i < myvector.size(); ++i) {
                newData->setColumnLabel(static_cast<int>(i)+k, myvector[i]);
            }
            mass = myvector.size();
        }
        myvector.clear();
    }

    std::vector<int> tester(0);
    std::vector<PlotCellValue> pCellVector_;
    std::vector<int>::iterator rit;
    double num;
    int counter = 0;
    std::stringstream ssLine;
    ssLine.imbue(std::locale::classic());
    while(inFile.good()&& std::getline(inFile, line)){
        csvline_populate(myvector, line);
        position += line.size()*1.f/(size*1.f);
        if (mass == 0) {
            mass = myvector.size();
            if (!constantOrder_.get()){
                newData->reset(countKeyColumn_.get(),
                    static_cast<int>(myvector.size()) - countKeyColumn_.get());
                k = 0;
            }
            else {
                newData->reset(1, static_cast<int>(myvector.size()));
                newData->setColumnLabel(0,"Index");
                k = 1;
            }
            for (int i = k; i < newData->getColumnCount(); ++i) {
                std::stringstream Str;
                Str << i;
                newData->setColumnLabel(i,Str.str());
            }
        }
        if (constantOrder_.get()){
            PlotCellValue cellValue(counter);
            pCellVector_.push_back(cellValue);
        }
        size_t found = 0;
        for(size_t i=0; i<mass; ++i){
            if (i < myvector.size()){
                found = myvector[i].find(',');
                ssLine.str("");
                ssLine.clear();
                if (found != std::string::npos) {
                    ssLine <<  myvector[i].replace(found,1,".");
                }
                else {
                    ssLine << myvector[i];
                }
                //the cast should be work for numbers x.y
                if (ssLine >> num) {
                    if(counter==0){
                        tester.push_back(2);
                    }
                    else {
                        if ((tester.at(i) == 1)){
                            pCellVector_.push_back(PlotCellValue(myvector[i]));
                            continue;
                        }
                        else if (tester.at(i) == 0) {
                            rit = tester.begin() + i;
                            *rit = 2;
                        }
                    }
                    pCellVector_.push_back(PlotCellValue(num));
                }
                else if (myvector[i].length() == 0) {
                    if (counter == 0) {
                        tester.push_back(0);
                    }
                    pCellVector_.push_back(PlotCellValue());
                }
                else {
                    if (counter==0){
                        tester.push_back(1);
                    }
                    else {
                        if ((tester.at(i) == 2)){
                            pCellVector_.push_back(PlotCellValue());
                            continue;
                        }
                        else if (tester.at(i) == 0) {
                            rit = tester.begin() + i;
                            *rit = 1;
                        }
                    }
                    pCellVector_.push_back(PlotCellValue(myvector[i]));
                }
            }
            else {
                if (i < tester.size()){
                    pCellVector_.push_back(PlotCellValue());
                }
                else {
                    pCellVector_.push_back(PlotCellValue());
                    tester.push_back(0);
                }
            }
        }
        newData->insert(pCellVector_);
        pCellVector_.clear();
        myvector.clear();
        ++counter;
        Processor::setProgress(position);
    }
    inFile.close();
//    delete[] buffer;
    return newData;
}

namespace {
    const std::string whitespaces(" \t");
}

std::string PlotDataSource::trimString(const std::string& oldString) {
    size_t start = oldString.find_first_not_of(whitespaces);
    if (start == std::string::npos) // oldString contains only whitespaces
        return "";
    size_t end = oldString.find_last_not_of(whitespaces);
    return oldString.substr(start, (end - start + 1));
}

void PlotDataSource::csvline_populate(std::vector<std::string> &record, const std::string& line){
    size_t linepos = 0;
    size_t endpos = 0;
    // safety first
    if (separator_.get().empty())
        return;
    char delimiter = separator_.get()[0];

    // we are at the beginning of an entry, skip whitespaces and check if not already reached end of line
    while ((endpos != std::string::npos) && (linepos = line.find_first_not_of(whitespaces, linepos)) != std::string::npos) {
        // now linepos points to the first non blank character, here starts the entry
        // check whether there are quotes
        if (line[linepos] == '"') {
            // find positon of closing quotes
            endpos = line.find_first_of('"', linepos + 1);

            std::string toPush = line.substr(linepos + 1 , endpos - linepos - 1);
            // ensure we haven't found double quotes ("") which shall be resolved to one double quote in resulting string
            while ((endpos != std::string::npos) && (endpos + 1 < line.length()) && (line[endpos + 1] == '"')) {
                linepos = endpos + 1;
                endpos = line.find_first_of('"', endpos + 2);
                toPush.append(line.substr(linepos, endpos - linepos));
            }

            // push string in quotes onto record
            record.push_back(trimString(toPush));

            // ignore everything until next delimiter
            endpos = line.find_first_of(delimiter, endpos);
        }
        // ok, this entry is not in quotes - just push everything until next delimiter onto record
        else {
            endpos = line.find_first_of(delimiter, linepos);
            record.push_back(trimString(line.substr(linepos, endpos - linepos)));
        }
        linepos = endpos + 1;
    }
}

}
