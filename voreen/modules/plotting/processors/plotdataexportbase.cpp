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

#include "plotdataexportbase.h"

#include "../datastructures/plotdata.h"
#include "../datastructures/plotcell.h"
#include "../datastructures/plotrow.h"
#include "../datastructures/plotfunction.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/tgt_math.h"
#include <iostream>
#include <fstream>

namespace voreen {

const std::string PlotDataExportBase::loggerCat_("voreen.PlotDataExportBase");

PlotDataExportBase::PlotDataExportBase():
    Processor()
    , inPort_(Port::INPORT,"PlotData.InPort")
    , separator_("Separator","CSV Separator",";",Processor::VALID)
    , valueBracket_("use inverted comma","Use Inverted Comma (txt and csv)",false,Processor::VALID)
    , pData_(0)
    , pPlotFunction_(0)
{
    addProperty(&separator_);
    addProperty(&valueBracket_);
    addPort(&inPort_);

    separator_.setGroupID("Fileparams");
    valueBracket_.setGroupID("Fileparams");

    setPropertyGroupGuiName("FileParams", "File-Parameter");
}

PlotDataExportBase::~PlotDataExportBase() {
    // nothing to do here
}

Processor* PlotDataExportBase::create() const {
    return new PlotDataExportBase();
}

bool PlotDataExportBase::isEndProcessor() const {
    return true;
}

void PlotDataExportBase::process() {
    pData_ = 0;
    pPlotFunction_ = 0;
    if (inPort_.hasData()) {
        if (dynamic_cast<const PlotData*>(inPort_.getData())) {
            pData_ = dynamic_cast<const PlotData*>(inPort_.getData());
            pPlotFunction_ = 0;
        }
        else if (dynamic_cast<const PlotFunction*>(inPort_.getData())) {
            pPlotFunction_ = dynamic_cast<const PlotFunction*>(inPort_.getData());
            pData_ = 0;
        }
        else {
            LWARNING("PlotDataExportBase can only handle with PlotData- and PlotFunction-Objects");
        }
    }
}

void PlotDataExportBase::deinitialize() {
    pData_ = 0;
    pPlotFunction_ = 0;
    Processor::deinitialize();
}

void PlotDataExportBase::exportFile(std::string filename) {
    if (pData_ || pPlotFunction_) {
        size_t pos = filename.find_last_of(".");
        std::string extension = (pos == std::string::npos) ? "" : filename.substr(pos);
        std::string separator = separator_.get();
        std::ofstream exportfile (filename.c_str());
        bool invertedComma = valueBracket_.get();
        if (extension == ".csv" || extension == ".txt") {
            if (pData_) {
                for (int i = 0; i < pData_->getColumnCount(); ++i) {
                    if (invertedComma)
                        exportfile << "\"" << pData_->getColumnLabel(i) << "\"";
                    else exportfile << pData_->getColumnLabel(i);
                    if (i < pData_->getColumnCount()-1)
                        exportfile << separator;
                    else
                        exportfile << "\n";
                }
                std::vector<PlotCellValue> row;
                for (int i = 0; i < pData_->getRowsCount(); ++i) {
                    row = pData_->getRow(i).getCells();
                    for (size_t j = 0; j < row.size(); ++j) {
                        if (invertedComma)
                            exportfile << "\"";
                        if (row[j].isValue())
                            exportfile << row.at(j).getValue();
                        else
                            exportfile << row.at(j).getTag();
                        if (invertedComma)
                            exportfile << "\"";
                        if (j < row.size()-1)
                            exportfile << separator;
                        else if (i < pData_->getRowsCount()-1)
                            exportfile << "\n";
                    }
                }
            }
            else {
                exportfile << pPlotFunction_->getPlotExpression().getExpressionString();
            }
        }
        else if (extension == ".xml") {
            exportfile << "<?xml version=\"1.0\" ?>\n";
            if (pData_) {
                exportfile << "<plotdata>\n<plotrowvalue>\n";
                for (int i = 0; i < pData_->getColumnCount(); ++i) {
                    exportfile << "<plotcellvalue>" << pData_->getColumnLabel(i);
                    exportfile << "</plotcellvalue>";
                }
                exportfile << "\n</plotrowvalue>\n";
                std::vector<PlotCellValue> row;
                for (int i = 0; i < pData_->getRowsCount(); ++i) {
                    row = pData_->getRow(i).getCells();
                    exportfile << "<plotrowvalue>";
                    for (size_t j = 0; j < row.size(); ++j) {
                        exportfile << "<plotcellvalue>";
                        if (row[j].isValue())
                            exportfile << row.at(j).getValue();
                        else
                            exportfile << row.at(j).getTag();
                        exportfile << "</plotcellvalue>";
                    }
                    exportfile << "</plotrowvalue>\n";
                }
                exportfile << "</plotdata>";
                std::vector<PlotCellImplicit> rowIm;
                if (pData_->getImplicitRowsCount() > 0) {
                    exportfile << "<plotrowimplicit>";
                    for (int i = 0; i < pData_->getColumnCount(); ++i) {
                        exportfile << "<plotcellimplicit>&nbsp;</plotcellimplicit>";
                    }
                    exportfile << "</plotrowimplicit>\n";
                    for (int i = 0; i < pData_->getImplicitRowsCount(); ++i) {
                        rowIm = pData_->getImplicitRow(i).getCells();
                        exportfile << "<plotrowimplicit>\n";
                        for (size_t j = 0; j < 2*rowIm.size(); ++j) {
                            if (j == rowIm.size())
                                exportfile << "</plotrowimplicit>\n<plotrowimplicit>\n";
                            if (rowIm.at(j%rowIm.size()).getAggregationFunction()) {
                                if (j/rowIm.size() == 0)
                                    exportfile << "<plotcellimplicit>" << rowIm.at(j%rowIm.size()).getAggregationFunction()->toString() << "</plotcellimplicit>";
                                else
                                    exportfile << "<plotcellimplicit>" << rowIm.at(j%rowIm.size()).getValue() << "</plotcellimplicit>";
                            }
                            else
                                exportfile << "<plotcellimplicit></plotcellimplicit>";

                        }
                        exportfile << "\n</plotrowimplicit>\n";
                    }
                }
            }
            else {
                exportfile << "<plotfunction>\n<expression>\n" << pPlotFunction_->getPlotExpression().getExpressionString()
                    << "\n</expression>\n<numberofvariables>\n" << pPlotFunction_->getPlotExpression().numberOfVariables()
                    << "\n</numberofvariables>\n</plotfunction>";
            }
        }
        else if (extension == ".htm" || extension == ".html") {
            exportfile << "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org"
                << "/TR/html4/strict.dtd\"><html>\n<head><title>";
            if (pData_)
                exportfile << "PlotData";
            else
                exportfile << "PlotFunction";
            exportfile << "</title></head><body>\n<table border=\"2px\">";
            exportfile << "<thead><tr>";
            if (pData_) {
                for (int i = 0; i < pData_->getColumnCount(); ++i) {
                    exportfile << "<td>" << pData_->getColumnLabel(i);
                    exportfile << "</td>";
                }
                exportfile << "</tr></thead>\n";
                std::vector<PlotCellValue> row;
                for (int i = 0; i < pData_->getRowsCount(); ++i) {
                    row = pData_->getRow(i).getCells();
                    exportfile << "<tr>";
                    for (size_t j = 0; j < row.size(); ++j) {
                        exportfile << "<td>";
                        if (row[j].isValue())
                            exportfile << row.at(j).getValue();
                        else
                            exportfile << row.at(j).getTag();
                        exportfile << "</td>";
                    }
                    exportfile << "</tr>\n";
                }
                std::vector<PlotCellImplicit> rowIm;
                if (pData_->getImplicitRowsCount() > 0) {
                    exportfile << "<tr>";
                    for (int i = 0; i < pData_->getColumnCount(); ++i) {
                        exportfile << "<td>&nbsp;</td>";
                    }
                    exportfile << "</tr>\n";
                    for (int i = 0; i < pData_->getImplicitRowsCount(); ++i) {
                        rowIm = pData_->getImplicitRow(i).getCells();
                        exportfile << "<tr>";
                        for (size_t j = 0; j < 2*rowIm.size(); ++j) {
                            if (j == rowIm.size())
                                exportfile << "</tr><tr>";
                            if (rowIm.at(j%rowIm.size()).getAggregationFunction()) {
                                if (j/rowIm.size() == 0)
                                    exportfile << "<td>" << rowIm.at(j%rowIm.size()).getAggregationFunction()->toString() << "</td>";
                                else
                                    exportfile << "<td>" << rowIm.at(j%rowIm.size()).getValue() << "</td>";
                            }
                            else
                                exportfile << "<td></td>";

                        }
                        exportfile << "</tr>\n";
                    }
                }
            }
            else {
                exportfile << "<td>Expression</td><td>Number of Variables</td></tr></head>\n<tr>\n<td>"
                    << pPlotFunction_->getPlotExpression().getExpressionString() << "</td><td>"
                    << pPlotFunction_->getPlotExpression().numberOfVariables() << "</td></tr>\n";
            }
            exportfile << "</table>\n</body>\n</html>";
        }
        exportfile.close();
    }
}


}
