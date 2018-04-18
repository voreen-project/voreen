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

#ifndef VRN_PLOTDATASELECT_H
#define VRN_PLOTDATASELECT_H

#include "../ports/plotport.h"
#include "../utils/functionlibrary.h"
#include "../properties/plotpredicateproperty.h"

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/vectorproperty.h"

namespace voreen {

class Serializable;
class PlotBase;
class PlotData;
class PlotPredicate;

/**
  * This class describes the Processor who implements the Select-Function
  * on PlotData-Objectes similar to the Select-Function of SQL.
  */
class VRN_CORE_API PlotDataSelect : public Processor {

public:
    PlotDataSelect();
    ~PlotDataSelect();
    virtual Processor* create() const;

    virtual std::string getCategory() const { return "Plotting"; }
    virtual std::string getClassName() const { return "PlotDataSelect"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_TESTING; }
    virtual bool isEndProcessor() const;
    virtual bool isReady() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);
    /// Provides the link to the object for the widget.
    const PlotData* getPlotData() const;
    /// deletes the specified column from the resulting PlotData.
    void disSelectColumn(int column);
    /// deletes the specified columns from the resulting PlotData.
    void disSelectColumns(const std::vector<int>& columns);
    /// specifies the columns that are remained included their sequence.
    void selectColumns(const std::vector<int>& columns);
    /// expected two points, the area between the points is taken ,the first point is positioned on the left above, and the second one at right bottom of the rectangle.
    void selectIncludeArea(const std::pair<int,int>& point1, const std::pair<int,int>& point2);
    /// expected a vector of columns and a vector of rows. the cross-selection of this vectors is in the Resultplotdata
    void selectColumnsRows(const std::vector<int>& columns, const std::vector<int>& rows);
    /// delete resetlist.
    void resetAllData();
    /// delete tha last dataset.
    void resetLastData();
    /// changes the Functionality between Select/Column Order/Filter
    void selectFunctionalityType(FunctionLibrary::ProcessorFunctionalityType type);
    /// move the old column before specified column
    void slideColumnBefore(int oldcolumn, int beforecolumn);
    /// move the old column after specified column
    void slideColumnAfter(int oldcolumn, int aftercolumn);
    /// exchange the old column with specified column
    void slideColumnSwitch(int firstcolumn, int secondcolumn);
    /// cells with empty values will be removed.
    void filterEmptyCells();
    /// cells with non alphanumeric values will be removed.
    void filterAlphaNumeric();
    /// execute the predicates of PlotPredicateProperties for objects of PloData .
    void selectPredicate();
    /// copy the column and set it at the expected position.
    void copyColumn(int copycolumn, int targetcolumn);
    /// specified colunmn can be renamed
    void renameColumn(int column, std::string caption);
    ///gives back the size of the resetlist
    size_t getResetListSize() const;

protected:
    virtual void setDescriptions() {
        setDescription("Select, Sort or Filter the Data on the PlotData.inPort.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

private:

    typedef std::vector< std::pair< int, PlotPredicate*> > PredicateVector;
    typedef std::vector< std::pair< int, std::string> > LabelVector;

    struct SelectValues : public Serializable {
    public:
        std::vector<int> projectionVector;
        std::vector<int> selectedRows;
        PredicateVector predicateVector;
        LabelVector labelVector;
        int keyColumnCount;
        int dataColumnCount;

        /// @see Serializer::serialize
        virtual void serialize(Serializer& s) const;

        /// @see Deserializer::deserialize
        virtual void deserialize(Deserializer& s);

    };

    PlotData* selectData(PlotData* plotData);
    void filterData(PlotData* plotData);
    void setOutPortData();
    void readData();
    void checkBoxChange();
    void calculate();
    void deleteList();
    void keyColumnChange();
    void updateView();

    void fillResetList(std::vector<int> columns, const std::vector<int>& rows,
        const PredicateVector& predVector, const LabelVector& columntexts,
        int keyColumnCount, int dataColumnCount);

    PlotPort inPort_;
    PlotPort outPort_;

    const PlotData* pData_;
    PlotData* pDataOut_;

    std::vector<SelectValues> resetList_;

    /// properties
    BoolProperty select_;
    IntProperty newKeyColumnCount_;
    PlotPredicateProperty plotPredicateVector_;

    BoolProperty columnorder_;

    BoolProperty filter_;

    IntProperty maxTableColumnWidth_;

    ColorProperty keyColumnColor_;
    ColorProperty dataColumnColor_;

    static const std::string loggerCat_;
};


} //namespace

#endif // VRN_PLOTDATASELECT_H
