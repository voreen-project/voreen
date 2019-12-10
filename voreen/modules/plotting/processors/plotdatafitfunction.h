/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_PLOTDATAFITFUNCTION_H
#define VRN_PLOTDATAFITFUNCTION_H

#include "../ports/plotport.h"

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/colorproperty.h"

#include "../utils/functionlibrary.h"

namespace voreen {

class PlotBase;
class PlotFunction;
class PlotData;


class VRN_CORE_API PlotDataFitFunction : public Processor {
public:
    /**
     * type to transport the values who are used to fit a PlotData-Table
     **/
    struct FittingValues : public Serializable {
    public:
        int column;
        FunctionLibrary::RegressionType regressionType;
        plot_t mse;
        int dimension;

        /// @see Serializer::serialize
        virtual void serialize(Serializer& s) const;

        /// @see Deserializer::deserialize
        virtual void deserialize(Deserializer& s);
    };

    PlotDataFitFunction();
    virtual Processor* create() const;
    virtual std::string getCategory() const  { return "Plotting"; }
    virtual std::string getClassName() const { return "PlotDataFitFunction"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_TESTING; }
    //virtual bool usesExpensiveComputation() const { return true; }
    virtual bool isEndProcessor() const;
    virtual bool isReady() const;

    /// execute linear regression on specified data
    void linearRegression(int column);
    /// execute sqrt regression on specified data
    void sqrtRegression(int column);
    /// execute square regression on specified data
    void squareRegression(int column);
    /// execute cubic regression on specified data
    void cubicRegression(int column);
    /// execute multi regression on specified data
    void multiRegression(int dimension,int column);
    /// execute logarithmic regression on specified data
    void logarithmicRegression(int column);
    /// execute power regression on specified data
    void powerRegression(int column);
    /// execute exponential regression on specified data
    void exponentialRegression(int column);
    /// execute sinus regression on specified data
    void sinRegression(int column);
    /// execute cosine regression on specified data
    void cosRegression(int column);
    /// execute interpolation regression on specified data
    void interpolRegression(int column);

    void constantSpline(int column);
    /// execute simple apline regression on specified data
    void simpleSplineRegression(int column);
    /// execute b spline regression on specified data
    void squareSplineRegression(int column);
    /// execute b spline regression on specified data
    void bSplineRegression(int column);
    /// execute cubic spline regression on specified data
    void cubicSplineRegression(int column);

    void ignoreFalseValues();

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

    /// returns the PlotDataObject which is on the inport
    const PlotData* getPlotData() const;
    /// returns the PlotFunction that is set on the outport
    const PlotFunction* getPlotFunction() const;
    /// return the aktuell fittingvalues
    const FittingValues getFittingValues() const;


    /// default destructor
    virtual ~PlotDataFitFunction();
protected:
    virtual void setDescriptions() {
        setDescription("interpolate the Data on the PlotData.inPort.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

private:


    void setOutPortData();
    void readData();
    void calculate();
    void recalculate();
    void deleteList();
    void updateView();
    void changeText();

    void fitData();

    PlotPort outPort_;
    PlotPort inPort_;

    const PlotData* pData_;
    PlotFunction* pDataOut_;

    FittingValues fittingValues_;

    BoolProperty ignoreFalseValues_;

    StringProperty expressionNameInput_;
    IntOptionProperty expressionText_;

    StringProperty selfDescription_;
    IntProperty maxLength_;

    ButtonProperty recalculate_;

    IntProperty maxTableColumnWidth_;


    ColorProperty keyColumnColor_;
    ColorProperty dataColumnColor_;
    ColorProperty fittingColumnColor_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_PLOTDATAFITFUNCTION_H
