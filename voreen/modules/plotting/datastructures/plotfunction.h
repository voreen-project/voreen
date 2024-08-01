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

#ifndef VRN_PLOTFUNCTION_H
#define VRN_PLOTFUNCTION_H

#include "plotexpression.h"
#include "interval.h"
#include "plotbase.h"

namespace voreen {

class PlotData;

/**
 * \brief This class stores an expression and can evaluate on a specified domain.
 *
 */
class VRN_CORE_API PlotFunction : public PlotBase {
public:
    /// type description how to handle with long Expressions
    enum ExpressionDescriptionLengthType {
        FULL = 0,
        MAXLENGTH = 1,
        CUSTOM = 2,
        ONLYNAME = 3
    };

    /// default constructor
    PlotFunction();

    /// initialises keyColumnCount, dataColumnCount by the dimensions of the expression
    PlotFunction(const std::string& expressionString, const std::string& expressionName = "f");

    /// initialises keyColumnCount, dataColumnCount by the dimensions of the expression
    PlotFunction(const PlotExpression& expr);

    /// copy constructor
    PlotFunction(const PlotFunction& rhs);

    /// default destructor
    virtual ~PlotFunction();

    PlotFunction& operator=(const PlotFunction& rhs);

    /// Evaluates the function
    plot_t evaluateAt(const std::vector<plot_t>& value) const;

    /// Evaluates the function in the interval
    std::vector<std::vector<plot_t> > evaluateAt(const std::vector<Interval<plot_t> >& interval,const std::vector<plot_t>& step) const;

    /// Evaluates the function in the specified interval with step size step  and returns a table of the type PlotData.
    /// The length of the arrays have to be same as the number of parameters in expr_.
    bool select(const std::vector<Interval<plot_t> >& interval, const std::vector<plot_t>& step, PlotData& target) const;

    /// select for functions with only one variable
    bool select(const Interval<plot_t>& interval, const plot_t& step, PlotData& target) const;

    /// select for functions with many variables
    bool select(const std::vector<std::pair<plot_t,plot_t> >& interval, const std::vector<plot_t>& step, PlotData& target) const;

    /// select for functions with many variables
    bool select(Interval<plot_t>* interval, plot_t* step, int count, PlotData& target) const;
    /// return the Plotexpression
    const PlotExpression& getPlotExpression() const;
    PlotExpression& getPlotExpression();
    /// setting the length of the output column string how represents the Expression
    void setExpressionLength(ExpressionDescriptionLengthType expressionType, int maxLength = 100, const std::string& customText = "function");

private:

    //! An expression with variables
    PlotExpression expr_;
    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_PLOTFUNCTION_H
