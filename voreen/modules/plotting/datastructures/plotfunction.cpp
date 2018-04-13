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

#include "plotfunction.h"
#include "../utils/functionlibrary.h"
#include "plotdata.h"

#include "tgt/logmanager.h"

namespace voreen {

const std::string PlotFunction::loggerCat_("voreen.plotting.PlotFunction");

PlotFunction::PlotFunction()
    : PlotBase(0,0)
    , expr_(PlotExpression())
{
}

PlotFunction::PlotFunction(const std::string& expressionString, const std::string& expressionName)
    : PlotBase(0,0)
    , expr_(PlotExpression(expressionString,expressionName))
{
    PlotBase::reset(std::max<int>(expr_.numberOfVariables(),1),1);
    setColumnLabel(0,"x-axis");
    for (int i = 0; i < expr_.numberOfVariables(); ++i) {
        setColumnLabel(i,expr_.getVariable(i));
    }
    setColumnLabel(std::max<int>(expr_.numberOfVariables(),1),expr_.toString());
}


PlotFunction::PlotFunction(const PlotExpression& expr)
    :PlotBase(std::max<int>(expr.numberOfVariables(),1),1)
    , expr_(expr)
{
    setColumnLabel(0,"x-axis");
    for (int i = 0; i < expr_.numberOfVariables(); ++i) {
        setColumnLabel(i,expr_.getVariable(i));
    }
    setColumnLabel(std::max<int>(expr_.numberOfVariables(),1),expr_.toString());
}

PlotFunction::PlotFunction(const PlotFunction& rhs)
    : PlotBase(rhs)
    , expr_(rhs.expr_)
{
}

PlotFunction::~PlotFunction() {
}

PlotFunction& PlotFunction::operator=(const PlotFunction& rhs) {
    PlotBase::operator=(rhs);
    PlotExpression cpy(rhs.expr_);
    std::swap(expr_, cpy);
    return *this;
}

plot_t PlotFunction::evaluateAt(const std::vector<plot_t>& value) const {
    return expr_.evaluateAt(value);
}

std::vector<std::vector<plot_t> > PlotFunction::evaluateAt(const std::vector<Interval<plot_t> >& interval,const std::vector<plot_t>& step) const {
    std::vector<std::vector<plot_t> > result;
    result.clear();
    if (interval.size() != step.size() || (expr_.numberOfVariables() > static_cast<int>(step.size())))
        return result;
    std::vector<plot_t> x;
    std::vector<plot_t> start;
    for (size_t z = 0; z < interval.size(); ++z) {
        if (interval.at(z).getLeftOpen()) {
            start.push_back(interval.at(z).getLeft() + step.at(z));
        }
        else {
            start.push_back(interval.at(z).getLeft());
        }
    }
    x = start;
    plot_t y;
    std::vector<plot_t> partVector;
    size_t j = interval.size()-1;
    while (interval.at(0).contains(x.at(0))) {
        partVector.clear();
        y = evaluateAt(x);
        for (size_t i = 0; i < step.size(); ++i) {
            partVector.push_back(x[i]);
        }
        partVector.push_back(y);
        result.push_back(partVector);
        x[j] += step.at(j);
        for (size_t k = 1; k < interval.size(); ++k) {
            if (!interval.at(k).contains(x[k])) {
                x[k] = start[k];
                x[k-1] += step[k-1];
            }
        }
    }
    return result;
}


bool PlotFunction::select(const std::vector<Interval<plot_t> >& interval, const std::vector<plot_t>& step, PlotData& target) const {
    if (interval.size() != step.size() || (expr_.numberOfVariables() > static_cast<int>(step.size())) ||
        keyColumnCount_ + dataColumnCount_ == 0)
        return false;
    std::vector<std::vector<plot_t> > result = evaluateAt(interval,step);
    target.reset(keyColumnCount_,1);
    for (size_t i = 0; i < result.size(); ++i) {
        target.insert(result.at(i));
    }
    for (int i = 0; i < target.getColumnCount(); ++i) {
        target.setColumnLabel(i,getColumnLabel(i));
    }
    target.sorted_ = true;
    return true;
}

bool PlotFunction::select(const Interval<plot_t>& interval, const plot_t& step, PlotData& target) const {
    std::vector<Interval<plot_t> > intervalVector;
    intervalVector.push_back(interval);
    std::vector<plot_t> stepVector;
    stepVector.push_back(step);
    return select(intervalVector, stepVector,target);
}

bool PlotFunction::select(const std::vector<std::pair<plot_t,plot_t> >& interval, const std::vector<plot_t>& step, PlotData& target) const {
    std::vector<Interval<plot_t> > intervals;
    Interval<plot_t> intV;
    for (size_t i = 0; i < interval.size(); ++i) {
        if (interval.at(i).first <= interval.at(i).second) {
            intV = Interval<plot_t>(interval.at(i).first,interval.at(i).second);
        }
        else {
            intV = Interval<plot_t>(interval.at(i).second,interval.at(i).first);
        }
        intervals.push_back(intV);
    }
    return select(intervals,step,target);
}

bool PlotFunction::select(Interval<plot_t>* interval, plot_t* step, int count, PlotData& target) const {
    std::vector<Interval<plot_t> > intervals;
    std::vector<plot_t> steps;
    for (int i = 0; i < count; ++i) {
        intervals.push_back(interval[i]);
        steps.push_back(step[i]);
    }
    return select(intervals,steps,target);
}

const PlotExpression& PlotFunction::getPlotExpression() const {
    return expr_;
}

PlotExpression& PlotFunction::getPlotExpression() {
    return expr_;
}

void PlotFunction::setExpressionLength(ExpressionDescriptionLengthType expressionType, int maxLength, const std::string& customText) {
    std::string text;
    if (getKeyColumnCount() <= 0 || getDataColumnCount() <= 0)
        return;
    for (size_t k = 0; k < customText.size(); ++k) {
        if (customText[k] > 0)
            text += customText[k];
    }
    if (expressionType == CUSTOM && text.size() > 0) {
        setColumnLabel(std::max<int>(1,expr_.numberOfVariables()),text);
    }
    else if (expressionType == MAXLENGTH && static_cast<int>(expr_.toString().size()) > maxLength) {
        setColumnLabel(std::max<int>(1,expr_.numberOfVariables()),expr_.toString().substr(0,maxLength)+"...");
    }
    else if (expressionType == ONLYNAME) {
        setColumnLabel(std::max<int>(1,expr_.numberOfVariables()),
            expr_.toString().substr(0,expr_.toString().size()-expr_.getExpressionString().size()-1));
    }
    else
        setColumnLabel(std::max<int>(1,expr_.numberOfVariables()),expr_.toString());
}

}

