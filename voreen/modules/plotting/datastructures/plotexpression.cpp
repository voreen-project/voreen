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

//activate this option to get the debug messages of the parser
//#define VRN_PLOTEXPRESSION_DEBUG

#include "plotexpression.h"
#include "../utils/parser/plotfunctionlexer.h"
#include "../utils/parser/plotfunctionparser.h"
#include "../utils/parser/plotfunctionnode.h"
#include "../utils/parser/plotfunctionvisitor.h"

#include <limits>
#include "tgt/logmanager.h"
#include <ctype.h>
#include <math.h>


namespace voreen {


const std::string PlotExpression::loggerCat_("voreen.plotting.PlotExpression");

const int PlotExpression::charOffset_(97);

PlotExpression::PlotExpression() {
    initialize();
}

PlotExpression::PlotExpression(const std::string& value, const std::string& expressionName) {
    initialize();
    for (size_t k = 0; k < value.size(); ++k) {
        if (value[k] > 0)
            representation_ += tolower(value[k]);
        else
            throw VoreenException("Not Validate String");
    }
    representation_.erase(remove_if(representation_.begin(), representation_.end(), isspace), representation_.end());
    std::istringstream is(representation_);
    if (expressionName.size() > 0)
        expressionName_ = expressionName;
#ifdef VRN_PLOTEXPRESSION_DEBUG
    LINFO("Try to parse \"" << value << "\"");
#endif
    glslparser::PlotFunctionParser parser(&is);
#ifdef VRN_PLOTEXPRESSION_DEBUG
    parser.setDebugging(true);
#endif
    node_ =  dynamic_cast<glslparser::PlotFunctionNode*>(parser.parse());
#ifdef VRN_PLOTEXPRESSION_DEBUG
    std::string parserLog = parser.getLog().str();
    LINFO("FunctionParser Debug Messages:\n" << parserLog);
#endif
    if (node_ == 0)
        throw VoreenException("Not Validate String");
#ifdef VRN_PLOTEXPRESSION_DEBUG
    std::ostringstream elog_;
    LINFO("Expression accepted.");
#endif
    glslparser::PlotFunctionVisitor plotfunctionvisitor;
    std::vector<glslparser::PlotFunctionVisitor::TokenVector> functionTokenVector = plotfunctionvisitor.getPlotFunctionToken(node_);
    functionVector_.resize(functionTokenVector.size());
#ifdef VRN_PLOTEXPRESSION_DEBUG
    elog_ << "Number of Partialfunctions: " << functionTokenVector.size() << "\n";
#endif
    for (size_t j = 0; j < functionTokenVector.size(); ++j) {
#ifdef VRN_PLOTEXPRESSION_DEBUG
        elog_ << j+1 << ". PartialFunction\n";
#endif
        functionVector_.at(j).function = functionTokenVector.at(j).function;
        functionVector_.at(j).interval = functionTokenVector.at(j).interval;
        for (size_t i = 0; i < functionTokenVector.at(j).function.size(); ++i) {
            if (functionTokenVector.at(j).function.at(i)->getTokenID() == glslparser::PlotFunctionTerminals::ID_VARIABLE){
                std::string var = dynamic_cast<glslparser::IdentifierToken* const>(functionTokenVector.at(j).function.at(i))->getValue();
                ++(variables_[var[0]-PlotExpression::charOffset_].count);
            }
#ifdef VRN_PLOTEXPRESSION_DEBUG
            elog_ << functionTokenVector.at(j).function.at(i)->toString() << "\n";
#endif
        }
    }
    numberOfVariables_ = calculateNumberOfVariables();
    calculateDomain();
#ifdef VRN_PLOTEXPRESSION_DEBUG
    elog_ << "\nDomain:";
    elog_ << "\nNumber of Variables: " << numberOfVariables() << "\n";
    for (size_t i = 0; i < domain_.size(); ++i) {
        elog_ << i+1 << ". PartialFunction\n";
        for (size_t j = 0; j < domain_[i].size(); ++j) {
            elog_ << domain_[i].at(j).toString() << "\n";
        }
    }
    std::string visitorLog = elog_.str();
    LINFO("FunctionVisitor Debug Messages:\n" << visitorLog);
#endif
    LINFO("Expression parsed correctly.");
}

PlotExpression::~PlotExpression(){
    variables_.clear();
    //for (size_t i = 0; i < functionVector_.size(); ++i) {
    //    for (size_t j = 0; j < functionVector_[i].function.size(); ++j) {
    //        delete functionVector_[i].function[j];
    //    }
    //    for (size_t k = 0; k < functionVector_[i].interval.size(); ++k) {
    //        delete functionVector_[i].interval[k];
    //    }
    //}
    functionVector_.clear();
    domain_.clear();
}

void PlotExpression::initialize() {
    representation_ = "";
    expressionName_ = "f";
    numberOfVariables_ = 0;
    domain_.clear();
    variables_.resize(26);
    node_ = 0;
}

int PlotExpression::calculateNumberOfVariables() {
    int z = 0;
    for (size_t i = 0; i < variables_.size(); ++i) {
        if (variables_[i].count > 0) {
            ++z;
            variables_[i].numberOfVariable = z;
        }
    }
    return z;
}

void PlotExpression::evaluateDomain(std::stack<glslparser::Token*>& tokens, int index, int& number) {
    if (tokens.size() == 0)
        return;

    glslparser::Token* token = tokens.top();
    tokens.pop();
    if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_COLON) { // ":"
    }
    else if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_COMMA) {  // ","
        glslparser::BracketToken* const leftbracket = dynamic_cast<glslparser::BracketToken* const>(tokens.top());
        tokens.pop();
        plot_t leftrange = evaluate(tokens,std::vector<plot_t>());
        glslparser::BracketToken* const rightbracket = dynamic_cast<glslparser::BracketToken* const>(tokens.top());
        tokens.pop();
        plot_t rightrange = evaluate(tokens,std::vector<plot_t>());
        bool leftopen = leftbracket->getTokenID() == glslparser::PlotFunctionTerminals::ID_LPAREN;
        bool rightopen = rightbracket->getTokenID() == glslparser::PlotFunctionTerminals::ID_RPAREN;
        Interval<plot_t> interval = Interval<plot_t>(leftrange,rightrange,leftopen,rightopen);
        functionVector_[index].domain[number] = interval;
        ++number;
    }
    else if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_VERTICAL_BAR) { // "|"
//        ++number;
    }
    else if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_DASH) { // "-"
//        ++number;
    }
}

void PlotExpression::calculateDomain() {
    int count = std::max<int>(numberOfVariables_,1);
    for (size_t i = 0; i < functionVector_.size(); ++i) {
        functionVector_[i].domain.resize(count);
        for (int j = 0; j < count; ++j) {
            functionVector_[i].domain[j] = Interval<plot_t>(-std::numeric_limits<plot_t>::infinity(),
                    std::numeric_limits<plot_t>::infinity(),true,true);
        }
        std::stack<glslparser::Token*> token;
        for (size_t j = functionVector_.at(i).interval.size(); j > 0; --j) {
            token.push(functionVector_.at(i).interval[j-1]);
        }
        int number = 0;
        while (token.size() > 0 && number < count) {
            evaluateDomain(token, static_cast<int>(i), number);
        }
    }
    domain_.clear();
    for (size_t j = 0; j < functionVector_.size(); ++j) {
        domain_.push_back(functionVector_.at(j).domain);
    }
}

std::string PlotExpression::getVariable(int number) const {
    std::string result = "";
    char ch = 0;
    int z = 0;
    for (size_t i = 0; i < variables_.size(); ++i) {
        if (variables_[i].count > 0)
            ++z;
        if (z == number+1) {
            ch = static_cast<char>(i) + PlotExpression::charOffset_;
            break;
        }
    }
    result += ch;
    return result;
}

plot_t PlotExpression::evaluateAt(const std::vector<plot_t>& value) const {
    if (static_cast<int>(value.size()) < numberOfVariables())
        return std::numeric_limits<plot_t>::quiet_NaN();
    else {
        int k = -1;
        for (size_t i = 0; i < functionVector_.size(); ++i) {
            k = static_cast<int>(i);
            for (size_t j = 0; j < functionVector_[i].domain.size(); ++j) {
                if (!functionVector_[i].domain[j].contains(value[j]))
                    k = -1;
            }
            if (k != -1)
                break;
        }
        if (k >= 0) {
            std::stack<glslparser::Token*> token;
            for (size_t i = functionVector_.at(k).function.size(); i > 0; --i)
                token.push(functionVector_.at(k).function[i-1]);
            return evaluate(token,value);
        }
    }
    return std::numeric_limits<plot_t>::quiet_NaN();
}

plot_t PlotExpression::evaluate(std::stack<glslparser::Token*>& tokens, const std::vector<plot_t>& value) const {
    if (tokens.size() == 0)
        return std::numeric_limits<plot_t>::quiet_NaN();
    glslparser::Token* token = tokens.top();
    tokens.pop();
    if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_INTCONST ||
        token->getTokenID() == glslparser::PlotFunctionTerminals::ID_FLOATCONST) {
        return dynamic_cast<glslparser::ConstantToken* const>(token)->convert<plot_t>();
    }
    else if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_VARIABLE) {
        char var = dynamic_cast<glslparser::IdentifierToken* const>(token)->getValue()[0];
        int number = variables_.at(var-PlotExpression::charOffset_).numberOfVariable-1;
        return value.at(number);
    }
    else if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_FUNCTION_TERM) {
        return evaluate(tokens,value);
    }
    else if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_FUNCTION) {
        glslparser::FunctionToken* const function = dynamic_cast<glslparser::FunctionToken* const>(token);
        std::string stringFunction = function->getValue();
        plot_t plotvalue = evaluate(tokens,value);
        if (stringFunction == "abs") {
            return abs(plotvalue);
        }
        else if (stringFunction == "sqrt") {
            return std::sqrt(plotvalue);
        }
        else if (stringFunction == "sin") {
            return std::sin(plotvalue);
        }
        else if (stringFunction == "cos") {
            return std::cos(plotvalue);
        }
        else if (stringFunction == "tan") {
            return std::tan(plotvalue);
        }
        else if (stringFunction == "arcsin") {
            return std::asin(plotvalue);
        }
        else if (stringFunction == "arccos") {
            return std::acos(plotvalue);
        }
        else if (stringFunction == "arctan") {
            return std::atan(plotvalue);
        }
        else if (stringFunction == "sinh") {
            return std::sinh(plotvalue);
        }
        else if (stringFunction == "cosh") {
            return std::cosh(plotvalue);
        }
        else if (stringFunction == "tanh") {
            return std::tanh(plotvalue);
        }
        else if (stringFunction == "ln") {
            return std::log(plotvalue);
        }
        else if (stringFunction == "exp") {
            return std::exp(plotvalue);
        }
        else if (stringFunction == "log") {
            return std::log10(plotvalue);
        }
        else if (stringFunction == "fac") {
            plot_t result = plotvalue;
            while (plotvalue -1 > 0) {
                plotvalue -= 1;
                result *= plotvalue;
            }
            return result;
        }
        else if (stringFunction == "int") {
            return int(plotvalue);
        }
        else if (stringFunction == "floor") {
            return std::floor(plotvalue);
        }
        else if (stringFunction == "ceil") {
            return std::ceil(plotvalue);
        }
        else if (stringFunction == "rnd") {
            return std::floor(plotvalue+0.5);
        }
        else if (stringFunction == "sgn") {
            return plotvalue > 0 ? 1 : (plotvalue == 0 ? 0 : -1);
        }
        else if (stringFunction == "sgx") {
            return plotvalue >= 0 ? 1 : 0;
        }
    }
    else if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_PLUS) {
        glslparser::OperatorToken* const op = dynamic_cast<glslparser::OperatorToken* const>(token);
        if (op->getParameter() == 1)
            return evaluate(tokens,value);
        else if (op->getParameter() == 2) {
            plot_t value1 = evaluate(tokens,value);
            plot_t value2 = evaluate(tokens,value);
            return value1 + value2;
        }
    }
    else if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_DASH) {
        glslparser::OperatorToken* const op  = dynamic_cast<glslparser::OperatorToken* const>(token);
        if (op->getParameter() == 1)
            return -evaluate(tokens,value);
        else if (op->getParameter() == 2) {
            plot_t value1 = evaluate(tokens,value);
            plot_t value2 = evaluate(tokens,value);
            return value1 - value2;
        }
    }
    else if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_STAR) {
        glslparser::OperatorToken* const op = dynamic_cast<glslparser::OperatorToken* const>(token);
        if (op->getParameter() == 2) {
            plot_t value1 = evaluate(tokens,value);
            plot_t value2 = evaluate(tokens,value);
            return value1 * value2;
        }
    }
    else if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_SLASH) {
        glslparser::OperatorToken* const op = dynamic_cast<glslparser::OperatorToken* const>(token);
        if (op->getParameter() == 2) {
            plot_t value1 = evaluate(tokens,value);
            plot_t value2 = evaluate(tokens,value);
            return value1 / value2;
        }
    }
    else if (token->getTokenID() == glslparser::PlotFunctionTerminals::ID_CARET) {
        glslparser::OperatorToken* const op = dynamic_cast<glslparser::OperatorToken* const>(token);
        if (op->getParameter() == 2) {
            plot_t value1 = evaluate(tokens,value);
            plot_t value2 = evaluate(tokens,value);
            return pow(value1,value2);
        }
    }
    return std::numeric_limits<plot_t>::quiet_NaN();
}

int PlotExpression::numberOfVariables() const {
    return numberOfVariables_;
}

std::string PlotExpression::getExpressionString() const {
    return representation_;
}

void PlotExpression::setExpressionName(const std::string& expressionName) {
    if (expressionName.size() > 0)
        expressionName_ = expressionName;
}

const std::string& PlotExpression::getExpressionName() const {
    return expressionName_;
}

std::string PlotExpression::toString() const {
    std::string result = "";
    std::string partresult = "";
    for (int i = 0; i < numberOfVariables_; ++i) {
        partresult += getVariable(i);
        if (i < numberOfVariables_-1)
            partresult += ",";
    }
    if (numberOfVariables_ == 0)
        partresult = "x";
    result = expressionName_ + "("+ partresult + ")="+representation_;
    return result;
}

const std::vector<std::vector<Interval<plot_t> > >& PlotExpression::getDomain() const {
    return domain_;
}

void PlotExpression::deletePlotExpressionNodes() {
    delete node_;
}

} // namespace voreen
