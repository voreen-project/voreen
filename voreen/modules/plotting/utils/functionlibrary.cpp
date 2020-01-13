/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "functionlibrary.h"
#include "../datastructures/plotexpression.h"
#include "tgt/logmanager.h"

#include <limits>

namespace voreen {

const std::string FunctionLibrary::loggerCat_("voreen.plotting.FunctionLibrary");

FunctionLibrary::FunctionLibrary()
{
}

FunctionLibrary::~FunctionLibrary(){
}

std::pair<plot_t,plot_t> FunctionLibrary::linearRegression(const std::vector<std::pair<plot_t,plot_t> >& points) {
    std::vector<plot_t> result;
    if (points.size() <= 0)
        return std::pair<plot_t,plot_t>(std::numeric_limits<plot_t>::quiet_NaN(),std::numeric_limits<plot_t>::quiet_NaN());

    plot_t sumx = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        if (points[i].first == points[i].first)
            sumx += points[i].first;
    }
    plot_t xMid = sumx/points.size();
    plot_t sumy = 0;

    for (size_t i = 0; i < points.size(); ++i) {
        if (points[i].second == points[i].second)
            sumy += points[i].second;
    }
    plot_t yMid = sumy/points.size();
    plot_t b = 0;
    plot_t numerator = 0;
    plot_t denominator = 0;
    plot_t diff1;
    plot_t diff2;
    for (size_t i = 0; i< points.size(); ++i) {
        diff1 = points[i].first - xMid;
        diff2 = points[i].second - yMid ;
        numerator += (diff1 * diff2);
        denominator += (diff1 * diff1);
    }
    b = numerator / denominator;
    plot_t a = yMid - ( b * xMid);

    return std::pair<plot_t,plot_t>(a,b);
}

std::pair<PlotExpression,plot_t> FunctionLibrary::fittingLinear(const std::vector<std::pair<plot_t,plot_t> >& points) {
    std::pair<plot_t,plot_t> result =  FunctionLibrary::linearRegression(points);
    if (result.first != result.first && result.second != result.second)
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    std::stringstream Str("");
    Str.imbue(std::locale::classic());
    Str << result.first << " + " << result.second << "*x";
    plot_t MSE = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        MSE += std::pow(result.first + result.second*points.at(i).first-points.at(i).second,2);
    }
    MSE = MSE/points.size();
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(Str.str()),MSE);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::fittingSqrt(const std::vector<std::pair<plot_t,plot_t> > &points) {
    std::vector<std::pair<plot_t,plot_t> > logPoints;
    for (size_t i = 0; i < points.size(); ++i) {
        logPoints.push_back(std::pair<plot_t,plot_t>(std::sqrt(points.at(i).first),points.at(i).second));
    }
    std::pair<plot_t,plot_t> result = FunctionLibrary::linearRegression(logPoints);
    if (result.first != result.first && result.second != result.second)
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    std::stringstream Str("");
    Str.imbue(std::locale::classic());
    Str << result.first << " + " << result.second << "*sqrt(x)";
    plot_t MSE = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        MSE += std::pow(result.first + result.second*std::sqrt(points.at(i).first)-points.at(i).second,2);
    }
    MSE = MSE/points.size();
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(Str.str()),MSE);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::fittingSquare(const std::vector<std::pair<plot_t,plot_t> >& points) {
    std::vector<std::pair<plot_t,plot_t> > logPoints;
    for (size_t i = 0; i < points.size(); ++i) {
        logPoints.push_back(std::pair<plot_t,plot_t>(std::pow(points.at(i).first,2),points.at(i).second));
    }
    std::pair<plot_t,plot_t> result = FunctionLibrary::linearRegression(logPoints);
    if (result.first != result.first && result.second != result.second)
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    std::stringstream Str("");
    Str.imbue(std::locale::classic());
    Str << result.first << " + " << result.second << "*x^2";
    plot_t MSE = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        MSE += std::pow(result.first + result.second*std::pow(points.at(i).first,2)-points.at(i).second,2);
    }
    MSE = MSE/points.size();
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(Str.str()),MSE);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::fittingCubic(const std::vector<std::pair<plot_t,plot_t> >& points) {
    std::vector<std::pair<plot_t,plot_t> > logPoints;
    for (size_t i = 0; i < points.size(); ++i) {
        logPoints.push_back(std::pair<plot_t,plot_t>(std::pow(points.at(i).first,3),points.at(i).second));
    }
    std::pair<plot_t,plot_t> result = FunctionLibrary::linearRegression(logPoints);
    if (result.first != result.first && result.second != result.second)
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    std::stringstream Str("");
    Str.imbue(std::locale::classic());
    Str << result.first << " + " << result.second << "*x^3";
    plot_t MSE = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        MSE += std::pow(result.first + result.second*std::pow(points.at(i).first,3)-points.at(i).second,2);
    }
    MSE = MSE/points.size();
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(Str.str()),MSE);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::fittingMulti(int dimension,const std::vector<std::pair<plot_t,plot_t> >& points) {
    std::vector<std::pair<plot_t,plot_t> > logPoints;
    for (size_t i = 0; i < points.size(); ++i) {
        logPoints.push_back(std::pair<plot_t,plot_t>(std::pow(points.at(i).first,dimension),points.at(i).second));
    }
    std::pair<plot_t,plot_t> result = FunctionLibrary::linearRegression(logPoints);
    if (result.first != result.first && result.second != result.second)
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    std::stringstream Str("");
    Str.imbue(std::locale::classic());
    Str << result.first << " + " << result.second << "*x^" << dimension;
    plot_t MSE = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        MSE += std::pow(result.first + result.second*std::pow(points.at(i).first,dimension)-points.at(i).second,2);
    }
    MSE = MSE/points.size();
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(Str.str()),MSE);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::constantSpline(const std::vector<std::pair<plot_t,plot_t> >& points) {
    std::stringstream str("");
    str.imbue(std::locale::classic());
    for (size_t i = 0; i < points.size()-1; ++i) {
        str << points[i].second << ":[" << points[i].first << "," << points[i+1].first << "]";
        if (i != points.size()-2)
            str << ";";
    }
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(str.str()),-1);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::simpleSpline(const std::vector<std::pair<plot_t,plot_t> >& points) {

    std::stringstream str("");
    str.imbue(std::locale::classic());
    plot_t m;
    plot_t b;
    for (size_t i=0; i<points.size()-1; ++i) {
        m = (points.at(i+1).second - points.at(i).second) / (points.at(i+1).first - points.at(i).first);
        b = points.at(i).second - m*points.at(i).first;
        str << m << " * x + " << b << ":[" << points.at(i).first << "," << points.at(i+1).first << "]";
        if (i != points.size()-2)
            str << ";";
    }
    std::string expr = str.str();
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(expr),-1);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}


std::string FunctionLibrary::bSplineN(int i, int p, const std::vector<std::pair<plot_t,plot_t> >& points) {
    std::stringstream str;
    std::string value1;
    std::string value2;
    if (p == 1) {
        str << "sgx(x-" << points[i].first << ")*sgx(-x+" << points[i+1].first << ")";
        return str.str();
    }
    else {
        value1 = bSplineN(i,p-1,points);
        value2 = bSplineN(i+1,p-1,points);
        plot_t a = points[i+p-1].first - points[i].first;
        plot_t b = points[i+p].first - points[i+1].first;
        if (points[i].first == 0.0 && a != 1)
            str << "(x/" << a << ")*";
        else if (points[i].first == 0.0 && a == 1)
            str << "x*";
        else if (points[i].first != 0.0 && a == 1)
            str << "(x-" << points[i].first << ")*";
        else
            str << "((x-" << points[i].first << ")/" << a << ")*";
        str << value1 << "+";
        if (points[i+p].first == 0 && b != 1)
            str << "(-x/" << b << ")*";
        else if (points[i+p].first == 0 && b == 1)
            str << "(-x)*";
        else if (points[i+p].first != 0 && b == 1)
            str << "(" << points[i+p].first << "-x)*";
        else
            str << "((" << points[i+p].first << "-x)/" << b << ")*";
        str << value2;
        return str.str();
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::bSpline(const std::vector<std::pair<plot_t,plot_t> >& points) {

    std::stringstream str("");
    str.imbue(std::locale::classic());
    plot_t k;
    plot_t m = 1;
    std::vector<plot_t> z;
    z.resize(points.size()+1);
    z.front() = 0.5;
    z.back() = 0.5;
    for (size_t i =0; i < points.size()-1; ++i) {
        m = (points.at(i+1).second - points.at(i).second) / (points.at(i+1).first - points.at(i).first);
        z[i+1] = -z[i] + 2 * (m);
    }

    for (size_t i=0; i<points.size()-1; ++i) {
        k = (z[i+1] - z[i]) / (2*(points.at(i+1).first - points.at(i).first));
        str << points.at(i).second << " + " << z[i] << " * (x - " << points.at(i).first << ") + " << k << " * ( x-" << points.at(i).first << ")^2";
        str << ":[" << points.at(i).first << "," << points.at(i+1).first << "]";
        if (i != points.size()-2)
            str << ";";
    }
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(str.str()),-1);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::squareSpline(const std::vector<std::pair<plot_t,plot_t> >& points) {

    std::stringstream str("");
    str.imbue(std::locale::classic());
    plot_t k;
    plot_t m = 1;
    std::vector<plot_t> z;
    z.resize(points.size()+1);
    z.front() = 0.5;
    z.back() = 0.5;
    for (size_t i =0; i < points.size()-1; ++i) {
        m = (points.at(i+1).second - points.at(i).second) / (points.at(i+1).first - points.at(i).first);
        z[i+1] = -z[i] + 2 * (m);
    }
    for (size_t i=0; i<points.size()-1; ++i) {
        k = (z[i+1] - z[i]) / (2*(points.at(i+1).first - points.at(i).first));
        str << points.at(i).second << " + " << z[i] << " * (x - " << points.at(i).first << ") + " << k << " * ( x-" << points.at(i).first << ")^2";
        str << ":[" << points.at(i).first << "," << points.at(i+1).first << "]";
        if (i != points.size()-2)
            str << ";";
    }
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(str.str()),-1);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::cubicSpline(const std::vector<std::pair<plot_t,plot_t> >& points) {

    std::stringstream str("");
    str.imbue(std::locale::classic());
    std::vector<plot_t> xi;
    std::vector<plot_t> z;
    std::vector<plot_t> h;
    std::vector<plot_t> y;
    std::vector<plot_t> a;
    std::vector<plot_t> b;
    std::vector<plot_t> c;
    xi.resize(points.size());
    y.resize(points.size());
    a.resize(points.size()-1);
    b.resize(points.size()-1);
    c.resize(points.size()-1);
    z.resize(points.size());
    h.resize(points.size());
    h.back() = 0;
    for (size_t i=0; i < points.size(); ++i) {
        if (i < points.size()-1)
            h[i] = points[i+1].first - points[i].first;
        xi[i] = points[i].first;
    }
    for (size_t i = 1; i < a.size(); ++i) {
        y[i] = 6*((points[i+1].second-points[i].second)/h[i] - (points[i].second- points[i-1].second)/h[i-1]);
        a[i] = h[i-1];
        b[i] = 2*(h[i-1]+h[i]);
        c[i] = h[i];
    }
    c[0] = h[0];
    b[0] = 2*h[0];
    b.back() = 2*(h[a.size()-2]+h[a.size()-1]);
    a.back() = h[a.size()-1];
    TridiagonalSolve(a,b,c,y,xi);
    z = xi;
    plot_t k1;
    plot_t k2;
    for (size_t i=0; i < points.size()-1; ++i) {
        k1 = points[i+1].second/h[i] - h[i]*z[i+1]/6;
        k2 = points[i].second/h[i] - h[i]*z[i]/6;
        str << "(" << z[i+1] << "*(x-" << points[i].first << ")^3+" << z[i] << "*("
            << points[i+1].first << "-x)^3)/" << 6*h[i];
        if (k1-k2 <0)
            str << "-x*" << k2-k1;
        else
            str << "+x*" << k1-k2;
        str << " + " << -points[i].first*k1 + points[i+1].first*k2 << "";
        str << ":[" << points.at(i).first << "," << points.at(i+1).first << "]";
        if (i < points.size()-2)
            str << ";";
    }
    std::string expr = str.str();
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(expr),-1);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::fittingLogarithmus(const std::vector<std::pair<plot_t,plot_t> >& points) {
    std::vector<std::pair<plot_t,plot_t> > logPoints;
    for (size_t i = 0; i < points.size(); ++i) {
        logPoints.push_back(std::pair<plot_t,plot_t>(std::log10(points.at(i).first),points.at(i).second));
    }
    std::pair<plot_t,plot_t> result = FunctionLibrary::linearRegression(logPoints);
    if (result.first != result.first && result.second != result.second)
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    std::stringstream Str("");
    Str.imbue(std::locale::classic());
    Str << result.first << " + " << result.second << "*log(x)";
    plot_t MSE = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        MSE += std::pow(result.first + result.second*std::log10(points.at(i).first)-points.at(i).second,2);
    }
    MSE = MSE/points.size();
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(Str.str()),MSE);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::fittingPower(const std::vector<std::pair<plot_t,plot_t> >& points) {
    std::vector<std::pair<plot_t,plot_t> > logPoints;
    for (size_t i = 0; i < points.size(); ++i) {
        logPoints.push_back(std::pair<plot_t,plot_t>(std::log10(points.at(i).first),std::log10(points.at(i).second)));
    }
    std::pair<plot_t,plot_t> result = FunctionLibrary::linearRegression(logPoints);
    if (result.first != result.first && result.second != result.second)
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    std::stringstream Str("");
    Str.imbue(std::locale::classic());
    Str << exp(result.first) << " * x^";
    if (result.second >= 0)
        Str << result.second;
    else
        Str << "(" << result.second << ")";
    plot_t MSE = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        MSE += std::pow(std::exp(result.first) * std::pow(points.at(i).first,result.second)-points.at(i).second,2);
    }
    MSE = MSE/points.size();
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(Str.str()),MSE);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::fittingExponential(const std::vector<std::pair<plot_t,plot_t> >& points) {
    std::vector<std::pair<plot_t,plot_t> > logPoints;
    for (size_t i = 0; i < points.size(); ++i) {
        logPoints.push_back(std::pair<plot_t,plot_t>(points.at(i).first,std::log10(points.at(i).second)));
    }
    std::pair<plot_t,plot_t> result = FunctionLibrary::linearRegression(logPoints);
    if (result.first != result.first && result.second != result.second)
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    std::stringstream Str("");
    Str.imbue(std::locale::classic());
    Str << exp(result.first) << " * (" << result.second << ")^x";
    plot_t MSE = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        MSE += std::pow(std::exp(result.first) * std::pow(result.second,points.at(i).first)-points.at(i).second,2);
    }
    MSE = MSE/points.size();
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(Str.str()),MSE);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::fittingSin(const std::vector<std::pair<plot_t,plot_t> >& points) {
    std::vector<std::pair<plot_t,plot_t> > logPoints;
    for (size_t i = 0; i < points.size(); ++i) {
        logPoints.push_back(std::pair<plot_t,plot_t>(points.at(i).first,std::asin(points.at(i).second)));
    }
    std::pair<plot_t,plot_t> result = FunctionLibrary::linearRegression(logPoints);
    if (result.first != result.first && result.second != result.second)
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    std::stringstream Str("");
    Str.imbue(std::locale::classic());
    Str << result.first << " + " << std::asin(result.second) << "*sin(x)";
    plot_t MSE = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        MSE += std::pow(result.first + std::asin(result.second)*std::sin(points.at(i).first)-points.at(i).second,2);
    }
    MSE = MSE/points.size();
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(Str.str()),MSE);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::fittingCos(const std::vector<std::pair<plot_t,plot_t> >& points) {
    std::vector<std::pair<plot_t,plot_t> > logPoints;
    for (size_t i = 0; i < points.size(); ++i) {
        logPoints.push_back(std::pair<plot_t,plot_t>(points.at(i).first,points.at(i).second));
    }
    std::pair<plot_t,plot_t> result = FunctionLibrary::linearRegression(logPoints);
    if (result.first != result.first && result.second != result.second)
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    std::stringstream Str("");
    Str.imbue(std::locale::classic());
    Str << result.first << " + " << std::acos(result.second) << "*cos(x)";
    plot_t MSE = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        MSE += std::pow(result.first + std::acos(result.second)*std::cos(points.at(i).first)-points.at(i).second,2);
    }
    MSE = MSE/points.size();
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(Str.str()),MSE);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }
}

std::pair<PlotExpression,plot_t> FunctionLibrary::fittingInterpol(const std::vector<std::pair<plot_t,plot_t> >& interpolPoints) {
    std::stringstream result;
    std::stringstream partresult;
    result.imbue(std::locale::classic());
    partresult.imbue(std::locale::classic());
    for ( size_t i=0; i<interpolPoints.size(); ++i) {
        partresult.clear();
        partresult.str("");
        for( size_t k=0; k<interpolPoints.size(); ++k) {
            if (k != i) {
                partresult << "((x-" << interpolPoints.at(k).first << ") / (" << interpolPoints.at(i).first
                    << " - " << interpolPoints.at(k).first << "))*";
            }
        }
        if (partresult.str().size() > 0) {
            if (result.str().size() > 0)
                result << "+";
            result << partresult.str()  << "(" << interpolPoints.at(i).second << ")";
        }
    }
    try {
        return std::pair<PlotExpression,plot_t>(PlotExpression(result.str()),-1);
    }
    catch (VoreenException& ve) {
        LERROR("Exception: " << ve.what());
        return std::pair<PlotExpression,plot_t>(PlotExpression(),-2);
    }


}

void FunctionLibrary::TridiagonalSolve(const std::vector<plot_t>& a1, const std::vector<plot_t>& a2,
                                       const std::vector<plot_t>& a3, const std::vector<plot_t>& y,
                                       std::vector<plot_t>& z)
{
    std::vector<plot_t> a;
    std::vector<plot_t> b;
    std::vector<plot_t> c;
    std::vector<plot_t> d;
    for (size_t j = 0; j < a1.size(); ++j) {
        a.push_back(a1[j]);
        b.push_back(a2[j]);
        c.push_back(a3[j]);
        d.push_back(y[j]);
    }
    if (b[0] == 0)
        return;
    /* Modify the coefficients. */
    c[0] /= b[0];    /* Division by zero risk. */
    d[0] /= b[0];    /* Division by zero would imply a singular matrix. */
    for (size_t i = 1; i < a.size(); ++i){
        plot_t id = 1.0 / (b[i] - c[i-1] * a[i]);  /* Division by zero risk. */
        c[i] *= id;                             /* Last value calculated is redundant. */
        d[i] = (d[i] - d[i-1] * a[i]) * id;
    }
    /* Now back substitute. */
    z[a.size() - 1] = d[a.size() - 1];
    for (int k = static_cast<int>(a.size()) - 2; k >= 0; --k) {
        z[k] = d[k] - c[k] * z[k + 1];
    }
}




//PlotFittingFunctionChomsky:: PlotFittingFunctionChomsky const {
// int i, j, k, n, m;
//
//
//double A[n][m+1];
//double L[n][m+1];
//
// // Einträge der Matrix A eingeben:
//    for (i=0;i<n;i++)
//      for (j=0;j<m+1;j++)
//        {
//         cout <<"Bitte Eintrag A[" <<i+1 <<"," <<j+1 <<"] eingeben:";
//         cin >>A[i][j];
//        }
//
//       for (i=0;i<n;i++)
//        for (j=0;j<m+1;j++)
//          L[i][j]=A[i][j];
//
//    // Jetzt kommt der Algorithmus für die Cholesky-Zerlegung:
//     int p;
//   for ( k = 0; k < n; k++ )
//    {
//            for ( j = 0; j < m+1; j++ )
//            L[k][m+1] = A[k][m+1] - L[k][j-k+m+1]*L[k][j-k+m+1];
//        if ( L[k][m+1] <= 0 ) break;
//        L[k][m+1] = sqrt( L[k][m+1] );
//
//        p=min(k+m, n);
//        for ( i = k+1; i < p; i++ )
//        {
//            for ( j = 0; j < k; j++ )
//                L[i][k-i+m+1] =A[i][k-i+m+1]- L[i][j]*L[k-i+m+1][j];
//            L[i][k-i+m+1] =L[i][m+1]/ L[k][m+1];
//           L[k][i] = L[i][k];
//        }
//    }
//
//    cout <<"\n\n Die Cholesky-Zerlegung der Matrix ist: \n";
//
//     for (i=0; i<n; i++)
//     {
//        {
//         for (j=0; j<m+1; j++)
//          cout << L[i][j] <<"\t";
//        }
//        cout <<"\n" ;
//     }
//
//   system("pause");
//    return 0;
//
//}

//float cubic_interpolate( float y0, float y1, float y2, float y3, float mu ) {

//double a, b, diff;
//        unsigned int j = 1;
//
//        while ((m_x[j] < input) && (j < count))
//            j++;
//
//        diff = m_x[j] - m_x[j-1];
//
//        a = (m_x[j] - input) / diff;    //div should not be 0
//        b = (input - m_x[j-1]) / diff;
//
//        return (a*m_y[j-1] + b*m_y[j] + ((a*a*a - a)*m_der[j - 1] +
//               (b*b*b - b)*m_der[j])*(diff*diff) / 6);

//std::vector<std::vector<plot_t>> matrix;
//for( size_t x=0; x<matrix.size(); ++x){
//    for(size_t y=0; y<matrix.at(0).size(); ++y) {
//        for( size_t k=0; k<matrix.size(); ++k){
//
//        }
//
//    }
//}


} // namespace voreen
