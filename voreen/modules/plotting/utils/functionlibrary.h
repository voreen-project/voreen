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

#ifndef VRN_FUNCTIONLIBRARY_H
#define VRN_FUNCTIONLIBRARY_H

#include "../datastructures/plotbase.h"

#include <vector>

namespace voreen {

class PlotExpression;

    /**
    * \brief this class provides  amount of function to fit Points to a function and solving linear systems.
    **/
class VRN_CORE_API FunctionLibrary {
public:

    /**
    * \brief enumeration for the type setting of plotdatasource and plotdatagroup processor and widget
    *    so they can use the same enumeration.
    **/
    enum ProcessorFunctionalityType {
        NONE = 0,
        SELECT = 1,
        COLUMNORDER = 2,
        FILTER = 3,
        GROUPBY = 4,
    };

    /**
    * \brief enumeration indetifiers the communication between PlotDatafitfunction processor and widget.
    **/
    enum RegressionType {
        NOREGRESSION = -1,
        LINEAR = 10, SQRT = 15,
        POLYNOMIAL = 20,
        SQUARE = 22,
        CUBIC = 23,
        SIMPLESPLINE = 30, CONSTANTSPLINE = 31,
        SQUARESPLINE = 32, BSPLINE = 33,
        CUBICSPLINE = 34,
        LOGARITHMIC = 40,
        POWER = 50,
        EXPONENTIAL = 60,
        SIN = 70,
        COS = 80,
        INTERPOLATION = 90
    };


    ~FunctionLibrary();

    /**
     * \brief calculate the lineare Regression about and gives back the parameter a and b
     * for the function y = a + b*x
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return both parameter of a linear function
     * \note  if the points are not regressionable the reasult are nan
     **/
    static std::pair<plot_t,plot_t> linearRegression(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the lineare Regression and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> fittingLinear(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the Squareroot Regression and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> fittingSqrt(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the Multidimension Regression and gives back the calculated expression und the MSE
     * y = a + b*x^dimension
     *
     * \param dimension the dimension of the Regression
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> fittingMulti(int dimension,const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the square Regression and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> fittingSquare(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the Cubic Regression and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> fittingCubic(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the constant Spline and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> constantSpline(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the simple Spline and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> simpleSpline(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the Square Spline and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> squareSpline(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the B-Spline and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> bSpline(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the cubic-spline and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> cubicSpline(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the logarithm Regression and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> fittingLogarithmus(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the power Regression and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> fittingPower(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the exponential Regression and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> fittingExponential(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the sin Regression and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> fittingSin(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the cos Regression and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> fittingCos(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the Interpolationsfunction and gives back the calculated expression und the MSE
     *
     * \param points  except a vector of pairs the first value is the x and the second is y
     * \return Pair of the Expression and MSE
     * \note  if the points are not regressionable the expression ist empty and the MSE is -2.
     *        if the function has no MSE the result of the MSE is -1;
     **/
    static std::pair<PlotExpression,plot_t> fittingInterpol(const std::vector<std::pair<plot_t,plot_t> >& points);

    /**
     * \brief calculate the result of a equation of a three diagonal matrix
     *
     * \param a the left diagonal
     * \param b the central diagonal
     * \param c the right diagonal
     * \param d the y values of the Points
     * \param x the x values of the Points and the factors for the Splinefunction
     **/
    static void TridiagonalSolve(const std::vector<plot_t>& a, const std::vector<plot_t>& b,
        const std::vector<plot_t>& c, const std::vector<plot_t>& d, std::vector<plot_t>& z);

private:

    FunctionLibrary();

    static std::string bSplineN(int i, int p, const std::vector<std::pair<plot_t,plot_t> >& points);

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_FUNCTIONLIBRARY_H
