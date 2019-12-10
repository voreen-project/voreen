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

#ifndef VRN_LP_EIGENFUNCTIONFIT_H
#define VRN_LP_EIGENFUNCTIONFIT_H

//plotdata
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"

//eigen
#include "Eigen/Dense"
#include "unsupported/Eigen/NonLinearOptimization"

//utils
#include "lp_plotmodifier.h"
#include <sstream>

namespace voreen {
//**********************************************************************************
// Structs for Fit Functions
//**********************************************************************************
//Basic struct
struct FF_Functor {
    FF_Functor(const PlotData* pData, int pColumn, int input, int value) : data_(pData), column_(pColumn), inputs_(input), values_(value) {}

    const PlotData* data_;
    int column_;
    int inputs_; //< number of parameters
    int values_; //< number of values i.e. rows in data_

    int inputs() const { return inputs_; }
    int values() const { return values_; }
};

struct FFS_Specification{};
//**********************************************************************************
//gaussian sum
struct FF_Gaussian_Sum : FF_Functor {
    FF_Gaussian_Sum(const PlotData* pData, int pColumn, int pN) : FF_Functor(pData, pColumn,3*pN+1, pData->getRowsCount()), n_(pN){}

    int n_; //<number of summands

    int operator()(const Eigen::VectorXf &x, Eigen::VectorXf &fvec) const
    {
        float value, helpV;
        for(int i = 0; i < values_; i++){
            value = 0.0f;
            for(int j = 0; j < n_; j++){
                helpV = (static_cast<float>(data_->getRow(i).getValueAt(0)) - x(3*j+1))/x(3*j+2);
                helpV = -1.0f*helpV*helpV;
                value += x(3*j)*exp(helpV);
            }
            fvec(i) = value+x(3*n_)-static_cast<float>(data_->getRow(i).getValueAt(column_));
        }
        return 0;
    }

    int df(const Eigen::VectorXf &x, Eigen::MatrixXf &fjac) const
    {
        float X, helpV;
        for(int i = 0; i < values_ ; i++){
            X = static_cast<float>(data_->getRow(i).getValueAt(0));
            for(int j = 0; j < n_ ; j++)
            {
                helpV = (X - x(3*j+1))/x(3*j+2);
                helpV *= -1.0f*helpV;
                fjac(i,3*j) = exp(helpV);
                fjac(i,3*j+1) = exp(helpV)*x(3*j)*(2*X-2*x(3*j+1))/(x(3*j+2)*x(3*j+2));
                fjac(i,3*j+2) = exp(helpV)*x(3*j)*(2*(X-x(3*j+1))*(X-x(3*j+1)))/(x(3*j+2)*x(3*j+2)*x(3*j+2));
            }
            fjac(i,3*n_) = 1;
        }
        return 0;
    }
};

struct FFS_GAUSSIAN_SUM : FFS_Specification{
    FFS_GAUSSIAN_SUM() : n(0),single(false){}
    FFS_GAUSSIAN_SUM(int pN, bool pSingle) : n(pN),single(pSingle){}
    int n;
    bool single;
};
//**********************************************************************************
//guassian single
struct FF_Gaussian_Single : FF_Functor{
    //constructor
    FF_Gaussian_Single(const PlotData* pData, int pColumn, int start, int end) : FF_Functor(pData, pColumn, 4, end - start + 1), start_(start), end_(end){}

    int start_;
    int end_;

    int operator()(const Eigen::VectorXf &x, Eigen::VectorXf &fvec) const
    {
        float helpV;
        for(int i = 0; i < values_; i++){
            helpV = (static_cast<float>(data_->getRow(start_+i).getValueAt(0)) - x(1))/x(2);
            helpV = -1.0f*helpV*helpV;
            fvec(i) = x(0)*exp(helpV)+x(3)-static_cast<float>(data_->getRow(start_+i).getValueAt(column_));
        }
        return 0;
    }

    int df(const Eigen::VectorXf &x, Eigen::MatrixXf &fjac) const
    {
        float X, helpV;
        for(int i = 0; i < values_ ; i++){
            X = static_cast<float>(data_->getRow(start_+i).getValueAt(0));
            helpV = (X - x(1))/x(2);
            helpV *= -1.0f*helpV;

            fjac(i,0) = exp(helpV);
            fjac(i,1) = exp(helpV)*x(0)*(2*X-2*x(1))/(x(2)*x(2));
            fjac(i,2) = exp(helpV)*x(0)*(2*(X-x(1))*(X-x(1)))/(x(2)*x(2)*x(2));
            fjac(i,3) = 1;
        }
        return 0;
    }
};

struct FFS_GAUSSIAN_SINGLE : FFS_Specification{
    FFS_GAUSSIAN_SINGLE() : n(0),single(false),intervals(){}
    FFS_GAUSSIAN_SINGLE(int pN, bool pSingle) : n(pN),single(pSingle),intervals(){}
    int n;
    bool single;
    std::vector< std::pair<int,int> > intervals;
};


//********************************************************************************************
//            Main Class
//********************************************************************************************
/**
 * Utility class beeing used in lineprofile.h
 * This class encapsulate the functions for the function fit.
 */
class LP_EigenFunctionFit {
public:
    /** enum for the types of function fit */
    enum fitFunctionMode {
        FFM_GAUSSIAN_SUM,        //<fits a gaussian sum
        FFM_GAUSSIAN_SINGLE      //<fits single gaussians
    };

    /**
     * Sets the members of this class
     * @param ff                the type of fit function beeing used
     * @param specification     additional specification for "functionFit"
     */
    static void setFitFunction(fitFunctionMode ffm, FFS_Specification* specification = 0);
    /**
     * Main function for function fit.
     * This function adds the fit function specified in "setFitFunction" from pInput into pOutput.
     * @param pInput    input plotdata
     * @param pColumn   used column in pInput for function fit
     * @param pOutput   plotdat, which will be modified
     */
    static void functionFit(const PlotData* pInput, int pColumn, PlotData* pOutput);
private:
    /**
     * Calculate start paramerter depending on FFM_.
     * @param x         result start parameter vector
     * @param pData     the plotdata
     * @param column    the column
     */
    static void getStartParameter(Eigen::VectorXf &x, const PlotData* pData, const int column);
    /** fit function mode beeing used*/
    static fitFunctionMode FFM_;
    /** Specifies the FFM */
    static FFS_Specification* functionSpecification_;

};



}   //namespace

#endif // VRN_LP_EIGENFUNCTIONFIT_H
