/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "lp_eigenfunctionfit.h"
#include "lp_plotmodifier.h"

namespace voreen {

//deklaration of static member
FFS_Specification* LP_EigenFunctionFit::functionSpecification_ = 0;
LP_EigenFunctionFit::fitFunctionMode LP_EigenFunctionFit::FFM_ = LP_EigenFunctionFit::FFM_GAUSSIAN_SUM;


void LP_EigenFunctionFit::setFitFunction(fitFunctionMode ffm, FFS_Specification* specification){
    LP_EigenFunctionFit::FFM_ = ffm;
    LP_EigenFunctionFit::functionSpecification_ = specification;
}

void LP_EigenFunctionFit::functionFit(const PlotData* pInput, int pColumn, PlotData* pOutput){
    if(!pInput || !pOutput)
        return;
    //initialize variables
    int count = pInput->getRowsCount();         //number of rows
    int labelCount = pInput->getColumnCount();  //number of column
    std::stringstream ss;                       //label stream
    //switch fit function
    //****************************************************************************************
    if(FFM_ == FFM_GAUSSIAN_SUM){
        //initialize variables
        FFS_GAUSSIAN_SUM* ffs = static_cast<FFS_GAUSSIAN_SUM*>(LP_EigenFunctionFit::functionSpecification_);
        Eigen::VectorXf x(3*ffs->n+1);
        Eigen::VectorXf y(count);
        getStartParameter(x, pInput, pColumn);
        //make function fit
        FF_Gaussian_Sum functor(pInput,pColumn,ffs->n);
        Eigen::LevenbergMarquardt<FF_Gaussian_Sum, float> lm(functor);
        lm.parameters.ftol = 1e-6;
        lm.parameters.gtol = 1e-6;
        lm.parameters.xtol = 1e-6;
        lm.minimize(x);
        //test, how many column must be added
        if(!ffs->single){
            std::vector<float> values(count);
            //calculate sum
            functor(x,y);
            for(int i = 0; i < count; i++)
                values[i] = y(i)+static_cast<float>(pInput->getRow(i).getValueAt(pColumn));
            LP_PlotModifier::addColumn(pOutput,values,"Gaussian Sum of " + pInput->getColumnLabel(pColumn));
        }
        else {
            std::vector< std::vector<PlotCellValue> > cells(count);
        //calculate sum
            functor(x,y);
            for(int i = 0; i < count; i++){
                std::vector<PlotCellValue> vec(ffs->n+1);
                vec[0] = y(i)+static_cast<float>(pInput->getRow(i).getValueAt(pColumn));
                cells[i] = vec;
            }
        //calculate singles
            FF_Gaussian_Single functorSingle(pInput,pColumn,0,count-1);
            Eigen::VectorXf x4(4);
            x4(3) = x(3*ffs->n); //last value is const for all gaussians
            for(int i = 0; i < ffs->n; i++){
                //set x
                x4(0) = x(3*i); x4(1) = x(3*i+1); x4(2) = x(3*i+2);
                //set y
                functorSingle(x4,y);
                for(int j = 0; j < count; j++)
                    cells[j][i+1] = y(j)+static_cast<float>(pInput->getRow(j).getValueAt(pColumn));
            }
            LP_PlotModifier::addRowVector(pOutput,cells);
            //set column labels
                pOutput->setColumnLabel(labelCount,"Gaussian Sum of " + pInput->getColumnLabel(pColumn));
            for(int i = 0; i < ffs->n; i++){
                ss.clear();
                ss.str("Gaussian "); ss << i+1 << ". Single of " << pInput->getColumnLabel(pColumn);
                pOutput->setColumnLabel(labelCount+i+1,ss.str());
            }
        }

    } else
    //****************************************************************************************
    if(FFM_ == FFM_GAUSSIAN_SINGLE){
        FFS_GAUSSIAN_SINGLE* ffs = static_cast<FFS_GAUSSIAN_SINGLE*>(LP_EigenFunctionFit::functionSpecification_);
        Eigen::VectorXf x(4*ffs->n);
        Eigen::VectorXf x4(4);
        getStartParameter(x,pInput,pColumn);
        std::vector< std::vector<PlotCellValue > > cells(count, std::vector<PlotCellValue>(ffs->n+1)); //needed for single representation
        std::vector<float> values(count); //needed for non single representation
        //calculation begins
        for(int i = 0; i < ffs->n; i++) {
            x4(0) = x(4*i); x4(1) = x(4*i+1); x4(2) = x(4*i+2); x4(3) = x(4*i+3);
            Eigen::VectorXf y(ffs->intervals[i].second - ffs->intervals[i].first +1);
            FF_Gaussian_Single functor(pInput, pColumn, ffs->intervals[i].first, ffs->intervals[i].second);
            Eigen::LevenbergMarquardt<FF_Gaussian_Single, float> lm(functor);
            lm.parameters.ftol = 1e-6;
            lm.parameters.gtol = 1e-6;
            lm.parameters.xtol = 1e-6;
            lm.minimize(x4);
            //switch between singele and non single
            //calculate y
            if(!ffs->single){
                functor(x4,y);
                for(int j = 0; j < y.size(); j++)
                    values[ffs->intervals[i].first+j] = y(j)+static_cast<float>(pInput->getRow(ffs->intervals[i].first+j).getValueAt(pColumn));
            } else {
                //add sum
                functor(x4,y);
                for(int j = 0; j < y.size(); j++)
                    cells[ffs->intervals[i].first+j][0] = y(j)+static_cast<float>(pInput->getRow(ffs->intervals[i].first+j).getValueAt(pColumn));
                //add singles
                FF_Gaussian_Single functorSingleAll(pInput,pColumn,0,count-1);
                Eigen::VectorXf y_all(count);
                //set y
                functorSingleAll(x4,y_all);
                for(int j = 0; j < count; j++)
                    cells[j][i+1] = y_all(j)+static_cast<float>(pInput->getRow(j).getValueAt(pColumn));
            }
        }
        //fill output
        if(!ffs->single){
            LP_PlotModifier::addColumn(pOutput,values,"Gaussian Singles of " + pInput->getColumnLabel(pColumn));
        } else {
            LP_PlotModifier::addRowVector(pOutput,cells);
            //set column labels
                pOutput->setColumnLabel(labelCount,"Gaussian Sum of " + pInput->getColumnLabel(pColumn));
            for(int i = 0; i < ffs->n; i++){
                ss.clear();
                ss.str("Gaussian "); ss << i+1 << ". Single of " << pInput->getColumnLabel(pColumn);
                pOutput->setColumnLabel(labelCount+i+1,ss.str());
            }
        }
    } else return;

}

void LP_EigenFunctionFit::getStartParameter(Eigen::VectorXf &x, const PlotData* pData, const int column){
    if(FFM_ == FFM_GAUSSIAN_SUM){
        FFS_GAUSSIAN_SUM* ffs = static_cast<FFS_GAUSSIAN_SUM*>(LP_EigenFunctionFit::functionSpecification_);
        for(int i = 0; i < ffs->n ; i++){
            x(i*3) = LP_PlotModifier::findMax(pData,column);
            x(i*3+1) = LP_PlotModifier::findXValue(pData, static_cast<float>(i+1)/static_cast<float>(ffs->n+1));
            x(i*3+2) = 0.4f;
        }
        x(3*ffs->n) = LP_PlotModifier::findMin(pData,column);
    } else if (FFM_ == FFM_GAUSSIAN_SINGLE){
        FFS_GAUSSIAN_SINGLE* ffs = static_cast<FFS_GAUSSIAN_SINGLE*>(LP_EigenFunctionFit::functionSpecification_);
        int length = pData->getRowsCount()/ffs->n;
        int rest = pData->getRowsCount() % ffs->n;
        ffs->intervals.resize(ffs->n);
            //set parameter and intervals
            for(int i = 0; i < ffs->n ; i++){
                x(i*4) = LP_PlotModifier::findMax(pData,column);
                x(i*4+1) = LP_PlotModifier::findXValue(pData, static_cast<float>(i+1)/static_cast<float>(ffs->n+1));
                x(i*4+2) = 0.4f;
                x(i*4+3) = LP_PlotModifier::findMin(pData,column);
                if(i == 0)
                    ffs->intervals[i].first = 0;
                else
                    ffs->intervals[i].first = ffs->intervals[i-1].second +1;
                if(i < rest)
                    ffs->intervals[i].second = ffs->intervals[i].first + length ;
                else
                    ffs->intervals[i].second = ffs->intervals[i].first + length -1;
            }
    }
}


} // namespace voreen
