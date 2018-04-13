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

#ifndef VRN_VOLUMEOPERATORDISTANCETRANSFORM_H
#define VRN_VOLUMEOPERATORDISTANCETRANSFORM_H

#include "voreen/core/datastructures/volume/volumeoperator.h"


/**
 * For Reference: The Operators defined here realize the algorithm of Saito and
 * Toriwaki from '94 with the optimizations proposed by hirata in '96 and the
 * special cases from Meijster et al from '02
 */


//This Namespace contains helper functions for special handling for infinite distances
namespace DistanceFunction{
    static const float INF = FLT_MAX;

    inline float Add(float a, float b){
        if ((a==INF) || (b==INF))
            return INF;
        else
            return a+b;
    }

    inline float Mult(float a, float b){
        if ((a==INF) || (b==INF))
            return INF;
        else
            return a*b;
    }

    inline float Neg(float a){
        if (a == INF) {
            return INF;
        }
        else {
            return -a;
        }
    }

    inline float Div(float divid, float divis) {
        if (divis == 0 || divid == INF){
            return  INF;
        }else{
            return  divid / divis;
        }
    }

    /**
     * Evaluates a parabola of the form
     * (x-i)^2 + gi2 at position x
     *
     */
    inline float EvalParabola(float x, float i, float gi2){
        return Add((x-i)*(x-i), gi2);
    }

    /**
     * Computes the abscissa of the intersection point of two parabolas
     * i.e. solves (i-x)^2+gi2 = (u-x)^2+gu2 at integer coordinates
     */
    inline float PointOfParabolaIntersection(float i, float u, float gi2, float gu2) {
        return Div(Add(Add((u*u - i*i),gu2), Neg(gi2) ), 2*(u-i));
    }

    /**
     * Evaluates a function of the form |x-i|+gi for the Manhattan metric
     */
    inline float EvalVShape(float x, float i, float gi){
        return Add(std::fabs(x - i), gi);
    }

    /**
     * Solves |x-i|+gi = |x-u|+gu. The Intersection of two Manhattan metric functions
     */
    inline float PointOfVShapeIntersection(float i, float u, float gi, float gu){
        if(gu >= gi + u - i){
            return INF;
        }
        if(gi > gu + u - i){
            return -INF;
        }
        return Div(Add(gu, Neg(gi)) + u + i, 2);
    }

    /**
     * Evaluates a function for the Chessboard metric, i.e. max(|x-i|,gi)
     */
    inline float EvalChebyFunction(float x, float i, float gi){
        return std::max(std::fabs(x-i), gi);
    }

    /**
     * Computes the intersection of two functions that describe a chessboard metric
     * i.e. max(|x-i|,gi) = max(|x-u|,gu)
     */
    inline float PointOfChebyFunctionIntersection(float i, float u, float gi, float gu){
        if(gi <= gu){
            return std::max(Add(i, gu), Div(i+u, 2));
        }else{
            return std::min(u-gi, Div(i+u, 2));
        }
    }
};

namespace voreen {

// Generic implementation Squared Euclidean distance transform:
class VRN_CORE_API VolumeOperatorSquaredEuclideanDistanceTransformBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressReporter* progressReporter = 0) const = 0;
};

//specific implementation squared euclidean distance transform
template<typename T>
class VolumeOperatorSquaredEuclideanDistanceTransformGeneric : public VolumeOperatorSquaredEuclideanDistanceTransformBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressReporter* progressReporter = 0) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

// Generic implementation Euclidean distance transform:
class VRN_CORE_API VolumeOperatorEuclideanDistanceTransformBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressReporter* progressReporter = 0) const = 0;
};

//specific implementation euclidean distance transform
template<typename T>
class VolumeOperatorEuclideanDistanceTransformGeneric : public VolumeOperatorEuclideanDistanceTransformBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressReporter* progressReporter = 0) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

// Generic implementation Manhattan distance transform:
class VRN_CORE_API VolumeOperatorManhattanDistanceTransformBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressReporter* progressReporter = 0) const = 0;
};

//specific implementation Manhattan distance transform
template<typename T>
class VolumeOperatorManhattanDistanceTransformGeneric : public VolumeOperatorManhattanDistanceTransformBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressReporter* progressReporter = 0) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

// Generic implementation Chebychev distance transform:
class VRN_CORE_API VolumeOperatorChebychevDistanceTransformBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressReporter* progressReporter = 0) const = 0;
};

//specific implementation Chebychev distance transform
template<typename T>
class VolumeOperatorChebychevDistanceTransformGeneric : public VolumeOperatorChebychevDistanceTransformBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressReporter* progressReporter = 0) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorSquaredEuclideanDistanceTransformGeneric<T>::apply(const VolumeBase* vh, ProgressReporter* pR) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if (!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if (!va)
        return 0;

    tgt::svec3 dim = va->getDimensions();

    VolumeAtomic<float>* tmp = new VolumeAtomic<float>(dim, true);
    VolumeAtomic<float>* output = new VolumeAtomic<float>(dim, true);

    const VolumeAtomic<T>* inputSource = va;
    const VolumeAtomic<float>* currentSource = 0;
    VolumeAtomic<float>* currentDestination = output;

    //First Compute the 1 dimensional distance function
    for(size_t z = 0; z < dim.z; ++z){
        for(size_t y = 0; y < dim.y; ++y){
            //begin forward scan
            if(inputSource->voxel(0,y,z) == T(0)){
                currentDestination->voxel(0,y,z) = 0.0f;
            }else{
                currentDestination->voxel(0,y,z) = DistanceFunction::INF;
            }
            //increment the distance for each foreground voxel
            for(size_t x = 1; x < dim.x; ++x){
                if(inputSource->voxel(x,y,z) == T(0)){
                    currentDestination->voxel(x,y,z) = 0.0f;
                }else{
                    currentDestination->voxel(x,y,z) = DistanceFunction::Add(1.0f, currentDestination->voxel(x-1,y,z));
                }
            }

            //do the same backwards
            for(size_t x = dim.x - 2; x != SIZE_MAX; --x){
                if(currentDestination->voxel(x+1,y,z) < currentDestination->voxel(x,y,z)){
                    currentDestination->voxel(x,y,z) = DistanceFunction::Add(1.0f, currentDestination->voxel(x+1,y,z));
                }
            }
        }

        if(pR){
            //report progress!
            pR->setProgress(static_cast<float>(z)/static_cast<float>(dim.z) / 3.f);
        }
    }

    //swap Buffers
    currentSource = currentDestination;
    currentDestination = tmp;

    std::vector<tgt::svec2> LowerEnvelope; //"stack" to save the lower envelope of parabolas
    LowerEnvelope.reserve(dim.y);

    size_t w;

    //solve Distance Transform slice-wise
    for (size_t z = 0; z< dim.z; ++z){
        for (size_t x = 0; x < dim.x; ++x){
            //clear stack
            LowerEnvelope.clear();
            LowerEnvelope.push_back(tgt::svec2(0,0));

            //Forward Scan
            for(size_t y = 1; y < dim.y ; ++y){
                //determine contribution of topmost parabola to the lowerenvelope repeatedly
                //if it doesnt contribute delete it
                while (!LowerEnvelope.empty() &&
                    (DistanceFunction::EvalParabola(static_cast<float>(LowerEnvelope.back().y),static_cast<float>(LowerEnvelope.back().x),
                    DistanceFunction::Mult(currentSource->voxel(x,LowerEnvelope.back().x,z),currentSource->voxel(x,LowerEnvelope.back().x,z))) >
                    DistanceFunction::EvalParabola(static_cast<float>(LowerEnvelope.back().y),static_cast<float>(y),
                    DistanceFunction::Mult(currentSource->voxel(x,y,z),currentSource->voxel(x,y,z))))){
                    LowerEnvelope.pop_back();
                }

                //lowerenvelope can't be empty of course
                if(LowerEnvelope.empty()){
                    LowerEnvelope.push_back(tgt::svec2(y,0));
                }
                else{
                    //compute intersection and save parabola to the lowerenvelope
                    w = 1 + static_cast<size_t>(
                        DistanceFunction::PointOfParabolaIntersection(static_cast<float>(LowerEnvelope.back().x),
                        static_cast<float>(y),
                        DistanceFunction::Mult(currentSource->voxel(x,LowerEnvelope.back().x,z),currentSource->voxel(x,LowerEnvelope.back().x,z)),
                        DistanceFunction::Mult(currentSource->voxel(x,y,z),currentSource->voxel(x,y,z))));

                    if(w < dim.y){
                        LowerEnvelope.push_back(tgt::svec2(y,w));
                    }
                }
            }

            //Backward Scan
            //Evaluate the Lowerenvelope function at integer coordinates to retrieve the exact squared euclidean distance
            for(size_t y = dim.y-1; y != SIZE_MAX; --y)
            {
                currentDestination->voxel(x,y,z) = DistanceFunction::EvalParabola(static_cast<float>(y),
                                                                static_cast<float>(LowerEnvelope.back().x),
                                                                DistanceFunction::Mult(currentSource->voxel(x,LowerEnvelope.back().x,z),currentSource->voxel(x,LowerEnvelope.back().x,z)));
                if(y == LowerEnvelope.back().y){
                    LowerEnvelope.pop_back();
                }
            }
        }

        if(pR){
            //report progress!
            pR->setProgress(1.f/3.f + static_cast<float>(z)/static_cast<float>(dim.z) / 3.f);
        }
    }

    currentSource = currentDestination;
    currentDestination = output;
    LowerEnvelope.reserve(dim.z);

    //repeat the same procedure as above but now analyze interslice distance by swapping y and z
    for(size_t y = 0; y< dim.y; ++y){
        for(size_t x = 0; x < dim.x; ++x){
            LowerEnvelope.clear();
            LowerEnvelope.push_back(tgt::svec2(0,0));

            //Forward Scan
            for(size_t z = 1; z < dim.z ; ++z){
                while(!LowerEnvelope.empty() &&
                    (DistanceFunction::EvalParabola(static_cast<float>(LowerEnvelope.back().y),static_cast<float>(LowerEnvelope.back().x), currentSource->voxel(x, y, LowerEnvelope.back().x))) >
                    DistanceFunction::EvalParabola(static_cast<float>(LowerEnvelope.back().y),static_cast<float>(z),currentSource->voxel(x,y,z))){
                    LowerEnvelope.pop_back();
                }

                if(LowerEnvelope.empty()){
                    LowerEnvelope.push_back(tgt::svec2(z,0));
                }
                else{
                    w = 1 + static_cast<size_t>(
                        DistanceFunction::PointOfParabolaIntersection(static_cast<float>(LowerEnvelope.back().x),
                        static_cast<float>(z),
                        currentSource->voxel(x,y,LowerEnvelope.back().x),
                        currentSource->voxel(x,y,z)));

                    if(w < dim.y){
                        LowerEnvelope.push_back(tgt::svec2(z,w));
                    }
                }
            }

            //Backward Scan
            for(size_t z = dim.z-1; z != SIZE_MAX; --z){
                currentDestination->voxel(x,y,z) = DistanceFunction::EvalParabola(static_cast<float>(z),static_cast<float>(LowerEnvelope.back().x),currentSource->voxel(x,y,LowerEnvelope.back().x));
                if(z == LowerEnvelope.back().y){
                    LowerEnvelope.pop_back();
                }
            }
        }

        if(pR){
            //report progress!
            pR->setProgress(2.f/3.f + static_cast<float>(y)/static_cast<float>(dim.y) / 3.f);
        }
    }

    delete tmp;

    if (pR)
        pR->setProgress(1.f);

    return new Volume(currentDestination, vh);
}

template<typename T>
Volume* VolumeOperatorEuclideanDistanceTransformGeneric<T>::apply(const VolumeBase* vb, ProgressReporter* pR) const {
    if (pR)
        pR->setProgressRange(tgt::vec2(0.f, 0.75f));
    Volume* v = UniversalUnaryVolumeOperatorGeneric<VolumeOperatorSquaredEuclideanDistanceTransformBase>::APPLY_OP(vb, pR);


    VolumeAtomic<float>* va = v->getWritableRepresentation<VolumeAtomic<float> >();

    if(pR){
        pR->setProgressRange(tgt::vec2(0.f,1.f));
        pR->setProgress(0.75f);
    }

    tgt::svec3 dim = va->getDimensions();
    VRN_FOR_EACH_VOXEL_WITH_PROGRESS_SUB_TASK(idx, tgt::svec3(0,0,0), dim, pR, 0.75f, 0.25f){
        va->voxel(idx) = sqrtf(va->voxel(idx));
    }

    if (pR){
        pR->setProgress(1.f);
    }

    return v;
}

template<typename T>
Volume* VolumeOperatorManhattanDistanceTransformGeneric<T>::apply(const VolumeBase* vb, ProgressReporter* pR) const{
    const VolumeRAM* v = vb->getRepresentation<VolumeRAM>();
    if (!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if (!va)
        return 0;

    tgt::svec3 dim = va->getDimensions();

    VolumeAtomic<float>* tmp = new VolumeAtomic<float>(dim, true);
    VolumeAtomic<float>* output = new VolumeAtomic<float>(dim, true);

    const VolumeAtomic<T>* inputSource = va;
    const VolumeAtomic<float>* currentSource = tmp;
    VolumeAtomic<float>* currentDestination = output;

    for(size_t z = 0; z < dim.z; ++z){
        for(size_t y = 0; y < dim.y; ++y){
            if(inputSource->voxel(0,y,z) == T(0)){
                currentDestination->voxel(0,y,z) = 0.0f;
            }else{
                currentDestination->voxel(0,y,z) = DistanceFunction::INF;
            }
            for(size_t x = 1; x < dim.x; ++x){
                if(inputSource->voxel(x,y,z) == T(0)){
                    currentDestination->voxel(x,y,z) = 0.0f;
                }else{
                    currentDestination->voxel(x,y,z) = DistanceFunction::Add(1.0f, currentDestination->voxel(x-1,y,z));
                }
            }

            for(size_t x = dim.x - 2; x != SIZE_MAX; --x){
                if(currentDestination->voxel(x+1,y,z) < currentDestination->voxel(x,y,z)){
                    currentDestination->voxel(x,y,z) = DistanceFunction::Add(1.0f, currentDestination->voxel(x+1,y,z));
                }
            }
        }

        if(pR){
            //report progress!
            pR->setProgress(static_cast<float>(z)/static_cast<float>(dim.z) / 3.f);
        }
    }

    currentSource = currentDestination;
    currentDestination = tmp;

    std::vector<tgt::svec2> LowerEnvelope;
    LowerEnvelope.reserve(dim.y);

    size_t w;

    for(size_t z = 0; z< dim.z; ++z){
        for(size_t x = 0; x < dim.x; ++x){
            LowerEnvelope.clear();
            LowerEnvelope.push_back(tgt::svec2(0,0));

            //Forward Scan
            for(size_t y = 1; y < dim.y ; ++y){
                while (!LowerEnvelope.empty() &&
                    DistanceFunction::EvalVShape(static_cast<float>(LowerEnvelope.back().y),static_cast<float>(LowerEnvelope.back().x),currentSource->voxel(x,LowerEnvelope.back().x,z)) >
                    DistanceFunction::EvalVShape(static_cast<float>(LowerEnvelope.back().y),static_cast<float>(y),currentSource->voxel(x,y,z))){
                    LowerEnvelope.pop_back();
                }

                if(LowerEnvelope.empty()){
                    LowerEnvelope.push_back(tgt::svec2(y,0));
                }
                else{
                    w = 1 + static_cast<size_t>(
                        DistanceFunction::PointOfVShapeIntersection(static_cast<float>(LowerEnvelope.back().x),
                        static_cast<float>(y),
                        currentSource->voxel(x,LowerEnvelope.back().x,z),
                        currentSource->voxel(x,y,z)));

                    if (w < dim.y){
                        LowerEnvelope.push_back(tgt::svec2(y,w));
                    }
                }
            }

            //Backward Scan
            for(size_t y = dim.y-1; y != SIZE_MAX; --y)
            {
                currentDestination->voxel(x,y,z) = DistanceFunction::EvalVShape(static_cast<float>(y),
                                                                static_cast<float>(LowerEnvelope.back().x),
                                                                currentSource->voxel(x,LowerEnvelope.back().x,z));
                if(y == LowerEnvelope.back().y){
                    LowerEnvelope.pop_back();
                }
            }
        }

        if(pR){
            //report progress!
            pR->setProgress(1.f/3.f + static_cast<float>(z)/static_cast<float>(dim.z) / 3.f);
        }
    }

    currentSource = currentDestination;
    currentDestination = output;
    LowerEnvelope.reserve(dim.z);

    for (size_t y = 0; y< dim.y; ++y){
        for (size_t x = 0; x < dim.x; ++x){
            LowerEnvelope.clear();
            LowerEnvelope.push_back(tgt::svec2(0,0));

            //Forward Scan
            for (size_t z = 1; z < dim.z ; ++z){
                while (!LowerEnvelope.empty() &&
                    DistanceFunction::EvalVShape(static_cast<float>(LowerEnvelope.back().y),static_cast<float>(LowerEnvelope.back().x),currentSource->voxel(x,y,LowerEnvelope.back().x)) >
                    DistanceFunction::EvalVShape(static_cast<float>(LowerEnvelope.back().y),static_cast<float>(z),currentSource->voxel(x,y,z))){
                    LowerEnvelope.pop_back();
                }

                if(LowerEnvelope.empty()){
                    LowerEnvelope.push_back(tgt::svec2(z,0));
                }
                else{
                    w = 1 + static_cast<size_t>(
                        DistanceFunction::PointOfVShapeIntersection(static_cast<float>(LowerEnvelope.back().x),
                        static_cast<float>(z),
                        currentSource->voxel(x,y,LowerEnvelope.back().x),
                        currentSource->voxel(x,y,z)));

                    if (w < dim.y){
                        LowerEnvelope.push_back(tgt::svec2(z,w));
                    }
                }
            }

            //Backward Scan
            for(size_t z = dim.z-1; z != SIZE_MAX; --z){
                currentDestination->voxel(x,y,z) = DistanceFunction::EvalVShape(static_cast<float>(z),
                                                                static_cast<float>(LowerEnvelope.back().x),
                                                                currentSource->voxel(x,y,LowerEnvelope.back().x));
                if(z == LowerEnvelope.back().y){
                    LowerEnvelope.pop_back();
                }
            }
        }

        if(pR){
            //report progress!
            pR->setProgress(2.f/3.f + static_cast<float>(y)/static_cast<float>(dim.y) / 3.f);
        }
    }

    delete tmp;

    if (pR)
        pR->setProgress(1.f);

    return new Volume(currentDestination, vb);
}

template<typename T>
Volume* VolumeOperatorChebychevDistanceTransformGeneric<T>::apply(const VolumeBase* vb, ProgressReporter* pR) const{
    const VolumeRAM* v = vb->getRepresentation<VolumeRAM>();
    if (!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if (!va)
        return 0;

    tgt::svec3 dim = va->getDimensions();

    VolumeAtomic<float>* tmp = new VolumeAtomic<float>(dim, true);
    VolumeAtomic<float>* output = new VolumeAtomic<float>(dim, true);

    const VolumeAtomic<T>* inputSource = va;
    const VolumeAtomic<float>* currentSource = tmp;
    VolumeAtomic<float>* currentDestination = output;

    for(size_t z = 0; z < dim.z; ++z){
        for(size_t y = 0; y < dim.y; ++y){
            if(inputSource->voxel(0,y,z) == T(0)){
                currentDestination->voxel(0,y,z) = 0.0f;
            }else{
                currentDestination->voxel(0,y,z) = DistanceFunction::INF;
            }
            for(size_t x = 1; x < dim.x; ++x){
                if(inputSource->voxel(x,y,z) == T(0)){
                    currentDestination->voxel(x,y,z) = 0.0f;
                }else{
                    currentDestination->voxel(x,y,z) = DistanceFunction::Add(1.0f, currentDestination->voxel(x-1,y,z));
                }
            }

            for(size_t x = dim.x - 2; x != SIZE_MAX; --x){
                if(currentDestination->voxel(x+1,y,z) < currentDestination->voxel(x,y,z)){
                    currentDestination->voxel(x,y,z) = DistanceFunction::Add(1.0f, currentDestination->voxel(x+1,y,z));
                }
            }
        }

        if(pR){
            //report progress!
            pR->setProgress(static_cast<float>(z)/static_cast<float>(dim.z) / 3.f);
        }
    }

    currentSource = currentDestination;
    currentDestination = tmp;

    std::vector<tgt::svec2> LowerEnvelope;
    LowerEnvelope.reserve(dim.y);

    size_t w;

    for(size_t z = 0; z< dim.z; ++z){
        for(size_t x = 0; x < dim.x; ++x){
            LowerEnvelope.clear();
            LowerEnvelope.push_back(tgt::svec2(0,0));

            //Forward Scan
            for(size_t y = 1; y < dim.y ; ++y){
                while (!LowerEnvelope.empty() &&
                    DistanceFunction::EvalChebyFunction(static_cast<float>(LowerEnvelope.back().y),static_cast<float>(LowerEnvelope.back().x),currentSource->voxel(x,LowerEnvelope.back().x,z)) >
                    DistanceFunction::EvalChebyFunction(static_cast<float>(LowerEnvelope.back().y),static_cast<float>(y),currentSource->voxel(x,y,z))){
                    LowerEnvelope.pop_back();
                }

                if(LowerEnvelope.empty()){
                    LowerEnvelope.push_back(tgt::svec2(y,0));
                }
                else{
                    w = 1 + static_cast<size_t>(
                        DistanceFunction::PointOfChebyFunctionIntersection(static_cast<float>(LowerEnvelope.back().x),
                        static_cast<float>(y),
                        currentSource->voxel(x,LowerEnvelope.back().x,z),
                        currentSource->voxel(x,y,z)));

                    if (w < dim.y){
                        LowerEnvelope.push_back(tgt::svec2(y,w));
                    }
                }
            }

            //Backward Scan
            for(size_t y = dim.y-1; y != SIZE_MAX; --y)
            {
                currentDestination->voxel(x,y,z) = DistanceFunction::EvalChebyFunction(static_cast<float>(y),
                                                                static_cast<float>(LowerEnvelope.back().x),
                                                                currentSource->voxel(x,LowerEnvelope.back().x,z));
                if(y == LowerEnvelope.back().y){
                    LowerEnvelope.pop_back();
                }
            }
        }

        if(pR){
            //report progress!
            pR->setProgress(1.f/3.f + static_cast<float>(z)/static_cast<float>(dim.z) / 3.f);
        }
    }

    currentSource = currentDestination;
    currentDestination = output;
    LowerEnvelope.reserve(dim.z);

    for (size_t y = 0; y< dim.y; ++y){
        for (size_t x = 0; x < dim.x; ++x){
            LowerEnvelope.clear();
            LowerEnvelope.push_back(tgt::svec2(0,0));

            //Forward Scan
            for (size_t z = 1; z < dim.z ; ++z){
                while (!LowerEnvelope.empty() &&
                    DistanceFunction::EvalChebyFunction(static_cast<float>(LowerEnvelope.back().y),static_cast<float>(LowerEnvelope.back().x),currentSource->voxel(x,y,LowerEnvelope.back().x)) >
                    DistanceFunction::EvalChebyFunction(static_cast<float>(LowerEnvelope.back().y),static_cast<float>(z),currentSource->voxel(x,y,z))){
                    LowerEnvelope.pop_back();
                }

                if(LowerEnvelope.empty()){
                    LowerEnvelope.push_back(tgt::svec2(z,0));
                }
                else{
                    w = 1 + static_cast<size_t>(
                        DistanceFunction::PointOfChebyFunctionIntersection(static_cast<float>(LowerEnvelope.back().x),
                        static_cast<float>(z),
                        currentSource->voxel(x,y,LowerEnvelope.back().x),
                        currentSource->voxel(x,y,z)));

                    if (w < dim.y){
                        LowerEnvelope.push_back(tgt::svec2(z,w));
                    }
                }
            }

            //Backward Scan
            for(size_t z = dim.z-1; z != SIZE_MAX; --z){
                currentDestination->voxel(x,y,z) = DistanceFunction::EvalChebyFunction(static_cast<float>(z),
                                                                static_cast<float>(LowerEnvelope.back().x),
                                                                currentSource->voxel(x,y,LowerEnvelope.back().x));
                if(z == LowerEnvelope.back().y){
                    LowerEnvelope.pop_back();
                }
            }
        }

        if(pR){
            //report progress!
            pR->setProgress(2.f/3.f + static_cast<float>(y)/static_cast<float>(dim.y) / 3.f);
        }
    }

    delete tmp;

    if (pR)
        pR->setProgress(1.f);

    return new Volume(currentDestination, vb);
}


typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorSquaredEuclideanDistanceTransformBase> VolumeOperatorSquaredEuclideanDistanceTransform;
typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorEuclideanDistanceTransformBase> VolumeOperatorEuclideanDistanceTransform;
typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorManhattanDistanceTransformBase> VolumeOperatorManhattanDistanceTransform;
typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorChebychevDistanceTransformBase> VolumeOperatorChebychevDistanceTransform;

} // namespace


#endif
