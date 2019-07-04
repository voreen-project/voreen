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

#include "foldVessel.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include <boost/math/constants/constants.hpp>

using tgt::vec3;

namespace voreen {

#define	TwoPi  6.28318530717958648
const double eps=1e-14;

const std::string foldVessel::loggerCat_("voreen.vesselnetworkanalysis.foldVessel");
foldVessel::foldVessel()
    : VolumeProcessor()

    , inport_(Port::INPORT, "volumethinning.inport", "Volume Input")
    , outport_(Port::OUTPORT, "volumethinning.outport", "Volume Output",Processor::VALID)
    ,enableProcessing_("enableProcessing","enable processing",true)
    ,bezierCurve_("createBezierCurve","create Bezier curve",false)
    ,threshold_("threshold","min threshold value for folding",0.0001,0,1.0)
    ,bezierX0_("bezX0","P0:X",0.5,0.1,1.0)
    ,bezierX1_("bezX1","P1:X",0.25,0,1.0)
    ,bezierY1_("bezY1","P1:Y",0.5,0,1.0)
    ,bezierX2_("bezX2","P2:X",0.5,0.1,1.0)



{
    addPort(outport_);
    addPort(inport_);

    addProperty(enableProcessing_);
    addProperty(bezierCurve_);
    addProperty(threshold_);
    addProperty(bezierX0_);
    addProperty(bezierX1_);
    addProperty(bezierY1_);
    addProperty(bezierX2_);


}

foldVessel::~foldVessel() {
}



static double _root3 ( double x )
{
    double s = 1.;
    while ( x < 1. )
    {
        x *= 8.;
        s *= 0.5;
    }
    while ( x > 8. )
    {
        x *= 0.125;
        s *= 2.;
    }
    double r = 1.5;
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    return r * s;
}

double root3 ( double x )
{
    if ( x > 0 ) return _root3 ( x ); else
    if ( x < 0 ) return-_root3 (-x ); else
    return 0.;
}



// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] � i*x[2], return 1
int SolveP3(double *x,double a,double b,double c) {	// solve cubic equation x^3 + a*x^2 + b*x + c = 0
        double a2 = a*a;
    double q  = (a2 - 3*b)/9;
        double r  = (a*(2*a2-9*b) + 27*c)/54;
        // equation x^3 + q*x + r = 0
    double r2 = r*r;
        double q3 = q*q*q;
        double A,B;
        if (r2 <= (q3 + eps)) {//<<-- FIXED!
                double t=r/sqrt(q3);
                if( t<-1) t=-1;
                if( t> 1) t= 1;
        t=acos(t);
        a/=3; q=-2*sqrt(q);
        x[0]=q*cos(t/3)-a;
        x[1]=q*cos((t+TwoPi)/3)-a;
        x[2]=q*cos((t-TwoPi)/3)-a;
        return(3);
    } else {
        //A =-pow(fabs(r)+sqrt(r2-q3),1./3);
        A =-root3(fabs(r)+sqrt(r2-q3));
                if( r<0 ) A=-A;
                B = A==0? 0 : B=q/A;

                a/=3;
                x[0] =(A+B)-a;
        x[1] =-0.5*(A+B)-a;
        x[2] = 0.5*sqrt(3.)*(A-B);
                if(fabs(x[2])<eps) { x[2]=x[1]; return(2); }
        return(1);
    }
}// SolveP3(double *x,double a,double b,double c) {

void findT(double *x,tgt::vec2 m,tgt::vec2 p0, tgt::vec2 p1, tgt::vec2 p2){

    tgt::vec2 A,B,pos;
    double a,b,c,d;

    A.x = p1.x - p0.x;
    A.y = p1.y - p0.y;
    B.x = p0.x - 2 * p1.x + p2.x;
    B.y = p0.y - 2 * p1.y + p2.y;

    pos.x = p0.x - m.x;
    pos.y = p0.y - m.y;

    a = tgt::dot(B,B);
    b = 3 * tgt::dot(A,B) / a;
    c = (2 * tgt::dot(A,A) + tgt::dot(pos,B) ) / a;
    d = tgt::dot(pos,A) / a;

    SolveP3(x,b,c,d);


}

void getPos(tgt::vec2 &pos, tgt::vec2 p0,tgt::vec2 p1,tgt::vec2 p2 ,double t){
        double a = (1 - t) * (1 - t);
        double b = 2 * t * (1 - t);
        double c = t * t;

        pos.x = a * p0.x + b * p1.x + c * p2.x;
        pos.y = a * p0.y + b * p1.y + c * p2.y;
}

/*
 * map values from input(radius,angle,z) to target(x,y,z)
 * sample from target to input with calculation of polar coordinates
 */
void foldVessel::mapToQuader(const VolumeRAM* input, VolumeRAM_UInt8* target){


    //convert angle coordinates from 360 degree
    double factor_roh = (double)input->getDimensions().y / 360.0;

    //center position
    tgt::vec2 center(target->getDimensions().x/2,target->getDimensions().y/2);

    for (size_t voxel_z=0; voxel_z<(target->getDimensions().z); voxel_z++) {
        for (size_t voxel_y=0; voxel_y<(target->getDimensions().y); voxel_y++) {
            for (size_t voxel_x=0; voxel_x<(target->getDimensions().x); voxel_x++) {

                tgt::vec3 point(voxel_x,voxel_y,voxel_z);

                //2D point -> z is unimportant
                tgt::vec2 xyPoint(point.x,point.y);
                //get point from center position
                xyPoint = xyPoint - center;


                double radius = 2*std::sqrt(xyPoint.x*xyPoint.x + xyPoint.y*xyPoint.y); //2 as factor relative size of cylinder
                double radian = std::atan2(xyPoint.x ,xyPoint.y );
                if(radian <= 0)
                    radian+=2*boost::math::constants::pi<double>();
                double roh = radian *(180.0/boost::math::constants::pi<double>());

                //convert angle coordinates from 360 degree
                roh *= factor_roh;

                //if angle is 360 its equal to getDimension and neeed to be decreased by one
                if(roh >= target->getDimensions().y)
                    roh--;

                tgt::vec3 newPoint(radius,roh,voxel_z);

                float value = input->getVoxelNormalizedLinear(newPoint);

                //check if values are over the given convert angle coordinates from 360 degreethreshold
                if(value > threshold_.get()/100){
                    /*
                    point.x *= 0.5;
                    point.x += target->getDimensions().x/4;
                    point.y *= 0.5;
                    point.y += target->getDimensions().y/4;
                    point.z *= 0.5;
                    point.z += target->getDimensions().z/4;
                    */
                    target->setVoxelNormalized(value,point);
                }



            }
        }
    }

}



int getDist(int x1,int x2 ){
    return std::sqrt((x1 - x2) * (x1 - x2));
}

bool checkInterval(double t){
    //set t to 0 or in 1 in corner cases
    if(t<1 && t > -0.01)
        t=0;
    if(t>0 && t < 1.01)
        t=1;
    return t>=0 && t<=1;
}


void findBezier(tgt::vec2 &finalPos,double x[3], const VolumeRAM* input, tgt::vec2 m,float bezierX0,float bezierX1, float bezierY1,float bezierX2){
    tgt::vec3 dimensions = input->getDimensions();
    tgt::vec2 p0(dimensions.x*bezierX0,0);
    tgt::vec2 p1(dimensions.x*bezierX1,dimensions.y*bezierY1);
    tgt::vec2 p2(dimensions.x*bezierX2,dimensions.y-1);

    tgt::vec2 pos1(0,0);
    tgt::vec2 pos2(0,0);
    tgt::vec2 pos3(0,0);
    int dist1 = std::numeric_limits<int>::max();
    int dist2 = std::numeric_limits<int>::max();
    int dist3 = std::numeric_limits<int>::max();

    findT(x,m,p0,p1,p2);
    getPos(pos1,p0,p1,p2,x[0]);
    getPos(pos2,p0,p1,p2,x[1]);
    getPos(pos3,p0,p1,p2,x[2]);

    if(checkInterval(x[0]))
        dist1 = getDist(m.x,pos1.x);
    if(checkInterval(x[1]))
        dist2 = getDist(m.x,pos2.x);
    if(checkInterval(x[2]))
        dist3 = getDist(m.x,pos3.x);

    if(dist1 <= dist2){
        finalPos=pos1;
        if(dist1<dist3){}
        else
            finalPos=pos3;
    }else{
        finalPos=pos2;
        if(dist2<=dist3){}
        else{
            finalPos=pos3;
        }
    }

}

void createLine(VolumeRAM_UInt8* target, std::vector<tgt::vec3>& line){
    int x_pos = target->getDimensions().x/2;
    int y_pos = target->getDimensions().y/2;

    for(size_t zi = 0; zi < target->getDimensions().z; zi++){
        tgt::ivec3 middle(x_pos,y_pos,zi);
        line.push_back(middle);
    }

}

void foldVessel::process() {

    const VolumeBase* invol = inport_.getData();
    const VolumeRAM* input = invol->getRepresentation<VolumeRAM>();

    tgt::vec3 dimensions = input->getDimensions();

    VolumeRAM_UInt8* target = new VolumeRAM_UInt8(dimensions);
    VolumeRAM_UInt8* targetBezier = new VolumeRAM_UInt8(dimensions);

    if (!enableProcessing_.get()) {
        outport_.setData(invol, false);
        return;
    }

    VolumeRAM* outputVolume = target;
    outputVolume->clear();
    targetBezier->clear();

    //create centered cylinder
    mapToQuader(input,target);

    //create Bezier curve
    if(bezierCurve_.get()){
        const int maxDiff=1; //max allowed difference between 2 iterations
        std::vector<tgt::vec3> centerline;
        createLine(target,centerline); //creates centerline interpretation
        tgt::vec3 m_3(0,0,0); //current point 3d
        tgt::vec2 m(0,0); //curent point in 2d
        tgt::vec2 finalPos(0.f,0.f); //position on bezier curve
        double x[3]; //solution for 3 degree solution
        int lastDistance=0;  //last saved distance
        int newDistance=0; //current distance from centerline to bezier
        float value=0; //normalized voxel value
        //iterate from Z:0-getDimension.z & get distance from centerline to bezier curve
        for(auto &p : centerline){
            m.x=p.x;
            m.y=p.z;
            findBezier(finalPos,x,input,m,bezierX0_.get(),bezierX1_.get(),bezierY1_.get(),bezierX2_.get());
            m_3.x = finalPos.x;
            m_3.y = p.y;
            m_3.z = p.z;

            newDistance = finalPos.x-m.x;

            //make sure that there are no jumps between the iteration (max maxDiff diferrence between 2 iterations)
            if(p.z>0){
                if(newDistance>0 && lastDistance>0){
                    if((lastDistance+maxDiff) < newDistance)
                        newDistance=lastDistance+maxDiff;
                    else if((newDistance+maxDiff) < lastDistance){
                        newDistance=lastDistance-maxDiff;
                    }
                }else if(newDistance<=0 && lastDistance>0 ){
                    newDistance=lastDistance-maxDiff;
                }else if(newDistance>0 && lastDistance<=0){
                    newDistance=lastDistance+maxDiff;
                }else if(newDistance<=0 && lastDistance<=0){
                    if((lastDistance+maxDiff) > newDistance)
                        newDistance=lastDistance-maxDiff;
                    else if((lastDistance+maxDiff) < newDistance){
                        newDistance=lastDistance+maxDiff;
                    }
                }
            }
            lastDistance = newDistance;

            //move each voxel in x direction with calculated distance
            for (size_t voxel_y=0; voxel_y<(target->getDimensions().y); voxel_y++) {
                for (size_t voxel_x=0; voxel_x<(target->getDimensions().x); voxel_x++) {
                    value = target->getVoxelNormalized(voxel_x,voxel_y,p.z);
                    int newX = voxel_x + newDistance;
                    if(newX < 0){
                        // ignore x position smaller than 0
                    }
                    else if((static_cast<unsigned int>(newX)) < target->getDimensions().x){
                            targetBezier->setVoxelNormalized(value,newX,voxel_y,p.z);
                    }else{
                        //ignore values greater than max x
                        voxel_x=target->getDimensions().x;
                    }
                }
            }
        }

    }

    if(bezierCurve_.get()){
        Volume* vh = new Volume(targetBezier, invol->getSpacing(), invol->getOffset() );
        outport_.setData(vh);
    }else{
        Volume* vh = new Volume(target, invol->getSpacing(), invol->getOffset() );
        outport_.setData(vh);
    }
}
VoreenSerializableObject* foldVessel::create() const {
    return new foldVessel();
}

} // namespace voreen
