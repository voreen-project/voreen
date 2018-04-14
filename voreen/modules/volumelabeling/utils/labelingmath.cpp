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

#include "labelingmath.h"
#include "tgt/logmanager.h"

// for least-squares-fitting
#include "jama/jama_qr.h"

using namespace labeling;
using namespace tgt;
using namespace std;

inline float faculty(int a) {
    if (a > 1)
        return a*faculty(a-1);
    else
        return 1.f;
}

// functions for evaluating bernstein polynomials
// formula taken from: Angel, Edward. Interactive Computer Graphics, 2006

// evaluates a bernstein polynomial
// @param u parameter in range [0;1]
// @param comp component of vector (starting with 0)
// @param deg degree of polynomial
float bernstein(float u, int comp, int deg) {
    float factor = faculty(deg) / ( faculty(comp)*faculty(deg - comp) );
    return factor*pow(u, static_cast<float>(comp))*pow(1-u, static_cast<float>(deg-comp));
}

// evaluates a bernstein polynomial's first derivative
// @param u parameter in range [0;1]
// @param comp component of vector (starting with 0)
// @param deg degree of polynomial
float bernsteinDerivative(float u, int comp, int deg) {
    float factor = faculty(deg) / ( faculty(comp)*faculty(deg - comp) );
    if (comp == 0)
        return -factor*deg*pow(1-u, static_cast<float>(deg-1));
    else if (deg == comp)
        return factor*comp*pow(u,static_cast<float>(comp-1));
    else
        return     factor*( comp*pow(u,static_cast<float>(comp-1))*pow(1-u, static_cast<float>(deg-comp))
        - pow(u, static_cast<float>(comp))*(deg-comp)*pow(1-u, static_cast<float>(deg-comp-1)) );
}

//----------------------------------------------------------------------------

// several functions for evaluating polynomials

float evalPolynomial(float t, void* params) {
    int degree = (static_cast<Polynomial*>(params))->degree;
    float* coeff = (static_cast<Polynomial*>(params))->coeff;
    float result = 0;
    for (int i=0; i<=degree; i++)
        result += coeff[i]*pow(t,i);
    return result;
}

float evalPolynomialDerivative(float t, void* params) {
    int degree = (static_cast<Polynomial*>(params))->degree;
    float* coeff = (static_cast<Polynomial*>(params))->coeff;
    float result = 0;
    for (int i=0; i<degree; i++)
        result += coeff[i+1]*(i+1)*pow(t,i);
    return result;
}

float evalPolynomialCurve2DTangentMagnitude(float t, void* params) {
    Polynomial xpoly = (static_cast<Polynomial*>(params))[0];
    Polynomial ypoly = (static_cast<Polynomial*>(params))[1];
    return sqrt( pow(evalPolynomialDerivative(t, &xpoly), 2) +
        pow(evalPolynomialDerivative(t, &ypoly), 2) );
}

float evalPolynomialCurve3DTangentMagnitude(float t, void* params) {
    Polynomial xpoly = (static_cast<Polynomial*>(params))[0];
    Polynomial ypoly = (static_cast<Polynomial*>(params))[1];
    Polynomial zpoly = (static_cast<Polynomial*>(params))[2];
    return sqrt( pow(evalPolynomialDerivative(t, &xpoly), 2) +
        pow(evalPolynomialDerivative(t, &ypoly), 2) +
        pow(evalPolynomialDerivative(t, &zpoly), 2) );
}

//----------------------------------------------------------------------------

Curve3D::Curve3D() {
}

bool Curve3D::setCtrlPoints(vector<tgt::vec3> &ctrlPoints, float curveLength) {
    ctrlPoints_ = ctrlPoints;
    return calcFunction(curveLength);
}

//----------------------------------------------------------------------------

Curve2D::Curve2D() {
}

bool Curve2D::setCtrlPoints(vector<tgt::vec2> &ctrlPoints, float curveLength) {
    ctrlPoints_ = ctrlPoints;
    return calcFunction(curveLength);
}

//----------------------------------------------------------------------------

Curve3DPolynomial::Curve3DPolynomial(int degree) {
    degree_ = degree;
    numCoeff_ = degree+1;
    xfunc_.degree = degree;
    xfunc_.coeff = 0;
    yfunc_.degree = degree;
    yfunc_.coeff = 0;
    zfunc_.degree = degree;
    zfunc_.coeff = 0;
}

bool Curve3DPolynomial::calcFunction(float curveLength) {
    if (ctrlPoints_.size() < (size_t)numCoeff_)
        return false;

    if (!xfunc_.coeff) {
        xfunc_.coeff = new float[numCoeff_];
        yfunc_.coeff = new float[numCoeff_];
        zfunc_.coeff = new float[numCoeff_];
    }

    /*
    Calculate a least-squares-fit of the polynomial function to the control points.
    For each coordinate (x,y,z) there is one polynomial function. The least-squares-fit
    is performed by solving the linear equation system X*c=y,
    where X(i,j)=i^j and y holds the control points coordinates.
    */

    TNT::Array1D<float> y(static_cast<int>(ctrlPoints_.size()));
    TNT::Array2D<float> X(static_cast<int>(ctrlPoints_.size()), static_cast<int>(numCoeff_));
    TNT::Array1D<float> c(numCoeff_);

    // Create matrix X and the corresponding QR-matrix
    for (int i=0; i<static_cast<int>(ctrlPoints_.size()); ++i) {
        for (int j=0; j<numCoeff_; ++j) {
            X[i][j] = pow(static_cast<float>(i), static_cast<float>(j));
        }
    }
    JAMA::QR<float> QR_Matrix(X);

    // x(t)
    for (int i=0; i<static_cast<int>(ctrlPoints_.size()); ++i)
        y[i] = ctrlPoints_[i].x;

    c = QR_Matrix.solve(y);
    for (int i=0; i<numCoeff_; ++i)
        xfunc_.coeff[i] = c[i];

    // y(t)
    for (int i=0; i<static_cast<int>(ctrlPoints_.size()); ++i)
        y[i] = ctrlPoints_[i].y;

    c = QR_Matrix.solve(y);
    for (int i=0; i<numCoeff_; ++i)
        yfunc_.coeff[i] = c[i];

    // z(t)
    for (int i=0; i<static_cast<int>(ctrlPoints_.size()); ++i)
        y[i] = ctrlPoints_[i].z;
    c = QR_Matrix.solve(y);
    for (int i=0; i<numCoeff_; ++i)
        zfunc_.coeff[i] = c[i];

    // map the curve segment's parameter range to [0;1]
    paramScale_ = static_cast<float>(ctrlPoints_.size()-1);
    paramShift_ = 0.f;
    if (curveLength > 0.f) {
        float length = getSegmentLength(0.f, 1.f);
        paramScale_ = (curveLength / length) * (ctrlPoints_.size()-1);
        paramShift_ = -(paramScale_-(ctrlPoints_.size()-1)) / 2.f;
    }

    return true;
}

tgt::vec3 Curve3DPolynomial::getCurvePoint(float t) {
    float p = t*paramScale_+paramShift_;
    tgt::vec3 result;
    result.x = evalPolynomial(p, &xfunc_);
    result.y = evalPolynomial(p, &yfunc_);
    result.z = evalPolynomial(p, &zfunc_);
    return result;
}

tgt::vec3 Curve3DPolynomial::getTangent(float t) {
    float p = t*paramScale_+paramShift_;
    tgt::vec3 result = vec3(0.f);
    result.x = evalPolynomialDerivative(p, &xfunc_);
    result.y = evalPolynomialDerivative(p, &yfunc_);
    result.z = evalPolynomialDerivative(p, &zfunc_);
    return result*paramScale_;
}

float Curve3DPolynomial::getTangentMagnitude(float t) {
    return length(getTangent(t));
}


float Curve3DPolynomial::getSegmentLength(float t1, float t2) {
    float avgTangentMagnitude = ( getTangentMagnitude(t1) +
        getTangentMagnitude((t1+t2)/2.f) +
        getTangentMagnitude(t2) ) / 3.f;
    return fabs(t2-t1)*avgTangentMagnitude;
}

tgt::vec3 Curve3DPolynomial::getNextPoint(float &curParam, float offset) {
    float pixelsPerParam = getTangentMagnitude(curParam);
    curParam += offset / pixelsPerParam;
    return getCurvePoint(curParam);
}

void Curve3DPolynomial::shift(vec3 shiftVector) {
    xfunc_.coeff[0] += shiftVector.x;
    yfunc_.coeff[0] += shiftVector.y;
    zfunc_.coeff[0] += shiftVector.z;

    for (size_t i=0; i<ctrlPoints_.size(); ++i)
        ctrlPoints_[i] = ctrlPoints_[i] + shiftVector;
}


//----------------------------------------------------------------------------

Curve2DPolynomial::Curve2DPolynomial(int degree) {
    degree_ = degree;
    numCoeff_ = degree+1;
    xfunc_.degree = degree;
    xfunc_.coeff = 0;
    yfunc_.degree = degree;
    yfunc_.coeff = 0;
}

bool Curve2DPolynomial::calcFunction(float curveLength) {
    if (ctrlPoints_.size() < (size_t)numCoeff_)
        return false;

    if (!xfunc_.coeff){
        xfunc_.coeff = new float[numCoeff_];
        yfunc_.coeff = new float[numCoeff_];
    }

    /*
    Calculate a least-squares-fit of the polynomial function to the control points.
    For each coordinate (x,y,z) there is one polynomial function. The least-squares-fit
    is performed by solving the linear equation system X*c=y,
    where X(i,j)=i^j and y holds the control points coordinates.
    */

    TNT::Array1D<float> y(static_cast<int>(ctrlPoints_.size()));
    TNT::Array2D<float> X(static_cast<int>(ctrlPoints_.size()), static_cast<int>(numCoeff_));
    TNT::Array1D<float> c(numCoeff_);

    // Create matrix X and the corresponding QR-matrix
    for (int i=0; i<static_cast<int>(ctrlPoints_.size()); i++)
        for (int j=0; j<numCoeff_; j++)
            X[i][j] = pow((float)i, (float)j);

    JAMA::QR<float> QR_Matrix(X);

    // x(t)
    for (int i=0; i<static_cast<int>(ctrlPoints_.size()); ++i)
        y[i] = ctrlPoints_[i].x;
    c = QR_Matrix.solve(y);
    for (int i=0; i<numCoeff_; ++i)
        xfunc_.coeff[i] = c[i];

    // y(t)
    for (int i=0; i<static_cast<int>(ctrlPoints_.size()); ++i)
        y[i] = ctrlPoints_[i].y;
    c = QR_Matrix.solve(y);
    for (int i=0; i<numCoeff_; ++i)
        yfunc_.coeff[i] = c[i];

    // map the curve segment's parameter range to [0;1]
    paramScale_ = static_cast<float>(ctrlPoints_.size()-1);
    paramShift_ = 0.f;
    if (curveLength > 0.f) {
        float length = getSegmentLength(0.f, 1.f);
        paramScale_ = (curveLength / length) * (ctrlPoints_.size()-1);
        paramShift_ = -(paramScale_-(ctrlPoints_.size()-1)) / 2.f;
    }
    return true;
}

tgt::vec2 Curve2DPolynomial::getCurvePoint(float t) {
    float p = t*paramScale_+paramShift_;
    tgt::vec2 result;
    result.x = evalPolynomial(p, &xfunc_);
    result.y = evalPolynomial(p, &yfunc_);
    return result;
}

tgt::vec2 Curve2DPolynomial::getTangent(float t) {
    float p = t*paramScale_+paramShift_;
    tgt::vec2 result = vec2(0.f);
    result.x = evalPolynomialDerivative(p, &xfunc_);
    result.y = evalPolynomialDerivative(p, &yfunc_);
    return result*paramScale_;
}

float Curve2DPolynomial::getTangentMagnitude(float t) {
    return length(getTangent(t));
}


float Curve2DPolynomial::getSegmentLength(float t1, float t2) {
    float avgTangentMagnitude = ( getTangentMagnitude(t1) +
        getTangentMagnitude((t1+t2)/2.f) +
        getTangentMagnitude(t2) ) / 3.f;
    return fabs(t2-t1)*avgTangentMagnitude;
}

tgt::vec2 Curve2DPolynomial::getNextPoint(float &curParam, float offset) {
    float pixelsPerParam = getTangentMagnitude(curParam);
    curParam += offset / pixelsPerParam;
    return getCurvePoint(curParam);
}

void Curve2DPolynomial::shift(vec2 shiftVector) {
    xfunc_.coeff[0] += shiftVector.x;
    yfunc_.coeff[0] += shiftVector.y;

    for (size_t i=0; i<ctrlPoints_.size(); i++)
        ctrlPoints_[i] = ctrlPoints_[i] + shiftVector;
}

//----------------------------------------------------------------------------

BezierPatch::BezierPatch(bool useDisplayList) {
    degreeS_ = -1;
    degreeT_ = -1;
    ctrlPoints_ = NULL;
    displayList_ = 0;
    useDisplayList_ = useDisplayList;
}

BezierPatch::~BezierPatch() {
    if (displayList_ > 0) {
        glDeleteLists(displayList_, 1);
        displayList_ = 0;
    }
}

int BezierPatch::getDegreeS(){
    return degreeS_;
}

int BezierPatch::getDegreeT(){
    return degreeT_;
}

void BezierPatch::setCtrlPoints(vec3* ctrlPoints, int degreeS, int degreeT) {
    delete ctrlPoints_;

    ctrlPoints_ = ctrlPoints;
    degreeS_ = degreeS;
    degreeT_ = degreeT;
    if (displayList_ > 0) {
        glDeleteLists(displayList_, 1);
        displayList_ = 0;
    }
}

vec3* BezierPatch::getCtrlPoints(int &degreeS, int &degreeT) {
    degreeS = degreeS_;
    degreeT = degreeT_;
    return ctrlPoints_;
}

vec3 BezierPatch::getPoint(float s, float t) {
    vec3 result = vec3(0.f);
    for (int i=0; i<=degreeS_; ++i) {
        for (int j=0; j<=degreeT_; ++j) {
            result += bernstein(s,i,degreeS_)*bernstein(t,j,degreeT_)*
                ctrlPoints_[j*(degreeS_+1)+i];
        }
    }
    return result;
}

vec3 BezierPatch::getTangentS(float s, float t) {
    vec3 result = vec3(0.f);
    for (int i=0; i<=degreeS_; ++i){
        for (int j=0; j<=degreeT_; ++j) {
            result += bernsteinDerivative(s,i,degreeS_)*bernstein(t,j,degreeT_)*
                ctrlPoints_[j*(degreeS_+1)+i];
        }
    }
    return result;
}

vec3 BezierPatch::getTangentT(float s, float t) {
    vec3 result = vec3(0.f);
    for (int i=0; i<=degreeS_; ++i) {
        for (int j=0; j<=degreeT_; ++j) {
            result += bernstein(s,i,degreeS_)*bernsteinDerivative(t,j,degreeT_)*
                ctrlPoints_[j*(degreeS_+1)+i];
        }
    }
    return result;
}

vec3 BezierPatch::getNormal(float s, float t){
    vec3 tangentS = getTangentS(s,t);
    vec3 tangentT = getTangentT(s,t);
    vec3 result = normalize(cross(tangentS,tangentT));

    return result;
}

void BezierPatch::render(int s_steps, int t_steps, bool genTexCoords, GLuint texUnit) {
    bool creatingDisplayList = false;
    if (displayList_ == 0 && useDisplayList_) {
        displayList_ = glGenLists(1);
        glNewList(displayList_, GL_COMPILE);
        creatingDisplayList = true;
    }

    if (displayList_ == 0 || creatingDisplayList) {
        float s_offset = 1/(s_steps+1.f);
        float t_offset = 1/(t_steps+1.f);
        for (float s=0.f; s<0.99; s += s_offset ) {
            vec3 last0 = getPoint(s,0.f);
            vec3 last1 = getPoint(s+s_offset,0.f);
            glBegin(GL_QUAD_STRIP);
            if (genTexCoords)
                glMultiTexCoord2f( texUnit, s, 0.f );
            glVertex3fv(last0.elem);
            if (genTexCoords)
                glMultiTexCoord2f( texUnit, s+s_offset, 0.f );
            glVertex3fv(last1.elem);
            for (float t=t_offset; t<1.01f; t += t_offset) {
                vec3 cur0 = getPoint(s, t);
                vec3 cur1 = getPoint(s+s_offset, t);
                vec3 normal0 = getNormal(s, t);
                vec3 normal1 = getNormal(s+s_offset, t);
                if (genTexCoords)
                    glMultiTexCoord2f(texUnit, s,t);
                glNormal3fv(normal0.elem);
                glVertex3fv(cur0.elem);
                if (genTexCoords)
                    glMultiTexCoord2f(texUnit, s+s_offset,t);
                glNormal3fv(normal1.elem);
                glVertex3fv(cur1.elem);

                last0 = cur0;
                last1 = cur1;
            }
            glEnd();
        }
    }

    if (creatingDisplayList)
        glEndList();

    if (displayList_ != 0)
        glCallList(displayList_);

    // tried to use opengl evaluators for rendering the patch. problem is
    // that I did not manage to generate texture coords for a texture unit
    // different from the first.

    /*
    GLfloat* ctrlPoints = new GLfloat[(degreeS_+1)*(degreeT_+1)*3];
    for (int i=0; i<(degreeS_+1)*(degreeT_+1); i++){
    vec3 point = ctrlPoints_[i];
    ctrlPoints[3*i] = point.x;
    ctrlPoints[3*i+1] = point.y;
    ctrlPoints[3*i+2] = point.z;
    }
    glMap2f(GL_MAP2_VERTEX_3, 0, 1, 3, degreeS_+1, 0, 1,
        (degreeS_+1)*3, degreeT_+1, ctrlPoints);

    GLfloat texpts[2][2][2] = {{{0.0, 0.0}, {0.0, 1.0}},
    {{1.0, 0.0}, {1.0, 1.0}}};
    glMap2f(GL_MAP2_TEXTURE_COORD_2, 0, 1, 2, 2,
    0, 1, 4, 2, &texpts[0][0][0]);

    glEnable(GL_MAP2_VERTEX_3);
    glEnable(GL_MAP2_TEXTURE_COORD_2);

    glMapGrid2f(10, 0, 1, 10, 0, 1);
    glEvalMesh2(GL_FILL, 0, 10, 0, 10);
    glDisable(GL_MAP2_VERTEX_3);
    glDisable(GL_MAP2_TEXTURE_COORD_2);

    delete ctrlPoints; */
}

//----------------------------------------------------------------------------
