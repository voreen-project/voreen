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

#include "plotlibraryfilebase.h"

namespace voreen {

const std::string PlotLibraryFileBase::loggerCat_("voreen.plotting.PlotLibraryFileBase");

PlotLibraryFileBase::PlotLibraryFileBase()
    : PlotLibrary()
    , lineColor_(0.f,0.f,0.f,1.f)
    , projectionsMatrix_(tgt::Matrix4<plot_t>::createIdentity())
    , modelviewMatrix_(tgt::Matrix4<plot_t>::createIdentity())
    , invertedmodelviewMatrix_(tgt::Matrix4<plot_t>::zero)
    , left_(0)
    , right_(0)
    , bottom_(0)
    , top_(0)
    , near_(0)
    , far_(0)
    , clippingPlanes_()
{}


void PlotLibraryFileBase::setOutputFile(const std::string& file) {
    outputFile_ = file;
}

void PlotLibraryFileBase::setWireColor(tgt::Color color) {
    lineColor_ = color;
}

tgt::dvec2 PlotLibraryFileBase::convertPlotCoordinatesToViewport(const tgt::dvec3& plotCoordinates) const {
    tgt::dvec3 copy = logScale3dtoLogCoordinates(plotCoordinates);

    double x, y;
    tgt::Vector3<plot_t> point = projectionsMatrix_ * modelviewMatrix_ * copy;
    x = windowSize_.x * (point.x + 1)/2.0;
    y = windowSize_.y * (point.y + 1)/2.0;
    return tgt::dvec2(x, y);
}

tgt::dvec2 PlotLibraryFileBase::convertPlotCoordinatesToViewport(const plot_t x, const plot_t y, const plot_t z) const {
    return convertPlotCoordinatesToViewport(tgt::dvec3(x,y,z));
}

tgt::dvec3 PlotLibraryFileBase::convertPlotCoordinatesToViewport3(const tgt::dvec3& plotCoordinates) const {
    tgt::dvec3 copy = logScale3dtoLogCoordinates(plotCoordinates);

    double x, y, z;
    tgt::Vector3<plot_t> point = projectionsMatrix_ * modelviewMatrix_ * copy;
    x = windowSize_.x * (point.x + 1)/2.0;
    y = windowSize_.y * (point.y + 1)/2.0;
    z = (point.z + 1)/2.0;
    return tgt::dvec3(x, y, z);
}

tgt::dvec3 PlotLibraryFileBase::convertPlotCoordinatesToViewport3(const plot_t x, const plot_t y, const plot_t z) const {
    return convertPlotCoordinatesToViewport3(tgt::dvec3(x,y,z));
}

tgt::dvec2 PlotLibraryFileBase::convertViewportToPlotCoordinates(tgt::ivec2 viewCoord) const {
    tgt::Vector4<plot_t> coord = tgt::Vector4<plot_t>(2*static_cast<double>(viewCoord.x)/static_cast<double>(windowSize_.x)-1,2*static_cast<double>(viewCoord.y)/static_cast<double>(windowSize_.y)-1,-1,1);
    tgt::Vector4<plot_t> plot_point;
    tgt::Matrix4<plot_t> matrix1 = (projectionsMatrix_ * modelviewMatrix_);
    tgt::Matrix4<plot_t> invertedmatrix;
    matrix1.invert(invertedmatrix);
    plot_point = invertedmatrix * coord;
    //tgt::Matrix4<plot_t> matrix2 = matrix1 * invertedmatrix;
    double x = plot_point.x;
    double y = plot_point.y;
    if (logarithmicAxisFlags_[X_AXIS])
        x = convertFromLogCoordinates(plot_point.x, X_AXIS);
    if (logarithmicAxisFlags_[Y_AXIS])
        y = convertFromLogCoordinates(plot_point.y, Y_AXIS);
    return tgt::dvec2(x, y);
}


PlotLibraryFileBase::Projection_Coordinates PlotLibraryFileBase::convertPlotCoordinatesToViewport3Projection(const tgt::dvec3& plotCoordinates) const {
    tgt::dvec3 copy = logScale3dtoLogCoordinates(plotCoordinates);

    double x, y, z;

    tgt::Vector3<plot_t> point2 = modelviewMatrix_ * copy;
    tgt::Vector3<plot_t> point = projectionsMatrix_ * point2;
    x = windowSize_.x * (point.x + 1)/2.0;
    y = windowSize_.y * (point.y + 1)/2.0;
    z = (point.z + 1)/2.0;
    return Projection_Coordinates(copy,point2,tgt::Vector3<plot_t>(x,y,z));
}

PlotLibraryFileBase::Projection_Coordinates PlotLibraryFileBase::convertPlotCoordinatesToViewport3Projection(const plot_t x, const plot_t y, const plot_t z) const {
    return convertPlotCoordinatesToViewport3Projection(tgt::dvec3(x,y,z));
}

PlotLibraryFileBase::Projection_Coordinates PlotLibraryFileBase::projection_neu(const tgt::Vector3<plot_t>& point) {
    plot_t x,y,z;

    tgt::Vector3<plot_t> point0 = modelviewMatrix_ * point;
    tgt::Vector3<plot_t> point1 = projectionsMatrix_ * point0;
    x = windowSize_.x * (point1.x + 1)/2.0;
    y = windowSize_.y * (point1.y + 1)/2.0;
    z = (point1.z + 1)/2.0;
    return Projection_Coordinates(point,point0,tgt::Vector3<plot_t>(x, y, z));

}

PlotLibraryFileBase::Projection_Coordinates PlotLibraryFileBase::projection_neu(plot_t x, plot_t y, plot_t z) {
    return projection_neu(tgt::Vector3<plot_t>(x,y,z));
}


int PlotLibraryFileBase::intersect_Triangle_Ray(const std::vector< tgt::Vector3<plot_t> >& ray, const std::vector< tgt::Vector3<plot_t> >& triangle, tgt::Vector3<plot_t>* result_point) {
    if ((ray.size() != 2) || (triangle.size() < 3))
        return -1;
    tgt::Vector3<plot_t> u,v,n, ray_dir, w0,w;
    double r,a,b;
    u = triangle[1] - triangle[0];
    v = triangle[2] - triangle[0];
    n = tgt::cross(u,v);
    if (n == tgt::Vector3<plot_t>::zero)
        return -1;
    ray_dir = ray[1] - ray[0];
    w0 = ray[0] - triangle[0];
    a = -tgt::dot(n,w0);
    b = tgt::dot(n,ray_dir);
    if (std::abs(b) < 0.0000001) {
        if (a == 0)
            return 2;
        else
            return 0;
    }

    r = a / b;
    if ((r < 0.0) || (r > 1.0))
        return 0;
    *result_point = ray[0] + r*ray_dir;
    double uu, uv, vv, wu, wv, D;
    uu = tgt::dot(u,u);
    uv = tgt::dot(u,v);
    vv = tgt::dot(v,v);
    w = *result_point - triangle[0];
    wu = tgt::dot(w,u);
    wv = tgt::dot(w,v);
    D = uv*uv - uu*vv;

    double s,t;
    s = (uv*wv-vv*wu)/D;
    if ((s < 0.0) || (s > 1.0))
        return 0;
    t = (uv*wu-uu*wv)/D;
    if ((t < 0.0) || ((s+t)>1.0))
        return 0;
    return 1;
}

int PlotLibraryFileBase::intersect_Ray_Ray(const std::vector< tgt::Vector2<plot_t> >& ray1, const std::vector< tgt::Vector2<plot_t> >& ray2, tgt::Vector2<plot_t>* result_values) {
    tgt::Vector2<plot_t> p0,p1,d0,d1, u;
    p0 = ray1[0];
    p1 = ray2[0];
    d0 = ray1[1] - ray1[0];
    d1 = ray2[1] - ray2[0];
    u = p0-p1;
    if (det(d1,-d0) == 0)
        return 0;
    double s,t;
    t = det(u,-d0)/det(d1,-d0);
    s = det(d1,u)/det(d1,-d0);
    *result_values = tgt::Vector2<plot_t>(t,s);
    return 1;
}

PlotLibraryFileBase::Node* PlotLibraryFileBase::searchNode(PlotLibraryFileBase::Node* startNode, const tgt::Vector3<plot_t>& point1,const tgt::Vector3<plot_t>& point2, double alpha) {
    Node* runnode = startNode;
    Node* runnode1 = 0;
    Node* runnode2 = 0;
    while (1) {
        if (runnode->x_ == point1.x && runnode->y_ == point1.y && runnode->z_ == point1.z) {
            runnode1 = runnode->next_;
            while (1) {
                if (runnode1->intersect_) {
                    if (runnode1->alpha_ > alpha)
                        runnode2 = runnode1->prev_;
                    runnode1 = runnode1->next_;
                }
                else {
                    if (runnode1->x_ == point2.x && runnode1->y_ == point2.y && runnode1->z_ == point2.z) {
                        if (runnode2 == 0)
                            return runnode;
                        else
                            return runnode2;
                    }
                    else
                        break;
                }
            }
        }
        runnode = runnode->next_;
        if (runnode == startNode)
            return 0;
    }
}

void PlotLibraryFileBase::clippolygons(const std::vector< tgt::Vector3<plot_t> >& subject, const std::vector< tgt::Vector3<plot_t> >& clipp, std::vector< std::vector< tgt::Vector3<plot_t> > >* resultPolygons) {
    double a,b;
    int subject_second, clipp_second;
    Node* subject_node = 0;
    Node* runnode = 0;
    Node* runnode1 = 0;
    Node* clipp_node = 0;
    for (size_t si = 0; si < subject.size(); ++si) {
        if (runnode == 0) {
            runnode = new Node(subject[si].x,subject[si].y,subject[si].z);
            subject_node = runnode;
        }
        else {
            runnode1 = runnode;
            runnode = new Node(subject[si].x,subject[si].y,subject[si].z);
            runnode->prev_ = runnode1;
            runnode1->next_ = runnode;
        }
    }
    runnode->next_ = subject_node;
    subject_node->prev_ = runnode;
    for (size_t cj = 0; cj < clipp.size(); ++cj) {
        if (runnode == 0) {
            runnode = new Node(clipp[cj].x,clipp[cj].y,clipp[cj].z);
            clipp_node = runnode;
        }
        else {
            runnode1 = runnode;
            runnode = new Node(clipp[cj].x,clipp[cj].y,clipp[cj].z);
            runnode->prev_ = runnode1;
            runnode1->next_ = runnode;
        }
    }
    runnode->next_ = clipp_node;
    clipp_node->prev_ = runnode;
    for (size_t si = 0; si < subject.size(); ++si) {
        subject_second = static_cast<int>((si + 1) % subject.size());
        for (size_t cj = 0; cj < clipp.size(); ++cj) {
            clipp_second = static_cast<int>((cj + 1) % clipp.size());
            if (intersect(subject[si].xy(),subject[subject_second].xy(),clipp[cj].xy(),clipp[clipp_second].xy(),&a,&b)) {
                Node* i1 = createVertex(subject[si],subject[subject_second],a);
                i1->alpha_ = a;
                Node* i2 = createVertex(clipp[cj],clipp[clipp_second],b);
                i2->alpha_ = b;
                i1->neighbor_ = i2;
                i2->neighbor_ = i1;
                runnode1 = searchNode(subject_node,subject[si],subject[subject_second],a);
                if (runnode == 0)
                    return;
                i1->next_ = runnode1->next_;
                i1->prev_ = runnode1;
                runnode1->next_ = i1;
                i1->neighbor_->prev_ = i1;
                runnode1 = searchNode(clipp_node,clipp[cj],clipp[clipp_second],b);
                if (runnode == 0)
                    return;
                i2->next_ = runnode1->next_;
                i2->prev_ = runnode1;
                runnode1->next_ = i2;
                i2->neighbor_->prev_ = i2;
            }
        }
    }
    double w1 = windingnumber(subject,clipp[0]);
    double w2 = windingnumber(clipp,subject[0]);
    bool entry_status;
    if (w2 <=-1 || w2 >=1)
        entry_status = false;
    else
        entry_status = true;

    runnode = subject_node;
    int i = 0;
    while (i == 0 || runnode != subject_node) {
        ++i;
        if (runnode->intersect_) {
            runnode->entry_ = entry_status;
            if (entry_status)
                entry_status = false;
            else
                entry_status = true;
        }
        runnode = runnode->next_;
    }

    if (w1 <=-1 || w1 >=1)
        entry_status = false;
    else
        entry_status = true;
    runnode = clipp_node;
    i = 0;
    while (i == 0 || runnode != subject_node) {
        ++i;
        if (runnode->intersect_) {
            runnode->entry_ = entry_status;
            if (entry_status)
                entry_status = false;
            else
                entry_status = true;
        }
        runnode = runnode->next_;
    }
    resultPolygons->clear();
    bool unprocessed = true;
    while (unprocessed) {
        runnode = subject_node;
        while (runnode->intersect_ == false || runnode->visited_ == true) {
            runnode = runnode->next_;
        }
        if (runnode == subject_node) {
            unprocessed = false;
        }
        else {
            resultPolygons->resize(resultPolygons->size()+1);
            runnode1 = runnode;
            resultPolygons->at(resultPolygons->size()-1).push_back(runnode->getPoint());
            do {
                runnode->visited_ = true;
                if (runnode->entry_) {
                    do {
                        runnode = runnode->next_;
                        resultPolygons->at(resultPolygons->size()-1).push_back(runnode->getPoint());
                    } while (runnode->intersect_);
                }
                else {
                    do {
                        runnode = runnode->prev_;
                        resultPolygons->at(resultPolygons->size()-1).push_back(runnode->getPoint());
                    } while (runnode->intersect_);
                }
                runnode = runnode->neighbor_;
            } while (runnode1 != runnode);
        }
    }

    runnode = subject_node->next_;
    runnode1 = clipp_node->next_;
    while (subject_node != 0 || clipp_node != 0) {
        if (runnode->next_ == subject_node) {
            delete runnode;
            delete subject_node;
            runnode = 0;
            subject_node = 0;
        }
        else {
            runnode = runnode->next_;
            delete runnode->prev_;
        }
        if (runnode1->next_ == clipp_node) {
            delete runnode;
            delete subject_node;
            runnode = 0;
            subject_node = 0;
        }
        else {
            runnode1 = runnode1->next_;
            delete runnode1->prev_;
        }
    }
}

PlotLibraryFileBase::Node* PlotLibraryFileBase::createVertex(const tgt::Vector3d& startPoint, const tgt::Vector3d& endPoint, double alpha) {
    const tgt::Vector3d point = (endPoint - startPoint)*alpha;
    Node* schnitt = new Node(point.x,point.y,point.z);
    schnitt->intersect_ = true;
    return schnitt;
}

bool PlotLibraryFileBase::intersect(const tgt::Vector2<plot_t>& startPoint1, const tgt::Vector2<plot_t>& endPoint1, const tgt::Vector2<plot_t>& startPoint2, const tgt::Vector2<plot_t>& endPoint2, double* alpha1, double* alpha2) {
    double wec_p1 = tgt::dot(startPoint1-startPoint2,perpendicularVector(endPoint2-startPoint2));
    double wec_p2 = tgt::dot(endPoint1-startPoint2,perpendicularVector(endPoint2-startPoint2));
    if (wec_p1*wec_p2 <= 0) {
        double wec_q1 = tgt::dot(startPoint2-startPoint1,perpendicularVector(endPoint1-startPoint1));
        double wec_q2 = tgt::dot(endPoint2-startPoint1,perpendicularVector(endPoint1-startPoint1));
        if (wec_q1*wec_q2 <= 0) {
            *alpha1 = wec_p1/(wec_p1-wec_p2);
            *alpha2 = wec_q1/(wec_q1-wec_q2);
            return true;
        }
    }
    return false;
}

bool PlotLibraryFileBase::clippOnPlane(const tgt::Vector4<plot_t>& clippplane, const tgt::Vector4<plot_t>& eyepoint) {
    if (invertedmodelviewMatrix_ == tgt::Matrix4<plot_t>::zero)
        modelviewMatrix_.invert(invertedmodelviewMatrix_);
    tgt::Vector4<plot_t> eyeplane = clippplane*invertedmodelviewMatrix_;
    if (eyeplane.x*eyepoint.x + eyeplane.y*eyepoint.y + eyeplane.z*eyepoint.z + eyeplane.w*eyepoint.w >= 0)
        return true;
    else
        return false;
}

bool PlotLibraryFileBase::inClippRegion(const tgt::Vector3<plot_t>& point) {
    tgt::Vector4<plot_t> eyePoint;
    for (size_t i = 0; i < clippingPlanes_.size(); ++i) {
        eyePoint = tgt::Vector4<plot_t>(point,1);
        if (clippOnPlane(clippingPlanes_[i],eyePoint) == false)
            return false;
    }
    return true;
}

bool PlotLibraryFileBase::inClippRegion(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points) {
    tgt::Vector3<plot_t> modelpoint;
    for (size_t i = 0; i < points.size(); ++i) {
        modelpoint = points[i].afterModelview_;
        if (inClippRegion(modelpoint))
            return true;
    }
    return false;
}

double PlotLibraryFileBase::windingnumber(const std::vector< tgt::Vector3<plot_t> >& polygon, const tgt::Vector3<plot_t>& point) {
    double w = 0;
    std::vector< tgt::Vector3<plot_t> > newpolygon;
    for (size_t i = 0; i < polygon.size(); ++i) {
        newpolygon.push_back(polygon[i]-point);
    }
    int second;
    for (size_t i = 0; i < newpolygon.size(); ++i) {
        second = static_cast<int>((i + 1) % newpolygon.size());
        if (newpolygon[i].y*newpolygon[second].y < 0) {
            double r = newpolygon[i].x+ (newpolygon[i].y*(newpolygon[second].x-newpolygon[i].x))/(newpolygon[i].y-newpolygon[second].y);
            if (r > 0) {
                if (newpolygon[i].y < 0)
                    w += 1;
                else
                    w -= 1;
            }
        }
        else if (newpolygon[i].y == 0 && newpolygon[i].x > 0) {
            if (newpolygon[i].y < 0)
                w += 0.5;
            else
                w -= 0.5;
        }
        else if (newpolygon[second].y == 0 && newpolygon[second].x > 0) {
            if (newpolygon[second].y > 0)
                w += 0.5;
            else
                w -= 0.5;
        }
    }

    return w;
}

int tri_tri_intersect3D(const tgt::Vector3<plot_t>& c1, const tgt::Vector3<plot_t>& c2, const tgt::Vector3<plot_t>& c3,const tgt::Vector3<plot_t>& d1, const tgt::Vector3<plot_t>& d2, const tgt::Vector3<plot_t>& d3) {
    tgt::Vector3<plot_t> P1,P2,Q1,Q2;
    tgt::Vector3<plot_t> p11, p21, r,r4;
    tgt::Vector3<plot_t> t = tgt::Vector3<plot_t>::zero;

    double beta1, beta2, beta3;
    double gama1, gama2, gama3;
    double det1, det2, det3;
    double dp0, dp1, dp2;
    double dq1,dq2,dq3,dr, dr3;
    double alpha1, alpha2;
    bool alpha1_legal, alpha2_legal;
    double SF = 0.0;
    bool beta1_legal, beta2_legal;
    P1 = c2-c1;
    P2 = c3-c1;
    Q1 = d2-d1;
    Q2 = d3-d1;
    r = d1 - c1;
    // determinant computation
    dp0 = P1[1]*P2[2]-P2[1]*P1[2];
    dp1 = P1[0]*P2[2]-P2[0]*P1[2];
    dp2 = P1[0]*P2[1]-P2[0]*P1[1];
    dq1 = Q1[0]*dp0 - Q1[1]*dp1 + Q1[2]*dp2;
    dq2 = Q2[0]*dp0 - Q2[1]*dp1 + Q2[2]*dp2;
    dr  = -r[0]*dp0  + r[1]*dp1  - r[2]*dp2;



    beta1 = dr*dq2;  // beta1, beta2 are scaled so that beta_i=beta_i*dq1*dq2
    beta2 = dr*dq1;
    beta1_legal = (beta2>=0) && (beta2 <=dq1*dq1) && (dq1 != 0);
    beta2_legal = (beta1>=0) && (beta1 <=dq2*dq2) && (dq2 != 0);

    dq3=dq2-dq1;
    dr3=+dr-dq1;   // actually this is -dr3


    if ((dq1 == 0) && (dq2 == 0))
    {
        if (dr!=0) return 0;  // triangles are on parallel planes
        else
        {                        // triangles are on the same plane
            tgt::Vector3<plot_t> C2,C3,D2,D3, N1;
            // We use the coplanar test of Moller which takes the 6 vertices and 2 normals
            //as input.
            C2 = c1 + P1;
            C3 = c1 + P2;
            D2 = d1 + Q1;
            D3 = d1 + Q2;
            N1 = P1 + P2;
            return 0;//coplanar_tri_tri(N1,C1, C2,C3,D1,D2,D3);
        }
    }

    else if (!beta2_legal && !beta1_legal) return 0;// fast reject-all vertices are on
                                                    // the same side of the triangle plane

    else if (beta2_legal && beta1_legal)    //beta1, beta2
    {
        SF = dq1*dq2;
        t = tgt::Vector3<plot_t>(beta2 * Q2.xy() +(-beta1) *Q1.xy(),0);
    }

    else if (beta1_legal && !beta2_legal)   //beta1, beta3
    {
        SF = dq1*dq3;
        beta1 =beta1-beta2;   // all betas are multiplied by a positive SF
        beta3 =dr3*dq1;
        t = tgt::Vector3<plot_t>((SF-beta3-beta1)*Q1.xy() + beta3*Q2.xy(),0);
    }

    else if (beta2_legal && !beta1_legal) //beta2, beta3
    {
        SF = dq2*dq3;
        beta2 =beta1-beta2;   // all betas are multiplied by a positive SF
        beta3 =dr3*dq2;
        t = tgt::Vector3<plot_t>((SF-beta3)*Q1.xy() + (beta3-beta2)*Q2.xy(),0);
        Q1=Q2;
        beta1=beta2;
    }
    r4 = tgt::Vector3<plot_t>(SF*r.xy() + beta1*Q1.xy(),0);
    //seg_collide3(t,r4);  // calculates the 2D intersection
    tgt::Vector2<plot_t> p1,p2;
    p1=SF*P1.xy();
    p2=SF*P2.xy();
    det1 = p1[0]*t[1]-t[0]*p1[1];
    gama1 = (p1[0]*t[1]-t[0]*p1[1])*det1;
    alpha1 = (t[0]*r4[1] - r4[0]*t[1])*det1;
    alpha1_legal = (alpha1>=0) && (alpha1<=(det1*det1)  && (det1!=0));
    det2 = p2[0]*r4[1] - r4[0]*p2[1];
    alpha2 = (t[0]*r4[1] - r4[0]*t[1]) *det2;
    gama2 = (p2[0]*t[1] - t[0]*p2[1]) * det2;
    alpha2_legal = (alpha2>=0) && (alpha2<=(det2*det2) && (det2 !=0));
    det3=det2-det1;
    gama3=((p2[0]-p1[0])*(r[1]-p1[1]) - (r[0]-p1[0])*(p2[1]-p1[1]))*det3;
    if (alpha1_legal)
    {
        if (alpha2_legal)
        {
            if ( ((gama1<=0) && (gama1>=-(det1*det1))) || ((gama2<=0) && (gama2>=-(det2*det2))) || (gama1*gama2<0)) return 12;
        }
        else
        {
            if ( ((gama1<=0) && (gama1>=-(det1*det1))) || ((gama3<=0) && (gama3>=-(det3*det3))) || (gama1*gama3<0)) return 13;
        }
    }
    else
    if (alpha2_legal)
    {
        if ( ((gama2<=0) && (gama2>=-(det2*det2))) || ((gama3<=0) && (gama3>=-(det3*det3))) || (gama2*gama3<0)) return 23;
    }
    return 0;
}


double PlotLibraryFileBase::det(const tgt::Vector2<plot_t>& column1, const tgt::Vector2<plot_t>& column2) {
    return column1.x*column2.y-column1.y*column2.x;
}

tgt::Vector2<plot_t> PlotLibraryFileBase::perpendicularVector(const tgt::Vector2<plot_t>& vector) {
    return tgt::Vector2<plot_t>(vector.y,-vector.x);
}

tgt::Vector3<plot_t> PlotLibraryFileBase::Projection_Coordinates::get(Coordinate_Type type) {
    if (type == PlotLibraryFileBase::ORIGIN)
        return origin_;
    else if (type == PlotLibraryFileBase::MODELVIEW)
        return afterModelview_;
    else
        return afterProjection_;
}

} // namespace
