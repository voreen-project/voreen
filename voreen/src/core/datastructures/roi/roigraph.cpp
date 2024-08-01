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

#include "voreen/core/datastructures/roi/roigraph.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

using tgt::vec3;
using tgt::vec4;
using tgt::mat4;

//-----------------------------------------------------------------------------

tgt::vec3 GraphNode::getTangent() const {
    vec3 n(0.0f);
    if(getNumNeighbors() == 1)
        n = normalize(getNeighbor(0)->getPosition() - getPosition());
    else if(getNumNeighbors() == 1)
        n = normalize(getNeighbor(0)->getPosition() - getNeighbor(1)->getPosition());

    return n;
}

void GraphNode::serialize(Serializer& s) const {
    s.serialize("pos", pos_);
    s.serialize("neighbors", nbs_);
}

void GraphNode::deserialize(Deserializer& s) {
    s.deserialize("pos", pos_);
    s.deserialize("neighbors", nbs_);
}

void GraphNode::deleteZeroLinks() {
    if(nbs_.size() == 0)
        return;

    for(int j=static_cast<int>(nbs_.size()); j>=0; j--) {
        if(nbs_.size() == 0)
            return;

        if(nbs_[j] == 0) {
            if(j < (int)nbs_.size())
                nbs_.erase(nbs_.begin()+j);
            //LWARNINGC("voreen.graphnode", "deleted zerolink");
        }
    }
}

//-----------------------------------------------------------------------------

vec3 catmullRom(const GraphNode* n1, const GraphNode* n2, float t) {
    const GraphNode* n0 = n1->getNextLinear(n2);
    vec3 m0 = vec3(0.0f);
    if(n0) {
        m0 = ( (n1->getPosition() - n0->getPosition()) + (n2->getPosition() - n1->getPosition()) ) * 0.5f;
    }
    const GraphNode* n3 = n2->getNextLinear(n1);
    vec3 m1 = vec3(0.0f);
    if(n3) {
        m1 = ( (n2->getPosition() - n1->getPosition()) + (n3->getPosition() - n2->getPosition()) ) * 0.5f;
    }
    return hermite(n1->getPosition(), n2->getPosition(), m0, m1, t);
}

//-----------------------------------------------------------------------------

void Graph::serialize(Serializer& s) const {
    s.serialize("nodes", nodes_);
}

void Graph::deserialize(Deserializer& s) {
    s.deserialize("nodes", nodes_);

    //Fix broken (null) links
    for(size_t i=0; i<nodes_.size(); i++) {
        nodes_[i]->deleteZeroLinks();

        for(int j=(static_cast<int>(nodes_[i]->getNumNeighbors())-1); j>=0; j--) {
            GraphNode* n = nodes_[i]->getNeighbor(j);

            bool found = false;
            for(size_t k=0; k<n->getNumNeighbors(); k++) {
                if(n->getNeighbor(k) == nodes_[i]) {
                    found = true;
                    break;
                }
            }
            if(!found) {
                n->nbs_.push_back(nodes_[i]);
                //LWARNINGC("voreen.graph", "Fixed link!");
            }
        }
    }

}

//-----------------------------------------------------------------------------
const std::string ROIGraph::loggerCat_ = "voreen.ROIGraph";

ROIGraph::ROIGraph() : ROISingle()
{
}

ROIGraph::ROIGraph(Grid grid, Graph g) : ROISingle(grid), g_(g)
{
}

bool ROIGraph::inROI(tgt::vec3 /*p*/) const {
    return false;
}

const Graph& ROIGraph::getGraph() const {
    return g_;
}

Graph& ROIGraph::getGraph() {
    return g_;
}

Geometry* ROIGraph::generateMesh() const {
    MeshListGeometry* mlg = new MeshListGeometry();

    MeshGeometry mg;
    for(size_t i=0; i<g_.getNumNodes(); i++) {
        const GraphNode* n = g_.getNode(i);
        vec3 p = n->getPosition();

        for(size_t j=0; j<n->getNumNeighbors(); j++) {
            if(n->getNeighbor(j) < n) { //prevent double edges
                vec3 last = p;
                for(int i=1; i<11; i++) {
                    FaceGeometry fg;
                    fg.addVertex(VertexGeometry(last, vec3(0.0f), getColor()));
                    last = catmullRom(n, n->getNeighbor(j), (float)i/10.0f);
                    fg.addVertex(VertexGeometry(last, vec3(0.0f), getColor()));
                    mg.addFace(fg);
                }
            }
        }
    }

    mlg->addMesh(mg);
    mlg->transform(getGrid().getPhysicalToWorldMatrix());
    return mlg;
}

Geometry* ROIGraph::generateMesh(tgt::plane pl) const {
    // Transform plane to grid coordinates:
    mat4 m = getGrid().getWorldToPhysicalMatrix();
    const tgt::plane physicalPlane = pl.transform(m);

    MeshListGeometry* mlg = new MeshListGeometry();

    size_t maxComp = tgt::maxElem(physicalPlane.n);
    float d = getGrid().getSpacing()[maxComp] * 0.5f;

    tgt::plane clipA = physicalPlane;
    clipA.d += d;
    tgt::plane clipB = physicalPlane;
    clipB.n *= -1.0f;
    clipB.d *= -1.0f;
    clipB.d += d;

    MeshGeometry mg;
    for(size_t i=0; i<g_.getNumNodes(); i++) {
        const GraphNode* n = g_.getNode(i);

        for(size_t j=0; j<n->getNumNeighbors(); j++) {
            if(n->getNeighbor(j) < n) { //prevent double edges
                vec3 p1 = n->getPosition();
                vec3 p2 = n->getNeighbor(j)->getPosition();

                float d1 = dot(p1, physicalPlane.n) - physicalPlane.d;
                float d2 = dot(p2, physicalPlane.n) - physicalPlane.d;
                if( (d1 * d2 < 0.0f)  || (fabs(d1) < d) || (fabs(d2) < d) ){
                    FaceGeometry fg;
                    fg.addVertex(VertexGeometry(p1, vec3(0.0f), getColor()));
                    fg.addVertex(VertexGeometry(p2, vec3(0.0f), getColor()));

                    fg.clip(clipA);
                    fg.clip(clipB);

                    mg.addFace(fg);
                }
            }
        }
    }

    mlg->addMesh(mg);
    mlg->transform(getGrid().getPhysicalToWorldMatrix());
    return mlg;
}

Geometry* ROIGraph::generateRasterMesh(const tgt::plane& /*pl*/, Grid /*g*/) const {
    return 0;
}

tgt::Bounds ROIGraph::getBoundingBox() const {
    return tgt::Bounds();
}

std::vector<const ControlPoint*> ROIGraph::getControlPoints() const {
    //return std::vector<const ControlPoint*>();
    std::vector<const ControlPoint*> cps;

    for(size_t i=0; i<g_.getNumNodes(); i++) {
        const GraphNode* n = g_.getNode(i);
        vec3 p = n->getPosition();

        cps.push_back(new ControlPoint(this, (int)i, getGrid().getPhysicalToWorldMatrix() * p));
    }

    return cps;
}

std::vector<const ControlPoint*> ROIGraph::getControlPoints(tgt::plane pl) const {
    std::vector<const ControlPoint*> cps;

    // Transform plane to grid coordinates:
    mat4 m = getGrid().getWorldToPhysicalMatrix();
    const tgt::plane physicalPlane = pl.transform(m);

    size_t maxComp = tgt::maxElem(physicalPlane.n);
    float d = getGrid().getSpacing()[maxComp] * 0.5f;

    for(size_t i=0; i<g_.getNumNodes(); i++) {
        const GraphNode* n = g_.getNode(i);
        vec3 p = n->getPosition();

        float d1 = dot(p, physicalPlane.n) - physicalPlane.d;
        if(fabs(d1) < d)
            cps.push_back(new ControlPoint(this, (int)i, getGrid().getPhysicalToWorldMatrix() * p));
    }

    return cps;
}

bool ROIGraph::moveControlPoint(const ControlPoint* cp, tgt::vec3 offset) {
    if(cp->getROI() != this)
        return false;

    offset = getGrid().getWorldToPhysicalMatrix().getRotationalPart() * offset;

    if((size_t)cp->getID() < g_.getNumNodes()) {
        GraphNode* n = g_.getNode(cp->getID());

        n->setPosition(n->getPosition() + offset);
        invalidate();
        return true;
    }

    return false;
}

void ROIGraph::serialize(Serializer& s) const {
    s.serialize("graph", g_);

    ROISingle::serialize(s);
}

void ROIGraph::deserialize(Deserializer& s) {
    s.deserialize("graph", g_);

    ROISingle::deserialize(s);
}

} // namespace
