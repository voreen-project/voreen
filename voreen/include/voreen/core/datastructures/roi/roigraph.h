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

#ifndef VRN_ROIGRAPH_H
#define VRN_ROIGRAPH_H

#include "voreen/core/datastructures/roi/roicube.h"

namespace voreen {

class Graph;

class GraphNode : public Serializable {
    friend class Graph;
protected:
    tgt::vec3 pos_;
    std::vector<GraphNode*> nbs_;

public:
    GraphNode() : pos_(0.0f) {}
    GraphNode(tgt::vec3 pos) : pos_(pos) {
    }

    tgt::vec3 getPosition() const { return pos_; }
    void setPosition(tgt::vec3 pos) { pos_ = pos; }

    tgt::vec3 getTangent() const;

    GraphNode* getNeighbor(size_t i) const { return nbs_[i]; }
    size_t getNumNeighbors() const { return nbs_.size(); }

    void removeLinks() {
        for(size_t i=0; i<nbs_.size(); ++i) {
            GraphNode* n = nbs_[i];

            for(size_t j=0; j<n->nbs_.size(); ++j) {
                if(n->nbs_[j] == this) {
                    n->nbs_.erase(n->nbs_.begin()+j);
                    break;
                }
            }
        }
        nbs_.clear();
    }

    void removeNode() {
        nbs_.front()->nbs_.push_back(nbs_.back());
        nbs_.back()->nbs_.push_back(nbs_.front());
        removeLinks();
    }

    void connect(GraphNode* nb) {
        nbs_.push_back(nb);
        nb->nbs_.push_back(this);
    }

    GraphNode* getNextCritical(const GraphNode* from) {
        if(getNumNeighbors() != 2)
            return this;

        for(size_t i=0; i<getNumNeighbors(); ++i) {
            GraphNode* n = nbs_[i];
            if(n != from)
                return n->getNextCritical(this);
        }

        return 0;
    }

    const GraphNode* getNextLinear(const GraphNode* from) const {
        if(getNumNeighbors() != 2)
            return this;

        if(nbs_[0] != from)
            return nbs_[0];
        else
            return nbs_[1];
    }

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    void deleteZeroLinks();
};

tgt::vec3 catmullRom(const GraphNode* n1, const GraphNode* n2, float t);

class Graph : public Serializable {
public:
    Graph() {}
    void addNode(GraphNode* gn) { nodes_.push_back(gn); }
    size_t getNumNodes() const { return nodes_.size(); }
    const GraphNode* getNode(size_t i) const { return nodes_[i]; }
    GraphNode* getNode(size_t i) { return nodes_[i]; }

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
protected:
    std::vector<GraphNode*> nodes_;
};

class ROIGraph : public ROISingle {
public:
    ROIGraph(Grid grid, Graph g);
    ROIGraph();
    virtual std::string getClassName() const { return "ROIGraph"; }
    virtual ROIBase* create() const { return new ROIGraph(); }

    virtual bool inROI(tgt::vec3 p) const;
    const Graph& getGraph() const;
    Graph& getGraph();

    virtual Geometry* generateMesh() const;
    virtual Geometry* generateMesh(tgt::plane pl) const;
    virtual Geometry* generateRasterMesh(const tgt::plane& pl, Grid g) const;
    virtual tgt::Bounds getBoundingBox() const;

    virtual std::vector<const ControlPoint*> getControlPoints() const;
    virtual std::vector<const ControlPoint*> getControlPoints(tgt::plane pl) const;
    bool moveControlPoint(const ControlPoint* cp, tgt::vec3 offset);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
private:
    Graph g_;

    static const std::string loggerCat_;
};

} //namespace

#endif
