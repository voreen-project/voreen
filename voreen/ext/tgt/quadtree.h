/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2019 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#ifndef TGT_QUADTREE_H
#define TGT_QUADTREE_H

#include <vector>

#include "tgt/assert.h"
#include "tgt/bounds.h"


/*
    NO 3DSP DEPENDENCIES IN THIS FILE !!!
    COMMIT ONLY WORKING FILES -- DO NOT COMMIT FILES HAVING CONFLICTS !!!
*/

namespace tgt {

/*
    TODO remove unnecessary(?) null-pointer checks (in destructor...)
    TODO remove leavecounter
    TODO add const_iterator
    TODO add const_leaf_iterator
*/

/**
    This is a generic quadtree.
    T must inherit from HasBounds.
    Use T directly and NOT as pointer like this:<br>
    <br>
    tgt::Quadtree\<FooBar*\> quadtree; <br>
    Instead do this: <br>
    tgt::Quadtree\<FooBar\> quadtree; <br>
    <br>
    because FooBar is already stored as pointer internally.<br>
    The Quadtree lies in the X-Z-plane. <br>
    The Quadtree DOES NOT delete any objects that has been inserted. <br>
*/
template<class T>
class Quadtree {
private:

    // forward declaration
    class QuadNode;

    /**
     * tgt::Quadtree::QuadNodeBase class
    */
    class QuadNodeBase {
    public:

        QuadNode*       parent_; // NULL if root

        Bounds  bounds_;
        Bounds  extendedBounds_;

        float   threshold_;
        bool    isLeaf_;

        size_t  finishedChildren_;

    public:

        QuadNodeBase(const Bounds& bounds, float threshold, bool isLeaf, QuadNode* parent)
          : parent_(parent),
            bounds_(bounds),
            extendedBounds_(bounds),
            threshold_(threshold),
            isLeaf_(isLeaf),
            finishedChildren_(0)
        {}

        virtual ~QuadNodeBase() {}

        virtual tgt::Bounds* insert(T* obj, float weight) = 0;
    };

    /**
     * tgt::Quadtree::QuadNode class
    */
    class QuadNode : public QuadNodeBase {
    private:
#if (defined(_MSC_VER) && _MSC_VER < 1400) || defined(__ICC)
        friend class iterator;
        friend class leaf_iterator;
#else
        friend class Quadtree<T>::iterator;
        friend class Quadtree<T>::leaf_iterator;
#endif
        friend class Quadtree<T>;

        QuadNodeBase*   children_[4];

    public:

        QuadNode(const Bounds& bounds, float threshold, QuadNode* parent)
          : QuadNodeBase(bounds, threshold, false, parent) {
            for (size_t i = 0; i < 4; ++i) {
                children_[i] = new QuadLeaf(
                    calculateChildBounds(bounds, i),
                    threshold,
                    this);
            }
        }

        virtual ~QuadNode() {
            for (size_t i = 0; i < 4; ++i)
                if (children_[i])
                    delete children_[i];
        }

        virtual tgt::Bounds* insert(T* obj, float weight);

    };

public:

    /**
     * tgt::Quadtree::QuadLeaf class
    */
    class QuadLeaf : public QuadNodeBase {
    public:

        std::vector<T*>     objs_;
        std::vector<float>  weights_;
        float               weightSum_;

    public:

        QuadLeaf(const Bounds& bounds, float threshold, QuadNode* parent)
          : QuadNodeBase(bounds, threshold, true, parent),
            weightSum_(0)
        {}
        virtual ~QuadLeaf() {}

        virtual tgt::Bounds* insert(T* obj, float weight);
    };

    /**
     * Leaf-Iterator class with a small "i" like in the STL.
     * Iterates over all leaves and only leaves.
    */
    class leaf_iterator {
    private:

        QuadLeaf* leaf_;

        static QuadLeaf* iterate(QuadNodeBase* child);

    public:

        leaf_iterator(QuadLeaf* leaf)
          : leaf_(leaf)
        {}

        ///  preincrement operator
        leaf_iterator& operator ++ () {
            leaf_ = iterate(leaf_);
            return *this;
        }

        /// postincrement operator
        leaf_iterator  operator ++(int) {
            leaf_iterator iter = *this;
            leaf_ = iterate(leaf_);
            return iter;
        }

        bool operator == (const leaf_iterator& iter) {
            return leaf_ == iter.leaf_;
        }

        bool operator != (const leaf_iterator& iter) {
            return leaf_ != iter.leaf_;
        }

        /// adress operator
        QuadLeaf* operator&() {
            return leaf_;
        }

        QuadLeaf* operator->() {
            return leaf_;
        }

        QuadLeaf& operator*() {
            return *leaf_;
        }
    };

    /**
    * Iterator class with a small "i" like in the STL, this is supposed
    * to be used when it is necessary to traverse the Quadtree and skip
    * certain branches.
    */
    class iterator {
    private:

        QuadNodeBase* node_;
        const QuadNodeBase* root_;

        QuadNodeBase* iterate(QuadNodeBase* child);

    public:

        iterator(QuadNodeBase* node)
          : node_(node), root_(node)
        {}

        ///  preincrement operator
        iterator& operator ++ () {
            node_ = iterate(node_);
            return *this;
        }

        /// postincrement operator
        iterator  operator ++(int) {
            iterator iter = *this;
            node_ = iterate(node_);
            return iter;
        }

        bool operator == (const iterator& iter) {
            return node_ == iter.node_;
        }

        bool operator != (const iterator& iter) {
            return node_ != iter.node_;
        }

        /// similar to the increment-operators, but skips an entire branch
        void skip();

        /// adress operator
        QuadNodeBase* operator&() {
            return node_;
        }

        QuadNodeBase* operator->() {
            return node_;
        }

        QuadNodeBase& operator*() {
            return *node_;
        }
    };

    enum {
        TOP_LEFT = 0,
        TOP_RIGHT,
        BOTTOM_LEFT,
        BOTTOM_RIGHT
    };

/*
    Quadtree declarations
*/
    Quadtree(const Bounds& bounds, float threshold)
        : rootNode_(new QuadNode(bounds, threshold, NULL))
    {}

    ~Quadtree() {
        delete rootNode_;
    }

    void insert(T* obj, float weight) {
        rootNode_->insert(obj, weight);
    }

    T* find(float x, float y);

    static Bounds calculateChildBounds(const Bounds& bounds, size_t childIndex);

    leaf_iterator leaf_begin() {

        QuadNode* tmpnode;
        tmpnode = rootNode_;
        while (!(tmpnode->children_[0]->isLeaf_))
            tmpnode = (QuadNode*)(tmpnode->children_[0]);

        return leaf_iterator((QuadLeaf*)(tmpnode->children_[0]));
    };

    leaf_iterator leaf_end() {
        return leaf_iterator(NULL);
    };

    iterator begin() {
        return iterator(rootNode_);
    }

    iterator end() {
        return iterator(NULL);
    };


    protected:
        QuadNode*    rootNode_;

    public:
        QuadNode* getRootNode(){ return rootNode_;}

};

/*
    Non-inline-implementation
*/

/*
    helper
*/

template<class T>
Bounds Quadtree<T>::calculateChildBounds(const Bounds& bounds, size_t childIndex) {
    // calculate new bounds
    vec3 llf = bounds.getLLF();
    vec3 urb = bounds.getURB();

    float xdelta = urb.x - llf.x;
    float zdelta = urb.z - llf.z;

    xdelta *= 0.5f;
    zdelta *= 0.5f;

    if (childIndex == 0 || childIndex == 1)
        urb.z -= zdelta;
    if (childIndex == 2 || childIndex == 3)
        llf.z += zdelta;
    if (childIndex == 0 || childIndex == 2)
        urb.x -= xdelta;
    if (childIndex == 1 || childIndex == 3)
        llf.x += xdelta;

    return Bounds(llf, urb);
};

/*
    insert
*/

template<class T>
Bounds* Quadtree<T>::QuadNode::insert(T* obj, float weight) {
    tgtAssert(obj->getBounds().isDefined() || obj->getBounds().onlyPoint(), "Bounding Box of inserted object not defined");

    for (size_t i = 0; i < 4; ++i) {
        // Must be inserted here?
        if ( children_[i]->bounds_.insideXZ(obj->getBounds().center()) ) {
            // Is the child a leaf?
            if (children_[i]->isLeaf_) {
                // Is the threshold exceeded?
                if (((QuadLeaf*) children_[i])->weightSum_ + weight > QuadNodeBase::threshold_) {
                    // Yes? - Change to QuadNode and insert contents there
                    // Store old values
                    std::vector<T*>    objs    = ((QuadLeaf*) children_[i])->objs_;
                    std::vector<float> weights = ((QuadLeaf*) children_[i])->weights_;
                    // destroy node
                    delete ((QuadLeaf*) children_[i]);
                    // create node
                    children_[i] = new QuadNode(
                                        calculateChildBounds(QuadNodeBase::bounds_, i),
                                        QuadNodeBase::threshold_,
                                        this); // this is the parent

                    // insert object
                    tgtAssert(objs.size() == weights.size(), "The size must be equal.");
                    typename std::vector<T*>::iterator iterObjs    = objs.begin();
                           std::vector<float>::iterator iterWeights = weights.begin();
                    while (iterObjs != objs.end()) {
                        QuadNodeBase::extendedBounds_.addVolume(
                            *children_[i]->insert( *iterObjs, *iterWeights) );
                        ++iterObjs;
                        ++iterWeights;
                    }
                }
            }
            // insert the actual element at all events
            QuadNodeBase::extendedBounds_.addVolume(
                *children_[i]->insert(obj, weight) );
        }
    }

// GCC 3.3 seems to have problems with this
#if (__GNUC__ <= 3) && (__GNUC_MINOR__ <= 3)
    Bounds* result = &extendedBounds_;
#else
    Bounds& bounds = QuadNodeBase::extendedBounds_;
    Bounds* result = &bounds;
#endif

    return result;
}

template<class T>
Bounds* Quadtree<T>::QuadLeaf::insert(T* obj, float weight) {
    objs_.   push_back(obj);
    weights_.push_back(weight);

    weightSum_ += weight;

    tgtAssert(weightSum_ <= QuadNodeBase::threshold_, "threshold_ greater than sum of weights");

    QuadNodeBase::extendedBounds_.addVolume( obj->getBounds() );

    // GCC 3.3 seems to have problems with this
#if (__GNUC__ <= 3) && (__GNUC_MINOR__ <= 3)
    Bounds* result = &extendedBounds_;
#else
    Bounds& bounds = QuadNodeBase::extendedBounds_;
    Bounds* result = &bounds;
#endif

    return result;
}

// FIXME remove this method: Use iterators
template<class T>
T* Quadtree<T>::find(float x, float y) {
    tgt::vec3 objectCoordinates(x,1.f,y);
    QuadNode* currentNode = rootNode_;
    bool done=false;
    int i=0;
    while (!done) {
        for (i=0;i<4;i++) {
            if (currentNode->children_[i]->bounds_.insideXZ(objectCoordinates)) {
                if (!currentNode->children_[i]->isLeaf_) {
                    currentNode =dynamic_cast<QuadNode*>(currentNode->children_[i]);
                    i=4;
                } else {
                    QuadLeaf* currentLeaf = (QuadLeaf*) currentNode->children_[i];
                    std::vector<T*> items = currentLeaf->objs_;
                    for (size_t k = 0; k < items.size(); k++) {
                        // wtf ist getWorldBoundingBox()? where does this come from?
                        if ((items[k])->getWorldBoundingBox().insideXZ(objectCoordinates)) {
                            done=true;
                            return(items[k]);
                        }
                    }
                    return(0);
                }
            }
        }
        if (i==4) done=true;
    }
    return(0);
}

/*
    leaf_iterator
*/
template<class T>
typename Quadtree<T>::QuadLeaf* Quadtree<T>::leaf_iterator::iterate(QuadNodeBase* child) {

    while (child->parent_) {

        size_t i = 0;
        // find childIdx
        while (child->parent_->children_[i] != child)
            ++i;

        if (i == 3) {
            child = child->parent_;
        }
        else {
            // go to sibling
            ++i;
            child = child->parent_->children_[i];

            while (!child->isLeaf_)
                child = ((QuadNode*)child)->children_[0];

            return (QuadLeaf*) child;
        }
    }

    // OK, finished
    return NULL;
}

/*
    iterator
*/
template<class T>
typename Quadtree<T>::QuadNodeBase* Quadtree<T>::iterator::iterate(QuadNodeBase* child) {

    if (child->isLeaf_ || child->finishedChildren_ == 4) { // we are finished with this node

        if (child->isLeaf_)
            child = child->parent_;

        while (child->finishedChildren_ == 4) { // find the first unfinished parent-node
            child->finishedChildren_ = 0;
            if (child == root_)
                return NULL; // we are done
            child = child->parent_;
        }
    }

    return ((QuadNode*)child)->children_[child->finishedChildren_++];
}

template<class T>
void Quadtree<T>::iterator::skip() {

    node_->finishedChildren_ = 0;

    if (node_->parent_)
        node_ = iterate(node_->parent_);

    else
        node_ = NULL; // the node we just skipped is the root-node
}

}

#endif //TGT_QUADTREE_H
