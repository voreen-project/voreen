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

#pragma once
#include <memory>
#include <vector>
#include "tgt/bounds.h"
#include "tgt/memory.h"

#include <boost/optional.hpp>


#define Bounds tgt::TemplateBounds<T>
namespace voreen {

enum class BoundsHierarchyNodeType {
    INNER,
    LEAF,
};

template<typename T, typename V>
union BoundsHierarchyNode;

template<typename T, typename V>
struct BoundsHierarchyInnerNode {
    BoundsHierarchyNodeType type;
    Bounds bounds;
    BoundsHierarchyNode<T,V>* left; //owning
    BoundsHierarchyNode<T,V>* right; //owning


    ~BoundsHierarchyInnerNode() {
        delete left;
        delete right;
    }
};

template<typename T, typename V>
struct BoundsHierarchyLeaf {
    BoundsHierarchyNodeType type;
    Bounds bounds;
    V val;
};


template<typename T, typename V>
union BoundsHierarchyNode {
public:

    BoundsHierarchyNode(BoundsHierarchyNode<T,V>&& left, BoundsHierarchyNode<T,V>&& right);
    BoundsHierarchyNode(V val, Bounds bounds);
    BoundsHierarchyNode(BoundsHierarchyNode<T,V>&& other);
    BoundsHierarchyNode& operator=(BoundsHierarchyNode<T,V>&& other);
    ~BoundsHierarchyNode();
    const Bounds& getBounds() const;

    void findBounds(const tgt::Vector3<T>& p, std::vector<V>& res) const;

    void checkValidity() const {
        switch(inner.type) {
            case BoundsHierarchyNodeType::INNER:
                tgtAssert(inner.bounds.isDefined() && inner.left && inner.right, "Invalid node");
                inner.left->checkValidity();
                inner.right->checkValidity();

                break;
            case BoundsHierarchyNodeType::LEAF:
                tgtAssert(leaf.bounds.isDefined(), "Invalid node");
                break;
        }
    }

//private:
    BoundsHierarchyInnerNode<T,V> inner;
    BoundsHierarchyLeaf<T,V> leaf;
};

template<typename T, typename V>
class BoundsHierarchy {
public:
    static boost::optional<BoundsHierarchy> buildBottomUp(std::vector<std::pair<V, Bounds>>&&);
    static boost::optional<BoundsHierarchy> buildTopDown(std::vector<std::pair<V, Bounds>>&&);

    BoundsHierarchy(BoundsHierarchyNode<T,V>&&) noexcept;
    BoundsHierarchy(BoundsHierarchy&&) noexcept;
    BoundsHierarchy& operator=(BoundsHierarchy&& other);

    std::vector<V> findBounds(const tgt::Vector3<T>& p) const;

private:
    BoundsHierarchyNode<T,V> root;
};

/// Implementation -------------------------------------------------------------


template<typename T>
static Bounds combine(const Bounds& left, const Bounds& right) {
    Bounds res;
    res.addVolume(left);
    res.addVolume(right);
    return res;
}

template<typename T, typename V>
BoundsHierarchyNode<T,V>::BoundsHierarchyNode(BoundsHierarchyNode<T,V>&& left, BoundsHierarchyNode<T,V>&& right) {
    inner.type = BoundsHierarchyNodeType::INNER;
    inner.bounds = combine(left.getBounds(), right.getBounds());
    inner.left = new BoundsHierarchyNode(std::move(left));
    inner.right = new BoundsHierarchyNode(std::move(right));
}

template<typename T, typename V>
BoundsHierarchyNode<T,V>::BoundsHierarchyNode(V val, Bounds bounds) {
    leaf.type = BoundsHierarchyNodeType::LEAF;
    leaf.val = val;
    leaf.bounds = bounds;
}

template<typename T, typename V>
BoundsHierarchyNode<T,V>::BoundsHierarchyNode(BoundsHierarchyNode<T,V>&& other) {
    switch(other.inner.type) {
        case BoundsHierarchyNodeType::INNER:
            inner.type = BoundsHierarchyNodeType::INNER;
            inner.bounds = other.inner.bounds;
            inner.left = other.inner.left;
            inner.right = other.inner.right;

            other.inner.left = nullptr;
            other.inner.right = nullptr;
            other.inner.bounds = Bounds();
            break;
        case BoundsHierarchyNodeType::LEAF:
            leaf.type = BoundsHierarchyNodeType::LEAF;
            leaf.val = other.leaf.val;
            leaf.bounds = other.leaf.bounds;

            other.leaf.bounds = Bounds();
            break;
    }
}

template<typename T, typename V>
BoundsHierarchyNode<T,V>& BoundsHierarchyNode<T,V>::operator=(BoundsHierarchyNode<T,V>&& other) {
    if(this != &other) {
        this->~BoundsHierarchyNode();
        new(this) BoundsHierarchyNode(std::move(other));
    }
    return *this;
}

template<typename T, typename V>
BoundsHierarchyNode<T,V>::~BoundsHierarchyNode() {
    switch(inner.type) {
        case BoundsHierarchyNodeType::INNER:
            inner.~BoundsHierarchyInnerNode();
            break;
        case BoundsHierarchyNodeType::LEAF:
            leaf.~BoundsHierarchyLeaf();
            break;
    }
}

template<typename T, typename V>
const Bounds& BoundsHierarchyNode<T,V>::getBounds() const {
    switch(inner.type) {
        case BoundsHierarchyNodeType::INNER:
            return inner.bounds;
        case BoundsHierarchyNodeType::LEAF:
            return leaf.bounds;
        default:
            tgtAssert(false, "Unknown type");
            return leaf.bounds;
    }
}

template<typename T, typename V>
void BoundsHierarchyNode<T,V>::findBounds(const tgt::Vector3<T>& p, std::vector<V>& res) const {
    if(!getBounds().containsPoint(p)) {
        return;
    }
    switch(inner.type) {
        case BoundsHierarchyNodeType::INNER:
            inner.left->findBounds(p, res);
            inner.right->findBounds(p, res);
            break;
        case BoundsHierarchyNodeType::LEAF:
            res.push_back(leaf.val);
            break;
    }
}

template<typename T, typename V>
std::vector<V> BoundsHierarchy<T,V>::findBounds(const tgt::Vector3<T>& p) const {
    std::vector<V> res;
    root.findBounds(p, res);
    return res;
}

template<typename T, typename V>
static BoundsHierarchyNode<T,V> buildHierarchyBottomUp(const std::vector<std::pair<V, Bounds>>&& bounds) {
    tgtAssert(bounds.size() > 0, "No bounds");

    std::vector<BoundsHierarchyNode<T,V>> nodes;
    for(auto& pair : bounds) {
        nodes.emplace_back(pair.first, pair.second);
    }

    while(nodes.size() > 1) {
        std::pair<size_t, size_t> best;
        size_t best_volume = -1;

        for(int i=0; i<nodes.size(); ++i) {
            for(int j=0; j<nodes.size(); ++j) {
                if(i==j) {
                    continue;
                }
                size_t vol = combine(nodes[i].getBounds(), nodes[j].getBounds()).volume();
                if(vol < best_volume) {
                    best = std::make_pair(i, j);
                    best_volume = vol;
                }
            }
        }

        BoundsHierarchyNode<T,V> new_node(std::move(nodes[best.first]), std::move(nodes[best.second]));

        // Remove empty nodes

        std::swap(nodes[best.first], nodes[nodes.size()-1]); // Restore element at best.first from end of the vec
        if(best.second != nodes.size()-1) {
            std::swap(nodes[best.second], nodes[nodes.size()-2]);
            //But!: The above does not make sense if best.second == last element, so...
        } else {
            // ...in that case we have to restore best.first
            std::swap(nodes[best.first], nodes[nodes.size()-2]);
        }
        nodes.pop_back();
        nodes.pop_back();

        nodes.push_back(std::move(new_node));
        tgtAssert(nodes.back().inner.left, "Nullptr ref");
        tgtAssert(nodes.back().inner.right, "Nullptr ref");


        for(auto& node : nodes) {
            node.checkValidity();
        }
    }
    tgtAssert(nodes.size() == 1, "Invalid nodes size");

    return std::move(nodes[0]);
}

template<typename T, typename V>
static BoundsHierarchyNode<T,V> buildHierarchyTopDown(std::vector<std::pair<V, Bounds>>&& bounds, int dim) {
    tgtAssert(bounds.size() > 0, "No bounds");

    if(bounds.size() == 1) {
        auto p = bounds.front();
        return BoundsHierarchyNode<T,V>(p.first, p.second);
    }

    auto begin = bounds.begin();
    auto end = bounds.end();

    std::sort(bounds.begin(), bounds.end(), [dim] (const std::pair<V, Bounds>& p1, const std::pair<V, Bounds>& p2) {
        return p1.second.getLLF()[dim] < p2.second.getLLF()[dim];
    });

    auto mid = bounds.begin() + bounds.size()/2;

    int nextDim = dim == 2 ? 0 : dim+1;
    auto l = buildHierarchyTopDown(std::vector<std::pair<V, Bounds>>(begin, mid), nextDim);
    auto r = buildHierarchyTopDown(std::vector<std::pair<V, Bounds>>(mid, end), nextDim);

    return BoundsHierarchyNode<T,V>(std::move(l), std::move(r));
}

template<typename T, typename V>
BoundsHierarchy<T,V>::BoundsHierarchy(BoundsHierarchyNode<T,V>&& root) noexcept
    : root(std::move(root))
{
}

template<typename T, typename V>
BoundsHierarchy<T,V>::BoundsHierarchy(BoundsHierarchy&& other) noexcept
    : root(std::move(other.root))
{
}
template<typename T, typename V>
BoundsHierarchy<T,V>& BoundsHierarchy<T,V>::operator=(BoundsHierarchy<T,V>&& other) {
    this->~BoundsHierarchyNode();
    new(this) BoundsHierarchy(std::move(other));
    return *this;
}
template<typename T, typename V>
boost::optional<BoundsHierarchy<T,V>> BoundsHierarchy<T,V>::buildBottomUp(std::vector<std::pair<V, Bounds>>&& bounds) {
    if(bounds.empty()) {
        return boost::none;
    }
    return BoundsHierarchy<T,V>(buildHierarchyBottomUp(std::move(bounds)));
}
template<typename T, typename V>
boost::optional<BoundsHierarchy<T,V>> BoundsHierarchy<T,V>::buildTopDown(std::vector<std::pair<V, Bounds>>&& bounds) {
    if(bounds.empty()) {
        return boost::none;
    }
    return BoundsHierarchy<T,V>(buildHierarchyTopDown(std::move(bounds), 0));
}

#undef Bounds

}
