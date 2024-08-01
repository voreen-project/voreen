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

#ifndef VRN_BRICKPOOLMANAGERQUEUE_H
#define VRN_BRICKPOOLMANAGERQUEUE_H

#include "tgt/types.h"
#include "tgt/assert.h"

namespace voreen {

template<class T>
class BrickPoolManagerQueueNode {
public:
    BrickPoolManagerQueueNode();
    BrickPoolManagerQueueNode(BrickPoolManagerQueueNode* previous, BrickPoolManagerQueueNode* next, T data);
    ~BrickPoolManagerQueueNode();

    BrickPoolManagerQueueNode* previous_;
    BrickPoolManagerQueueNode* next_;
    T data_;
};

template<class T>
class BrickPoolManagerQueue {
public:
    BrickPoolManagerQueue();
    ~BrickPoolManagerQueue();

    BrickPoolManagerQueueNode<T>* insertToFront(T data);
    void remove(BrickPoolManagerQueueNode<T>* node);
    T removeLast();
    void pushToFront(BrickPoolManagerQueueNode<T>* node);
    void clear();

    size_t size_;
    BrickPoolManagerQueueNode<T>* first_;
    BrickPoolManagerQueueNode<T>* last_;
};

//------------------------------------------------
//  Implementation
//------------------------------------------------
//  Node
template<class T>
BrickPoolManagerQueueNode<T>::BrickPoolManagerQueueNode()
    : previous_(0)
    , next_(0)
    //, data_(0)
{}

template<class T>
BrickPoolManagerQueueNode<T>::BrickPoolManagerQueueNode( BrickPoolManagerQueueNode* previous, BrickPoolManagerQueueNode* next, T data)
    : previous_(previous)
    , next_(next)
    , data_(data)
{}

template<class T>
BrickPoolManagerQueueNode<T>::~BrickPoolManagerQueueNode() {
}

//Queue
template<class T>
BrickPoolManagerQueue<T>::BrickPoolManagerQueue()
    : size_(0)
{
    first_ = new BrickPoolManagerQueueNode<T>();
    last_ = new BrickPoolManagerQueueNode<T>(first_,0,0);
    first_->next_ = last_;
}

template<class T>
BrickPoolManagerQueue<T>::~BrickPoolManagerQueue() {
    clear();
    delete first_;
    delete last_;
}

template<class T>
void BrickPoolManagerQueue<T>::remove(BrickPoolManagerQueueNode<T>* node) {
    node->previous_->next_ = node->next_;
    node->next_->previous_ = node->previous_;
    delete node;
}

template<class T>
BrickPoolManagerQueueNode<T>* BrickPoolManagerQueue<T>::insertToFront(T data) {
    BrickPoolManagerQueueNode<T>* newNode = new BrickPoolManagerQueueNode<T>(first_,first_->next_,data);
    first_->next_->previous_ = newNode;
    first_->next_ = newNode;
    size_++;
    return newNode;
}

template<class T>
T BrickPoolManagerQueue<T>::removeLast() {
    tgtAssert(size_ > 0, "Remove called on empty queue!");
    BrickPoolManagerQueueNode<T>* tmp = last_->previous_;
    tmp->previous_->next_ = last_;
    last_->previous_ = tmp->previous_;
    T data = tmp->data_;
    delete tmp;
    size_--;
    return data;
}

template<class T>
void BrickPoolManagerQueue<T>::pushToFront(BrickPoolManagerQueueNode<T>* node) {
    tgtAssert(node,"No node passed!");
    //remove old node
    node->previous_->next_ = node->next_;
    node->next_->previous_ = node->previous_;
    //insert into fromt
    node->previous_ = first_;
    node->next_ = first_->next_;
    first_->next_->previous_ = node;
    first_->next_ = node;
}

template<class T>
void BrickPoolManagerQueue<T>::clear() {
    BrickPoolManagerQueueNode<T>* next = first_->next_;
    while(next != last_) {
        BrickPoolManagerQueueNode<T>* tmp = next->next_;
        delete next;
        next = tmp;
    }
    first_->next_ = last_;
    last_->previous_ = first_;
}

} // voreen namespace

#endif //end VRN_BRICKPOOLMANAGERQUEUE_H
