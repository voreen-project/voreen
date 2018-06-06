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

#ifndef VRN_LRITEM_H
#define VRN_LRITEM_H

#include "voreen/core/utils/GLSLparser/production.h"

#include <list>
#include <stdexcept>

namespace voreen {

namespace glslparser {

class GrammarSymbol;

/**
 * [A -> B.C]
 */
class LR0Item {
public:
    LR0Item(const std::string& name, const Production& production, int position);
    virtual ~LR0Item() {}

    bool operator==(const LR0Item& rhs) const;
    bool operator!=(const LR0Item& rhs) const;

    const std::string& getName() const { return name_; }
    int getPosition() const { return position_; }
    const Production& getProduction() const { return production_; }

    GrammarSymbol* getLeftSymbol() const;
    Production getLeftString() const;
    GrammarSymbol* getRightSymbol() const;
    Production getRightString() const;

    virtual std::string toString() const;

protected:
    std::string name_;
    Production production_;
    int position_;
};

// ============================================================================

/**
 * [A -> B.C, x]
 */
class LR1Item : public LR0Item {
public:
    LR1Item(const std::string& name, const Production& production, int position,
        GrammarSymbol* const lookahead);

    bool operator==(const LR1Item& rhs) const;
    bool operator!=(const LR1Item& rhs) const;

    GrammarSymbol* getLookAhead() const { return lookahead_; }

    virtual std::string toString() const;

protected:
    GrammarSymbol* lookahead_;
};

// ============================================================================

template<typename T>
class ItemSet {
public:
    ItemSet();
    ItemSet(const std::vector<T>& items);
    ItemSet(const ItemSet<T>& rhs);
    virtual ~ItemSet();

    const ItemSet<T>& operator=(const ItemSet<T>& rhs);
    bool operator==(const ItemSet<T>& rhs) const;
    bool operator!=(const ItemSet<T>& rhs) const;
    T& operator[](size_t index);
    const T& operator[](size_t index) const;

    bool addItem(const T& item);
    bool contained(const std::list<ItemSet<T>* >& collection) const;
    bool containsItem(const T& ref) const;
    typename std::list<ItemSet<T>* >::const_iterator findIn(const std::list<ItemSet<T>* >& collection) const;
    int findPositionIn(const std::list<ItemSet<T>* >& collection) const;
    int findKernelsPositionIn(const std::list<ItemSet<T>* >& collection) const;
    bool hasNextSet() const;
    int merge(const ItemSet<T>& other);

    const std::vector<T>& getItems() const { return items_; }
    const ItemSet<LR0Item>* getKernel() const;

    size_t getNumItems() const { return items_.size(); }
    bool isEmpty() const { return items_.empty(); }

    std::string toString() const;

private:
    void clone(const ItemSet<T>& rhs);
    void buildKernel();

protected:
    std::vector<T> items_;
    ItemSet<LR0Item>* kernel_;
};

typedef ItemSet<LR0Item> LR0ItemSet;
typedef ItemSet<LR1Item> LR1ItemSet;

}   // namespace glslparser

}   // namespace voreen

#endif
