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

#include "voreen/core/utils/GLSLparser/generator/lritem.h"

#include "voreen/core/utils/GLSLparser/grammarsymbol.h"
#include <sstream>

namespace voreen {

namespace glslparser {

LR0Item::LR0Item(const std::string& name, const Production& production, int position)
    : name_(name),
    production_(production),
    position_(position)
{
}

bool LR0Item::operator!=(const LR0Item& rhs) const {
    return ((name_ != rhs.name_) || (production_ != rhs.production_)
        || (position_ != rhs.position_));
}

bool LR0Item::operator==(const LR0Item& rhs) const {
    return ((name_ == rhs.name_) && (production_ == rhs.production_)
        && (position_ == rhs.position_));
}

GrammarSymbol* LR0Item::getLeftSymbol() const {
    if ( (position_ > 0) && (position_ <= static_cast<int>(production_.getBody().size())) )
        return production_.getBody()[position_ - 1];

    return 0;
}

Production LR0Item::getLeftString() const {
    return production_.getLeftPart(position_ - 1);
}

GrammarSymbol* LR0Item::getRightSymbol() const {
    if ( (position_ >= 0) && (position_ < static_cast<int>(production_.getBody().size())) )
        return production_.getBody()[position_];

    return 0;
}

Production LR0Item::getRightString() const {
    return production_.getRightPart(position_);
}

std::string LR0Item::toString() const {
    std::ostringstream oss;
    oss << "[" << name_ << " ::= ";
    const std::vector<GrammarSymbol*>& symbols = production_.getBody();
    for (size_t i = 0; i < symbols.size(); ++i) {
        if (i != 0)
            oss << " ";

        if (static_cast<int>(i) == position_)
            oss << ".";

        oss << ((symbols[i]->isTerminal() == true) ? symbols[i]->getName()
            : ("<" + symbols[i]->getName() + ">"));
    }

    if (position_ == static_cast<int>(symbols.size()))
        oss << ".";

    oss << "]";
    return oss.str();
}

// ============================================================================

LR1Item::LR1Item(const std::string& name, const Production& production, int position,
                 GrammarSymbol* const lookahead)
    : LR0Item(name, production, position),
    lookahead_(lookahead)
{
}

bool LR1Item::operator!=(const LR1Item& rhs) const {
    if ((lookahead_ == 0) || (rhs.lookahead_ == 0))
        return true;

    return (LR0Item::operator!=(rhs) || (*lookahead_ != *(rhs.lookahead_)));
}

bool LR1Item::operator==(const LR1Item& rhs) const {
    if ((lookahead_ == 0) || (rhs.lookahead_ == 0))
        return false;

    return (LR0Item::operator==(rhs) && (*lookahead_ == *(rhs.lookahead_)));
}

std::string LR1Item::toString() const {
    std::ostringstream oss;
    oss << "[" << name_ << " ::= ";
    const std::vector<GrammarSymbol*>& symbols = production_.getBody();
    for (size_t i = 0; i < symbols.size(); ++i) {
        if (i != 0)
            oss << " ";

        if (static_cast<int>(i) == position_)
            oss << ".";

        oss << ((symbols[i]->isTerminal() == true) ? symbols[i]->getName()
            : ("<" + symbols[i]->getName() + ">"));
    }

    if (position_ == static_cast<int>(symbols.size()))
        oss << ".";

    oss << ", " << lookahead_->getName() << "]";
    return oss.str();
}

// ============================================================================


template<typename T>
ItemSet<T>::ItemSet()
    : items_(), kernel_(0)
{
}

template<typename T>
ItemSet<T>::ItemSet(const std::vector<T>& items)
: items_(items), kernel_(0)
{
}

template<typename T>
ItemSet<T>::ItemSet(const ItemSet<T>& rhs) {
    clone(rhs);
}

template<typename T>
ItemSet<T>::~ItemSet() {
    delete kernel_;
}

template<typename T>
const ItemSet<LR0Item>* ItemSet<T>::getKernel() const {
    if (kernel_ == 0)
        const_cast<ItemSet<T>*>(this)->buildKernel();

    return kernel_;
}

template<>
const ItemSet<LR0Item>* ItemSet<LR0Item>::getKernel() const {
    return this;
}

template<typename T>
const ItemSet<T>& ItemSet<T>::operator=(const ItemSet<T>& rhs) {
    clone(rhs);
    return *this;
}

template<typename T>
bool ItemSet<T>::operator!=(const ItemSet<T>& rhs) const {
    if (rhs.getNumItems() != items_.size())
        return true;

    bool result = false;
    for (size_t i = 0; ((result == false) && (i < items_.size())); ++i)
        result = (items_[i] != rhs.items_[i]);

    return result;
}

template<typename T>
bool ItemSet<T>::operator==(const ItemSet<T>& rhs) const {
    if (rhs.getNumItems() != items_.size())
        return false;

    bool result = true;
    for (size_t i = 0; ((result == true) && (i < items_.size())); ++i)
        result = (items_[i] == rhs.items_[i]);

    return result;
}

template<typename T>
T& ItemSet<T>::operator[](size_t index) {
    if (index >= items_.size())
        throw std::runtime_error("ItemSet::operator[]: index out of bounds!");

    return items_[index];
}

template<typename T>
const T& ItemSet<T>::operator[](size_t index) const {
    if (index >= items_.size())
        throw std::runtime_error("ItemSet::operator[]: index out of bounds!");

    return items_[index];
}

template<typename T>
bool ItemSet<T>::addItem(const T& item) {
    if (containsItem(item) == false) {
        items_.push_back(item);
        return true;
    }
    return false;
}

template<typename T>
bool ItemSet<T>::contained(const std::list<ItemSet<T>* >& collection) const {
    return (findIn(collection) != collection.end());
}

template<typename T> typename std::list<ItemSet<T>* >::const_iterator
    ItemSet<T>::findIn(const std::list<ItemSet<T>* >& collection) const
{
    typename std::list<ItemSet<T>* >::const_iterator it = collection.begin();
    for ( ; (it != collection.end()); ++it) {
        if (*(*it) == *this)
            return it;
    }
    return collection.end();
}

template<typename T>
int ItemSet<T>::findPositionIn(const std::list<ItemSet<T>* >& collection) const
{
    typename std::list<ItemSet<T>* >::const_iterator it = collection.begin();
    for (int i = 0; it != collection.end(); ++it) {
        if (*(*it) == *this)
            return i;
        ++i;
    }
    return -1;
}

template<typename T>
int ItemSet<T>::findKernelsPositionIn(const std::list<ItemSet<T>* >& collection) const
{
    typename std::list<ItemSet<T>* >::const_iterator it = collection.begin();
    for (int i = 0; it != collection.end(); ++it) {
        if (*(getKernel()) == *((*it)->getKernel()))
            return i;
        ++i;
    }
    return -1;
}

template<typename T>
bool ItemSet<T>::hasNextSet() const {
    bool hasNextSet = false;
    for (size_t i = 0; ((! hasNextSet) && (i < items_.size())); ++i)
        hasNextSet = (items_[i].getRightSymbol() != 0);

    return hasNextSet;
}

template<typename T>
int ItemSet<T>::merge(const ItemSet<T>& other) {
    int numAdded = 0;

    // Insert all Items from the other ItemSet which are not
    // already in this one.
    //
    for (size_t i = 0; i < other.items_.size(); ++i) {
        if (! containsItem(other.items_[i])) {
            items_.push_back(other.items_[i]);
            ++numAdded;
        }
    }
    return numAdded;
}

template<typename T>
bool ItemSet<T>::containsItem(const T& ref) const {
    bool contained = false;
    for (size_t i = 0; ((contained == false) && (i < items_.size())); ++i)
        contained = (items_[i] == ref);
    return contained;
}

template<typename T>
std::string ItemSet<T>::toString() const {
    std::ostringstream oss;

    oss << "{";
    for (size_t i = 0; i < items_.size(); ++i) {
        if (i > 0)
            oss << ", \n";

        oss << items_[i].toString();
    }
    oss << "}";

    return oss.str();
}

// private methods
//

template<typename T>
void ItemSet<T>::clone(const ItemSet<T>& rhs) {
    items_ = rhs.items_;
    // Do not copy the pointer to the kernel or the kernel itself.
    // It will be re-created on demand by calling getKernel().
    //
    kernel_ = 0;
}

template<typename T>
void ItemSet<T>::buildKernel() {
    delete kernel_;
    kernel_ = new ItemSet<LR0Item>();

    for (size_t i = 0; i < items_.size(); ++i)
        kernel_->addItem(items_[i]);
}

template class ItemSet<LR0Item>;
template class ItemSet<LR1Item>;

}   // namespace glslparser

}   // namespace voreen
