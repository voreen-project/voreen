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

#ifndef VRN_PLOTPREDICATEPROPERTY_H
#define VRN_PLOTPREDICATEPROPERTY_H

#include "tgt/vector.h"
#include "tgt/logmanager.h"
#include "voreen/core/properties/templateproperty.h"

#include <vector>
#include <string>
#include <utility>

namespace voreen {

// forward declaration
class PlotPredicate;

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API TemplateProperty<std::vector<std::pair<int, PlotPredicate*> > >;
#endif

/**
 * \brief   Property encapsulating a buch of PlotPredicates together with the
 *          column index to which apply the predicate to.
 *
 *  \note   This property makes deep copies of all given PlotPredicates.
 **/
class VRN_CORE_API PlotPredicateProperty : public TemplateProperty<std::vector<std::pair<int, PlotPredicate*> > > {
public:
    /// Creates an empty PlotPredicateProperty
    PlotPredicateProperty(const std::string& id, const std::string& guiText,
             Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT);
    PlotPredicateProperty();
    virtual ~PlotPredicateProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "PlotPredicateProperty"; }
    virtual std::string getTypeDescription() const { return "PlotPredicate"; }

    /**
     * \brief   Adds predicate \a pred associated with column \a column to this property.
     *
     * \param   column  column index to which \a pred shall later be applied
     * \param   pred    PlotPredicate to add
     **/
    template <class T>
    void addPredicate(int column, const T& pred) {
        T* p = new T(pred);
        value_.push_back(std::make_pair(column, p));
        notifyAll();
    };

    /// clears this property
    void clear();

    /// returns all held PlotPredicates with its column indices
    const std::vector<std::pair<int, PlotPredicate*> >& get() const;

    /// returns a deep copy of all all held PlotPredicates with its column indices
    std::vector<std::pair<int,PlotPredicate*> > getCloned() const;

    /**
     * \brief   Sets \a value as value of this Property
     *
     * \note    Makes deep copies of all PlotPredicates in \a value, so they can safely be
     *          deleted at any time.
     **/
    void set(const std::vector<std::pair<int, PlotPredicate*> >& value);

    /// Returns the number of PlotPredicates in this property.
    int size() const;

    /**
     * \brief   change the Order of the Predicates and delete not selected predicates
     *
     * \param columns the values in the vector are the old column Position and the position
     *          in this vector is the new position.
     **/
    void select(const std::vector<int>& columns);

    /**
     * @see Property::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Property::deserialize
     */
    virtual void deserialize(Deserializer& s);

protected:
    void notifyAll();

    static const std::string loggerCat_;

};

}   // namespace

#endif //VRN_PLOTPREDICATEPROPERTY_H
