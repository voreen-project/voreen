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

#ifndef VRN_PLOTPREDICATE_H
#define VRN_PLOTPREDICATE_H

#include "voreen/core/io/serialization/serialization.h"
#include "interval.h"
#include "plotcell.h"

namespace voreen {

/**
 * Abstract super class for predicate classes used for selecting subsets out of PlotData tables.
 *
 * Subclasses have to override check function applying the predicate at all given PlotCells and
 * returning a boolean indicating whether all PlotCells match the predicate or not.
 */
class PlotPredicate : public AbstractSerializable  {
friend class PlotPredicateFactory;
public:
    /// Destructor
    virtual ~PlotPredicate() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const = 0;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds) = 0;
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const = 0;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const = 0;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    virtual bool check(const PlotCellValue& value) const = 0;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const = 0;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const = 0;

    /// returns a string representation of this predicate
    virtual std::string toString() const = 0;
};

// - PlotPredicateLess ----------------------------------------------------------------------------

/**
 * PlotPredicate checking PlotCells for being less than threshold given at instantiation.
 */
class PlotPredicateLess : public PlotPredicate {
friend class PlotPredicateFactory;

public:
    /**
     * Constructor: Creates a new PlotPredicate checking PlotCells for being less than \a threshold.
     *
     * \param   threshold   threshold of this PlotPredicate
     */
    PlotPredicateLess(plot_t threshold);

    /**
     * Constructor: Creates a new PlotPredicate checking PlotCells for being less than \a threshold.
     *
     * \param   threshold   threshold of this PlotPredicate
     */
    PlotPredicateLess(std::string threshold);

    /**
     * Constructor: Creates a new PlotPredicate checking PlotCells for being less than \a threshold.
     *
     * \param   threshold   threshold of this PlotPredicate
     */
    PlotPredicateLess(const PlotCellValue& threshold);

    /// Destructor
    virtual ~PlotPredicateLess() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds);
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    bool check(const PlotCellValue& value) const;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const;

    /// returns a string representation of this predicate
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

private:
    PlotPredicateLess();            ///< private default constructor only for PlotPredicateFactory
    PlotCellValue threshold_;       ///< threshold againts the given PlotCells shall be checked
};

// - PlotPredicateEqual --------------------------------------------------------------------------

/**
 * PlotPredicate checking PlotCells for being equal to threshold given at instantiation.
 */
class PlotPredicateEqual : public PlotPredicate {
friend class PlotPredicateFactory;

public:
    /**
     * Constructor: creates a new PlotPredicate checking PlotCells for being equal to threshold given at instantiation.
     *
     * \param   threshold   threshold of this PlotPredicate
     */
    PlotPredicateEqual(plot_t threshold);

    /**
     * Constructor: creates a new PlotPredicate checking PlotCells for being equal to threshold given at instantiation.
     *
     * \param   threshold   threshold of this PlotPredicate
     */
    PlotPredicateEqual(std::string threshold);

    /**
     * Constructor: creates a new PlotPredicate checking PlotCells for being equal to threshold given at instantiation.
     *
     * \param   threshold   threshold of this PlotPredicate
     */
    PlotPredicateEqual(const PlotCellValue& threshold);

    /// Destructor
    virtual ~PlotPredicateEqual() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds);
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    bool check(const PlotCellValue& value) const;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const;

    /// returns a string representation of this predicate
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;
    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

private:
    PlotPredicateEqual();       ///< private default constructor only for PlotPredicateFactory
    PlotCellValue threshold_;   ///< threshold againts the given PlotCells shall be checked
};

// - PlotPredicateGreater -------------------------------------------------------------------------

/**
 * PlotPredicate checking PlotCells for being greater than threshold given at instantiation.
 */
class PlotPredicateGreater : public PlotPredicate {
friend class PlotPredicateFactory;

public:
    /**
     * Constructor: creates a new PlotPredicate checking PlotCells for being greater than threshold given at instantiation.
     *
     * \param   threshold   threshold of this PlotPredicate
     */
    PlotPredicateGreater(plot_t threshold);

    /**
     * Constructor: creates a new PlotPredicate checking PlotCells for being greater than threshold given at instantiation.
     *
     * \param   threshold   threshold of this PlotPredicate
     */
    PlotPredicateGreater(std::string threshold);

    /**
     * Constructor: creates a new PlotPredicate checking PlotCells for being greater than threshold given at instantiation.
     *
     * \param   threshold   threshold of this PlotPredicate
     */
    PlotPredicateGreater(const PlotCellValue& threshold);

    /// Destructor
    virtual ~PlotPredicateGreater() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds);
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    bool check(const PlotCellValue& value) const;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const;

    /// returns a string representation of this predicate
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

private:
    PlotPredicateGreater();         ///< private default constructor only for PlotPredicateFactory
    PlotCellValue threshold_;       ///< threshold againts the given PlotCells shall be checked
};

// - PlotPredicateBetween -------------------------------------------------------------------------

/**
 * PlotPredicate checking PlotCells for being between the two thresholds given at instantiation.
 */
class PlotPredicateBetween : public PlotPredicate {
friend class PlotPredicateFactory;

public:
    /**
     * Constructor: creates a new PlotPredicate checking PlotCells for being between \a lowerThreshold and \a upperThreshold
     *
     * \param lowerThreshold    lower threshold againts the given PlotCells shall be checked
     * \param upperThreshold    upper threshold againts the given PlotCells shall be checked
     */
    PlotPredicateBetween(plot_t lowerThreshold, plot_t upperThreshold);

    /**
     * Constructor: creates a new PlotPredicate checking PlotCells for being between \a lowerThreshold and \a upperThreshold
     *
     * \param lowerThreshold    lower threshold againts the given PlotCells shall be checked
     * \param upperThreshold    upper threshold againts the given PlotCells shall be checked
     */
    PlotPredicateBetween(std::string lowerThreshold, std::string upperThreshold);

    /**
     * Constructor: creates a new PlotPredicate checking PlotCells for being between \a lowerThreshold and \a upperThreshold
     *
     * \param lowerThreshold    lower threshold againts the given PlotCells shall be checked
     * \param upperThreshold    upper threshold againts the given PlotCells shall be checked
     */
    PlotPredicateBetween(const PlotCellValue& lowerThreshold, const PlotCellValue& upperThreshold);

    /// Destructor
    virtual ~PlotPredicateBetween() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds);
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    bool check(const PlotCellValue& value) const;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const;

    /// returns a string representation of this predicate
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;
    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

private:
    PlotPredicateBetween();         ///< private default constructor only for PlotPredicateFactory
    PlotCellValue lowerThreshold_;  ///< lower threshold againts the given PlotCells shall be checked
    PlotCellValue upperThreshold_;  ///< upper threshold againts the given PlotCells shall be checked
};

// - PlotPredicateNotBetween -------------------------------------------------------------------------

/**
 * PlotPredicate checking PlotCells for being not between the two thresholds given at instantiation.
 */
class PlotPredicateNotBetween : public PlotPredicate {
friend class PlotPredicateFactory;

public:
    /**
     * Constructor: creates a new PlotPredicate checking PlotCells for being not between \a lowerThreshold and \a upperThreshold
     *
     * \param lowerThreshold    lower threshold againts the given PlotCells shall be checked
     * \param upperThreshold    upper threshold againts the given PlotCells shall be checked
     */
    PlotPredicateNotBetween(plot_t lowerThreshold, plot_t upperThreshold);

    /**
     * Constructor: creates a new PlotPredicate checking PlotCells for being not between \a lowerThreshold and \a upperThreshold
     *
     * \param lowerThreshold    lower threshold againts the given PlotCells shall be checked
     * \param upperThreshold    upper threshold againts the given PlotCells shall be checked
     */
    PlotPredicateNotBetween(std::string lowerThreshold, std::string upperThreshold);

    /**
     * Constructor: creates a new PlotPredicate checking PlotCells for being not between \a lowerThreshold and \a upperThreshold
     *
     * \param lowerThreshold    lower threshold againts the given PlotCells shall be checked
     * \param upperThreshold    upper threshold againts the given PlotCells shall be checked
     */
    PlotPredicateNotBetween(const PlotCellValue& lowerThreshold, const PlotCellValue& upperThreshold);

    /// Destructor
    virtual ~PlotPredicateNotBetween() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds);
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    bool check(const PlotCellValue& value) const;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const;

    /// returns a string representation of this predicate
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;
    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

private:
    PlotPredicateNotBetween();      ///< private default constructor only for PlotPredicateFactory
    PlotCellValue lowerThreshold_;  ///< lower threshold againts the given PlotCells shall be checked
    PlotCellValue upperThreshold_;  ///< upper threshold againts the given PlotCells shall be checked
};

// - PlotPredicateBetweenOrEqual ------------------------------------------------------------------

/**
 * PlotPredicate checking PlotCells for being between the two thresholds given at instantiation or equal to one of them.
 */
class PlotPredicateBetweenOrEqual : public PlotPredicate {
friend class PlotPredicateFactory;

public:
    /**
     * Constructor: creates a new PlotPredicate checking whether a given value is between \a lowerThreshold and \a upperThreshold or equals one of them.
     *
     * \param lowerThreshold    lower threshold againts the given PlotCells shall be checked
     * \param upperThreshold    upper threshold againts the given PlotCells shall be checked
     */
    PlotPredicateBetweenOrEqual(plot_t lowerThreshold, plot_t upperThreshold);

    /**
     * Constructor: creates a new PlotPredicate checking whether a given value is between \a lowerThreshold and \a upperThreshold or equals one of them.
     *
     * \param lowerThreshold    lower threshold againts the given PlotCells shall be checked
     * \param upperThreshold    upper threshold againts the given PlotCells shall be checked
     */
    PlotPredicateBetweenOrEqual(std::string lowerThreshold, std::string upperThreshold);

    /**
     * Constructor: creates a new PlotPredicate checking whether a given value is between \a lowerThreshold and \a upperThreshold or equals one of them.
     *
     * \param lowerThreshold    lower threshold againts the given PlotCells shall be checked
     * \param upperThreshold    upper threshold againts the given PlotCells shall be checked
     */
    PlotPredicateBetweenOrEqual(const PlotCellValue& lowerThreshold, const PlotCellValue& upperThreshold);

    /// Destructor
    virtual ~PlotPredicateBetweenOrEqual() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds);
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    virtual bool check(const PlotCellValue& value) const;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const;

    /// returns a string representation of this predicate
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;
    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

private:
    PlotPredicateBetweenOrEqual();  ///< private default constructor only for PlotPredicateFactory
    PlotCellValue lowerThreshold_;  ///< lower threshold againts the given PlotCells shall be checked
    PlotCellValue upperThreshold_;  ///< upper threshold againts the given PlotCells shall be checked
};

// - PlotPredicateNotBetweenOrEqual ------------------------------------------------------------------

/**
 * PlotPredicate checking PlotCells for being neither between the two thresholds given at instantiation nor equal to one of them.
 */
class PlotPredicateNotBetweenOrEqual : public PlotPredicate {
friend class PlotPredicateFactory;

public:
    /**
     * Constructor: creates a new PlotPredicate checking whether a given value is neither between \a lowerThreshold and \a upperThreshold nor equals one of them.
     *
     * \param lowerThreshold    lower threshold againts the given PlotCells shall be checked
     * \param upperThreshold    upper threshold againts the given PlotCells shall be checked
     */
    PlotPredicateNotBetweenOrEqual(plot_t lowerThreshold, plot_t upperThreshold);

    /**
     * Constructor: creates a new PlotPredicate checking whether a given value is neither between \a lowerThreshold and \a upperThreshold nor equals one of them.
     *
     * \param lowerThreshold    lower threshold againts the given PlotCells shall be checked
     * \param upperThreshold    upper threshold againts the given PlotCells shall be checked
     */
    PlotPredicateNotBetweenOrEqual(std::string lowerThreshold, std::string upperThreshold);

    /**
     * Constructor: creates a new PlotPredicate checking whether a given value is neither between \a lowerThreshold and \a upperThreshold nor equals one of them.
     *
     * \param lowerThreshold    lower threshold againts the given PlotCells shall be checked
     * \param upperThreshold    upper threshold againts the given PlotCells shall be checked
     */
    PlotPredicateNotBetweenOrEqual(const PlotCellValue& lowerThreshold, const PlotCellValue& upperThreshold);

    /// Destructor
    virtual ~PlotPredicateNotBetweenOrEqual() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds);
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    virtual bool check(const PlotCellValue& value) const;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const;

    /// returns a string representation of this predicate
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;
    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

private:
    PlotPredicateNotBetweenOrEqual();  ///< private default constructor only for PlotPredicateFactory
    PlotCellValue lowerThreshold_;  ///< lower threshold againts the given PlotCells shall be checked
    PlotCellValue upperThreshold_;  ///< upper threshold againts the given PlotCells shall be checked
};

// - PlotPredicateIsSubStr ------------------------------------------------------------------------

/**
 * PlotPredicate checking PlotCells for containing the string given at instantiation as substring.
 *
 * \note    Returns false for all PlotCells not containing a tag.
 */
class PlotPredicateIsSubStr : public PlotPredicate {
friend class PlotPredicateFactory;
public:
    /**
     * Constructor: Creates a new PlotPredicate checking PlotCells for containing \a substring as substring.
     *
     * \param   substring   string which shall be searched in given PlotCells
     */
    PlotPredicateIsSubStr(std::string threshold);

    /// Destructor
    virtual ~PlotPredicateIsSubStr() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds);
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    bool check(const PlotCellValue& value) const;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const;

    /// returns a string representation of this predicate
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;
    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

private:
    PlotPredicateIsSubStr();    ///< private default constructor only for PlotPredicateFactory
    std::string threshold_;     ///< threshold againts the given PlotCells shall be checked
};

// - PlotPredicateNotAlphaNumeric -----------------------------------------------------------------

/**
 * PlotPredicate checking PlotCells for being not alphanumeric.
 */
class PlotPredicateNotAlphaNumeric : public PlotPredicate {
friend class PlotPredicateFactory;
public:
    /// Destructor
    virtual ~PlotPredicateNotAlphaNumeric() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds);
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    bool check(const PlotCellValue& value) const;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const;

    /// returns a string representation of this predicate
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;
    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);
};

// - PlotPredicateAlphaNumeric --------------------------------------------------------------------

/**
 * PlotPredicate checking PlotCells for being alphanumeric.
 */
class PlotPredicateAlphaNumeric : public PlotPredicate {
friend class PlotPredicateFactory;
public:
    /// Destructor
    virtual ~PlotPredicateAlphaNumeric() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds);
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    bool check(const PlotCellValue& value) const;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const;

    /// returns a string representation of this predicate
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;
    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);
};

// - PlotPredicateEmpty ---------------------------------------------------------------------------

/**
 * PlotPredicate checking PlotCells for being Null.
 */
class PlotPredicateEmpty : public PlotPredicate {
friend class PlotPredicateFactory;
public:
    /// Destructor
    virtual ~PlotPredicateEmpty() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds);
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    bool check(const PlotCellValue& value) const;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const;

    /// returns a string representation of this predicate
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;
    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);
};

// - PlotPredicateNotEmpty ------------------------------------------------------------------------

/**
 * PlotPredicate checking PlotCells for being not Null.
 */
class PlotPredicateNotEmpty : public PlotPredicate {
friend class PlotPredicateFactory;
public:
    /// Destructor
    virtual ~PlotPredicateNotEmpty() {};

    /// Returns the number of stored thresholds.
    virtual int getNumberOfThresholdValues() const;
    /// Sets the stored threshold values to the values in \a thresholds
    virtual void setThresholdValues(const std::vector<PlotCellValue>& thresholds);
    /// Returns a vector of the stored thresholds.
    virtual std::vector<PlotCellValue> getThresholdValues() const;
    /// Returns a vector of descriptions for the stored thresholds
    virtual std::vector<std::string> getThresholdTitles() const;

    /// checks whether value stored in PlotCell \a value fulfills the predicate
    bool check(const PlotCellValue& value) const;

    /// Returns an interval representation of the PlotPredicate if possible, non numeric predicates return an empty interval.
    virtual Interval<plot_t> getIntervalRepresentation() const;

    /// creates a deep copy of the current PlotPredicate
    virtual PlotPredicate* clone() const;

    /// returns a string representation of this predicate
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;
    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};


} // namespace voreen

#endif // VRN_PLOTPREDICATE_H
