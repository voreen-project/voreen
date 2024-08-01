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

#ifndef VRN_AGGREGATIONFUNCTION_H
#define VRN_AGGREGATIONFUNCTION_H

#include "plotbase.h"
#include "voreen/core/voreenobject.h"

namespace voreen {

/**
 * \brief Abstract super class for aggregation functions which can be applied to the data.
 *
 * Subclasses have to override virtual evaluate function which applies the aggregation function
 * to a bunch of given plot_t.
 */
class VRN_CORE_API AggregationFunction : public VoreenSerializableObject {
friend class AggregationFunctionFactory;
public:

    enum AggregationFunctionEnum {
        NONE = 0,
        COUNT = 1,
        COUNTHISTOGRAM = 2,
        MODE = 3,
        MIN = 20,
        MINHISTOGRAM = 21,
        MAX = 30,
        MAXHISTOGRAM = 31,
        SUM = 40,
        SUMHISTOGRAM = 41,
        MEDIAN = 50,
        MEDIANHISTOGRAM = 51,
        AVERAGE = 60,
        GEOMETRICMEAN = 61,
        HARMONICMEAN = 62,
        STANDARDDEVIATION = 70,
        VARIANCE = 75
    };

    /// Destructor
    virtual ~AggregationFunction() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const = 0;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const = 0;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const = 0;

};

/**
 * \brief Aggregation function which finds the count of given plot_t values in Histograms.
 */
class VRN_CORE_API AggregationFunctionCountHistogram : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionCountHistogram";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionCountHistogram(); }

    /// Destructor
    virtual ~AggregationFunctionCountHistogram() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which finds the count of given plot_t values.
 */
class VRN_CORE_API AggregationFunctionCount : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionCount";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionCount(); }

    /// Destructor
    virtual ~AggregationFunctionCount() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which finds the minimum value in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionMinHistogram : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionMinHistogram";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionMinHistogram(); }

    /// Destructor
    virtual ~AggregationFunctionMinHistogram() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which finds the minimum value in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionMin : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionMin";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionMin(); }

    /// Destructor
    virtual ~AggregationFunctionMin() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which finds the maximum value in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionMaxHistogram : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionMaxHistogram";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionMaxHistogram(); }

    /// Destructor
    virtual ~AggregationFunctionMaxHistogram() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which finds the maximum value in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionMax : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionMax";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionMax(); }

    /// Destructor
    virtual ~AggregationFunctionMax() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which calculates the sum of all values in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionSumHistogram : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionSumHistogram";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionSumHistogram(); }

    /// Destructor
    virtual ~AggregationFunctionSumHistogram() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which calculates the sum of all values in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionSum : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionSum";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionSum(); }

    /// Destructor
    virtual ~AggregationFunctionSum() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which calculates the arithmetic average (mean) of all values in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionAverage : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionAverage";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionAverage(); }

    /// Destructor
    virtual ~AggregationFunctionAverage() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};
/**
 * \brief Aggregation function which calculates the geomatric average (mean) of all values in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionGeometricAverage : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionGeometricAverage";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionGeometricAverage(); }

    /// Destructor
    virtual ~AggregationFunctionGeometricAverage() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};
/**
 * \brief Aggregation function which calculates the harmonic average (mean) of all values in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionHarmonicAverage : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionHarmonicAverage";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionHarmonicAverage(); }

    /// Destructor
    virtual ~AggregationFunctionHarmonicAverage() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which calculates the Median of all values in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionMedian : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionMedian";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionMedian(); }

    /// Destructor
    virtual ~AggregationFunctionMedian() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which calculates the Median of all values in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionMedianHistogram : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionMedianHistogram";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionMedianHistogram(); }

    /// Destructor
    virtual ~AggregationFunctionMedianHistogram() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which calculates the standard deviation of all values in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionStandardDeviation : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionStandardDeviation";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionStandardDeviation(); }

    /// Destructor
    virtual ~AggregationFunctionStandardDeviation() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which calculates the Variance of all values in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionVariance : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionVariance";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionVariance(); }

    /// Destructor
    virtual ~AggregationFunctionVariance() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

/**
 * \brief Aggregation function which calculates the Mode of all values in given vector<plot_t>.
 */
class VRN_CORE_API AggregationFunctionMode : public AggregationFunction {
friend class AggregationFunctionFactory;
public:
    virtual std::string getClassName() const    { return "AggregationFunctionMode";     }
    virtual AggregationFunction* create() const { return new AggregationFunctionMode(); }

    /// Destructor
    virtual ~AggregationFunctionMode() {};

    /// Applies the aggregation function to values \a values.
    virtual plot_t evaluate(std::vector<plot_t>& values) const;

    /// Create a copy of the actuell AggregationFunction.
    virtual AggregationFunction* clone() const;

    /// Returns a string representation of this AggregationFunction.
    virtual std::string toString() const;

    /// @see Serializer::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Deserializer::deserialize
    virtual void deserialize(Deserializer& s);

};

} // namespace voreen

#endif // VRN_AGGREGATIONFUNCTION_H
