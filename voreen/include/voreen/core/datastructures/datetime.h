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

#ifndef VRN_DATETIME_H
#define VRN_DATETIME_H

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/voreencoreapi.h"

#include <string>
#include <ostream>
#include <ctime>

namespace voreen {

/**
 * Stores a date value, a time value, or a combined date and time value.
 */
class VRN_CORE_API DateTime : public Serializable {

public:
    /// Specifies the data type stored in the DateTime object.
    enum Type {
        DATETIME,
        DATE,
        TIME
    };

    /// Creates a DateTime object of the passed type and initializes it to all zero values.
    DateTime(Type type = DATETIME);

    /// Constructs a DateTime object of type DATETIME and initializes it with the passed parameters.
    DateTime(int year, int month, int day, int hour, int min, int sec, int milliSecond);

    /// Constructs a DateTime object of type DATETIME and initializes it from the passed timestamp.
    DateTime(time_t timestamp);

    /// Creates a DateTime object of type DATE and initializes it with the passed parameters.
    static DateTime createDate(int year, int month, int day);

    /// Creates a DateTime object of type TIME and initializes it with the passed parameters.
    static DateTime createTime(int hour, int minute, int second, int milliSecond);

    /// Creates a DateTime object of type TIME from milliseconds.
    static DateTime createTimeFromMilliseconds(int numMilliseconds);

    friend std::ostream& operator<< (std::ostream& out, const DateTime &dateTime);
    std::string toString(bool withMillis = true) const;

    /// Returns the type of the DateTime object.
    Type getType() const;

    void setDateTime(int year, int month, int day, int hour, int min, int sec, int milliSecond);
    void setDate(int year, int month, int day);
    void setTime(int hour, int min, int sec, int milliSecond);

    int getYear() const;
    int getMonth() const;
    int getDay() const;
    int getHour() const;
    int getMinute() const;
    int getSecond() const;
    int getMillisecond() const;

    /**
     * Updates the date/time value from the passed timestamp.
     */
    void setTimestamp(time_t timestamp);

    /**
     * Constructs a timestamp from the date/time value.
     *
     * @note The timestamp does not contain milliseconds
     */
    time_t getTimestamp() const;

    /**
     * Returns the difference to dt in seconds (the return value may be negative)
     */
    double diffSeconds(const DateTime& dt) const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    /**
     * Returns a DateTime object storing the current calender time.
     */
    static DateTime now();

    bool operator==(const DateTime& dt) const;

    bool operator>(const DateTime& dt) const;

    bool operator<(const DateTime& dt) const;

protected:
    /// Checks all values and logs a warning, if out of range
    void checkValueRanges() const;

private:
    Type type_;

    int year_;
    int month_;
    int day_;
    int hour_;
    int minute_;
    int second_;
    int millisecond_;

    static const std::string loggerCat_;
};

} // namespace

#endif
