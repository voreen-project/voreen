/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/core/datastructures/datetime.h"

#include "voreen/core/io/serialization/serialization.h"

#include <sstream>
#include <cstdio>

namespace voreen {

const std::string DateTime::loggerCat_ = "voreen.DateTime";

DateTime::DateTime(Type type /*= DATETIME*/)
    : type_(type)
{
    setDateTime(0, 0, 0, 0, 0, 0, 0);
}

DateTime::DateTime(int year, int month, int day, int hour, int min, int sec, int milliSecond)
    : type_(DATETIME)
{
    setDateTime(year, month, day, hour, min, sec, milliSecond);
}

DateTime::DateTime(time_t timestamp)
    : type_(DATETIME)
{
    setTimestamp(timestamp);
}

DateTime DateTime::createDate(int year, int month, int day) {
    DateTime dateTime(DATE);
    dateTime.setDate(year, month, day);
    return dateTime;
}

DateTime DateTime::createTime(int hour, int minute, int second, int milliSecond) {
    DateTime dateTime(TIME);
    dateTime.setTime(hour, minute, second, milliSecond);
    return dateTime;
}

DateTime DateTime::createTimeFromMilliseconds(int numMilliseconds) {
    int temp = numMilliseconds;

    int millis = temp % 1000;
    temp /= 1000;
    int seconds = temp % 60;
    temp /= 60;
    int minutes = temp % 60;
    temp /= 60;
    int hours = temp % 24;

    if (temp > 23)
        LWARNING("createTimeFromMilliseconds(): number of milliseconds larger than one day: " << numMilliseconds);

    return createTime(hours, minutes, seconds, millis);
}

std::ostream& operator<<(std::ostream& os, const DateTime& dateTime) {
    return os << dateTime.toString();
}

std::string DateTime::toString(bool withMillis) const {
    char* charString = new char[30];
    if (type_ == DATETIME) {
        if (withMillis)
            sprintf(charString, "%04d-%02d-%02d %02d:%02d:%02d,%03d",
                getYear(), getMonth(), getDay(), getHour(), getMinute(), getSecond(), getMillisecond());
        else
            sprintf(charString, "%04d-%02d-%02d %02d:%02d:%02d",
                getYear(), getMonth(), getDay(), getHour(), getMinute(), getSecond());
    }
    else if (type_ == DATE) {
        sprintf(charString, "%04d-%02d-%02d", getYear(), getMonth(), getDay());
    }
    else if (type_ == TIME) {
        if (withMillis)
            sprintf(charString, "%02d:%02d:%02d,%03d", getHour(), getMinute(), getSecond(), getMillisecond());
        else
            sprintf(charString, "%02d:%02d:%02d", getHour(), getMinute(), getSecond());
    }
    else {
        LERRORC("voreen.DateTime", "invalid type");
    }

    std::string result(charString);
    delete[] charString;
    return result;
}

DateTime::Type DateTime::getType() const {
    return type_;
}

void DateTime::setDateTime(int year, int month, int day, int hour, int min, int sec, int milliSecond) {
    year_ = year;
    month_ = month;
    day_ = day;
    hour_ = hour;
    minute_ = min;
    second_ = sec;
    millisecond_ = milliSecond;

    checkValueRanges();
}

void DateTime::setDate(int year, int month, int day) {
    if (type_ == TIME) {
        LWARNING("setDate() used on instance of type TIME");
    }
    year_ = year;
    month_ = month;
    day_ = day;

    checkValueRanges();
}

void DateTime::setTime(int hour, int min, int sec, int milliSecond) {
    if (type_ == DATE) {
        LWARNING("setDate() used on instance of type DATE");
    }
    hour_ = hour;
    minute_ = min;
    second_ = sec;
    millisecond_ = milliSecond;

    checkValueRanges();
}

int DateTime::getYear() const {
    return year_;
}

int DateTime::getMonth() const  {
    return month_;
}

int DateTime::getDay() const {
    return day_;
}

int DateTime::getHour() const {
    return hour_;
}

int DateTime::getMinute() const {
    return minute_;
}

int DateTime::getSecond() const {
    return second_;
}

int DateTime::getMillisecond() const {
    return millisecond_;
}

void DateTime::setTimestamp(time_t timestamp) {
    tm* t = gmtime(&timestamp);
    setDateTime(
        t->tm_year + 1900,
        t->tm_mon + 1,
        t->tm_mday,
        t->tm_hour,
        t->tm_min,
        t->tm_sec,
        0);
}

time_t DateTime::getTimestamp() const {
    tm dateTime;
    dateTime.tm_year = year_ - 1900;
    dateTime.tm_mon = month_ - 1;
    dateTime.tm_mday = day_;
    dateTime.tm_hour = hour_;
    dateTime.tm_min = minute_;
    dateTime.tm_sec = second_;

    return mktime(&dateTime);
}

void DateTime::serialize(Serializer& s) const {
    s.serialize("datetype", type_); //< "type" is a reserved keyword

    s.serialize("year", year_);
    s.serialize("month", month_);
    s.serialize("day", day_);
    s.serialize("hour", hour_);
    s.serialize("minute", minute_);
    s.serialize("second", second_);
    s.serialize("millisecond", millisecond_);
}

void DateTime::deserialize(Deserializer& s) {
    int typeInt = 0;
    s.optionalDeserialize<int>("datetype", typeInt, DATETIME);
    type_ = Type(typeInt);

    s.deserialize("year", year_);
    s.deserialize("month", month_);
    s.deserialize("day", day_);
    s.deserialize("hour", hour_);
    s.deserialize("minute", minute_);
    s.deserialize("second", second_);
    s.deserialize("millisecond", millisecond_);
}

DateTime DateTime::now() {
    return DateTime(time(0));
}

void DateTime::checkValueRanges() const {
    if (year_ < 0 || year_ > 2100)
        LWARNING("year out of valid range [0;2100]: " << year_);
    if (month_ < 0 || month_ > 12)
        LWARNING("month out of valid range [0;12]: " << month_);
    if (day_ < 0 || day_ > 31)
        LWARNING("day out of valid range [0;31]: " << day_);
    if (hour_ < 0 || hour_ > 23)
        LWARNING("hour out of valid range [0;23]: " << hour_);
    if (minute_ < 0 || minute_ > 59)
        LWARNING("minute out of valid range [0;59]: " << minute_);
    if (second_ < 0 || second_ > 59)
        LWARNING("second out of valid range [0;59]: " << second_);
    if (millisecond_ < 0 || millisecond_ > 999)
        LWARNING("millisecond out of valid range [0;999]: " << millisecond_);
}

bool DateTime::operator==(const DateTime& dt) const {
    if ((year_) != dt.getYear())
        return false;
    if ((month_) != dt.getMonth())
        return false;
    if ((day_) != dt.getDay())
        return false;
    if ((hour_) != dt.getHour())
        return false;
    if ((minute_) != dt.getMinute())
        return false;
    if ((second_) != dt.getSecond())
        return false;
    if(millisecond_ != dt.getMillisecond())
        return false;
    return true;
}

bool DateTime::operator>(const DateTime& dt) const {
    if (year_ > dt.getYear())
        return true;
    else if (year_ < dt.getYear())
        return false;

    if (month_ > dt.getMonth())
        return true;
    else if (month_ < dt.getMonth())
        return false;

    if (day_ > dt.getDay())
        return true;
    else if (day_ < dt.getDay())
        return false;

    if (hour_ > dt.getHour())
        return true;
    else if (hour_ < dt.getHour())
        return false;

    if (minute_ > dt.getMinute())
        return true;
    else if (minute_ < dt.getMinute())
        return false;

    if (second_ > dt.getSecond())
        return true;
    else if (second_ < dt.getSecond())
        return false;

    if (millisecond_ > dt.getMillisecond())
        return true;
    else if (millisecond_ < dt.getMillisecond())
        return false;

    return false;
}

bool DateTime::operator<(const DateTime& dt) const {
    if (year_ > dt.getYear())
        return false;
    else if (year_ < dt.getYear())
        return true;

    if (month_ > dt.getMonth())
        return false;
    else if (month_ < dt.getMonth())
        return true;

    if (day_ > dt.getDay())
        return false;
    else if (day_ < dt.getDay())
        return true;

    if (hour_ > dt.getHour())
        return false;
    else if (hour_ < dt.getHour())
        return true;

    if (minute_ > dt.getMinute())
        return false;
    else if (minute_ < dt.getMinute())
        return true;

    if (second_ > dt.getSecond())
        return false;
    else if (second_ < dt.getSecond())
        return true;

    if (millisecond_ > dt.getMillisecond())
        return false;
    else if (millisecond_ < dt.getMillisecond())
        return true;

    return false;
}

double DateTime::diffSeconds(const DateTime& dt) const {
    return std::difftime(dt.getTimestamp(), getTimestamp());
}

} // namespace
