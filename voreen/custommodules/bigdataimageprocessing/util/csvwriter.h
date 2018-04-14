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

#ifndef VRN_CSVWRITER_H
#define VRN_CSVWRITER_H

#include <fstream>

#include <tgt/exception.h>

namespace voreen {

enum class CSVWriteMode {
    TRUNCATE,
    APPEND
};

template<typename... ColumnTypes>
class CSVWriter {
public:
    //CSVWriter(CSVWriter<ColumnTypes...>&&);
    //CSVWriter(std::ofstream&& file, const std::string& separator = ";");
    CSVWriter(const std::string& filePath, char separator = ';', CSVWriteMode writeMode = CSVWriteMode::TRUNCATE);

    template<typename... T, class=typename std::enable_if<sizeof...(ColumnTypes) == sizeof...(T)>::type>
    void writeHeader(T... names);

    void write(ColumnTypes... columns);

private:
    template<typename ColumnType>
    void writeRecursive(ColumnType currentColumn);

    template<typename ColumnType, typename... NextColumnTypes>
    void writeRecursive(ColumnType currentColumn, NextColumnTypes... nextColumns);

    template<typename T>
    void writeElement(T element);
    void writeElement(const std::string& element);
    void writeElement(const char* element);

    /*
    void writeHeaderRecursive(const std::string& currentColumn);
    template<typename... NextColumnTypes>
    void writeHeaderRecursive(const std::string& currentColumn, NextColumnTypes... nextColumns);
    */

    void writeString(const std::string&);

    std::ofstream file_;
    char separator_;
};

/*
template<typename... ColumnTypes>
CSVWriter<ColumnTypes...>::CSVWriter(CSVWriter<ColumnTypes...>&& other)
    : file_(std::move(other.file_))
    , separator_(other.separator_)
{
}
*/


/*
template<typename... ColumnTypes>
CSVWriter<ColumnTypes...>::CSVWriter(std::ofstream&& file, const std::string& separator)
    : file_(std::move(file))
    , separator_(separator)
{
}
template<typename... ColumnTypes>
CSVWriter<ColumnTypes...>::CSVWriter(const std::string& filePath, const std::string& separator)
    : CSVWriter<ColumnTypes...>(std::move(std::ofstream(filePath, std::ofstream::trunc)), separator)
{
}
*/
template<typename... ColumnTypes>
CSVWriter<ColumnTypes...>::CSVWriter(const std::string& filePath, char separator, CSVWriteMode writeMode) 
    : file_(filePath, writeMode == CSVWriteMode::TRUNCATE ? std::ofstream::trunc : std::ofstream::app)
    , separator_(separator)
{
    if(file_.fail()) {
        throw tgt::IOException("Could not create output file.");
    }
}

template<typename... ColumnTypes>
void CSVWriter<ColumnTypes...>::write(ColumnTypes... columns) {
    writeRecursive(columns...);
}


template<typename... ColumnTypes>
template<typename... T, class>
void CSVWriter<ColumnTypes...>::writeHeader(T... names) {
    writeRecursive(names...);
}


template<typename... ColumnTypes>
template<typename ColumnType>
void CSVWriter<ColumnTypes...>::writeRecursive(ColumnType currentColumn) {
    writeElement(currentColumn);
    file_ << "\n";
}


template<typename... ColumnTypes>
template<typename ColumnType, typename... NextColumnTypes>
void CSVWriter<ColumnTypes...>::writeRecursive(ColumnType currentColumn, NextColumnTypes... nextColumns) {
    writeElement(currentColumn);
    file_ << separator_;
    writeRecursive(nextColumns...);
}

template<typename... ColumnTypes>
template<typename T>
void CSVWriter<ColumnTypes...>::writeElement(T element) {
    file_ << element;
    // Can we assume that element does not need escaping? hmm...
    // For now the user is responsible to make sure that custom types
    // do not print characters that mess up the format
}

template<typename... ColumnTypes>
void CSVWriter<ColumnTypes...>::writeElement(const std::string& element) {
    writeString(element);
}

template<typename... ColumnTypes>
void CSVWriter<ColumnTypes...>::writeElement(const char* element) {
    writeString(element);
}

template<typename... ColumnTypes>
void CSVWriter<ColumnTypes...>::writeString(const std::string& str) {

    bool needsQuoting = false;
    for(char c : str) {
        if(c == '"' || c == separator_ || c == '\r' || c == '\n') {
            needsQuoting = true;
            break;
        }
    }

    if(needsQuoting) {
        // If the string contains character that need to be quoted, we surround it in quotes
        file_ << '"';
        for(char c : str) {
            if(c == '"') {
                file_ << '"'; //double quotes are escaped by preceded by another double quote
            }
            file_ << c;
        }
        file_ << '"';
    } else {
        file_ << str;
    }
}

} // namespace voreen

#endif // VRN_CSVWRITER_H
