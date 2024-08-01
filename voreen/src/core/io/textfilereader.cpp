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

#include "voreen/core/io/textfilereader.h"

#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

namespace voreen {

namespace {
const string DEFAULT_WHITESPACE = " \t\n\r";
const string DEFAULT_COMMENTCHARS = "#";
const string DEFAULT_SEPARATORS = DEFAULT_WHITESPACE;
}

TextFileReader::TextFileReader(const string& filename)
    : whitespace_(DEFAULT_WHITESPACE), commentChars_(DEFAULT_COMMENTCHARS),
      separators_(DEFAULT_SEPARATORS), isMyOwnStream_(true),
      file_(new ifstream(filename.c_str()))
{}

TextFileReader::TextFileReader(istream* stream)
    : whitespace_(DEFAULT_WHITESPACE), commentChars_(DEFAULT_COMMENTCHARS),
      separators_(DEFAULT_SEPARATORS), isMyOwnStream_(false),
      file_(stream)
{}

TextFileReader::~TextFileReader() {
    if (isMyOwnStream_)
        delete file_;
}

bool TextFileReader::getNextLine(string& type, string& args, bool toLowercase) {
    if (file_->eof())
        return false;

    do {
        string line;
        getline(*file_, line);

        size_t type_start = line.find_first_not_of(whitespace_);
        if (!line.empty() && type_start != string::npos && (commentChars_.find(line[type_start]) == string::npos)) {
            size_t sep_start = line.find_first_of(separators_, type_start);
            size_t type_size = string::npos;
            if (sep_start != string::npos) {
                type_size = sep_start - type_start;
                size_t args_start = line.find_first_not_of(separators_, sep_start);
                if (args_start != string::npos) {
                    // skip trailing whitespace
                    size_t args_end = line.find_last_not_of(whitespace_);
                    args = line.substr(args_start,
                                       (args_end != string::npos ? args_end - args_start + 1 : string::npos));
                }
                else
                    args = "";

            }
            type = line.substr(type_start, type_size);
            if (toLowercase)
                transform(type.begin(), type.end(), type.begin(), (int (*)(int))tolower);

            return true;
        }
    } while (!file_->eof());

    return false;
}

bool TextFileReader::getNextLine(string& type, istringstream& args, bool toLowercase/*=true*/) {
    string s;
    bool r = getNextLine(type, s, toLowercase);
    args.clear();
    args.rdbuf()->str(s);
    return r;
}

bool TextFileReader::getNextLinePlain(std::string& line, bool toLowercase) {
    if (file_->eof())
        return false;

    getline(*file_, line);
    if (toLowercase)
        transform(line.begin(), line.end(), line.begin(), (int (*)(int))tolower);
    return true;
}

bool TextFileReader::getNextLinePlain(std::istringstream& line, bool toLowercase) {
    string s;
    bool r = getNextLinePlain(s, toLowercase);
    line.clear();
    line.rdbuf()->str(s);
    return r;
}

bool TextFileReader::getVec3(std::istringstream& args, tgt::vec3& out) {
    args.ignore(1024, '('); // 1024 here is somewhat random - hopefully enough
    args >> out.x;
    args.ignore(1024, ',');
    args >> out.y;
    args.ignore(1024, ',');
    args >> out.z;
    return args.good();
}

bool TextFileReader::operator!() {
    return !(*file_);
}

void TextFileReader::setCommentChars(const std::string& commentChars) {
    commentChars_ = commentChars;
}

void TextFileReader::setSeparators(const std::string& separators) {
    separators_ = separators;
}

void TextFileReader::setWhitespace(const std::string& whitespace) {
    whitespace_ = whitespace;
}

string TextFileReader::getMagicNumber() {
    file_->seekg(0);

    if (file_->eof())
        return "";

    string line;
    getline(*file_, line);
    return line;
}

} // namespace voreen
