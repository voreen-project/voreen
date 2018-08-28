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

#include "stringutility.h"
#include <set>

std::vector<std::string> StringUtility::explode(const std::string& input, const char separator) {
    size_t lastPos = 0;
    std::vector<std::string> result;
    for (size_t i = 0; i < input.length(); ++i) {
        if (input[i] == separator) {
            result.push_back(input.substr(lastPos, (i - lastPos)));
            lastPos = i + 1;
        }
    }
    if (lastPos < input.length())
        result.push_back(input.substr(lastPos, (lastPos - input.length() - 1)));
    return result;
}

std::string StringUtility::getFileExtension(const std::string& filename) {
    size_t pos = filename.rfind('.');
    if (pos != std::string::npos)
        return filename.substr(pos + 1);
    return "";
}

bool StringUtility::isCharInSet(const char input, const std::string& options) {
    std::set<char> optionsSet;
    for (size_t i = 0; i < options.size(); ++i)
        optionsSet.insert(options[i]);
    return (optionsSet.find(input) != optionsSet.end());
}
