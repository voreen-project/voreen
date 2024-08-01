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

#include <string>
#include <fstream>

#include "template.h"
#include "voreen/core/utils/exception.h"

namespace voreen {

const std::string Template::loggerCat_("itk-wrapper.Template");

Template::Template(std::string filename) {
    filename_ = filename;
    if (!filename_.empty()) {
        templateString_ = readTemplate();
    }
}

void Template::replace(std::string find, std::string replace) {
    while (templateString_.find(find) != std::string::npos) {
        templateString_.replace(templateString_.find(find), find.length(), replace);
    }
}

std::string Template::getTemplateString() {
    return templateString_;
}

void Template::setTemplateString(std::string templateString) {
    templateString_ = templateString;
}

void Template::writeToFile(std::string file) {
    std::fstream stream;
    stream.open(file.c_str(), std::ios::out);
    stream.write(templateString_.c_str(), templateString_.length());
    stream.clear();
    stream.close();
}

std::string Template::readTemplate() {
    std::fstream stream;
    std::string temp;
    std::string result;

    //read templateFile into string
    stream.open(filename_.c_str(), std::ios::in);
    if (stream.is_open()) {
        std::getline(stream, temp);
        result += temp;
        while (stream.good()) {
            std::getline(stream, temp);
            result += "\n"+temp;
        }
        stream.clear();
        stream.close();
    }
    else {
        throw VoreenException("\nFailed to open template file: " + filename_);
    }
    return result;
}

} // namespace
