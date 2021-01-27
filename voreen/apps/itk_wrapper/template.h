/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef TEMPLATE_H
#define TEMPLATE_H

#include "tgt/logmanager.h"

namespace voreen {

/**
 * A Template provides methods to read a template-file (e.g. a textfile) into a string,
 * find and replace words in this string, and write the string to a new file.
*/
class Template {
public:
    /**
     * Constructor
     * Reads the template-file into a string
     *
     * @param filename, path of the template-file
    */
    Template(std::string filename);

    /**
     * Finds all words which are equal to "find" in the templateString_ and replaces
     * them by "replace".
     *
     * @param find, word that should be replaced
     * @param replace, word that replaces the word "find"
    */
    void replace(std::string find, std::string replace);

    /**
     * Returns the current template-string
     *
     * @return templateString_
    */
    std::string getTemplateString();

    /**
     * Sets the template-string
     *
     * @param templateString
    */
    void setTemplateString(std::string templateString);

    /**
     * Writes the current template-string to the file "file"
     *
     * @param file, path of the file to which the templateString_ should be written
    */
    void writeToFile(std::string file);

    static const std::string loggerCat_;

private:
    /**
     * Reads the template-file into a string.
     *
     * @return result, string containing the template-file
    */
    std::string readTemplate();

    std::string filename_;          ///< Path of the template-file.
    std::string templateString_;    ///< Template-string
};

} //namespace voreen

#endif // TEMPLATE_H
