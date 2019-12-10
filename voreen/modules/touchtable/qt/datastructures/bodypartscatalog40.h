/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_BODYPARTSCATALOG40_H
#define VRN_BODYPARTSCATALOG40_H

#include "bodypartscatalogbase.h"
#include "voreen/core/voreenapplication.h"
namespace voreen{

    /**
     * This Class realizes the construction of the base structure for BodyParts3D 4.0
     */
    class BodyParts3DCatalog40 : public BodyParts3DCatalogBase{
    public:
        BodyParts3DCatalog40(const std::string& dir, bool colorVersion);
        BodyParts3DCatalog40();
        ~BodyParts3DCatalog40();

        /**
         * Loads BodyParts from the file in the catalog (for files without color)
         *
         * @param filename name of the file
         */
        void readBodyParts(const std::string& filename);

        /**
         * Loads BodyParts with color from the file in the catalog
         *
         * @param filename name of the file
         */
        void readBodyPartsColored(const std::string& filename);

        /**
         * fills the files for Each BodyPart and also files_
         */
        void addPrimitives(const std::string& filename);

        /**
         * This function exports the catalog in the format filename (tab) descripton (tab) color (newline)
         *
         * @param filename path to save the file
         */
        void exportCatalog(const std::string& fileName);

        /**
         * Called after the catalog and files are all initialized, this function copies the
         * essential data from the BodyParts to the BodyPartfiles
         */
        void setFileAttributes();

        /**
         * returns a pointer to a file with name fileName
         */
        BodyPartfile* getFilePointer(const std::string& fileName);
    };
}
#endif // VRN_BODYPARTSCATALOG40_H
