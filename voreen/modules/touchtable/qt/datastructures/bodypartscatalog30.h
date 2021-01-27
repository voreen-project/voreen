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

#ifndef VRN_BODYPARTSCATALOG30_H
#define VRN_BODYPARTSCATALOG30_H

#include "bodypartscatalogbase.h"
#include "voreen/core/voreenapplication.h"

namespace voreen{

    /**
    * This Class realizes the construction of the base structure for BodyParts3D 3.0
    */
    class BodyParts3DCatalog30 : public BodyParts3DCatalogBase{
    public:
        BodyParts3DCatalog30(const std::string& dir, bool colorVersion);
        BodyParts3DCatalog30();
        ~BodyParts3DCatalog30();

        /**
         * Loads BodyParts from the file in the catalog (for files without color)
         */
        void readBodyParts(const std::string& filename);

        /**
         * Loads BodyParts with color from the file in the catalog
         */
        void readBodyPartsColored(const std::string& filename);

        /**
         * Connects the BodyParts. BodyParts Connected in this way have to load their "child" BodyParts
         * also when they are loaded.
         */
        void addPrimitives(const std::string& filename);

        /**
         * This function exports the catalog in the format filename (tab) descripton (tab) color (newline)
         *
         * @param filename path to save the file
         */
        void exportCatalog(const std::string& fileName);
    };

}
#endif // VRN_BODYPARTSCATALOG30_H
