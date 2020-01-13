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

#ifndef VRN_BODYPARTSCATALOGBASE_H
#define VRN_BODYPARTSCATALOGBASE_H

#include "voreen/core/voreenapplication.h"
namespace voreen{
    /**
     * This struct represents a BodyParts3D File that contains geometry data
     *
     * @member fileName name of the file
     * @member description description (usually the name) of the organ
     * @member color Color of the geometry in the file
     * @member render determines if the geometry should be loaded
     */
    struct BodyPartfile{
        BodyPartfile(const std::string& filename)
            : fileName(filename)
            , description("")
            , color(tgt::vec4(0.5f,0.5f,0.5f,1.0f))
            , render(false)
        {}
        std::string fileName;
        std::string description;
        tgt::vec4 color;
        bool render;
    };

    /**
     * This struct represents an organ from BodyParts3D
     *
     * @member id Organ id from FMA
     * @member checkState determines the state of the chackbox in the source processor widget
     * @member color Color of the organ
     * @member description description (usually the name) of the organ
     * @member partOf pointer to the parent Organ (relation depending on model)
     * @member parts pointer to child Organs (relation depending on model)
     * @member primitives pointer to Organs that must be loaded as well to load this Organ
     * @member files files that need to be loaded for this particular organ
     */
    struct BodyPart{
        int id;
        bool checkState;
        tgt::vec4 color;
        std::string description;
        std::string filename;
        BodyPart* partOf;
        std::vector<BodyPart*> parts;
        std::vector<BodyPart*> primitives;
        std::vector<BodyPartfile*> files;
        BodyPart(int iD, std::string fileName, std::string desc, bool r, tgt::vec4 c)
            : id(iD)
            ,filename(fileName)
            ,description(desc)
            ,checkState(r)
            ,color(c)
            ,partOf(0){
        };
    };

    /**
     * This Class provides the Base anatomy-atlas structure
     */
    class BodyParts3DCatalogBase{
    public:
        void setPath(const std::string& path){ path_ = path;};
        size_t getSize()const{return catalog_.size();};
        BodyPart* getBodyPart(size_t index);
        /**
         * Loads BodyParts from the file in the catalog (for files without color)
         *
         * @param filename name of the file
         */
        virtual void readBodyParts(const std::string& filename)=0;

        /**
         * Loads BodyParts with color from the file in the catalog
         *
         * @param filename name of the file
         */
        virtual void readBodyPartsColored(const std::string& filename)=0;

        /**
         * Connects BodyParts depending on a relation defined in the file filename
         */
        void connectBodyParts(const std::string& filename);

        /**
         * This function exports the catalog in the format filename (tab) descripton (tab) color (newline)
         *
         * @param filename path to save the file
         */
        virtual void exportCatalog(const std::string& fileName)=0;

        /**
         * Sorts the catalog using the id of the BodyParts
         */
        void sort();

        /**
         * This function gets called recursively and sets for each file wether it needs to be loaded or not
         */
        virtual void setRenderState(BodyPart* part);

        /**
         * Resets the loading state of all files
         */
        virtual void resetRenderState();

        size_t findBodyPart(int id);
        BodyPart* getHead();

        /**
         * Returns all Filenames of checked Organs
         */
        std::vector<std::string> fileNameOutput();

        /**
         * Returns all Colors of checked Organs
         */
        std::vector<tgt::vec4> colorOutput();

        /**
         * Returns all descriptions of checked Organs
         */
        std::vector<std::string> descriptionOutput();
    protected:
        size_t head_; ///< "root" of the model
        std::string path_; ///< path in which all files are
        std::vector<BodyPart*> catalog_; ///< list of all BodyParts
        std::vector<BodyPartfile*> files_; ///< list of all files
    };
}//namespace voreen
#endif // VRN_BODYPARTSCATALOGBASE_H
