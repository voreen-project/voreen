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

#ifndef VRN_BODYPARTS3DSOURCE_H
#define VRN_BODYPARTS3DSOURCE_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/interaction/idmanager.h"
#include "modules/touchtable/datastructures/geometry/trianglemeshgeometrysinglecolor.h"
#include "voreen/core/ports/geometryport.h"

namespace voreen{

    enum BodyPartsVersion{
        Version3_0,
        Version3_0HQ,
        Version4_0
    };

    /**
    * Reads Geometry files for BodyParts3D
    * currently supports versions 3.0 and 4.0
    */
    class VRN_CORE_API BodyParts3DSource : public Processor {
    public:
        BodyParts3DSource();
        ~BodyParts3DSource();

        virtual Processor* create() const { return new BodyParts3DSource();}

        virtual std::string getClassName() const { return "BodyParts3DSource";}
        virtual std::string getCategory() const  { return "Input";}
        virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;}
        std::string getPath() const { return folderPath_.get();}
        void setPath(const std::string& p) { folderPath_.set(p);}
        BodyPartsVersion getVersion() const { return version_;}
        virtual bool usesExpensiveComputation() const { return true; }

        bool isReady()const;

        /**
         * Sets version_ and adjusts the the selection of versionoption_ accordingly
         */
        void setVersion( BodyPartsVersion v);

        /**
         * Sets pathValid_ ...
         */
        void setPathValid( bool b );

        /**
         *  called when versionoption_ is changed. sets the version_ according to selection of versionoption_
         */
        void versionChanged();

        /**
         * calls updateFromProcessor for the widget
         */
        void updateWidget();

        /**
         * calls Processor::deserialize() and trys to initialize (load the Catalog) the Processorwidget
         */
        virtual void initialize();

        /**
        * Decides which version of BodyParts3D is currently selected and wether a Wavefront-obj or a bpvrn file should be loaded
        *
        * @param color the color for the geometry
        * @param filename the filename of the geometry
        * @return Geometry
        */
        TriangleMeshGeometryBodyParts3D* readPart(const tgt::vec4& color,const std::string& filename);

        /**
        * Reads an .obj-file, creates a .bpvrn-file in the same folder with the same name and
        * returns the geometry created from the file
        *
        * @param color the color of the geometry
        * @param filepath the path to the file to be read
        * @return the Geometry described by the file
        */
        TriangleMeshGeometryBodyParts3D* readObj(tgt::vec4 color, std::string& filepath);

        /**
        * Reads a .bpvrn-file and returns the geometry it describes
        *
        * @param color color of the geometry
        * @param filepath the path to the file to be read
        * @return the Geometry described by the file
        */
        TriangleMeshGeometryBodyParts3D* readBpvrn(tgt::vec4 color, std::string& filepath);

        /**
        * Handles the createMesh request from the widget. there is a 1 to 1 (to 1) correlation between the
        * vectors
        *
        * @param filenames the filenames of the geometries
        * @param colors the colors of the geometries
        * @param descriptions contains strings with the names of the bodyparts
        */
        void createMesh(const std::vector<std::string>& filenames, const std::vector<tgt::vec4>& colors, const std::vector<std::string>& descriptions);
    protected:
        virtual void setDescriptions(){
            setDescription("This is the Processor for parsing and loading BodyParts3D");
        }
        virtual void process();

    private:

        ButtonProperty loadBodyParts_;

        FileDialogProperty folderPath_;    ///< path to the folder containing the data for BodyParts3D 3.0

        StringOptionProperty versionoption_;        ///< dropdown menu to pick a version

        BodyPartsVersion version_;

        GeometryPort geometryPort_;

        BoolProperty pathValid_;

        static const std::string loggerCat_;
    };

}

#endif //VRN_BODYPARTS3DSOURCE_H
