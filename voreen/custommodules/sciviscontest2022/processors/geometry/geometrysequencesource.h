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

#ifndef VRN_GEOMETRYSEQUENCESOURCE_H
#define VRN_GEOMETRYSEQUENCESOURCE_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/buttonproperty.h"

#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"

namespace voreen {

/**
 * Reads a Voreen Geometry, point list or segmented list from file
 * and provides it as geometry through its outport.
 */
class VRN_CORE_API GeometrySequenceSource : public Processor {
public:
    GeometrySequenceSource();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "GeometrySequenceSource"; }
    virtual std::string getCategory() const { return "Input"; }
    virtual CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }
    virtual bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Loads a list of serialized Voreen Geometry (.vge), triangle meshes from .ply, .obj or .stl files, a point list or a segmented point list from a file.\
\
<p>Point lists and segmented point lists are read from ASCII files (.txt or loadable by selecting all file types). In point lists, each point is expected to consist of three components that are separated by white space.\
For a segmented point list, each point is expected to be followed by a numeric segment identifier. The segments have to be listed in ascending order. Additionally, the number of items that are to be skipped after each point can be specified.</p>");
    }

    virtual void process();
    virtual void initialize();

private:
    enum PointListType {
        PointList,          ///< files contain one point per line, elements separated by white space
        SegmentedPointList  ///< files contain one point per line, but with an additional segment id
    };

    /**
     * Delegates to readVoreenGeometry() or readPointList() depending on the selected geometry type
     * and assigns the returned geometry to the outport.
     */
    void readGeometry();

    /**
     * Deserializes a geometry from a Voreen Geometry file (*.vge).
     *
     * @return the read geometry, if deserialization succeeded
     *
     * @throw VoreenException if deserialization failed
     */
    Geometry* readVoreenGeometry(const std::string& filename);

    /**
     * Read a .PLY geometry.
     *
     * @return the read geometry
     * @throw VoreenException if reading failed
     */
    Geometry* readPLYGeometry(const std::string& filename);

    Geometry* readOBJGeometry(const std::string& filename);

    Geometry* readOBJGeometryWithColor(const std::string& filename, tgt::vec4 color);

    Geometry* readSTLGeometry(const std::string& filename);

    /**
     * Reads a point list or point segment list from a text file.
     *
     * @param filename the file to read
     * @param listType the type of the point list to read
     * @param skipItems the number of items to skip after each point
     *
     * @return the point list, if the file could be read successfully
     *
     * @throw VoreenException if the file could not be read
     */
    Geometry* readPointList(const std::string& filename, PointListType listType, int skipItems);

    /**
     * Removed the geometry from the output and clears the file property.
     */
    void clearGeometry();

    /**
     * Triggers reloading the geometry file.
     */
    void forceReload();

    /// Adjusts the visibility of the skipItemCount_ property.
    void updatePropertyVisibility();

    FileDialogProperty geometryDirectory_;   ///< filename of the file containing the point information
    StringOptionProperty geometryType_; ///< Determines whether the file is read as Voreen geometry, point or segment list.
    IntProperty skipItemCount_;         ///< data item to skip in the point list file after reading each point
    ButtonProperty loadGeometry_;
    ButtonProperty clearGeometry_;
    BoolProperty calculateNormals_;

    BoolProperty addColorToOBJ_;
    ColorProperty objColor_;

    GeometryPort outport_;

    bool forceReload_;
    static const std::string loggerCat_;
};

}   //namespace

#endif // VRN_GEOMETRYLISTSOURCE_H
