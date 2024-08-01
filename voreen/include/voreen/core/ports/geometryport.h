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

#ifndef VRN_GEOMETRYPORT_H
#define VRN_GEOMETRYPORT_H

#include "voreen/core/ports/genericport.h"
#include "voreen/core/datastructures/geometry/geometry.h"

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API GenericPort<Geometry>;
#endif

class VRN_CORE_API GeometryPort : public GenericPort<Geometry> {
public:
    GeometryPort(PortDirection direction, const std::string& id, const std::string& guiName = "",
                 bool allowMultipleConnections = false,
                 Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT);

    virtual std::string getClassName() const {return "GeometryPort";}
    virtual Port* create(PortDirection direction, const std::string& id, const std::string& guiName = "") const {return new GeometryPort(direction,id,guiName);}
    /// This port type supports caching.
    virtual bool supportsCaching() const;

    /**
     * Returns an hash of the stored Geometry retrieved
     * via Geometry::getHash.
     *
     * If multiple geometries are assigned to the port,
     * the MD5 hash of their concatenated hashes is returned.
     * For an empty port, an emptry string is returned.
     */
    virtual std::string getHash() const;

    /**
     * Saves the assigned Geometry to the given path.
     * If a filename without extension is passed,
     * ".xml" is appended to it.
     *
     * @throws VoreenException If saving failed or
     *      no Geometry is assigned.
     */
    virtual void saveData(const std::string& path) const;

    /**
     * Loads a Geometry from the given path and assigns it
     * to the port. If a filename without extension is passed,
     * ".xml" is appended to it.
     *
     * @throws VoreenException if loading failed.
     */
    virtual void loadData(const std::string& path);

    virtual tgt::col3 getColorHint() const;
    virtual std::string getContentDescription() const;
    virtual std::string getContentDescriptionHTML() const;
};

} // namespace

#endif // VRN_GEOMETRYPORT_H
