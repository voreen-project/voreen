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

#include "binarygeometrysource.h"
#include "binarygeometry.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometryindexed.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/meta/serializablevectormetadata.h"

#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "tgt/filesystem.h"
#include "tgt/exception.h"

#include <vector>
#include <fstream>

using tgt::vec3;
using tgt::ivec3;
using tgt::ivec2;
using std::vector;
using std::string;

namespace voreen {

const std::string BinaryGeometrySource::loggerCat_("voreen.sciviscontest2022.BinaryGeometrySource");

BinaryGeometrySource::BinaryGeometrySource()
    : Processor()
    , geometryFile_("geometryFile", "Geometry File", "Open Geometry File", VoreenApplication::app()->getUserDataPath(), "Geometry (*.vbge)")
    , loadGeometry_("loadGeometry", "Load Geometry")
    , clearGeometry_("clearGeometry", "Clear Geometry")
    , outport_(Port::OUTPORT, "geometry.pointlist", "PointList Output", false)
    , forceReload_(false)
{
    Processor::addPort(outport_);

    geometryFile_.onChange(MemberFunctionCallback<BinaryGeometrySource>(this, &BinaryGeometrySource::forceReload));
    loadGeometry_.onChange(MemberFunctionCallback<BinaryGeometrySource>(this, &BinaryGeometrySource::forceReload));
    clearGeometry_.onChange(MemberFunctionCallback<BinaryGeometrySource>(this, &BinaryGeometrySource::clearGeometry));

	Processor::addProperty(geometryFile_);
	Processor::addProperty(loadGeometry_);
	Processor::addProperty(clearGeometry_);
}

Processor* BinaryGeometrySource::create() const {
    return new BinaryGeometrySource();
}

void BinaryGeometrySource::process() {
    if (geometryFile_.get() != "" && forceReload_) {
        try {
            readGeometry();
        }
        catch (tgt::FileNotFoundException& f) {
            LERROR(f.what());
        }
        forceReload_ = false;
        updatePropertyVisibility();
    }
}

void BinaryGeometrySource::initialize() {
    Processor::initialize();
    forceReload_ = true;
}

void BinaryGeometrySource::readGeometry() {
    std::string filename = geometryFile_.get();
    setProgress(0.f);

    if (geometryFile_.get() == "")
        return;

    if (!tgt::FileSystem::fileExists(geometryFile_.get()))
        throw tgt::FileNotFoundException("File does not exists", geometryFile_.get());

	LINFO("Reading geometry file: " << geometryFile_.get());
	try {
		Geometry* geometry = readVoreenGeometry(geometryFile_.get());
		tgtAssert(geometry, "null pointer returned (exception expected)");
		if(TriangleMeshGeometryBase* tmgb = dynamic_cast<TriangleMeshGeometryBase*>(geometry))
			tmgb->loadTextureData();
		else if (GlMeshGeometryBase* glmgb = dynamic_cast<GlMeshGeometryBase*>(geometry))
			glmgb->loadTextureData();
		outport_.setData(geometry);
		setProgress(1.f);
	}
	catch (VoreenException& e) {
		LERROR(e.what());
		setProgress(0.f);
	}

    updatePropertyVisibility();
}

Geometry* BinaryGeometrySource::readVoreenGeometry(const std::string& filename) {
	// read Voreen geometry serialization (.vbge)
	auto geometry = BinaryFileParser::read_file(filename.c_str());
	if (!geometry)
		throw VoreenException("Failed to read file " + geometryFile_.get());

	return geometry.release();
}

void BinaryGeometrySource::clearGeometry() {
    outport_.setData(0);
    geometryFile_.set("");
    updatePropertyVisibility();
}

void BinaryGeometrySource::forceReload() {
    forceReload_ = true;
    invalidate();
}

void BinaryGeometrySource::updatePropertyVisibility() {
    loadGeometry_.setReadOnlyFlag(geometryFile_.get() == "");
    clearGeometry_.setReadOnlyFlag(!outport_.getData());
}

} // namespace voreen
