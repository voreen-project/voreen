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

#include "binarygeometrysequencesource.h"
#include "binarygeometry.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"
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

const std::string BinaryGeometrySequenceSource::loggerCat_("voreen.sciviscontest2022.BinaryGeometrySequenceSource");

BinaryGeometrySequenceSource::BinaryGeometrySequenceSource()
    : Processor()
    , geometryDirectory_("geometryDirectory", "Geometry Directory", "Open Geometry Directory",
						 VoreenApplication::app()->getUserDataPath(), "Geometry (*.vbge)",
						 FileDialogProperty::FileMode::DIRECTORY)
    , loadGeometry_("loadGeometry", "Load Geometry")
    , clearGeometry_("clearGeometry", "Clear Geometry")
    , outport_(Port::OUTPORT, "geometry.pointlist", "PointList Output", false)
    , forceReload_(false)
{
    Processor::addPort(outport_);

	geometryDirectory_.onChange(MemberFunctionCallback<BinaryGeometrySequenceSource>(this, &BinaryGeometrySequenceSource::forceReload));
    loadGeometry_.onChange(MemberFunctionCallback<BinaryGeometrySequenceSource>(this, &BinaryGeometrySequenceSource::forceReload));
    clearGeometry_.onChange(MemberFunctionCallback<BinaryGeometrySequenceSource>(this, &BinaryGeometrySequenceSource::clearGeometry));

	Processor::addProperty(geometryDirectory_);
	Processor::addProperty(loadGeometry_);
	Processor::addProperty(clearGeometry_);
}

Processor* BinaryGeometrySequenceSource::create() const {
    return new BinaryGeometrySequenceSource();
}

void BinaryGeometrySequenceSource::process() {
    if (!geometryDirectory_.get().empty() && forceReload_) {
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

void BinaryGeometrySequenceSource::initialize() {
    Processor::initialize();
    forceReload_ = true;
}

void BinaryGeometrySequenceSource::readGeometry() {
	// Iterate over list and find all vbge files
	const std::string& filedirectory = geometryDirectory_.get();
	std::vector<std::string> filenames = tgt::FileSystem::readDirectory(filedirectory, true, false);
	std::vector<std::string> geometryFilenames;
	for(const auto& filename : filenames) {
		if(tgt::FileSystem::fileExtension(filename, true) == "vbge") {
			geometryFilenames.push_back(filename);
		}
	}

	size_t idx = 0;
	std::unique_ptr<GeometrySequence> geometries { new GeometrySequence(true) }; // take ownership
	for (const auto& geometryFilename : geometryFilenames) {
		std::string filename = filedirectory + "/";
		filename += geometryFilename;

		LINFO("Load " + filename);
		setProgress(static_cast<float>(idx)/static_cast<float>(geometryFilenames.size()-1));

		if (filename.empty()) {
			++idx;
			continue;
		}

		if (!tgt::FileSystem::fileExists(filename)) {
			LERROR("File does not exists: " + filename);
			++idx;
			continue;
		}

		LINFO("Reading geometry file: " << filename);
		try {
			Geometry* geometry = readVoreenGeometry(filename);
			tgtAssert(geometry, "null pointer returned (exception expected)");
			if(auto tmgb = dynamic_cast<TriangleMeshGeometryBase*>(geometry))
				tmgb->loadTextureData();
			else if (auto glmgb = dynamic_cast<GlMeshGeometryBase*>(geometry))
				glmgb->loadTextureData();

			geometries->addGeometry(geometry);
		}
		catch (VoreenException& e) {
			LERROR(e.what());
		}

		++idx;
	}

	outport_.setData(geometries.release());
	updatePropertyVisibility();
}

Geometry* BinaryGeometrySequenceSource::readVoreenGeometry(const std::string& filename) {
	// read Voreen geometry serialization (.vbge)
	auto geometry = BinaryFileParser::read_file(filename.c_str());
	if (!geometry)
		throw VoreenException("Failed to open file " + geometryDirectory_.get() + " for reading");

	return geometry.release();
}

void BinaryGeometrySequenceSource::clearGeometry() {
    outport_.setData(nullptr);
	geometryDirectory_.set("");
    updatePropertyVisibility();
}

void BinaryGeometrySequenceSource::forceReload() {
    forceReload_ = true;
    invalidate();
}

void BinaryGeometrySequenceSource::updatePropertyVisibility() {
    loadGeometry_.setReadOnlyFlag(geometryDirectory_.get().empty());
    clearGeometry_.setReadOnlyFlag(!outport_.getData());
}

} // namespace voreen
