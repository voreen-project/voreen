/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "poimodule.h"
#include "processors/poipointsegmentgeometryimporter.h"
#include "processors/poipointsegmentgeometryexporter.h"
#include "processors/poicsvexport.h"
#include "processors/poicsvimport.h"
#include "processors/poisave.h"
#include "processors/poiselectionmanipulation.h"
#include "processors/poistorage.h"
#include "processors/poisource.h"
#include "processors/poirenderer2d.h"
#include "processors/poirenderer3d.h"
#include "processors/poitextinfo.h"
#include "datastructures/poilist.h"


namespace voreen {

POIModule::POIModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("POI Module");

    setGuiName("POI Module");
    addShaderPath(getModulePath("glsl"));

    registerProcessor(new POICSVImport());
    registerProcessor(new POICSVExport());
    registerProcessor(new POISave());
    registerProcessor(new POISelectionManipulation());
    registerProcessor(new POIStorage());
    registerProcessor(new POISource());
    registerProcessor(new POIRenderer2d());
    registerProcessor(new POIRenderer3d());
    registerProcessor(new POITextInfo());
    registerProcessor(new POIPointSegmentExporter());
    registerProcessor(new POIPointSegmentImporter());
}

std::string POIModule::getDescription() const {
    return "My first sample module.";
}

} // namespace

