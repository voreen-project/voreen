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

#include "gdcmmodule.h"

#include "modules/gdcm/io/gdcmvolumereader.h"

namespace voreen {

GdcmModule::GdcmModule(const std::string& modulePath)
    : VoreenModule(modulePath),
      aetProperty_("aeTitle", "DICOM AE Title", "Voreen-GDCM", Processor::INVALID_RESULT, Property::LOD_APPLICATION),
      incomingPortProperty_("incomingPort", "Incoming DICOM Port", 104, 1, 65535, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_APPLICATION),
      scpUrlProperty_("scpUrl", "Default SCP URL", "www.dicomserver.co.uk" , Processor::INVALID_RESULT, Property::LOD_APPLICATION),
      scpAetProperty_("scpAet", "Default SCP AE Title", "ANY-SCP" , Processor::INVALID_RESULT, Property::LOD_APPLICATION),
      scpPortProperty_("scpPort", "Default SCP Port", 104, 1, 65535 , Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_APPLICATION),
      ignoreSliceSpacing_("ignorespacing", "Ignore differences in slice spacing above 10%", false , Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT),
      useGdcmRescaling_("gdcmrescaling", "Use GDCM Rescaling", false , Processor::INVALID_RESULT, Property::LOD_DEVELOPMENT)
{
    //module name
    setID("GrassrootsDICOM");
    setGuiName("Grassroots DICOM");

    //register each processor
    registerVolumeReader(new GdcmVolumeReader());

    addProperty(aetProperty_);
    addProperty(incomingPortProperty_);
    addProperty(scpUrlProperty_);
    addProperty(scpAetProperty_);
    addProperty(scpPortProperty_);
    addProperty(ignoreSliceSpacing_);
    addProperty(useGdcmRescaling_);

}

} //end namespace

