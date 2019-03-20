/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "hdf5module.h"

#include "io/hdf5volumereader.h"
#include "io/hdf5volumewriter.h"

namespace voreen {

HDF5Module::HDF5Module(const std::string& modulePath)
    : VoreenModule(modulePath)
    , separateChannels_("separateChannels", "Separate Channels?", true)
{
    setID("HDF5");
    setGuiName("HDF5");

    addProperty(separateChannels_);
    ON_CHANGE_LAMBDA(separateChannels_, []{
        VoreenApplication::app()->showMessageBox("HDF5 VolumeReader changed", "In order for this change to take affect, the application has to be restarted.");
    });

    // Note: The order is chosen on purpose, such that the combined reader will be taken as default.
    // If incompatibilities were observed, swap order!
    if(separateChannels_.get()) {
        registerVolumeReader(new HDF5VolumeReader());
    }
    else {
        registerVolumeReader(new HDF5VolumeReaderCombinedChannels());
    }
    registerVolumeWriter(new HDF5VolumeWriter());

}

} // namespace
