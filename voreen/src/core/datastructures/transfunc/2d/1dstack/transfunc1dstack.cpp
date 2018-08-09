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

#include "voreen/core/datastructures/transfunc/2d/1dstack/transfunc1dstack.h"

#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"

namespace voreen {

const std::string TransFunc1DStack::loggerCat_("voreen.TransFunc1DStack");

TransFunc1DStack::TransFunc1DStack(int resolution)
    : TransFunc2D(resolution, resolution, TF_UBYTE, tgt::Texture::LINEAR)
    , tfResolution_(resolution)
{
    tfStack_.insert(std::make_pair<>(0.f, new TransFunc1DKeys(resolution)));
}

TransFunc1DStack::~TransFunc1DStack() {
    for (std::map<float, TransFunc1D*>::iterator it = tfStack_.begin(); it != tfStack_.end(); it++)
        delete it->second;
    tfStack_.clear();
}


} // namespace voreen
