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

#ifndef VRN_VORTICITYFILTER_H
#define VRN_VORTICITYFILTER_H

#include "parallelvolumefilter.h"
#include "gradientfilter.h"

namespace voreen {

class VorticityFilter : public ParallelVolumeFilter<ParallelFilterValue3D, ParallelFilterValue3D> {
public:
    VorticityFilter(GradientType gradientType, const tgt::vec3& spacing, const SamplingStrategy<ParallelFilterValue3D>& samplingStrategy);
    virtual ~VorticityFilter() {}

    ParallelFilterValue3D getValue(const Sample& sample, const tgt::ivec3& pos, const SliceReaderMetaData& inputMetaData, const SliceReaderMetaData& outputMetaData) const;

    GradientType getGradientType() const;

private:

    const GradientFilter gradientX_;
    const GradientFilter gradientY_;
    const GradientFilter gradientZ_;

    const GradientType gradientType_;
};

} // namespace voreen

#endif
