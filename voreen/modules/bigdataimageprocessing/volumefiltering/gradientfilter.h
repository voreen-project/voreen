/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_GRADIENTFILTER_H
#define VRN_GRADIENTFILTER_H

#include "parallelvolumefilter.h"

namespace voreen {

enum GradientType{
    VOG_CENTRAL_DIFFERENCE = 0,
    VOG_LINEAR_REGRESSION = 1,
    VOG_SOBEL = 2,
};

class GradientFilter : public ParallelVolumeFilter<ParallelFilterValue1D, ParallelFilterValue3D> {
public:
    GradientFilter(GradientType gradientType, const tgt::vec3& spacing, const SamplingStrategy<ParallelFilterValue1D>& samplingStrategy, const std::string& sliceBaseType);
    virtual ~GradientFilter() {}

    ParallelFilterValue3D getValue(const Sample& sample, const tgt::ivec3& pos) const;

    GradientType getGradientType() const;

private:

    ParallelFilterValue3D calcGradientCentralDifferences(const Sample& sample, const tgt::ivec3& pos);
    ParallelFilterValue3D calcGradientLinearRegression(const Sample& sample, const tgt::ivec3& pos);
    ParallelFilterValue3D calcGradientSobel(const Sample& sample, const tgt::ivec3& pos);

    std::function<ParallelFilterValue3D(const Sample&, const tgt::ivec3&)> gradientFunction_;

    const GradientType gradientType_;
    const tgt::vec3 spacing_;
};

} // namespace voreen

#endif
