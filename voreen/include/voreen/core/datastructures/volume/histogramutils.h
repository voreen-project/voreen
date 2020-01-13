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

#ifndef VRN_HISTOGRAMUTILS_H
#define VRN_HISTOGRAMUTILS_H

#include "voreen/core/datastructures/volume/histogram.h"

namespace voreen {

    struct GaussianCurve {
        float mean;
        float variance;
    };

    /** Creates a new histogram for the interval [minVal, maxVal] filled with the values from histogram.
     * Values from the old histogram that lay outside the new interval will be dismissed.
     *
     * @param histogram the base histogram
     * @param minVal the minimum sample value for the new histogram
     * @param maxVal the maximum sample value for the new histogram
     * @return the new histogram for [minVal, maxVal]
     **/
    Histogram1D* subHistogram(Histogram1D* histogram, float minVal, float maxVal);

    /** Divides the histogram data in two classes of gaussian distributions based on the
     * automatic tresholding method by Otsu (see Otsu 1979).
     *
     * @param histogram the histogram to be classified
     * @return an array with the two resulting gaussian distributions
     **/
    GaussianCurve* otsuTreshold(Histogram1D* histogram);

    /**
     * Computes the optimal threshold using Otsu's method for a given histogram.
     *
     * @param histogram the histogram for which the threshold should be computed
     * @return the threshold in normalized intensity values
     */
    float computeOptimalOtsuThreshold(const Histogram1D* histogram);

    /** Maps a gaussian curve form interval 'from' to interval 'to' i.e. changes
     * its mean and variance so that the curve has the same proportions over the new interval
     * like it had in the old interval.
     *
     * @param c the curve to be mapped proportionally
     * @param from the current interval as vec2
     * @param to the interval c is mapped to
     * @return a proportionally mapped curve
     */
    GaussianCurve mapCurve(GaussianCurve& c, tgt::vec2 from, tgt::vec2 to);


} // namespace voreen

#endif // VRN_HISTOGRAUTILS_H
