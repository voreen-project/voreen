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

#ifndef VRN_HISTOGRAMPAINTER_H
#define VRN_HISTOGRAMPAINTER_H

#include "tgt/vector.h"
#include <QWidget>

namespace voreen {

// forward declaration
class Histogram1D;

/**
 * This class draws a polygon that shows the histogram of a volume dataset.
 * To use this painter you have to create an instance of this class with the parent set
 * to the widget where the histogram should be painted. Propagate the histogram that
 * should be painted via the method setHistogram(). The histogram is painted in the method
 * paintEvent() which is called by Qt automatically.
 */
class HistogramPainter : public QWidget {
    Q_OBJECT
public:
    /**
     * Constructor - inits member variables
     *
     * @param parent the parent widget
     * @param xRange range in x direction (usually [0,1])
     * @param yRange range in y direction (usually [0,1])
     * @param padding the number of pixels that are border of the parent widget
     * @param arrowLength length of the arrows in the parent widget
     */
    HistogramPainter(QWidget* parent = 0, tgt::vec2 xRange = tgt::vec2(0.f, 1.f),
                     tgt::vec2 yRange = tgt::vec2(0.f, 1.f), int padding = 0, int arrowLength = 0);

    /// destructor - destroys the histogram object if one was set
    ~HistogramPainter();

    /**
     * Paints the histogram on the screen.
     *
     * @param event the paint event
     */
    void paintEvent(QPaintEvent* event);

    /**
     * Sets the range in x direction to given value.
     *
     * @param xRange new x range
     */
    void setXRange(const tgt::vec2& xRange);

    const Histogram1D* getHistogram() const;

public slots:
    /**
     * Sets the painted histogram to the given value.
     *
     * @param histogram the histogram that will be painted
     */
    void setHistogram(const Histogram1D* histogram);

    /**
     * Determines whether a logarithmic scale is to be used on the y-axis.
     */
    void setYAxisLogarithmic(bool logarithmic);

private:
    /**
     * Helper function for calculation of pixel coordinates from relative coordinates.
     *
     * @param p relative coordinates in the interval [0,1]
     * @return pixel coordinates
     */
    inline tgt::vec2 wtos(const tgt::vec2& p) const;

    tgt::vec2 xRange_;      ///< range in x direction
    tgt::vec2 yRange_;      ///< range in y direction
    bool yAxisLogarithmic_; ///< use logarithmic y-axis
    int padding_;           ///< additional padding at the border of the parent widget
    int arrowLength_;       ///< length of the arrows in the parent widget

    const Histogram1D* histogram_; ///< the histogram that is painted

    QPixmap* cache_;   ///< pixmap for caching the painted histogram
};

} // namespace voreen

#endif // VRN_HISTOGRAMPAINTER_H
