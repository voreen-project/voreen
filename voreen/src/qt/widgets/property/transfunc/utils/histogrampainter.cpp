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

#include "voreen/qt/widgets/property/transfunc/utils/histogrampainter.h"

#include "voreen/core/datastructures/volume/histogram.h"

#include <QPainter>
#include <QPaintEvent>
#include <QPoint>

namespace voreen {

HistogramPainter::HistogramPainter(QWidget* parent, tgt::vec2 xRange, tgt::vec2 yRange,
                                   int padding, int arrowLength)
    : QWidget(parent)
    , xRange_(xRange)
    , yRange_(yRange)
    , yAxisLogarithmic_(true)
    , padding_(padding)
    , arrowLength_(arrowLength)
    , histogram_(0)
    , cache_(0)
{
}

HistogramPainter::~HistogramPainter() {
    delete cache_;
}

void HistogramPainter::setHistogram(const Histogram1D* histogram) {
    if(histogram_ != histogram) {
        histogram_ = histogram;
        delete cache_;
        cache_ = 0;
        update();
    }
}

const Histogram1D* HistogramPainter::getHistogram() const {
    return histogram_;
}

void HistogramPainter::setYAxisLogarithmic(bool logarithmic) {
    yAxisLogarithmic_ = logarithmic;
    delete cache_;
    cache_ = 0;
}

void HistogramPainter::setXRange(const tgt::vec2& xRange) {
    if (xRange != xRange_) {
        xRange_ = xRange;
        delete cache_;
        cache_ = 0;
    }
}

void HistogramPainter::paintEvent(QPaintEvent* event) {
    event->accept();

    if (cache_ == 0 || cache_->rect() != rect()) {
        delete cache_;
        cache_ = new QPixmap(rect().size());
        cache_->fill(Qt::transparent);

        QPainter paint(cache_);

        // put origin in lower lefthand corner
        QMatrix m;
        m.translate(0.0, static_cast<float>(height())-1);
        m.scale(1.f, -1.f);
        paint.setMatrix(m);
        paint.setMatrixEnabled(true);

        if (histogram_) {
            // draw histogram
            paint.setPen(Qt::NoPen);
            paint.setBrush(QColor(255, 135, 135, 255)); //200 0 0 120
            paint.setRenderHint(QPainter::Antialiasing, true);

            int histogramWidth = static_cast<int>(histogram_->getNumBuckets());
            tgt::vec2 p;

            // Qt can't handle polygons that have more than 65536 points
            // so we have to split the polygon
            int maxSize = 65536; //max size of polygon
            std::vector<QPointF*> points;
            int vi = 0; //iterator in the points vector
            points.push_back(new QPointF[maxSize]);
            int count = 0;

            for (int x=0; x < histogramWidth; ++x) {
                float xpos = static_cast<float>(x) / histogramWidth;
                xpos = histogram_->getMinValue() + (histogram_->getMaxValue() - histogram_->getMinValue()) * xpos;
                // Do some simple clipping here, as the automatic clipping of drawPolygon()
                // gets very slow if lots of polygons have to be clipped away, e.g. when
                // zooming to small part of the histogram.
                if (xpos >= xRange_[0] && xpos <= xRange_[1]) {
                    //open new list, if old one is full
                    if( count == maxSize-2 ){
                        count = 0;
                        points.push_back(new QPointF[maxSize]);
                        vi++;
                        //copy last point to connect two polygons
                        points[vi][count].rx() = p.x;
                        points[vi][count].ry() = p.y;
                        count++;
                    }
                    float value = (yAxisLogarithmic_ ? histogram_->getBucketLogNormalized(x) : histogram_->getBucketNormalized(x));
                    p = wtos(tgt::vec2(xpos, value * (yRange_[1] - yRange_[0]) + yRange_[0]));

                    // optimization: if the y-coord has not changed from the two last points
                    // then just update the last point's x-coord to the current one
                    if( (count >= 2 ) && (points[vi][count - 2].ry() == p.y) && (points[vi][count - 1].ry() == p.y) && (count >= 2) ){
                        points[vi][count - 1].rx() = p.x;
                    } else {
                        points[vi][count].rx() = p.x;
                        points[vi][count].ry() = p.y;
                        count++;
                    }
                }
            }

            for(size_t i = 0; i < points.size(); ++i){
                if (count > 0) {
                    if (i == vi){
                        // needed for a closed polygon
                        p = wtos(tgt::vec2(0.f, yRange_[0]));
                        points[i][count].rx() = points[i][count -1].rx();
                        points[i][count].ry() = p.y;
                        count++;
                        points[i][count].rx() = points[i][0].rx();
                        points[i][count].ry() = p.y;
                        count++;

                        paint.drawPolygon(points[i], count);
                    } else {
                        // needed for a closed polygon
                        p = wtos(tgt::vec2(0.f, yRange_[0]));
                        points[i][maxSize - 2].rx() = points[i][maxSize - 3].rx();
                        points[i][maxSize - 2].ry() = p.y;
                        points[i][maxSize - 1].rx() = points[i][0].rx();
                        points[i][maxSize - 1].ry() = p.y;

                        paint.drawPolygon(points[i], maxSize);
                    }
                }
            }

            for (size_t i = 0; i < points.size(); ++i)
                delete[] points[i];
        }
    }

    QPainter paint(this);
    paint.drawPixmap(0, 0, *cache_);
}

tgt::vec2 HistogramPainter::wtos(const tgt::vec2& p) const {
    float sx = (p.x - xRange_[0]) / (xRange_[1] - xRange_[0]) * (static_cast<float>(width())  - 2 * padding_ - 1.5 * arrowLength_) + padding_;
    float sy = (p.y - yRange_[0]) / (yRange_[1] - yRange_[0]) * (static_cast<float>(height()) - 2 * padding_ - 1.5 * arrowLength_) + padding_;
    return tgt::vec2(sx, sy);
}

} // namespace
