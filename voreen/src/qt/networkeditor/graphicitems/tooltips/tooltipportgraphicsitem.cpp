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

#include "voreen/qt/networkeditor/graphicitems/tooltips/tooltipportgraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

#include "voreen/qt/networkeditor/graphicitems/core/portgraphicsitem.h"
//ports
#include "voreen/core/ports/port.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/volumeport.h"
    //volume
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumepreview.h"

#include <QGraphicsTextItem>

namespace voreen {

ToolTipPortGraphicsItem::ToolTipPortGraphicsItem(PortGraphicsItem* item)
    : ToolTipBaseGraphicsItem(item)
    , image_(0)
    , imageAspect_(0.f)
{
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
const QImage* ToolTipPortGraphicsItem::getImage() const{
    return image_;
}

const float ToolTipPortGraphicsItem::getImageAspect() const {
    return imageAspect_;
}

//---------------------------------------------------------------------------------------------------------------
//                  nwebasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
QRectF ToolTipPortGraphicsItem::boundingRect() const {
    return currentStyle()->ToolTipPortGI_boundingRect(this);
}

QPainterPath ToolTipPortGraphicsItem::shape() const {
    return currentStyle()->ToolTipPortGI_shape(this);
}

void ToolTipPortGraphicsItem::initializePaintSettings() {
    currentStyle()->ToolTipPortGI_initializePaintSettings(this);
}

void ToolTipPortGraphicsItem::mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting) {
    currentStyle()->ToolTipPortGI_paint(this,painter,option,widget,setting);
}

//---------------------------------------------------------------------------------------------------------------
//                  tooltipbasegraphicsitem functions
//---------------------------------------------------------------------------------------------------------------
bool ToolTipPortGraphicsItem::updateToolTip(qreal x, qreal y) {
    prepareGeometryChange();
    Port* port = dynamic_cast<PortGraphicsItem*>(parent())->getPort();

    // Check for port condition not being met.
    if (port->hasData()) {
        std::string failures;
        if (!port->areConditionsMet(&failures)) {
            // Delete image if set before.
            delete image_;
            image_ = 0;

            // Set proper error message as tooltip.
            toolTipText_->setHtml(QString::fromStdString(failures));
            return ToolTipBaseGraphicsItem::updateToolTip(x, y);
        }
    }
    
    // Set proper description.
    toolTipText_->setHtml(QString::fromStdString(port->getContentDescriptionHTML()));

    //test, if we have a image to show
    RenderPort* rp = 0;
    VolumePort* vp = 0;
    //try+catch to guarantee a call of TollTipBase::updateToolTip
    try {
        if (rp = dynamic_cast<RenderPort*>(port)) {
            delete image_;
            image_ = 0;

            if (!rp->hasData())
                throw "exception";

            tgt::ivec2 size = rp->getSize();
            tgt::Texture* tex = rp->getColorTexture();
            if (!tex)
                throw "exception";

            // download pixel data to buffer
            tgt::col4* pixels = 0;
            try {
                GLubyte* pixels_b = tex->downloadTextureToBuffer(GL_RGBA, GL_UNSIGNED_BYTE);
                LGL_ERROR;
                if (pixels_b)
                    pixels = reinterpret_cast<tgt::col4*>(pixels_b);
                else {
                    LERRORC("ToolTipPortGraphicsItem::updateToolTip", "failed to download texture");
                    throw "exception";
                }
            }
            catch (std::bad_alloc&) {
                LERRORC("ToolTipPortGraphicsItem::updateToolTip", "bad allocation");
                throw "exception";
            }

            image_ = new QImage(size.x, size.y, QImage::Format_ARGB32);

            // The pixels are stored row by row from bottom to top an in each row from left to right
            QColor color;
            for (int y=0; y < size.y; ++y) {
                for (int x=0; x < size.x; ++x) {
                    tgt::col4 col = pixels[x + y*size.x];
                    color.setRgb(col.r, col.g, col.b, col.a);
                    //(0,0) is top left
                    image_->setPixel(x, size.y - 1 - y, color.rgba());
                }
            }
            delete[] pixels;
            pixels = 0;

            imageAspect_ = static_cast<float>(size.x) / static_cast<float>(size.y);
        }
        else if(vp = dynamic_cast<VolumePort*>(port)) {
            delete image_;
            image_ = 0;
            if (!vp->hasData())
                throw "exception";

            VolumePreview* prev = vp->getData()->getDerivedData<VolumePreview>();
            //no preview present
            if(!prev) {
                return ToolTipBaseGraphicsItem::updateToolTip(x,y);
            }
            int internHeight = prev->getHeight();

            QImage origImg = QImage(internHeight, internHeight, QImage::Format_ARGB32);
            for (int y=0; y<internHeight; y++) {
                for (int x=0; x<internHeight; x++) {
                    int previewIndex = y * internHeight + x;
                    int greyVal = prev->getData()[previewIndex];
                    origImg.setPixel(x, y, qRgb(greyVal, greyVal, greyVal));
                }
            }
            image_ = new QImage(origImg.scaled(63, 63, Qt::IgnoreAspectRatio, Qt::SmoothTransformation));

            imageAspect_ = 1;
        }
        return ToolTipBaseGraphicsItem::updateToolTip(x,y);
    }
    catch (...){
        return ToolTipBaseGraphicsItem::updateToolTip(x,y);
    }
}


} // namespace voreen
