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

#ifndef VRN_TRANSFUNCPROPERTYWIDGETPAINTERSLIDER_H
#define VRN_TRANSFUNCPROPERTYWIDGETPAINTERSLIDER_H

#include "tgt/vector.h"

#include <QObject>

namespace voreen {

/***********************************************************
 * Base Slider Class                                       *
 * Provides sliders which can be moved inside a sub canvas *
 * of a given canvas.                                      *
 ***********************************************************/
class TransFuncPropertyWidgetPainterBaseSlider {
public:
    /**
     * Defines where in the canvas the slider should be located.
     */
    enum SliderOrientation {
        SO_LEFT,
        SO_RIGHT,
        SO_TOP,
        SO_BOTTOM
    };

    /** Constructor */
    TransFuncPropertyWidgetPainterBaseSlider(size_t associatedDimension, SliderOrientation orientation, int sliderSize);
    /** Paints the slider. Has to be implemented by sub classes. */
    virtual void paint() = 0;
    /** Checks, if the slider is licked by the mouse at pos. Has to be implemented by sub classes. */
    virtual bool isHit(tgt::ivec2 pos) = 0;
    //--------------------
    //  getter and setter
    //--------------------
    /** Sets the current canvas size. Calls "updateSliderPosition()". */
    void setCanvasSize(tgt::ivec2 size);
    /** Sets the current sub canvas offset in pixel. Calls "updateSliderPosition()". */
    void setSubCanvasOffset(tgt::ivec2 offset);
    /** Sets the current sub canvas width in pixel. Calls "updateSliderPosition()". */
    void setSubCanvasSize(tgt::ivec2 size);
    /** Sets the current value range of the sub canvas. Calls "updateSliderPosition()". */
    void setValueRange(tgt::vec2 range);
protected:
    /**
     * Calculates and sets the sider position according to valueRange_ and canvasSizeInPixel_.
     * Has to be implemented by sub classes.
     */
    virtual void updateSliderPosition() = 0;
    /**
     * Calculates and sets the stored value according to the slider position.
     * Has to be implemented by sub classes.
     */
    virtual void updateStoredValue() = 0;

    //--------------------
    //  member
    //--------------------
    tgt::ivec2 canvasSizeInPixel_;          ///< current canvas size
    tgt::ivec2 subCanvasOffsetInPixel_;     ///< left bottom sub canvas offset
    tgt::ivec2 subCanvasSizeInPixel_;       ///< sub canvas size

    tgt::vec2 valueRange_;          ///< value range of the sub canvas

    SliderOrientation orientation_;  ///< orientation of the slider
    size_t associatedDimension_;     ///< dimension represented by the slider of the slider
    int sliderSize_;                  ///< size of the slider in pixel
};

/***********************************************************
 * Gamma Slider                                            *
 ***********************************************************/
class TransFuncPropertyWidgetPainterGammaSlider : public QObject, public TransFuncPropertyWidgetPainterBaseSlider {
    Q_OBJECT
public:
    /** Constructor */
    TransFuncPropertyWidgetPainterGammaSlider(size_t associatedDimension, SliderOrientation orientation, int sliderSize);
    /** @see TransFuncPropertyWidgetPainterBaseSlider::paint */
    void paint();
    /** @see TransFuncPropertyWidgetPainterBaseSlider::isHit */
    bool isHit(tgt::ivec2 pos);
    //--------------------
    //  getter and setter
    //--------------------
    /** Sets the slider position. Calls updateStoredValue(). */
    void setPosition(tgt::ivec2 pos);
    /** Returns the relevant x/y position of the gamma slider in pixel coordinates. */
    int getPosition();
    /** Sets the gamma value. Calls updateSliderPosition(). */
    void setGammaValue(float gamma);
    /** Returns the stored gamma value. */
    float getGammaValue();
protected:
    /** @see TransFuncPropertyWidgetPainterBaseSlider::updateSliderPosition */
    void updateSliderPosition();
    /** @see TransFuncPropertyWidgetPainterBaseSlider::updateStoredValue */
    void updateStoredValue();
private:
    //--------------------
    //  member
    //--------------------
    int position_;  ///< position of the slider in pixel
    float gamma_;   ///< gamma value
signals:
    /** Signal emited by updateStoredValue. */
    void gammaChanged(float gamma, size_t dimension);
};

/***********************************************************
 * Domain Slider                                           *
 ***********************************************************/
class TransFuncPropertyWidgetPainterDomainSlider : public QObject, public TransFuncPropertyWidgetPainterBaseSlider {
    Q_OBJECT
public:
    /** Constructor */
    TransFuncPropertyWidgetPainterDomainSlider(size_t associatedDimension, SliderOrientation orientation, int sliderSize);
    /** @see TransFuncPropertyWidgetPainterBaseSlider::paint */
    void paint();
    /** @see TransFuncPropertyWidgetPainterBaseSlider::isHit */
    bool isHit(tgt::ivec2 pos);
    /** Checks, if left slider is hit. */
    bool isLeftHit(tgt::ivec2 pos);
    /** Checks, if right slider is hit. */
    bool isRightHit(tgt::ivec2 pos);
    /** Checks, if the interval between both sliders is hit. */
    bool isIntervalHit(tgt::ivec2 pos);
    //--------------------
    //  getter and setter
    //--------------------
    /** Sets the left slider position. Calls updateStoredValue(). */
    void setLeftPosition(tgt::ivec2 pos);
    /** Sets the right slider position. Calls updateStoredValue(). */
    void setRightPosition(tgt::ivec2 pos);
    /** Returns both relevant x/y slider postitions. */
    tgt::ivec2 getPosition();
    /** Sets the domainvalue. Calls updateSliderPosition(). */
    void setDomainValue(tgt::vec2 domain);
    /** returns the stored domain value. */
    tgt::vec2 getDomainValue();
protected:
    /** @see TransFuncPropertyWidgetPainterBaseSlider::updateSliderPosition */
    void updateSliderPosition();
    /** @see TransFuncPropertyWidgetPainterBaseSlider::updateStoredValue */
    void updateStoredValue();
private:
    //--------------------
    //  member
    //--------------------
    int leftPosition_;      ///< left slider position
    int rightPosition_;     ///< right slider position
    tgt::vec2 domain_;      ///< stored domain
signals:
    /** Signal emited by updateStoredValue() */
    void domainChanged(tgt::vec2 domain, size_t dimension);
};


} // namespace voreen

#endif // VRN_TRANSFUNCIOHELPERQT_H
