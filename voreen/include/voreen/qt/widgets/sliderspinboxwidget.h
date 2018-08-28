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

#ifndef VRN_SLIDERSPINBOXWIDGET_H
#define VRN_SLIDERSPINBOXWIDGET_H

#include <QWidget>

class QDoubleSpinBox;
class QHBoxLayout;
class QLineEdit;
class QSlider;
class QSpinBox;

#include "voreen/core/properties/property.h"
#include "voreen/qt/voreenqtapi.h"

namespace voreen {

class VRN_QT_API SliderSpinBoxWidget : public QWidget {
Q_OBJECT
Q_PROPERTY(int value READ getValue WRITE setValue)
Q_PROPERTY(int minValue READ getMinValue WRITE setMinValue)
Q_PROPERTY(int maxValue READ getMaxValue WRITE setMaxValue)

public:
    SliderSpinBoxWidget(QWidget* parent = 0);
    int getValue() const;
    int getMinValue() const;
    int getMaxValue() const;
    bool isSliderDown() const;
    void setFocusPolicy(Qt::FocusPolicy policy);

    /**
     * If slider tracking is disabled, valueChanged signals
     * are not fired during user interaction, i.e.,
     * while the user is dragging the slider.
     * Instead, valueChanged is emitted after
     * the slider has been released.
     * Slider tracking is enabled by default.
     *
     * @note The spin box is always kept in sync with
     *  the slider, regardless of the slider tracking mode.
     *
     * @see QAbstractSlider
     */
    void setSliderTracking(bool tracking);

    /// Returns whether slider tracking is enabled (default: true).
    bool hasSliderTracking() const;

    /**
     * If spin box tracking is disabled, valueChanged signals
     * are not fired during user interaction, i.e.,
     * while the user is typing into the spinbox.
     * Instead, valueChanged is emitted on return
     * or when the spin box looses focus.
     * Spinbox tracking is disabled by default.
     *
     * @see QAbstractSpinbox
     */
    virtual void setSpinboxTracking(bool tracking);

    /// Returns whether spin box tracking is enabled (default: false).
    bool hasSpinboxTracking() const;

signals:
    void valueChanged(int);
    void sliderPressedChanged(bool);
    void editingFinished();

public slots:
    void setValue(int value);
    void setMaxValue(int value);
    void setMinValue(int value);
    void setSingleStep(int value);

private slots:
    void spinEditingFinished();
    void sliderPressed();
    void sliderReleased();

protected:
    void changeEvent(QEvent*);

    QHBoxLayout* layout_;
    QSlider* slider_;
    QSpinBox* spinbox_;

    bool isSliderTrackingEnabled_;
    bool isSpinboxTrackingEnabled_;
    int value_;
};


// ---------------------------------------------------------------------------


class VRN_QT_API DoubleSliderSpinBoxWidget : public QWidget {
Q_OBJECT
Q_PROPERTY(double value READ getValue WRITE setValue)
Q_PROPERTY(double minValue READ getMinValue WRITE setMinValue)
Q_PROPERTY(double maxValue READ getMaxValue WRITE setMaxValue)

public:
    DoubleSliderSpinBoxWidget(QWidget* parent = 0);
    double getValue() const;
    double getMinValue() const;
    double getMaxValue() const;
    double getSingleStep() const;
    int getDecimals() const;
    bool isSliderDown() const;
    void setFocusPolicy(Qt::FocusPolicy policy);

    /**
     * If slider tracking is disabled, valueChanged signals
     * are not fired during user interaction, i.e.,
     * while the user is dragging the slider.
     * Instead, valueChanged is emitted after
     * the slider has been released.
     * Slider tracking is enabled by default.
     *
     * @note The spin box is always kept in sync with
     *  the slider, regardless of the slider tracking mode.
     *
     * @see QAbstractSlider
     */
    void setSliderTracking(bool tracking);

    /// Returns whether slider tracking is enabled (default: true).
    bool hasSliderTracking() const;

    /**
     * If spin box tracking is disabled, valueChanged signals
     * are not fired during user interaction, i.e.,
     * while the user is typing into the spinbox.
     * Instead, valueChanged is emitted on return
     * or when the spin box looses focus.
     * Spinbox tracking is disabled by default.
     *
     * @see QAbstractSpinbox
     */
    void setSpinboxTracking(bool tracking);

    /// Returns whether spin box tracking is enabled (default: false).
    bool hasSpinboxTracking() const;

signals:
    void valueChanged(double);
    void sliderPressedChanged(bool);
    void editingFinished();

public slots:
    void setValue(double value);
    void setMaxValue(double value);
    void setMinValue(double value);
    void setSingleStep(double step);
    void setDecimals(int decimals);
    void sliderPressed();
    void sliderReleased();

protected slots:
    void sliderValueChanged(int value);
    void spinEditingFinished();
    void adjustSliderScale();

protected:
    void changeEvent(QEvent*);

    QHBoxLayout* layout_;
    QSlider* slider_;
    QDoubleSpinBox* spinbox_;

    bool isSliderTrackingEnabled_;
    bool isSpinboxTrackingEnabled_;
    double value_;
    /// The actual min/max values. Min/max of widgets can be different because of step values.
    double min_;
    double max_;


    static const int SLIDER_STEPS;
};


} // namespace voreen;

#endif //VRN_SLIDERSPINBOXWIDGET_H
