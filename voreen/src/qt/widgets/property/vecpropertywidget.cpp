/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/property/vecpropertywidget.h"

#include "voreen/core/properties/vectorproperty.h"
#include "voreen/qt/widgets/sliderspinboxwidget.h"

#include <QMouseEvent>
#include <QMenu>

namespace voreen {

template<class WIDGETTYPE, class VECTORPROP, typename ELEMTYPE>
VecPropertyWidget<WIDGETTYPE, VECTORPROP, ELEMTYPE>::VecPropertyWidget(VECTORPROP* const prop,
                                                                       const size_t numComponents,
                                                                       QWidget* parent)
    : QPropertyWidget(prop, parent),
    numComponents_(numComponents),
    myLayout_(new QVBoxLayout()),
    widgets_(0),
    vectorProp_(prop)
{
    myLayout_->setSpacing(0);
    myLayout_->setMargin(1);
    myLayout_->setContentsMargins(0, 0, 0, 0);

    widgets_ = new WIDGETTYPE*[numComponents_];
    for (size_t i = 0; i < numComponents_; ++i) {
        widgets_[i] = new WIDGETTYPE(this);
        widgets_[i]->setSliderTracking(vectorProp_->hasTracking());
        myLayout_->addWidget(widgets_[i], 1);
        connect((const QObject*) widgets_[i], SIGNAL(sliderPressedChanged(bool)), this, SLOT(toggleInteractionMode(bool)));
    }

    QPropertyWidget::addLayout(myLayout_);
    updateFromPropertySlot();

    precisionMenu_ = new QMenu(this);
    highAction_ = precisionMenu_->addAction("High Precision");
    highAction_->setCheckable(true);
    precisionMenu_->addSeparator();
    instantValueChangeAction_ = precisionMenu_->addAction("Tracking Mode");
    instantValueChangeAction_->setCheckable(true);
    if(vectorProp_->hasTracking())
        instantValueChangeAction_->toggle();
}

template<class WIDGETTYPE, class VECTORPROP, typename ELEMTYPE>
VecPropertyWidget<WIDGETTYPE, VECTORPROP, ELEMTYPE>::~VecPropertyWidget() {
    for (size_t i = 0; i < numComponents_; ++i) {
        delete widgets_[i];
    }
    delete[] widgets_;
    delete myLayout_;
}

template<class WIDGETTYPE, class VECTORPROP, typename ELEMTYPE>
void VecPropertyWidget<WIDGETTYPE, VECTORPROP, ELEMTYPE>::updateFromPropertySlot() {
    const typename VECTORPROP::ElemType& values = vectorProp_->get();
    const typename VECTORPROP::ElemType& minValues = vectorProp_->getMinValue();
    const typename VECTORPROP::ElemType& maxValues = vectorProp_->getMaxValue();
    const typename VECTORPROP::ElemType& steppings = vectorProp_->getStepping();

    if(typeid(vectorProp_) == typeid(FloatVec4Property*) || typeid(vectorProp_) == typeid(FloatVec3Property*) || typeid(vectorProp_) == typeid(FloatVec2Property*)){
        int decimals = vectorProp_->getNumDecimals();
        for (size_t i = 0; i < numComponents_; ++i) {
            dynamic_cast<DoubleSliderSpinBoxWidget*>(widgets_[i])->setDecimals(decimals);
        }
    }

    for (size_t i = 0; i < numComponents_; ++i) {
        widgets_[i]->blockSignals(true);
        widgets_[i]->setMaxValue(maxValues[i]);
        widgets_[i]->setMinValue(minValues[i]);
        widgets_[i]->setValue(values[i]);
        widgets_[i]->setSingleStep(steppings[i]);
        widgets_[i]->setSliderTracking(vectorProp_->hasTracking());
        widgets_[i]->blockSignals(false);
    }
}

template<class WIDGETTYPE, class VECTORPROP, typename ELEMTYPE>
typename VECTORPROP::ElemType
VecPropertyWidget<WIDGETTYPE, VECTORPROP, ELEMTYPE>::setPropertyComponent(QObject* sender,
                                                                          ELEMTYPE value)
{
    typename VECTORPROP::ElemType newValue = vectorProp_->get();
    for (size_t i = 0; ((sender != 0) && (i < numComponents_)); ++i) {
        if (widgets_[i] == sender) {
            newValue[i] = value;
            vectorProp_->set(newValue);
            emit valueModifiedByUser();
            break;
        }
    }
    return newValue;
}

template<class WIDGETTYPE, class VECTORPROP, typename ELEMTYPE>
void VecPropertyWidget<WIDGETTYPE, VECTORPROP, ELEMTYPE>::mousePressEvent(QMouseEvent* event) {
    if(event->button() == Qt::RightButton) {
        QAction* prec = precisionMenu_->exec(QCursor::pos());
        if(prec == highAction_){
            if(highAction_->isChecked()) {
                if(std::is_same<typename VECTORPROP::ElemType::ElemType, float>()) {
                    vectorProp_->setStepping(typename VECTORPROP::ElemType(0.0001f));
                } else {
                    vectorProp_->setStepping(typename VECTORPROP::ElemType(1));
                }
                vectorProp_->setNumDecimals(4);
            }
            else {
                if(std::is_same<typename VECTORPROP::ElemType::ElemType, float>()) {
                    vectorProp_->setStepping(typename VECTORPROP::ElemType(0.05f));
                } else {
                    vectorProp_->setStepping(typename VECTORPROP::ElemType(1));
                }
                vectorProp_->setNumDecimals(2);
            }
        }
        else if(prec == instantValueChangeAction_) {
            vectorProp_->setTracking(!vectorProp_->hasTracking());
        }
        updateFromPropertySlot();
    }
    QWidget::mousePressEvent(event);
}

template class VecPropertyWidget<SliderSpinBoxWidget, IntVec2Property, int>;
template class VecPropertyWidget<SliderSpinBoxWidget, IntVec3Property, int>;
template class VecPropertyWidget<SliderSpinBoxWidget, IntVec4Property, int>;
template class VecPropertyWidget<DoubleSliderSpinBoxWidget, FloatVec2Property, float>;
template class VecPropertyWidget<DoubleSliderSpinBoxWidget, FloatVec3Property, float>;
template class VecPropertyWidget<DoubleSliderSpinBoxWidget, FloatVec4Property, float>;

}   // namespace
