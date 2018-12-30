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

#ifndef VRN_FLOATMAT4PROPTERYWIDGET_H
#define VRN_FLOATMAT4PROPTERYWIDGET_H

#include "voreen/core/properties/matrixproperty.h"
#include "voreen/qt/widgets/property/matrixpropertywidget.h"

#include <QLineEdit>
#include <QValidator>
#include <QWidget>

namespace voreen {

class FloatMat4PropertyWidget : public MatrixPropertyWidget {
Q_OBJECT
public:
    FloatMat4PropertyWidget(FloatMat4Property*, QWidget* = 0);

protected:
    virtual std::string getExtension() const { return "vm4"; }
    virtual void saveMatrix(const std::string& filename) const;
    virtual void loadMatrix(const std::string& filename);

    QDoubleValidator* doubleValidator_;
    FloatMat4Property* prop_;
    QLineEdit* x0_;
    QLineEdit* y0_;
    QLineEdit* z0_;
    QLineEdit* w0_;
    QLineEdit* x1_;
    QLineEdit* y1_;
    QLineEdit* z1_;
    QLineEdit* w1_;
    QLineEdit* x2_;
    QLineEdit* y2_;
    QLineEdit* z2_;
    QLineEdit* w2_;
    QLineEdit* x3_;
    QLineEdit* y3_;
    QLineEdit* z3_;
    QLineEdit* w3_;
    QLineEdit* x4_;
    QLineEdit* y4_;
    QLineEdit* z4_;
    QLineEdit* w4_;

protected slots:
    virtual void updateFromPropertySlot();
    void updateValue();
    void identity();
};

}// namespace voreen

#endif // FLOATMAT4PROPTERYWIDGET
