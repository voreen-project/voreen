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

#ifndef VRN_RAWVOLUMEWIDGET_H
#define VRN_RAWVOLUMEWIDGET_H

#include <QDialog>

#include "tgt/vector.h"
#include "tgt/matrix.h"

class QComboBox;
class QDoubleSpinBox;
class FloatMat4Property;
class QLabel;
class QPushButton;
class QSpinBox;

namespace voreen {

class RawVolumeWidget : public QDialog {
    Q_OBJECT
public:
    RawVolumeWidget(QWidget* parent, const QString& filename, std::string& objectModel, std::string& format, int& numFrames,
        tgt::ivec3& dim, tgt::vec3& spacing, int& headerSkip, bool& bigEndian, tgt::mat4& trafoMat,  int fixedZDim = -1);
    ~RawVolumeWidget();

protected:
    void resizeEvent(QResizeEvent *);

private:
    QPushButton* submit_;
    QComboBox* datatypeComboBox_;
    QComboBox* objectModelComboBox_;
    QSpinBox* numTimeFramesSpin_;
    QSpinBox* headerSkipSpin_;
    QComboBox* endiannessCombo_;
    QSpinBox* xDimension_;
    QSpinBox* yDimension_;
    QSpinBox* zDimension_;
    QDoubleSpinBox* xSpacing_;
    QDoubleSpinBox* ySpacing_;
    QDoubleSpinBox* zSpacing_;

    std::string& objectModel_;
    std::string& format_;
    int& numFrames_;
    tgt::ivec3& dim_;
    tgt::vec3& spacing_;
    tgt::mat4& trafoMat_;
    FloatMat4Property* trafoMatProp_;
    QWidget* trafoMatWidget_;
    int& headerSkip_;
    bool& bigEndian_;

private slots:
    void updateValues();
};

} // namespace voreen

#endif // VRN_RAWVOLUMEWIDGET_H
