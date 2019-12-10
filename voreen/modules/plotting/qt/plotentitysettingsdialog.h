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

#ifndef VRN_PLOTENTITYSETTINGSDIALOG_H
#define VRN_PLOTENTITYSETTINGSDIALOG_H

#include "../datastructures/plotentitysettings.h"

#include <QDialog>

class QCheckBox;
class QComboBox;
class QTabWidget;
class QDoubleSpinBox;
class QPushButton;
class QDialogButtonBox;

namespace voreen {

class ColorSelectorWidget;
class PlotEntitiesProperty;
class PlotData;

class PlotEntitySettingsDialog : public QDialog {
    Q_OBJECT
public:

    /**
     * \brief   Constructor for a new Dialog offering editing opportunities for a PlotEntitySettings.
     *
     * \param   entitySettings  PlotEntitySettings to edit (initial value), dialog will edit a copy
     * \param   plotData        PlotData reference for determining the column names
     * \param   prop            PlotEntitiesProperty
     * \param   parent          parent widget
     *
     * \note    This dialog will edit a copy of \a entry.
     **/
    PlotEntitySettingsDialog(const PlotEntitySettings& entitySettings, const PlotData* plotData,
            const PlotEntitiesProperty* prop, QWidget* parent = 0);

    // convert tgtcolor to qcolor and vice versa
    tgt::Color toTgtColor(QColor color);
    QColor toQColor(tgt::Color color);

    /**
     * Returns edited PlotEntitySettings.
     **/
    const PlotEntitySettings& getEntitySettings() const;

private slots:
    void firstColorChanged(QColor);
    void secondColorChanged(QColor);
    void cbMainColumnChanged(int index);
    void cbCandleTopColumnChanged(int index);
    void cbCandleBottomColumnChanged(int index);
    void cbStickTopColumnChanged(int index);
    void cbStickBottomColumnChanged(int index);
    void cbOptionalColumnChanged(int index);
    void cbSecondOptionalColumnChanged(int index);
    void cbStyleChanged(int index);
    void chbSplineChanged();
    void chbErrorbarChanged();
    void chbWireOnlyChanged();
    void chbHeightMapChanged();
    void cbColorMapChanged(int index);
    void tabChanged(int index);
    void sbMinGlyphSizeChanged(double value);
    void sbMaxGlyphSizeChanged(double value);
    void btOpenFileDialogClicked();

private:

    /// creates and layouts all widgets
    void createWidgets();
    /// updates widgets e.g. selections
    void updateWidgets();

    // fills the ComboBox cb with number columns of the associated plotdata
    // if empty entry is true, an empty entry will be inserted first (with value -1)
    void fillComboBox(QComboBox* cb, bool emptyEntry = false);

    /// sets the text of btOpenFileDialog
    void updateBtOpenFileDialogText(const std::string& text);

    // Qt widgets:
    QComboBox* cbMainColumn_;
    QComboBox* cbCandleTopColumn_;
    QComboBox* cbCandleBottomColumn_;
    QComboBox* cbStickTopColumn_;
    QComboBox* cbStickBottomColumn_;
    QComboBox* cbOptionalColumn_;
    QComboBox* cbSecondOptionalColumn_;
    QComboBox* cbStyle_;
    QCheckBox* chbSpline_;
    QCheckBox* chbErrorbar_;
    QCheckBox* chbWireOnly_;
    QCheckBox* chbHeightMap_;
    QComboBox* cbColorMap_;
    QTabWidget* tabWidget_;
    ColorSelectorWidget* cclrFirst_;
    ColorSelectorWidget* cclrSecond_;
    ColorSelectorWidget* cclrThird_;
    ColorSelectorWidget* cclrFourth_;
    QDoubleSpinBox* sbMinGlyphSize_;
    QDoubleSpinBox* sbMaxGlyphSize_;
    QPushButton* btOpenFileDialog_;
    QDialogButtonBox* buttonBox_;

    // private members
    PlotEntitySettings entitySettings_; ///< PlotEntitySettings to edit
    const PlotData* plotData_;          ///< reference to current PlotData (used for column names)
    const PlotEntitiesProperty* prop_;  ///< reference to PlotEntitiesProperty
};

} // namespace

#endif // VRN_PLOTENTITYSETTINGSDIALOG_H
