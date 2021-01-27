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

#ifndef VRN_TRANSFUNCPROPERTYWIDGETBASE_H
#define VRN_TRANSFUNCPROPERTYWIDGETBASE_H

#include "voreen/qt/widgets/property/qpropertywidgetwithtoolwindow.h"

#include "voreen/core/datastructures/transfunc/transfuncbase.h"
#include "voreen/core/datastructures/volume/volume.h"

#include <QPixmap>

    class QToolButton;
    class QMenu;

namespace tgt {
    class QtCanvas;
}

namespace voreen {

class Histogram;
class TransFuncPropertyEditorBase;
class TransFuncPropertyBase;
class TransFuncPropertyWidgetPainterBase;

/**
 * Base transfer function widget. It provides the base buttons, every widget should implement.
 */
class VRN_QT_API TransFuncPropertyWidgetBase : public QPropertyWidgetWithToolWindow {
    Q_OBJECT
protected:
    /** Defines the base layout of the widget. */
    enum WidgetBaseLayout {
        WBL_1D,
        WBL_2D
    };
    /** @note Implement in sub-class. */
    virtual WidgetBaseLayout getBaseLayout() const = 0;
    /** @note Implement in sub-class. */
    virtual TransFuncPropertyWidgetPainterBase* createPreviewCanvasPainter(tgt::QtCanvas* canvas) const = 0;
    /** @note Implement in sub-class. */
    virtual std::string getPathSubFolder() const = 0;
    /** @note Implement in sub-class. */
    virtual Histogram* getHistogramFromVolume(const VolumeBase* vb, const size_t channel) const = 0;
    /** @note Implement in sub-class. */
    virtual void restoreZoomMetaData() const = 0;
protected slots:
    /** @note Implement in sub-class. */
    virtual void storeZoomMetaData() const = 0;
public:
    /**
     * Constructor
     * @param prop the associated trans func property.
     */
    TransFuncPropertyWidgetBase(TransFuncPropertyBase* prop, QWidget* parent);

    /** Destructor */
    ~TransFuncPropertyWidgetBase();

    //--------------------
    //  update functions
    //--------------------
protected:
    /** Updates the property and sets the menu icon */
    void updateAlpha(TransFuncBase::AlphaMode mode);
    /**
     * Sets the window to floating.
     * @override QPropertyWidgetWithToolWindow
     */
    virtual void customizeToolWindow();
protected slots:
    /**
     * Updates the Widget from the property.
     * @override QPropertyWidget
     */
    virtual void updateFromPropertySlot();
    /**
     * Shows the tooltip.
     * Is connected to SIGNAL::TransFuncPropertyWidgetPainterBase::showInfoToolTip.
     */
    void showToolTipSlot(QPoint pos, QString tip);
    /**
     * Hides the tooltip.
     * Is connected to SIGNAL::TransFuncPropertyWidgetPainterBase::hideInfoToolTip.
     */
    void hideToolTipSlot();

    //--------------------
    //  gui functions
    //--------------------
protected:
    /**
     * Calls the layout functions.
     * @override PropertyWidget
     */
    virtual void initialize();

    /**
     * Returns null, since this widget does not need a separate label.
     * @override QPropertyWidget
     * @see PropertyOwnerWidget
     */
    virtual CustomLabel* getOrCreateNameLabel() const;

    /**
     * Creates and layouts all basic gui elements.
     */
    void initializeBaseLayout();
    /** @see initializeBaseLayout() */
    void initialize1DLayout(QLayout* mainLayout);
    /** @see initializeBaseLayout() */
    void initialize2DLayout(QLayout* mainLayout);

    //--------------------
    //  menu handling
    //--------------------
protected:
    /** Creates the icon shown in the preset menu. */
    virtual const QPixmap createPresetPreviewIcon(std::string& presetPath) const = 0;
protected slots:
    /** Generates the preset menu on click. */
    void createPresetMenuSlot();
    /** Loads the preset. Is called from the preset menu. */
    virtual void loadPresetSlot(QAction* action) = 0;
    /** Generates the domain menu on click. */
    virtual void createDomainMenuSlot();
    /** Loads the domain. Is called from the menu. */
    virtual void loadDomainSlot(QAction* action);
    /** Fit transfer function to data domain. */
    void fitToDataSlot();
    /** Creates the alpha menu. */
    void createAlphaMenuSlot();
    /** action being clled by the alphaMenu_. */
    void loadAlphaSlot(QAction* action);
    /** Called by the zoom in button. */
    void zoomInSlot();
    /** Called by the zoom out button. */
    void zoomOutSlot();
    /** Helper slot to reset the zoom. */
    void resetZoomSlot();
    /** Creates the advanced menu. */
    void createAdvancedMenuSlot();
    /** Action being clled by the advancedMenu_. */
    void doAdvancedActionSlot(QAction* action);

    //--------------------
    //  gui member
    //--------------------
protected:
    QToolButton* presetButton_;     ///< button to load tf presets
        QMenu* presetMenu_;         ///< tf menu
    QToolButton* domainButton_;     ///< button to load domain windows
        QMenu* domainMenu_;        ///< domain window menu
    QToolButton* alphaButton_;      ///< button to enable alpha
        QMenu* alphaMenu_;          ///< menu for alpha value

    QToolButton* zoomInButton_;     ///< button to zoom in
    QToolButton* zoomOutButton_;    ///< button to zoom out

    QToolButton* advancedButton_;   ///< button with advanced settings
        QMenu* advancedMenu_;       ///< advanced menu

    tgt::QtCanvas* previewCanvas_;  ///< canvas to visualize the tf texture
    TransFuncPropertyWidgetPainterBase* previewPainter_;  ///< painter of the previewCanvas

    //--------------------
    //  member
    //--------------------
    TransFuncPropertyBase* property_;       ///< property belonging to this widget
    TransFuncPropertyEditorBase* editor_;   ///< editor of this widget

};

} // namespace voreen

#endif // VRN_TRANSFUNCPROPERTYWIDGETBASE_H
