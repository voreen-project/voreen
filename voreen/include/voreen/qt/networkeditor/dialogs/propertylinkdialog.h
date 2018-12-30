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

#ifndef VRN_PROPERTYLINKDIALOG_H
#define VRN_PROPERTYLINKDIALOG_H

#include <QDialog>

//class QGraphicsView;
class QPushButton;
class QComboBox;
class QButtonGroup;
class QAbstractButton;

namespace voreen {

    class SplitGraphicsView;
    class PortOwnerGraphicsItem;
    class PropertyGraphicsItem;
    class PropertyOwnerGraphicsItem;
    class PropertyLinkArrowGraphicsItem;
    class PropertyListGraphicsItem;
    class LinkEvaluatorBase;

/**
 * This dialog is been used to link properties of one or two processors.
 */
class PropertyLinkDialog : public QDialog {
Q_OBJECT
protected:

    /**
     * This struct is used to memorize changes in the linking structure.
     * Since changes should only apply, if the ok button is pressed,
     * new, deleted and modified links must be stored temporal.
     * On ok all changes apply. On undo, all changes are droped.
     *
     */
    struct Modification {
        PropertyGraphicsItem* src_; ///< source property item
        PropertyGraphicsItem* dst_; ///< destination property item
        LinkEvaluatorBase* eval_;   ///< the link evaluator between both properties

        Modification(PropertyGraphicsItem* src, PropertyGraphicsItem* dst, LinkEvaluatorBase* eval)
            : src_(src), dst_(dst), eval_(eval) {};

        /** Compare operator, ignoring the evaluator */
        bool operator== (Modification& mod) {
            return (src_ == mod.src_ && dst_ == mod.dst_);
        }
    };
public:
    //------------------------------
    //  constructor and destructor
    //------------------------------
    PropertyLinkDialog(QWidget* parent, PortOwnerGraphicsItem* source, PortOwnerGraphicsItem* destination);
    ~PropertyLinkDialog();
    //------------------------------
    //  general members
    //------------------------------
protected:
    PortOwnerGraphicsItem* sourceItem_;         ///< pointer to source item (left)
    PortOwnerGraphicsItem* destinationItem_;    ///< pointer to destination item (right or sourceItem_ in only one view)

    //---------------------------------------------------------------------------------
    //      layout functions
    //---------------------------------------------------------------------------------
protected:
    /** Initializes the layout of the dialog. */
    void initializeLayout(bool hasOnlyOneView);

private:
    //------------------------------
    //  layout members
    //------------------------------
    SplitGraphicsView* sourceView_;                     ///< left QGraphicsView
    SplitGraphicsView* destinationView_;                ///< right QGraphicsView (can be null if only one view is active)
    QComboBox* evaluatorSrcDstBox_;                     ///< evaluator combobox from left to right
    QComboBox* evaluatorDstSrcBox_;                     ///< evaluator combobox from right to left
    QButtonGroup* arrowButtonGroup_;                    ///< group of all three direction buttons
        QAbstractButton* lastPressedDirectionButton_;   ///< stores the last pressed button to check, if selected button is reselected
    QPushButton* leftUpArrowButton_;                    ///< link from right to left or from down to up in single view
    QPushButton* bidirectionalArrowButton_;             ///< bidirectional link
    QPushButton* rightDownArrowButton_;                 ///< link from left to right or from up to down in single view
    QPushButton* chainButton_;                          ///< determines if two evaluators are chained
        bool areEvaluatorsChained_;                     ///< value of chainButton_

    //---------------------------------------------------------------------------------
    //      selection functions
    //---------------------------------------------------------------------------------
private:
    /** Returns the currently selected arrow in the source view. */
    PropertyLinkArrowGraphicsItem* getCurrentlySelectedArrow() const;
    /** Selects an arrow. Is used to auto-select new generated links. */
    void setCurrentlySelectedArrow(PropertyGraphicsItem* src, PropertyGraphicsItem* dst);

public slots:
    /** Is called if an arrow is selected or deselected by the user. */
    void sceneSelectionChanged();

    //---------------------------------------------------------------------------------
    //      evaluator functions
    //---------------------------------------------------------------------------------
protected:
    /** Returns the actual(temporal) evaluator between these two properties. */
    LinkEvaluatorBase* getEvaluator(PropertyGraphicsItem* src, PropertyGraphicsItem* dst);
    /**
     * Updates the evaluator boxes depending on arrow.
     * Calls updateSrcDstEvaluatorBox and updateDstSrcEvaluatorBox.
     */
    void updateEvaluatorBoxes(PropertyLinkArrowGraphicsItem* arrow);
    /** Fills the combo box with all possible evaluators. */
    void fillEvaluatorBox(QComboBox* box, PropertyGraphicsItem* src = 0,PropertyGraphicsItem* dst = 0);
    /** Called by the two changed slots. */
    void evaluatorBoxSelectionChanged(QComboBox* box, int index);
protected slots:
    /** Called if the selection has changed. */
    void srcDstEvaluatorBoxSelectionChanged(int index);
    /** Called if the selection has changed. */
    void dstSrcEvaluatorBoxSelectionChanged(int index);
    //chain button
protected:
    /** Changes the icon of the button. Is called by 'toggleChainButton'. */
    void toggleChainButton(bool b);
protected slots:
    /** Is triggered when clicked. */
    void toggleChainButton();

    //---------------------------------------------------------------------------------
    //      direction button functions
    //---------------------------------------------------------------------------------
protected:
    /** Updates the direction buttons depending on arrow. */
    void updateArrowButtonGroup(PropertyLinkArrowGraphicsItem* arrow);
protected slots:
    /** Triggered if a direction button is clicked. */
    void directionButtonClicked(QAbstractButton* button);

    //---------------------------------------------------------------------------------
    //      delete functions
    //---------------------------------------------------------------------------------
protected slots:
    /** Triggered by ckicking DEL. */
    void deleteArrowAction();
    /** Triggered by clicking delete all processor links button. */
    void deleteAllProcessorLinks();
    /** Triggered by clicking delete all port links button. */
    void deleteAllPortLinks();

    //---------------------------------------------------------------------------------
    //      auto linking functions
    //---------------------------------------------------------------------------------
protected slots:
    /** Triggered if mouse hovers the auto link button (processor) */
    void showProcessorAutoLinks();
    /** Triggered if mouse hovers the auto link button (port) */
    void showPortAutoLinks();
    /** Triggered if mouse leaves the auto link button (both) */
    void hideAutoLinks();
    /** Triggered if the auto link button is clicked (both) */
    void confirmAutoLinks();
protected:
    QList<PropertyLinkArrowGraphicsItem*> autoArrows_;  ///< list of temporal added auto links

    //---------------------------------------------------------------------------------
    //      settings functions
    //---------------------------------------------------------------------------------
protected:
    /** Saves all dialog window settings. Is called on closeEvent. */
    void saveSettings();
    /** Loads all dialog window settings. Is called in the constructor. */
    void loadSettings();

    //---------------------------------------------------------------------------------
    //      modification functions
    //---------------------------------------------------------------------------------
protected:
    /** Adds an add-modification to the add-list. */
    void addAddModification(Modification* mod);
    /** Adds an remove modification to the remove-list. */
    void addRemoveModification(Modification* mod);
    /** Adds an evaluator-modification to the evaluator-list. */
    void addEvaluatorModification(Modification* mod);
protected slots:
    /** Applys all temporal modifications. */
    void applyAllModifications();
    /** Undoes all temporal modifications. */
    void undoAllModifications();
    /** Links two properties temporal. Is called by connecting two properties in the editor. */
    void linkProperties(PropertyGraphicsItem* src, PropertyGraphicsItem* dst);

private:
    //member
    QList<Modification*> addModificationList_;          ///< list of all add modifications
    QList<Modification*> removeModificationList_;       ///< list of all remove modifications
    QList<Modification*> evaluatorModificationList_;    ///< list of all evaluator changed modifications

    //---------------------------------------------------------------------------------
    //      events
    //---------------------------------------------------------------------------------
protected:
    /** Handling all key shortcuts */
    void keyPressEvent(QKeyEvent* event);
    /** Used to save the window settings */
    void closeEvent(QCloseEvent* event);

};

} //namespace voreen

#endif // VRN_PROPERTYLINKDIALOG_H

