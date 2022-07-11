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

#ifndef VRN_NETWORKEDITOR_H
#define VRN_NETWORKEDITOR_H

#include "voreen/core/properties/propertyowner.h"
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/network/workspace.h"
#include "voreen/qt/voreenqtapi.h"

#include "voreen/qt/networkeditor/editor_settings.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"

#include <QGraphicsView>
#include <QMap>

class QToolButton;
class QAction;
class QTextEdit;

namespace voreen {
//voreen core
class NetworkEditor;
class NetworkEvaluator;
//style
class NWEStyle_Base;
//graph layouts
class NWEGL_Base;
//graphic items
class ProcessorGraphicsItem;
class PropertyGraphicsItem;
class PortGraphicsItem;
class PortOwnerGraphicsItem;
class PortArrowGraphicsItem;
class PropertyLinkArrowGraphicsItem;
class PortOwnerLinkArrowGraphicsItem;
class PortSizeLinkArrowGraphicsItem;
class TextBoxBaseGraphicsItem;
class WorkspaceDescriptionGraphicsItem;

/// Editor Layers
enum NetworkEditorLayer {
    NetworkEditorLayerUndefined,
    NetworkEditorLayerDataFlow,
    NetworkEditorLayerLinking,
    NetworkEditorLayerGeneralLinking,
    NetworkEditorLayerCameraLinking,
    NetworkEditorLayerPortSizeLinking
};
enum NetworkEditorCursorMode {
    NetworkEditorCursorSelectMode,
    NetworkEditorCursorMoveMode
};


class VRN_QT_API NetworkEditor : public QGraphicsView, public ProcessorNetworkObserver, public PropertyOwner {
    Q_OBJECT
        //friend classes to call protected functions
        friend class SaveWorkspaceMenuEntity;
        friend class SaveWorkspaceAsMenuEntity;
public:
    NetworkEditor(QWidget* parent = 0, NetworkEvaluator* evaluator = 0);
    ~NetworkEditor();

    virtual NetworkEditor* create() const {return new NetworkEditor();}
    virtual std::string getClassName() const {return "NetworkEditor";}

//---------------------------------------------------------------------------------------------------------------
//                  general members and functions
//---------------------------------------------------------------------------------------------------------------
protected:
    Workspace* workspace_;                  /// assigned workspace whose network is to be edited
    NetworkEvaluator* evaluator_;           /// evaluator

    /// connection between processors and their graphicitems
    QMap<Processor*,ProcessorGraphicsItem*> processorItemMap_;

    static const std::string loggerCat_;    /// loggarcat

    /// properties
        //style
    IntProperty scaleProcessorFontSizeProperty_;
    OptionProperty<NetworkEditorStyles> networkEditorStyleProperty_;
        //layout
    OptionProperty<NetworkEditorGraphLayouts> networkEditorGraphLayoutsProperty_;
    FloatProperty sugiShiftXProperty_;
    BoolProperty sugiOverlapProperty_;
    BoolProperty sugiMedianProperty_;
    BoolProperty sugiPortFlushProperty_;
        //settings
    BoolProperty showDocumentationProperty_;
public:
    void setWorkspace(Workspace* workspace);
    Workspace* getWorkspace() const;
    ProcessorNetwork* getProcessorNetwork() const;
    float getProcessorFontSizeScale() const;

//---------------------------------------------------------------------------------------------------------------
//                  ProcessorNetworkObserver functions
//---------------------------------------------------------------------------------------------------------------
public:
    void networkChanged();
    void processorAdded(const Processor* processor);
    void processorRemoved(const Processor* processor);
    void propertyLinkAdded(const PropertyLink* link);
    void propertyLinkRemoved(const PropertyLink* link);
    void processorRenamed(const Processor* processor, const std::string& prevName);
    void portConnectionAdded(const Port* outport, const Port* inport);
    void portConnectionRemoved(const Port* outport, const Port* inport);

public slots:
    void processorAdded(QString id, Processor* selectedProcessor);

//---------------------------------------------------------------------------------------------------------------
//                  scene transformations
//---------------------------------------------------------------------------------------------------------------
private:
    bool needsScale_; /// true if scene has to be scaled
public:
    void scale(qreal sx, qreal sy);
    QSize sizeHint() const;
protected:
    void scaleView();

//---------------------------------------------------------------------------------------------------------------
//                  create and handle graphicitems
//---------------------------------------------------------------------------------------------------------------
public:
    ProcessorGraphicsItem* getProcessorGraphicsItem(const Processor* processor) const;
    void selectPreviouslySelectedProcessors();

protected:
        //getter
    QList<PortOwnerGraphicsItem*> getSelectedPortOwnerGraphicsItems();
    QList<ProcessorGraphicsItem*> getSelectedProcessorGraphicsItems();
    PortGraphicsItem* getPortGraphicsItem(const Port* port) const;
    PropertyGraphicsItem* getPropertyGraphicsItem(const Property* prop) const;
        //create
    void generateGraphicsItems();
    ProcessorGraphicsItem* createProcessorGraphicsItem(Processor* processor);
    void createLinkArrowForPropertyLink(const PropertyLink* link);
        //delete
    void resetScene();
    void removeItems(QList<QGraphicsItem*> items);
    void removeItemsInDataFlow(QList<QGraphicsItem*> items);
    void removeItemsInGeneralLinking(QList<QGraphicsItem*> items);
    void removeItemsInCameraLinking(QList<QGraphicsItem*> items);
    void removeItemsInPortSizeLinking(QList<QGraphicsItem*> items);
    void removePortOwnerGraphicsItem(PortOwnerGraphicsItem* poItem);
    void removePortArrowGraphicsItem(PortArrowGraphicsItem* arrow);
    void removePropertyLinkArrowGraphicsItem(PropertyLinkArrowGraphicsItem* arrow);
    void removePortSizeLinkArrowGraphicsItem(PortSizeLinkArrowGraphicsItem* arrow);
    void removePortOwnerLinkArrowGraphicsItem(PortOwnerLinkArrowGraphicsItem* arrow);
    void removeTextBoxGraphicsItem(TextBoxBaseGraphicsItem* item);
private slots:
    void updateSelectedItems();

signals:
    void processorsSelected(const QList<Processor*>& processors);   ///< signal is emitted, if the selection has changed
    void networkEditor_visibleSceneSizeChanged_Signal();            ///< signal is emitted, if the visible scene size has changed

//---------------------------------------------------------------------------------------------------------------
//                  style, layer and cursor management
//---------------------------------------------------------------------------------------------------------------
private:
    NetworkEditorLayer currentLayer_;
    NetworkEditorCursorMode currentCursorMode_;
    NWEStyle_Base* currentStyle_;
    NWEGL_Base* currentGraphLayout_;
    bool currentToolTipMode_;

public:
    NetworkEditorLayer getCurrentLayer() const;
    NetworkEditorCursorMode getCurrentCursorMode() const;
    NWEStyle_Base* getCurrentStyle() const;
    bool getCurrentToolTipMode() const;
protected:
    void setLayer(NetworkEditorLayer layer);
    void setCursorMode(NetworkEditorCursorMode mode);
    void setStyle(NWEStyle_Base* style);
    void updateStyle();
    void updateGraphLayout();
    void processorFontOnChange();
    void setToolTipMode(bool mode);
public slots:
    void toggleToolTip();
//---------------------------------------------------------------------------------------------------------------
//                  button management
//---------------------------------------------------------------------------------------------------------------
private:
    QWidget* layerButtonContainer_;
        QToolButton* dataFlowLayerButton_;
        QToolButton* linkingLayerButton_;
    QWidget* linkingLayerButtonContainer_;
        QToolButton* generalLinkingButton_;
        QToolButton* cameraLinkingButton_;
        QToolButton* portSizeLinkingButton_;
    QWidget* generalLinkingLayerButtonContainer_;
        QToolButton* hideCameraLinksButton_;
        QToolButton* hidePortSizeLinksButton_;
        QToolButton* removeAllPropertyLinksButton_;
    QWidget* cameraLinkingLayerButtonContainer_;
        QToolButton* linkCamerasAutoButton_;
        QToolButton* linkCamerasButton_;
        QToolButton* removeAllCameraLinksButton_;
    QWidget* portSizeLinkingLayerButtonContainer_;
        QToolButton* linkPortSizeAutoButton_;
        QToolButton* linkPortSizeButton_;
        QToolButton* removeAllPortSizeLinksButton_;
    QWidget* stopButtonContainer_;
        QToolButton* stopNetworkEvaluatorButton_;
        bool networkEvaluatorIsLockedByButton_;
    QWidget* navigationButtonContainer_;
        QToolButton* selectCursorButton_;
        QToolButton* moveCursorButton_;
    QWidget* layoutButtonContainer_;
        QToolButton* centerViewButton_;
        QToolButton* graphLayoutButton_;

    QWidget* documentationButtonContainer_;
        QToolButton* documentationButton_;
public:
    bool cameraLinksHidden();
    bool portSizeLinksHidden();
protected:
    void initilizeEditorButtons();
    void layoutEditorButtons();
public slots:
    void setLayerToDataFlow();
    void setLayerToLinking();
    void setLayerToGeneralLinking();
    void setLayerToCameraLinking();
    void setLayerToPortSizeLinking();
private slots:
    void toggleNetworkEvaluator();
    void hideCameraLinks();
    void hidePortSizeLinks();
    void removeAllPropertyLinks();
    void linkCameras();
    void linkCamerasAutoChanged();
    void removeAllCameraLinks();
    void linkPortSize();
    void linkPortSizeAutoChanged();
    void removeAllPortSizeLinks();
    void sortEntireNetwork();
    void sortSubNetwork();
    void setViewCenter();
    void setCursorSelect();
    void setCursorMove();
    void linkCamerasOfProcessor(const Processor* processor);
    void toggleDocumentation();                     //toggles all documentations aka workspace description and text notes

//---------------------------------------------------------------------------------------------------------------
//                  documentation related functions
//---------------------------------------------------------------------------------------------------------------
private:
    //member
    QList<TextBoxBaseGraphicsItem*> documentationItems_; //< list containing all docu items

    QPoint nextTextBoxPos_; //< storage to determine, were a new text item should be created

    //serialize functions
protected:
    void serializeTextItems();   //is been called in SaveWorkspace(As)MenuEntity
private:
    void deserializeTextItems(); //< is been called in setWorkspace

//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
protected:
    //resize
    void resizeEvent(QResizeEvent* event);
    //wheel
    void wheelEvent(QWheelEvent* event);
    //mouse
    void mousePressEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
        //members
        bool translateScene_; /// true, if we move scene by mouse
        QPointF translateSceneVector_; // translation vector for scene movement
        QPointF lastTranslateCenter_;
        QPointF rightClickPosition_;
    //key
    void keyPressEvent(QKeyEvent* event);
    void dragEnterEvent(QDragEnterEvent* event);
    void dropEvent(QDropEvent* event);
    void dragMoveEvent(QDragMoveEvent* event);
    void dragLeaveEvent(QDragLeaveEvent* event);
        //members
        QGraphicsItem* selectedItem_;

//---------------------------------------------------------------------------------------------------------------
//                  action slots and context menu
//---------------------------------------------------------------------------------------------------------------
protected:
    void createContextMenuActions();
    void contextMenuEvent(QContextMenuEvent* event);

    QAction* copyAction_;
    QAction* pasteAction_;
    QAction* deleteAction_;
    QAction* editLinkAction_;
    QAction* deleteAllLinksAction_;
    QAction* deleteAllCameraLinksAction_;
    QAction* deleteAllPortSizeLinksAction_;
    QAction* createAllCameraLinksAction_;
    QAction* createAllPortSizeLinksAction_;
    QAction* deleteInnerLinksAction_;
    QAction* deleteInnerCameraLinksAction_;
    QAction* deleteInnerPortSizeLinksAction_;
    QAction* createInnerCameraLinksAction_;
    QAction* createInnerPortSizeLinksAction_;
    QAction* deletePortOwnerLinksAction_;
    QAction* deletePortOwnerCameraLinksAction_;
    QAction* deletePortOwnerPortSizeLinksAction_;
    QAction* createPortOwnerCameraLinksAction_;
    QAction* createPortOwnerPortSizeLinksAction_;
    QAction* sortSubNetworkAction_;
    QAction* createNewTextNoteAction_;
    QAction* createNewTextFrameAction_;

private slots:
    void copyActionSlot();
    void pasteActionSlot();
        bool clipboardHasValidContent();
    void deleteActionSlot();
    void editLinkActionSlot();
    void deleteLinksActionSlot();
    void deleteCameraLinksActionSlot();
    void deletePortSizeLinksActionSlot();
    void createCameraLinksActionSlot();
    void createPortSizeLinksActionSlot();

    void createNewTextNoteSlot();
    void createNewTextFrameSlot();
    void openPropertyLinkDialog(PortOwnerGraphicsItem* src, PortOwnerGraphicsItem* dest);
};

} // namespace voreen

#endif // VRN_NETWORKEDITOR_H
