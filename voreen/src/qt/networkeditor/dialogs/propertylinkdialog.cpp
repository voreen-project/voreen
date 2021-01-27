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

#include "voreen/qt/networkeditor/dialogs/propertylinkdialog.h"

#include "voreen/qt/networkeditor/networkeditor.h"
#include "voreen/qt/networkeditor/editor_settings.h"

#include "tgt/assert.h"

//core
#include "voreen/core/properties/property.h"
#include "voreen/core/properties/link/propertylink.h"
#include "voreen/core/properties/link/linkevaluatorbase.h"
#include "voreen/core/properties/link/linkevaluatorhelper.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/voreenapplication.h"

//gi
#include "voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/portownergraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/core/portgraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/connections/propertylinkarrowgraphicsitem.h"
#include "voreen/qt/networkeditor/dialogs/utils/splitgraphicsview.h"

#include "voreen/qt/widgets/enterexitpushbutton.h"

#include <QDialogButtonBox>
#include <QPushButton>
#include <QVBoxLayout>
#include <QGroupBox>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QButtonGroup>
#include <QLabel>
#include <QKeyEvent>
#include <QAbstractButton>
#include <QComboBox>
#include <QMessageBox>
#include <QSettings>

namespace voreen {

PropertyLinkDialog::PropertyLinkDialog(QWidget* parent, PortOwnerGraphicsItem* source, PortOwnerGraphicsItem* destination)
    : QDialog(parent)
    , sourceItem_(source), destinationItem_(destination)
    , sourceView_(0), destinationView_(0)
    , evaluatorSrcDstBox_(0), evaluatorDstSrcBox_(0)
    , arrowButtonGroup_(0), lastPressedDirectionButton_(0)
    , leftUpArrowButton_(0), bidirectionalArrowButton_(0), rightDownArrowButton_(0)
    , chainButton_(0), areEvaluatorsChained_(true)
{
    tgtAssert(source, "null pointer");
    tgtAssert(destination, "null pointer");

    //load dialog window settings
    loadSettings();

    //set dialog layout
    if (source != destination) {
        initializeLayout(false);
        sourceView_->setPropertyList(source->getPropertyList(),destination);
        destinationView_->setPropertyList(destination->getPropertyList(),source);
        connect(sourceView_,SIGNAL(linkProperties(PropertyGraphicsItem*, PropertyGraphicsItem*)),this,SLOT(linkProperties(PropertyGraphicsItem*, PropertyGraphicsItem*)));
        connect(destinationView_,SIGNAL(linkProperties(PropertyGraphicsItem*, PropertyGraphicsItem*)),this,SLOT(linkProperties(PropertyGraphicsItem*, PropertyGraphicsItem*)));
    }
    else {
        initializeLayout(true);
        sourceView_->setPropertyList(source->getPropertyList(),source);
        connect(sourceView_,SIGNAL(linkProperties(PropertyGraphicsItem*, PropertyGraphicsItem*)),this,SLOT(linkProperties(PropertyGraphicsItem*, PropertyGraphicsItem*)));
    }

    //clear selection
    sceneSelectionChanged();
    sourceView_->scene()->clearSelection();
    if(destinationView_) destinationView_->scene()->clearSelection();
}

PropertyLinkDialog::~PropertyLinkDialog() {
    //undo all temporal changes
    undoAllModifications();
    //move source PropertyList back to the NetworkEditor scene
    if(sourceItem_)
        sourceItem_->getPropertyList()->setScene(dynamic_cast<NetworkEditor*>(parent())->scene());
    //same for destination list
    if(sourceItem_ != destinationItem_)
        destinationItem_->getPropertyList()->setScene(dynamic_cast<NetworkEditor*>(parent())->scene());
}

//---------------------------------------------------------------------------------------------------------------
//                  layout
//---------------------------------------------------------------------------------------------------------------
void PropertyLinkDialog::initializeLayout(bool hasOnlyOneView) {

    setWindowTitle(tr("Property Link Editor"));

    QBoxLayout* layout = new QVBoxLayout(this);

    // GraphicsViews
    QWidget* viewContainer = new QWidget(this);
    QBoxLayout* viewLayout = new QHBoxLayout(viewContainer);
    viewLayout->setSpacing(0);
    if(hasOnlyOneView){
        sourceView_ = new SplitGraphicsView(this);
        sourceView_->setGlobalViewAlignment(SplitGraphicsView::GVA_CENTER);
        viewLayout->addWidget(sourceView_);
        connect(sourceView_->scene(), SIGNAL(selectionChanged()), this, SLOT(sceneSelectionChanged()));
    } else {
        sourceView_ = new SplitGraphicsView(this);
        sourceView_->setGlobalViewAlignment(SplitGraphicsView::GVA_LEFT);
        destinationView_ = new SplitGraphicsView(this);
        destinationView_->setGlobalViewAlignment(SplitGraphicsView::GVA_RIGHT);
        sourceView_->setBuddyView(destinationView_);
        destinationView_->setBuddyView(sourceView_);
        viewLayout->addWidget(sourceView_);
        viewLayout->addWidget(destinationView_);
        connect(sourceView_->scene(), SIGNAL(selectionChanged()), this, SLOT(sceneSelectionChanged()));
        connect(destinationView_->scene(), SIGNAL(selectionChanged()), this, SLOT(sceneSelectionChanged()));
    }
    layout->addWidget(viewContainer);


    //layout under the view
    QWidget* controlContainer = new QWidget(this);
    QBoxLayout* controlLayout = new QHBoxLayout(controlContainer);

    // left side of controls (all link specific)
    QWidget* selectedLinkContainer = new QWidget(this);
    QBoxLayout* selectedLinkLayout = new QVBoxLayout(selectedLinkContainer);

        //the evaluators
    QGroupBox* evaluatorContainer = new QGroupBox(tr("Selected Link Evaluator:"),this);
    QGridLayout* evaluatorLayout = new QGridLayout(evaluatorContainer);

    QLabel* evaluatorSrcDstLabel;
    if(hasOnlyOneView)
        evaluatorSrcDstLabel = new QLabel(QString("Evaluator ^ :"),this);
    else
        evaluatorSrcDstLabel = new QLabel(QString("Evaluator -> :"),this);
    evaluatorLayout->addWidget(evaluatorSrcDstLabel,0,0);

    evaluatorSrcDstBox_ = new QComboBox(this);
    evaluatorSrcDstBox_->setMinimumWidth(200);
    evaluatorSrcDstBox_->setEnabled(false);
    connect(evaluatorSrcDstBox_, SIGNAL(currentIndexChanged(int)), this, SLOT(srcDstEvaluatorBoxSelectionChanged(int)));
    evaluatorLayout->addWidget(evaluatorSrcDstBox_,0,1);

    QLabel* evaluatorDstSrcLabel;
    if(hasOnlyOneView)
        evaluatorDstSrcLabel = new QLabel(QString("Evaluator v :"),this);
    else
        evaluatorDstSrcLabel = new QLabel(QString("Evaluator <- :"),this);
    evaluatorLayout->addWidget(evaluatorDstSrcLabel,1,0);

    evaluatorDstSrcBox_ = new QComboBox(this);
    evaluatorDstSrcBox_->setMinimumWidth(200);
    evaluatorDstSrcBox_->setEnabled(false);
    connect(evaluatorDstSrcBox_, SIGNAL(currentIndexChanged(int)), this, SLOT(dstSrcEvaluatorBoxSelectionChanged(int)));
    evaluatorLayout->addWidget(evaluatorDstSrcBox_,1,1);

    chainButton_ = new QPushButton(QIcon(":/qt/icons/chain.png"), "");
    chainButton_->setIconSize(QSize(40,40));
    chainButton_->setEnabled(false);
    connect(chainButton_,SIGNAL(clicked()),this,SLOT(toggleChainButton()));
    chainButton_->setToolTip(tr("Evaluators are Connected"));
    evaluatorLayout->addWidget(chainButton_,0,2,2,1);

    selectedLinkLayout->addWidget(evaluatorContainer);

        //the direction buttons
    QGroupBox* directionButtonContainer = new QGroupBox(tr("Direction of Selected Link:"),this);
    QBoxLayout* directionButtonLayout = new QHBoxLayout(directionButtonContainer);

    if(hasOnlyOneView)
        leftUpArrowButton_ = new QPushButton(QIcon(":/qt/icons/arrow-up.png"), "");
    else
        leftUpArrowButton_ = new QPushButton(QIcon(":/qt/icons/arrow-left.png"), "");
    leftUpArrowButton_->setIconSize(QSize(40,40));
    leftUpArrowButton_->setCheckable(true);
    leftUpArrowButton_->setToolTip(tr("Direction: Right to Left"));
    directionButtonLayout->addWidget(leftUpArrowButton_);

    if(hasOnlyOneView)
        bidirectionalArrowButton_ = new QPushButton(QIcon(":/qt/icons/arrow-updown.png"), "");
    else
        bidirectionalArrowButton_ = new QPushButton(QIcon(":/qt/icons/arrow-leftright.png"), "");
    bidirectionalArrowButton_->setIconSize(QSize(40,40));
    bidirectionalArrowButton_->setCheckable(true);
    bidirectionalArrowButton_->setToolTip(tr("Direction: Bidirectional"));
    directionButtonLayout->addWidget(bidirectionalArrowButton_);

    if(hasOnlyOneView)
        rightDownArrowButton_ = new QPushButton(QIcon(":/qt/icons/arrow-down.png"), "");
    else
        rightDownArrowButton_ = new QPushButton(QIcon(":/qt/icons/arrow-right.png"), "");
    rightDownArrowButton_->setIconSize(QSize(40,40));
    rightDownArrowButton_->setCheckable(true);
    rightDownArrowButton_->setToolTip(tr("Direction: Left to Right"));
    directionButtonLayout->addWidget(rightDownArrowButton_);

    arrowButtonGroup_ = new QButtonGroup(controlContainer);
    arrowButtonGroup_->addButton(leftUpArrowButton_);
    arrowButtonGroup_->addButton(bidirectionalArrowButton_);
    arrowButtonGroup_->addButton(rightDownArrowButton_);
    connect(arrowButtonGroup_, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(directionButtonClicked(QAbstractButton*)));

    selectedLinkLayout->addWidget(directionButtonContainer);

    controlLayout->addWidget(selectedLinkContainer);

    controlLayout->addStretch(1);

    // right side of controls (target all links
    QWidget* allLinksContainer = new QWidget;
    QBoxLayout* allLinksLayout = new QVBoxLayout(allLinksContainer);

        //create all
    QGroupBox* createAllContainer = new QGroupBox(tr("Auto-Link by Name and Type:"),this);
    QBoxLayout* createAllLayout = new QHBoxLayout(createAllContainer);

    EnterExitPushButton* autolinkProcessor = new EnterExitPushButton(tr("Processor Links"), this);
    autolinkProcessor->setFocusPolicy(Qt::NoFocus);
    connect(autolinkProcessor, SIGNAL(clicked(bool)), this, SLOT(confirmAutoLinks()));
    connect(autolinkProcessor, SIGNAL(enterEventSignal()), this, SLOT(showProcessorAutoLinks()));
    connect(autolinkProcessor, SIGNAL(leaveEventSignal()), this, SLOT(hideAutoLinks()));
    createAllLayout->addWidget(autolinkProcessor);
    EnterExitPushButton* autolinkPort = new EnterExitPushButton(tr("Port Links"), this);
    autolinkPort->setFocusPolicy(Qt::NoFocus);
    connect(autolinkPort, SIGNAL(clicked(bool)), this, SLOT(confirmAutoLinks()));
    connect(autolinkPort, SIGNAL(enterEventSignal()), this, SLOT(showPortAutoLinks()));
    connect(autolinkPort, SIGNAL(leaveEventSignal()), this, SLOT(hideAutoLinks()));
    createAllLayout->addWidget(autolinkPort);

    allLinksLayout->addWidget(createAllContainer);

        //delete all
    QGroupBox* deleteAllContainer = new QGroupBox(tr("Delete all Links:"),this);
    QBoxLayout* deleteAllLayout = new QHBoxLayout(deleteAllContainer);

    QPushButton* deleteProcessorLinks = new QPushButton(tr("Processor Links"));
    deleteProcessorLinks->setFocusPolicy(Qt::NoFocus);
    connect(deleteProcessorLinks, SIGNAL(clicked(bool)), this, SLOT(deleteAllProcessorLinks()));
    deleteAllLayout->addWidget(deleteProcessorLinks);
    QPushButton* deletePortLinks = new QPushButton(tr("Port Links"));
    deletePortLinks->setFocusPolicy(Qt::NoFocus);
    connect(deletePortLinks, SIGNAL(clicked(bool)), this, SLOT(deleteAllPortLinks()));
    deleteAllLayout->addWidget(deletePortLinks);
    allLinksLayout->addWidget(deleteAllContainer);

    controlLayout->addWidget(allLinksContainer);

    layout->addWidget(controlContainer);

    // Buttons
    QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(buttonBox, SIGNAL(accepted()), this, SLOT(applyAllModifications()));
    connect(buttonBox, SIGNAL(rejected()), this, SLOT(close()));

    layout->addWidget(buttonBox);
}

//---------------------------------------------------------------------------------------------------------------
//                  view interaction
//---------------------------------------------------------------------------------------------------------------
void PropertyLinkDialog::sceneSelectionChanged() {
    //get currently selected arrow
    PropertyLinkArrowGraphicsItem* arrow = getCurrentlySelectedArrow();
    //if no arrow is selected:
    if (arrow == 0) {
        //block all signals
        arrowButtonGroup_->setExclusive(false);
        arrowButtonGroup_->blockSignals(true);
        evaluatorSrcDstBox_->blockSignals(true);
        evaluatorDstSrcBox_->blockSignals(true);
        evaluatorSrcDstBox_->setEnabled(false);
        evaluatorDstSrcBox_->setEnabled(false);
        chainButton_->blockSignals(true);
        chainButton_->setEnabled(false);
        //reset evaluator combo boxes
        evaluatorSrcDstBox_->clear();
        evaluatorDstSrcBox_->clear();
        //reset arrow buttons
        leftUpArrowButton_->setChecked(false);
        leftUpArrowButton_->setEnabled(false);
        bidirectionalArrowButton_->setChecked(true);    //bidirectinal default
        bidirectionalArrowButton_->setEnabled(false);
        rightDownArrowButton_->setChecked(false);
        rightDownArrowButton_->setEnabled(false);
        //set evaluators chained
        toggleChainButton(true);
        //enable all signals
        evaluatorSrcDstBox_->blockSignals(false);
        evaluatorDstSrcBox_->blockSignals(false);
        arrowButtonGroup_->blockSignals(false);
        arrowButtonGroup_->setExclusive(true);
        chainButton_->blockSignals(false);
        return;
    }
    else { // id we have an selected arrow update evaluator and arrow elements
        updateEvaluatorBoxes(arrow);
        updateArrowButtonGroup(arrow);
    }
}

PropertyLinkArrowGraphicsItem* PropertyLinkDialog::getCurrentlySelectedArrow() const {
    QList<QGraphicsItem*> items = sourceView_->scene()->selectedItems();
    //return null, if nothing is selected
    if (items.count() == 0)
        return 0;
    //if only one arrow is selected return it
    if(items.count() == 1 && items[0]->type() == UserTypesPropertyLinkArrowGraphicsItem)
        return qgraphicsitem_cast<PropertyLinkArrowGraphicsItem*>(items[0]);
    //all other cases return null
    return 0;
}

void PropertyLinkDialog::setCurrentlySelectedArrow(PropertyGraphicsItem* src, PropertyGraphicsItem* dst) {
    //src item has to be in sourceView_
    tgtAssert(sourceView_->getPropertyList()->getAllPropertyItems().contains(src),"src is not from sourceView");
    //clear selection
    sourceView_->scene()->blockSignals(true);
    sourceView_->scene()->clearSelection();
    sourceView_->scene()->blockSignals(false);
    //find arrow
    foreach(PropertyLinkArrowGraphicsItem* item, src->getDestinationLinkList())
        if(item->getSourceItem() == dst){
            item->setSelected(true);
            return;
    }
    foreach(PropertyLinkArrowGraphicsItem* item, src->getSourceLinkList())
        if(item->getDestinationItem() == dst){
            item->setSelected(true);
            return;
        }
    //should not get here
    tgtAssert(false,"No link between src and dst found!!!");
    sourceView_->scene()->clearSelection();
}

//---------------------------------------------------------------------------------------------------------------
//                  evaluator functions
//---------------------------------------------------------------------------------------------------------------
void PropertyLinkDialog::toggleChainButton() {
    toggleChainButton(!areEvaluatorsChained_);
}

void PropertyLinkDialog::toggleChainButton(bool b) {
    if(areEvaluatorsChained_ == b) return;
    areEvaluatorsChained_ = b;
    if(areEvaluatorsChained_) {
        chainButton_->setIcon(QIcon(":/qt/icons/chain.png"));
        chainButton_->setToolTip(tr("Evaluators are Connected"));
    } else {
        chainButton_->setIcon(QIcon(":/qt/icons/unchain.png"));
        chainButton_->setToolTip(tr("Evaluators are not Connected"));
    }
}

LinkEvaluatorBase* PropertyLinkDialog::getEvaluator(PropertyGraphicsItem* src, PropertyGraphicsItem* dst) {
    foreach(Modification* m, addModificationList_) {
        if(m->src_ == src && m->dst_ == dst)
            return m->eval_;
    }
    foreach(Modification* m, evaluatorModificationList_) {
        if(m->src_ == src && m->dst_ == dst)
            return m->eval_;
    }
    PropertyLink* link = src->getProperty()->getLink(dst->getProperty());
    if(link)
        return link->getLinkEvaluator();

    return 0;
}

void PropertyLinkDialog::updateEvaluatorBoxes(PropertyLinkArrowGraphicsItem* arrow) {
    tgtAssert(arrow,"PropertyLinkDialog::updateCombobox: null pointer passed");
    //pre settings
    evaluatorSrcDstBox_->blockSignals(true);
    evaluatorDstSrcBox_->blockSignals(true);
    evaluatorSrcDstBox_->setEnabled(true);
    evaluatorDstSrcBox_->setEnabled(true);
    evaluatorSrcDstBox_->clear();
    evaluatorDstSrcBox_->clear();

    //set left/up property of selected arrow as src
    PropertyGraphicsItem* src = 0;
    //set right/down property of selected arrow as dst
    PropertyGraphicsItem* dst = 0;

    //set src and dst
    bool isLeftUpArrow = false;
    //check if single or dual view
    if(sourceItem_ == destinationItem_) { //case single view
       if(arrow->getSourceItem()->pos().y() < arrow->getDestinationItem()->pos().y()) { //link from up to down
            src = arrow->getSourceItem();
            dst = arrow->getDestinationItem();
        } else { // link from down to up
            src = arrow->getDestinationItem();
            dst = arrow->getSourceItem();
            isLeftUpArrow = true;
        }
    } else { // case dual view
        if(sourceItem_->getPropertyList()->getAllPropertyItems().contains(arrow->getSourceItem())) { // link from left to right
            src = arrow->getSourceItem();
            dst = arrow->getDestinationItem();
        } else { //link from right ro left
            src = arrow->getDestinationItem();
            dst = arrow->getSourceItem();
            isLeftUpArrow = true;
        }
    }
    //update evaluators
    if(arrow->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {//bi
        //differ single/dual view
        if(sourceItem_ == destinationItem_) { //case single view
            fillEvaluatorBox(evaluatorSrcDstBox_,dst,src);
            fillEvaluatorBox(evaluatorDstSrcBox_,src,dst);
        } else { //case dual view
            fillEvaluatorBox(evaluatorSrcDstBox_,src,dst);
            fillEvaluatorBox(evaluatorDstSrcBox_,dst,src);
        }
        //set chain button
        chainButton_->setEnabled(true);
        if(getEvaluator(src,dst)->getClassName() == getEvaluator(dst,src)->getClassName()) {
            toggleChainButton(true);
            chainButton_->setEnabled(true);
        } else {
            toggleChainButton(false);
            chainButton_->setEnabled(false);
        }
    } else { //uni
        //differ single/dual view
        if(sourceItem_ == destinationItem_) { //case single view
            if(isLeftUpArrow) {
                fillEvaluatorBox(evaluatorSrcDstBox_,dst,src);
                fillEvaluatorBox(evaluatorDstSrcBox_,0,0);
            } else {
                fillEvaluatorBox(evaluatorSrcDstBox_,0,0);
                fillEvaluatorBox(evaluatorDstSrcBox_,src,dst);
            }
        } else { // case dual view
            if(isLeftUpArrow) {
                fillEvaluatorBox(evaluatorSrcDstBox_,0,0);
                fillEvaluatorBox(evaluatorDstSrcBox_,dst,src);
            } else {
                fillEvaluatorBox(evaluatorSrcDstBox_,src,dst);
                fillEvaluatorBox(evaluatorDstSrcBox_,0,0);
            }
        }
        //set chain button
        chainButton_->setEnabled(false);
        toggleChainButton(false);
    }
    //post settings
    evaluatorSrcDstBox_->blockSignals(false);
    evaluatorDstSrcBox_->blockSignals(false);
}

void PropertyLinkDialog::fillEvaluatorBox(QComboBox* box, PropertyGraphicsItem* src,PropertyGraphicsItem* dst) {
    //get other box
    QComboBox* otherBox = evaluatorDstSrcBox_;
    if(box == evaluatorDstSrcBox_)
        otherBox = evaluatorSrcDstBox_;

    //pre settings
    box->blockSignals(true);
    chainButton_->blockSignals(true);
    box->clear();
    box->setEnabled(true);

    //is empty?
    if(!src || !dst){
        box->setEnabled(false);
        toggleChainButton(false);
        chainButton_->setEnabled(false);
        chainButton_->blockSignals(false);
        box->blockSignals(false);
        return;
    }

    //get available evaluators
    std::vector<std::pair<std::string, std::string> > availableEvaluators =
            src->getProperty()->getCompatibleEvaluators(dst->getProperty());
    //is list empty?
    if(availableEvaluators.empty()){
        box->setEnabled(false);
        toggleChainButton(false);
        chainButton_->setEnabled(false);
        chainButton_->blockSignals(false);
        box->blockSignals(false);
        return;
    }

    //fill combo box
    for (std::vector<std::pair<std::string, std::string> >::iterator i = availableEvaluators.begin();
            i != availableEvaluators.end(); ++i) {
                box->addItem(QString::fromStdString(i->second),QVariant(QString::fromStdString(i->first)));
    }

    //set selected evaluator
    LinkEvaluatorBase* boxEval = 0;     //evaluator of 'box' direction
    LinkEvaluatorBase* oppEval = 0;   //evaluator of other direction
    if(!(boxEval = getEvaluator(src,dst))){ //case, link has no actual evaluator (new link)
        bool hasOppEvaluator = false;
        if((oppEval = getEvaluator(dst,src))) { //if link in other direction exists and is suitable, take that one
            for (std::vector<std::pair<std::string, std::string> >::iterator i = availableEvaluators.begin(); i != availableEvaluators.end(); ++i) {
                if(i->first == oppEval->getClassName()) {
                    box->setCurrentIndex(box->findText(QString::fromStdString(i->second)));
                    addEvaluatorModification(new Modification(src,dst, LinkEvaluatorHelper::createEvaluator(i->first)));
                    hasOppEvaluator = true;
                    break;
                }
            }
        }
        if(!hasOppEvaluator) { //link in other direction does not exist or is not suitable, take first found evaluator
            box->setCurrentIndex(box->findText(QString::fromStdString(availableEvaluators[0].second)));
            addEvaluatorModification(new Modification(src,dst, LinkEvaluatorHelper::createEvaluator(availableEvaluators[0].first)));
        }
    } else { //evaluator have been found
        box->setCurrentIndex(box->findText(QString::fromStdString(
            boxEval->getGuiName())));
    }

    //are evaluators chained?
    if(otherBox->currentIndex() == -1) {
        toggleChainButton(true);
        chainButton_->setEnabled(false);
    } else
    if(box->itemData(box->currentIndex()).value<QString>() ==
       otherBox->itemData(evaluatorDstSrcBox_->currentIndex()).value<QString>()) {
        toggleChainButton(true);
        chainButton_->setEnabled(true);
    } else {
        toggleChainButton(false);
        chainButton_->setEnabled(true);
    }
    chainButton_->blockSignals(false);
    box->blockSignals(false);
}

void PropertyLinkDialog::evaluatorBoxSelectionChanged(QComboBox* box, int index) {
    //return if empty or no arrow is selected
    if(index == -1) return;
    PropertyLinkArrowGraphicsItem* arrow = getCurrentlySelectedArrow();
    if(!arrow) return;

    //get other box
    QComboBox* otherBox = evaluatorDstSrcBox_;
    if(box == evaluatorDstSrcBox_)
        otherBox = evaluatorSrcDstBox_;

    //set source property of selected arrow according to the combobox
    PropertyGraphicsItem* src = 0;
    //set destination property of selected arrow according to the combobox
    PropertyGraphicsItem* dst = 0;

    //set src and dst
    //check if single or dual view
    if(sourceItem_ == destinationItem_) { //case single view
       if(arrow->getSourceItem()->pos().y() < arrow->getDestinationItem()->pos().y()) { //link from up to down
            src = arrow->getSourceItem();
            dst = arrow->getDestinationItem();
        } else { // link from down to up
            src = arrow->getDestinationItem();
            dst = arrow->getSourceItem();
        }
        //switch src and dst if
        if(box == evaluatorSrcDstBox_) {
            PropertyGraphicsItem* tmpItem = src;
            src = dst;
            dst = tmpItem;
        }
    } else { // case dual view
        if(sourceItem_->getPropertyList()->getAllPropertyItems().contains(arrow->getSourceItem())) { // link from left to right
            src = arrow->getSourceItem();
            dst = arrow->getDestinationItem();
        } else { //link from right ro left
            src = arrow->getDestinationItem();
            dst = arrow->getSourceItem();
        }
        //switch src and dst if
        if(box == evaluatorDstSrcBox_) {
            PropertyGraphicsItem* tmpItem = src;
            src = dst;
            dst = tmpItem;
        }
    }

    //better save than sorry
    if (arrow) {
        //generate link and modify evaluator
        LinkEvaluatorBase* leb = LinkEvaluatorHelper::createEvaluator(box->itemData(index).value<QString>().toStdString());
        addEvaluatorModification(new Modification(src,dst,leb));

        //set other evaluator, if evaluators are chained
        if(arrow->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone && areEvaluatorsChained_) {
            int otherIndex = otherBox->findData(box->itemData(index));
            if(otherIndex != -1) { //evaluator in other direction exists
                LinkEvaluatorBase* leb2 = LinkEvaluatorHelper::createEvaluator(otherBox->itemData(otherIndex).value<QString>().toStdString());
                addEvaluatorModification(new Modification(dst,src,leb2));
                otherBox->blockSignals(true);
                otherBox->setCurrentIndex(otherIndex);
                otherBox->blockSignals(false);
            } else { //evaluator in other direction does not exist
                toggleChainButton(false);
                chainButton_->setEnabled(false);
            }
        }  else { //test if evaluators are equal and set chain button
            if(getEvaluator(dst,src) && getEvaluator(src,dst)->getClassName() == getEvaluator(dst,src)->getClassName()) {
                toggleChainButton(true);
                chainButton_->setEnabled(true);
            } else {
                toggleChainButton(false);
                chainButton_->setEnabled(false);
            }
        }
        //update arrow buttons
        updateArrowButtonGroup(arrow);
    }
}

void PropertyLinkDialog::srcDstEvaluatorBoxSelectionChanged(int index) {
    evaluatorBoxSelectionChanged(evaluatorSrcDstBox_,index);
}

void PropertyLinkDialog::dstSrcEvaluatorBoxSelectionChanged(int index) {
    evaluatorBoxSelectionChanged(evaluatorDstSrcBox_,index);
}

//---------------------------------------------------------------------------------------------------------------
//                  direction button functions
//---------------------------------------------------------------------------------------------------------------
void PropertyLinkDialog::updateArrowButtonGroup(PropertyLinkArrowGraphicsItem* arrow) {
    tgtAssert(arrow,"PropertyLinkDialog::updateArrowButtonGroup: null pointer passed");
    //enable all buttons
    arrowButtonGroup_->blockSignals(true);
    leftUpArrowButton_->setEnabled(true);
    rightDownArrowButton_->setEnabled(true);
    bidirectionalArrowButton_->setEnabled(true);
    //if the sourde item has an arrow head, the link must be bidirectional
    if(arrow->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
        leftUpArrowButton_->setEnabled(true);
        rightDownArrowButton_->setEnabled(true);
        bidirectionalArrowButton_->setEnabled(true);
        bidirectionalArrowButton_->setChecked(true);
    } else { //differ between single view and dual view
        if(destinationItem_ == sourceItem_) {//case single view
            if(arrow->getSourceItem()->pos().y() < arrow->getDestinationItem()->pos().y()) { //case link from up to down
                rightDownArrowButton_->setEnabled(true);
                rightDownArrowButton_->setChecked(true);
                if(arrow->getDestinationItem()->getProperty()->isLinkableWith(arrow->getSourceItem()->getProperty())) {
                    leftUpArrowButton_->setEnabled(true);
                    bidirectionalArrowButton_->setEnabled(true);
                } else {
                    leftUpArrowButton_->setEnabled(false);
                    bidirectionalArrowButton_->setEnabled(false);
                }
            }
            else { //case link from down to up
                leftUpArrowButton_->setEnabled(true);
                leftUpArrowButton_->setChecked(true);
                if(arrow->getDestinationItem()->getProperty()->isLinkableWith(arrow->getSourceItem()->getProperty())) {
                    rightDownArrowButton_->setEnabled(true);
                    bidirectionalArrowButton_->setEnabled(true);
                } else {
                    rightDownArrowButton_->setEnabled(false);
                    bidirectionalArrowButton_->setEnabled(false);
                }
            }
        } else { //case dual view
            if(sourceItem_->getPropertyList()->getAllPropertyItems().contains(arrow->getSourceItem())) { //case link from left to right
                rightDownArrowButton_->setEnabled(true);
                rightDownArrowButton_->setChecked(true);
                if(arrow->getDestinationItem()->getProperty()->isLinkableWith(arrow->getSourceItem()->getProperty())) {
                    leftUpArrowButton_->setEnabled(true);
                    bidirectionalArrowButton_->setEnabled(true);
                } else {
                    leftUpArrowButton_->setEnabled(false);
                    bidirectionalArrowButton_->setEnabled(false);
                }
            }
            else { // case link from right to left
                leftUpArrowButton_->setEnabled(true);
                leftUpArrowButton_->setChecked(true);
                if(arrow->getDestinationItem()->getProperty()->isLinkableWith(arrow->getSourceItem()->getProperty())) {
                    rightDownArrowButton_->setEnabled(true);
                    bidirectionalArrowButton_->setEnabled(true);
                } else {
                    rightDownArrowButton_->setEnabled(false);
                    bidirectionalArrowButton_->setEnabled(false);
                }
            }
        }
    }
    //set last pressed button to stop reupdating if the same button is clicked in a row
    lastPressedDirectionButton_ = arrowButtonGroup_->checkedButton();
    //enable signal
    arrowButtonGroup_->blockSignals(false);
}

void PropertyLinkDialog::directionButtonClicked(QAbstractButton* button) {
    //return, if the selected button was clicked
    if(lastPressedDirectionButton_ == button) return;
    //set new last pressed button
    lastPressedDirectionButton_ = button;
    //get selected arrow.
    PropertyLinkArrowGraphicsItem* arrow = getCurrentlySelectedArrow();
    tgtAssert(arrow,"Arrow button was enabled and clicked while no arrow has been selected!!!");

    //set left/up property of selected arrow as src
    PropertyGraphicsItem* src = 0;
    //set right/down property of selected arrow as dst
    PropertyGraphicsItem* dst = 0;

    //better save than sorry
    if (arrow) {
        //set src and dst
        bool isLeftUpArrow = false;
        arrow->setZValue(ZValuesPropertyLinkArrowGraphicsItem); //prevents flickering
        arrow->getCorrespondingPropertyLinkArrowItem()->setZValue(ZValuesPropertyLinkArrowGraphicsItem);
        //check if single or dual view
        if(sourceItem_ == destinationItem_) { //case single view
            if(arrow->getSourceItem()->pos().y() < arrow->getDestinationItem()->pos().y()) { //link from up to down
                src = arrow->getSourceItem();
                dst = arrow->getDestinationItem();
            } else { // link from down to up
                src = arrow->getDestinationItem();
                dst = arrow->getSourceItem();
                isLeftUpArrow = true;
            }
        } else { // case dual view
            if(sourceItem_->getPropertyList()->getAllPropertyItems().contains(arrow->getSourceItem())) { // link from left to right
                src = arrow->getSourceItem();
                dst = arrow->getDestinationItem();
            } else { //link from right ro left
                src = arrow->getDestinationItem();
                dst = arrow->getSourceItem();
                isLeftUpArrow = true;
            }
        }
        //modify arrow depending on button had been clicked
        if (button == leftUpArrowButton_){
            if(arrow->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
                //from bi to left/up
                addRemoveModification(new Modification(src,dst,0));
                src->removeGraphicalLink(dst);
            } else {
                //from right/down to left/up
                fillEvaluatorBox(evaluatorDstSrcBox_,dst,src);
                addAddModification(new Modification(dst,src, LinkEvaluatorHelper::createEvaluator
                                  (getEvaluator(dst,src)->getClassName())));
                dst->addGraphicalLink(src);
                addRemoveModification(new Modification(src,dst,0));
                src->removeGraphicalLink(dst);
                src->setVisible(true); //to make arrow visible
            }
        } else if(button == rightDownArrowButton_) {
            if(arrow->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone) {
                //from bi to right/down
                addRemoveModification(new Modification(dst,src,0));
                dst->removeGraphicalLink(src);
            } else {
                //from left/up to right/down
                fillEvaluatorBox(evaluatorSrcDstBox_,src,dst);
                addAddModification(new Modification(src,dst, LinkEvaluatorHelper::createEvaluator
                                    (getEvaluator(src,dst)->getClassName())));
                src->addGraphicalLink(dst);
                addRemoveModification(new Modification(dst,src,0));
                dst->removeGraphicalLink(src);
                dst->setVisible(true); //to make arrow visible
            }
        } else if(button == bidirectionalArrowButton_) {
                if(isLeftUpArrow) { //left/up to bi
                    fillEvaluatorBox(evaluatorSrcDstBox_,src,dst);
                    addAddModification(new Modification(src,dst, LinkEvaluatorHelper::createEvaluator
                                        (getEvaluator(src,dst)->getClassName())));
                    src->addGraphicalLink(dst);
                    dst->setVisible(true);
                    src->setVisible(true);
                } else { //right/down to bi
                    fillEvaluatorBox(evaluatorDstSrcBox_,dst,src);
                    addAddModification(new Modification(dst,src, LinkEvaluatorHelper::createEvaluator
                                        (getEvaluator(dst,src)->getClassName())));
                    dst->addGraphicalLink(src);
                    dst->setVisible(true);
                    src->setVisible(true);
                }
        } else {
                tgtAssert(false,"unknown direction");
        }
        setCurrentlySelectedArrow(src,dst);
        updateEvaluatorBoxes(getCurrentlySelectedArrow());
        getCurrentlySelectedArrow()->setZValue(ZValuesSelectedPropertyLinkArrowGraphicsItem); //prevents flickering
    }
    //update views
    sourceView_->updateViewport();
}

//---------------------------------------------------------------------------------------------------------------
//                  delete functions
//---------------------------------------------------------------------------------------------------------------
void PropertyLinkDialog::deleteArrowAction() {
    PropertyLinkArrowGraphicsItem* arrow = getCurrentlySelectedArrow();

    if (arrow) {
        if(arrow->getSourceArrowHeadDirection() != ConnectionBaseGraphicsItem::ArrowHeadDirectionNone){
            arrow->getDestinationItem()->removeGraphicalLink(arrow->getSourceItem());
            addRemoveModification(new Modification(arrow->getDestinationItem(),arrow->getSourceItem(),0));
        }
        addRemoveModification(new Modification(arrow->getSourceItem(),arrow->getDestinationItem(),0));
        arrow->getSourceItem()->removeGraphicalLink(arrow->getDestinationItem());
        sourceView_->viewport()->update();
        if(destinationView_) destinationView_->viewport()->update();
    }
}

void PropertyLinkDialog::deleteAllProcessorLinks() {
    foreach (PropertyGraphicsItem* sProp, sourceItem_->getPropertyList()->getPropertyItems()){
        if(sProp->isVisible()){
            foreach (PropertyLinkArrowGraphicsItem* arrow, sProp->getSourceLinkList()){
                if(arrow->getDestinationItem()->isVisible() && arrow->getDestinationItem()->getPropertyOwnerItem() != sourceItem_){
                    addRemoveModification(new Modification(arrow->getSourceItem(), arrow->getDestinationItem(), 0));
                    arrow->getSourceItem()->removeGraphicalLink(arrow->getDestinationItem());
                }
            }
            foreach (PropertyLinkArrowGraphicsItem* arrow, sProp->getDestinationLinkList()){
                if(arrow->getDestinationItem()->isVisible() && arrow->getDestinationItem()->getPropertyOwnerItem() != sourceItem_){
                    addRemoveModification(new Modification(arrow->getSourceItem(), arrow->getDestinationItem(), 0));
                    arrow->getSourceItem()->removeGraphicalLink(arrow->getDestinationItem());
                }
            }
        }
    }
    foreach (PropertyGraphicsItem* dProp, destinationItem_->getPropertyList()->getPropertyItems()){
        if(dProp->isVisible()){
            foreach (PropertyLinkArrowGraphicsItem* arrow, dProp->getSourceLinkList()){
                if(arrow->getDestinationItem()->isVisible() && arrow->getDestinationItem()->getPropertyOwnerItem() != destinationItem_){
                    addRemoveModification(new Modification(arrow->getSourceItem(), arrow->getDestinationItem(), 0));
                    arrow->getSourceItem()->removeGraphicalLink(arrow->getDestinationItem());
                }
            }
            foreach (PropertyLinkArrowGraphicsItem* arrow, dProp->getDestinationLinkList()){
                if(arrow->getDestinationItem()->isVisible() && arrow->getDestinationItem()->getPropertyOwnerItem() != destinationItem_){
                    addRemoveModification(new Modification(arrow->getSourceItem(), arrow->getDestinationItem(), 0));
                    arrow->getSourceItem()->removeGraphicalLink(arrow->getDestinationItem());
                }
            }
        }
    }
    sourceView_->updateViewport();
}

void PropertyLinkDialog::deleteAllPortLinks() {
    foreach(PropertyListGraphicsItem* list, sourceItem_->getPropertyList()->getSubPropertyListItems()){
        if(dynamic_cast<PortGraphicsItem*>(list->getParentPropertyOwnerItem())) {
            foreach (PropertyGraphicsItem* sProp, list->getPropertyItems()){
                if(sProp->isVisible()){
                    foreach (PropertyLinkArrowGraphicsItem* arrow, sProp->getSourceLinkList()){
                        if(arrow->getDestinationItem()->isVisible()
                            && !sourceItem_->getPropertyList()->getAllPropertyItems().contains(arrow->getDestinationItem())
                            && getEvaluator(arrow->getSourceItem(),arrow->getDestinationItem())->getClassName() != "LinkEvaluatorRenderSize"){
                            addRemoveModification(new Modification(arrow->getSourceItem(), arrow->getDestinationItem(), 0));
                            arrow->getSourceItem()->removeGraphicalLink(arrow->getDestinationItem());
                        }
                    }
                    foreach (PropertyLinkArrowGraphicsItem* arrow, sProp->getDestinationLinkList()){
                        if(arrow->getDestinationItem()->isVisible()
                            && !sourceItem_->getPropertyList()->getAllPropertyItems().contains(arrow->getDestinationItem())
                            && getEvaluator(arrow->getSourceItem(),arrow->getDestinationItem())->getClassName() != "LinkEvaluatorRenderSize"){
                            addRemoveModification(new Modification(arrow->getSourceItem(), arrow->getDestinationItem(), 0));
                            arrow->getSourceItem()->removeGraphicalLink(arrow->getDestinationItem());
                        }
                    }
                }
            }
        }
    }
    foreach(PropertyListGraphicsItem* list, destinationItem_->getPropertyList()->getSubPropertyListItems()){
        if(dynamic_cast<PortGraphicsItem*>(list->getParentPropertyOwnerItem())) {
            foreach (PropertyGraphicsItem* dProp, list->getPropertyItems()){
                if(dProp->isVisible()){
                    foreach (PropertyLinkArrowGraphicsItem* arrow, dProp->getSourceLinkList()){
                        if(arrow->getDestinationItem()->isVisible()
                            && !destinationItem_->getPropertyList()->getAllPropertyItems().contains(arrow->getDestinationItem())
                            && getEvaluator(arrow->getSourceItem(),arrow->getDestinationItem())->getClassName() != "LinkEvaluatorRenderSize"){
                            addRemoveModification(new Modification(arrow->getSourceItem(), arrow->getDestinationItem(), 0));
                            arrow->getSourceItem()->removeGraphicalLink(arrow->getDestinationItem());
                        }
                    }
                    foreach (PropertyLinkArrowGraphicsItem* arrow, dProp->getDestinationLinkList()){
                        if(arrow->getDestinationItem()->isVisible()
                            && !destinationItem_->getPropertyList()->getAllPropertyItems().contains(arrow->getDestinationItem())
                            && getEvaluator(arrow->getSourceItem(),arrow->getDestinationItem())->getClassName() != "LinkEvaluatorRenderSize"){
                            addRemoveModification(new Modification(arrow->getSourceItem(), arrow->getDestinationItem(), 0));
                            arrow->getSourceItem()->removeGraphicalLink(arrow->getDestinationItem());
                        }
                    }
                }
            }
        }
    }
    sourceView_->updateViewport();
}

//---------------------------------------------------------------------------------------------------------------
//                  auto linking functions
//---------------------------------------------------------------------------------------------------------------
void PropertyLinkDialog::showProcessorAutoLinks() {
    autoArrows_.clear();
    foreach (PropertyGraphicsItem* leftProp, sourceItem_->getPropertyList()->getPropertyItems()) {
        foreach (PropertyGraphicsItem* rightProp, destinationItem_->getPropertyList()->getPropertyItems()) {
            bool namesCompatible =
                   leftProp->getProperty()->getGuiName().find(rightProp->getProperty()->getGuiName()) != std::string::npos
                || rightProp->getProperty()->getGuiName().find(leftProp->getProperty()->getGuiName()) != std::string::npos;

            if (namesCompatible
                    && (typeid(*leftProp->getProperty()) == typeid(*rightProp->getProperty()))
                    && (leftProp != rightProp)
                    && !(leftProp->getProperty()->isLinkedWith(rightProp->getProperty(),true))) {

                if(leftProp->getProperty()->isLinkableWith(rightProp->getProperty()) &&
                   rightProp->getProperty()->isLinkableWith(leftProp->getProperty())) {
                    PropertyLinkArrowGraphicsItem* arrow = new PropertyLinkArrowGraphicsItem(leftProp,rightProp);
                    arrow->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
                    arrow->setVisible(true); autoArrows_.append(arrow);
                    arrow = new PropertyLinkArrowGraphicsItem(rightProp,leftProp);
                    arrow->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
                    arrow->setVisible(true); autoArrows_.append(arrow);
                }
                else if(leftProp->getProperty()->isLinkableWith(rightProp->getProperty())) {
                    PropertyLinkArrowGraphicsItem* arrow = new PropertyLinkArrowGraphicsItem(leftProp,rightProp);
                    arrow->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
                    arrow->setVisible(true); autoArrows_.append(arrow);
                }
                else if(rightProp->getProperty()->isLinkableWith(leftProp->getProperty())) {
                    PropertyLinkArrowGraphicsItem* arrow = new PropertyLinkArrowGraphicsItem(rightProp,leftProp);
                    arrow->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
                    arrow->setVisible(true); autoArrows_.append(arrow);
                }
            }
        }
    }
    sourceView_->updateViewport();
}

void PropertyLinkDialog::showPortAutoLinks() {
    autoArrows_.clear();
    foreach (PropertyListGraphicsItem* leftList, sourceItem_->getPropertyList()->getSubPropertyListItems()) {
        if(dynamic_cast<PortGraphicsItem*>(leftList->getParentPropertyOwnerItem())) {
            foreach (PropertyListGraphicsItem* rightList, destinationItem_->getPropertyList()->getSubPropertyListItems()) {
                if(dynamic_cast<PortGraphicsItem*>(leftList->getParentPropertyOwnerItem())) {
                    foreach (PropertyGraphicsItem* leftProp, leftList->getPropertyItems()) {
                        foreach (PropertyGraphicsItem* rightProp, rightList->getPropertyItems()) {

                            if (leftProp->getProperty()->getGuiName() == rightProp->getProperty()->getGuiName()
                                    && (typeid(*leftProp->getProperty()) == typeid(*rightProp->getProperty()))
                                    && (leftProp != rightProp)
                                    && !(leftProp->getProperty()->isLinkedWith(rightProp->getProperty(),true))
                                    && !(dynamic_cast<RenderSizeOriginProperty*>(const_cast<Property*>(leftProp->getProperty())))
                                    && !(dynamic_cast<RenderSizeReceiveProperty*>(const_cast<Property*>(leftProp->getProperty())))
                                    && !(dynamic_cast<RenderSizeOriginProperty*>(const_cast<Property*>(rightProp->getProperty())))
                                    && !(dynamic_cast<RenderSizeReceiveProperty*>(const_cast<Property*>(rightProp->getProperty())))) {

                                if(leftProp->getProperty()->isLinkableWith(rightProp->getProperty()) &&
                                   rightProp->getProperty()->isLinkableWith(leftProp->getProperty())) {
                                    PropertyLinkArrowGraphicsItem* arrow = new PropertyLinkArrowGraphicsItem(leftProp,rightProp);
                                    arrow->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
                                    arrow->setVisible(true); autoArrows_.append(arrow);
                                    arrow = new PropertyLinkArrowGraphicsItem(rightProp,leftProp);
                                    arrow->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
                                    arrow->setVisible(true); autoArrows_.append(arrow);
                                }
                                else if(leftProp->getProperty()->isLinkableWith(rightProp->getProperty())) {
                                    PropertyLinkArrowGraphicsItem* arrow = new PropertyLinkArrowGraphicsItem(leftProp,rightProp);
                                    arrow->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
                                    arrow->setVisible(true); autoArrows_.append(arrow);
                                }
                                else if(rightProp->getProperty()->isLinkableWith(leftProp->getProperty())) {
                                    PropertyLinkArrowGraphicsItem* arrow = new PropertyLinkArrowGraphicsItem(rightProp,leftProp);
                                    arrow->setSourceArrowHeadDirection(ConnectionBaseGraphicsItem::ArrowHeadDirectionAny);
                                    arrow->setVisible(true); autoArrows_.append(arrow);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    sourceView_->updateViewport();
}

void PropertyLinkDialog::hideAutoLinks() {
    foreach (PropertyLinkArrowGraphicsItem* arrow, autoArrows_) {
        delete arrow;
    }
    autoArrows_.clear();
    sourceView_->updateViewport();
}

void PropertyLinkDialog::confirmAutoLinks() {
    foreach (PropertyLinkArrowGraphicsItem* arrow, autoArrows_) {
        if(sourceItem_->getPropertyList()->getPropertyItems().contains(arrow->getSourceItem()))
            linkProperties(arrow->getSourceItem(),arrow->getDestinationItem());
        delete arrow;
    }
    autoArrows_.clear();
    sourceView_->updateViewport();
}

//---------------------------------------------------------------------------------------------------------------
//                  settings functions
//---------------------------------------------------------------------------------------------------------------
void PropertyLinkDialog::saveSettings() {
    QSettings settings;
    settings.beginGroup("propertylinkdialog");
    settings.setValue("position",pos());
    settings.setValue("size",size());
    settings.endGroup();
}

void PropertyLinkDialog::loadSettings() {
    //set dialog flags
    setAttribute(Qt::WA_DeleteOnClose);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    //load settings
    QSettings settings;
    settings.beginGroup("propertylinkdialog");
    if(settings.contains("position"))
        move(settings.value("position",QPoint(0,0)).toPoint());
    resize(settings.value("size",QSize(600,480)).toSize());
    settings.endGroup();
}

//---------------------------------------------------------------------------------------------------------------
//                  modification functions
//---------------------------------------------------------------------------------------------------------------
void PropertyLinkDialog::addAddModification(Modification* mod){
    foreach(Modification* m, removeModificationList_) {
        if(*m == *mod) {
            removeModificationList_.removeOne(m);
            if(m->eval_ != mod->eval_)
                evaluatorModificationList_.append(mod);
            delete m;
            return;
        }
    }
    foreach(Modification* m, evaluatorModificationList_) {
        if(*m == *mod) {
            evaluatorModificationList_.removeOne(m);
            delete m;
            break;
        }
    }
    addModificationList_.append(mod);
}

void PropertyLinkDialog::addRemoveModification(Modification* mod) {
    foreach(Modification* m, addModificationList_) {
        if(*m == *mod) {
            addModificationList_.removeOne(m);
            delete m->eval_;
            delete m;
            return;
        }
    }
    removeModificationList_.append(mod);
    foreach(Modification* m, evaluatorModificationList_) {
        if(*m == *mod) {
            evaluatorModificationList_.removeOne(m);
            delete m->eval_;
            delete m;
            return;
        }
    }
}

void PropertyLinkDialog::addEvaluatorModification(Modification* mod) {
    foreach(Modification* m, addModificationList_) {
        if(*m == *mod) {
            delete m->eval_;
            m->eval_ = mod->eval_;
            return;
        }
    }
    foreach(Modification* m, evaluatorModificationList_) {
        if(*m == *mod) {
            delete m->eval_;
            m->eval_ = mod->eval_;
            return;
        }
    }
    evaluatorModificationList_.append(mod);
}

void PropertyLinkDialog::applyAllModifications() {
    foreach(Modification* m, evaluatorModificationList_) {
        PropertyLink* link = m->src_->getProperty()->getLink(m->dst_->getProperty());
        link->setLinkEvaluator(m->eval_);
        evaluatorModificationList_.removeOne(m);
        delete m;
    }

    foreach(Modification* m, addModificationList_) {
        dynamic_cast<NetworkEditor*>(parent())->getProcessorNetwork()->createPropertyLink(
                    const_cast<Property*>(m->src_->getProperty()),const_cast<Property*>(m->dst_->getProperty()),m->eval_);
        addModificationList_.removeOne(m);
        delete m;
    }

    foreach(Modification* m, removeModificationList_) {
        dynamic_cast<NetworkEditor*>(parent())->getProcessorNetwork()->removePropertyLink(
                    const_cast<Property*>(m->src_->getProperty()),const_cast<Property*>(m->dst_->getProperty()));
        removeModificationList_.removeOne(m);
        delete m;
    }
    close();
//    emit accept();
}

void PropertyLinkDialog::undoAllModifications() {
    foreach(Modification* m, evaluatorModificationList_) {
        evaluatorModificationList_.removeOne(m);
        delete m->eval_;
        delete m;
    }

    foreach(Modification* m, addModificationList_) {
        m->src_->removeGraphicalLink(m->dst_);
        addModificationList_.removeOne(m);
        delete m->eval_;
        delete m;
    }

    foreach(Modification* m, removeModificationList_) {
        m->src_->addGraphicalLink(m->dst_);
        removeModificationList_.removeOne(m);
        delete m;
    }
}

void PropertyLinkDialog::linkProperties(PropertyGraphicsItem* src, PropertyGraphicsItem* dst) {
    tgtAssert(src != dst, "Link between the same Property is not allowed");
    //switch source left
    if(!sourceItem_->getPropertyList()->getAllPropertyItems().contains(src)){
        PropertyGraphicsItem* help = src;
        src = dst; dst = help;
    }

    //block render size links
    if(dynamic_cast<RenderSizeOriginProperty*>(const_cast<Property*>(src->getProperty())) ||
       dynamic_cast<RenderSizeOriginProperty*>(const_cast<Property*>(dst->getProperty())) ||
       dynamic_cast<RenderSizeReceiveProperty*>(const_cast<Property*>(src->getProperty())) ||
       dynamic_cast<RenderSizeReceiveProperty*>(const_cast<Property*>(dst->getProperty()))) {
           if (QMessageBox::warning(this, tr("Warning"), QString::fromStdString("<b>You are trying to link RenderSize-Properties. </b><br> <br>"\
               "Please use the RenderSize-Layer in the NetworkEditor to link such properties. Linking RenderSize-Properties "\
               "outside this layer can cause nondeterministic behavior. <br><br> <b>Are you sure you still want to link these properties?</b>"), QMessageBox::Yes | QMessageBox::Cancel) == QMessageBox::Cancel) {
                return;
           }
    }

    if(bidirectionalArrowButton_->isChecked()) {
        //check, if one direction exists already
        foreach(PropertyLinkArrowGraphicsItem* arrow, dst->getSourceLinkList()) {
            if(arrow->getDestinationItem() == src)
                return;
        }
        foreach(PropertyLinkArrowGraphicsItem* arrow, src->getSourceLinkList()) {
            if(arrow->getDestinationItem() == dst)
                return;
        }

        //first direction, if possible
        if(src->getProperty()->getCompatibleEvaluators(dst->getProperty()).size()) {
            addAddModification(new Modification(src,dst, LinkEvaluatorHelper::createEvaluator(
                        src->getProperty()->getCompatibleEvaluators(dst->getProperty())[0].first)));
            src->addGraphicalLink(dst);
        }
        //second direction
        if(dst->getProperty()->getCompatibleEvaluators(src->getProperty()).size()) {
            addAddModification(new Modification(dst,src, LinkEvaluatorHelper::createEvaluator(
                        dst->getProperty()->getCompatibleEvaluators(src->getProperty())[0].first)));
            dst->addGraphicalLink(src);
        }
    } else if (leftUpArrowButton_->isChecked()) {
        //check if link exists already
        foreach(PropertyLinkArrowGraphicsItem* arrow, dst->getSourceLinkList()) {
            if(arrow->getDestinationItem() == src)
                return;
        }

        if(dst->getProperty()->getCompatibleEvaluators(src->getProperty()).size()) {
            addAddModification(new Modification(dst,src, LinkEvaluatorHelper::createEvaluator(
                                dst->getProperty()->getCompatibleEvaluators(src->getProperty())[0].first)));
            dst->addGraphicalLink(src);
        }

    } else if (rightDownArrowButton_->isChecked()) {
        //check, if link exists already
        foreach(PropertyLinkArrowGraphicsItem* arrow, src->getSourceLinkList()) {
            if(arrow->getDestinationItem() == dst)
                return;
        }

        if(src->getProperty()->getCompatibleEvaluators(dst->getProperty()).size()) {
            addAddModification(new Modification(src,dst, LinkEvaluatorHelper::createEvaluator(
                            src->getProperty()->getCompatibleEvaluators(dst->getProperty())[0].first)));
            src->addGraphicalLink(dst);
        }
    }
    src->setVisible(true);
    dst->setVisible(true);
    setCurrentlySelectedArrow(src,dst);
}

//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
void PropertyLinkDialog::keyPressEvent(QKeyEvent* event) {
    QDialog::keyPressEvent(event);

    if (!event->isAccepted() && (event->key() == Qt::Key_Delete || event->key() == Qt::Key_Backspace))
        deleteArrowAction();
}

void PropertyLinkDialog::closeEvent(QCloseEvent* event) {
    saveSettings();
    QDialog::closeEvent(event);
}

} //voreen namespace
