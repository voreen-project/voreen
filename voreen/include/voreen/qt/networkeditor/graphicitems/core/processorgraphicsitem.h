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

#ifndef VRN_PROCESSORGRAPHICSITEM_H
#define VRN_PROCESSORGRAPHICSITEM_H

#include "portownergraphicsitem.h"

namespace voreen {

class Processor;

/**
 * Qt graphic item class corresponding to voreen/core/processor
 */
class ProcessorGraphicsItem : public PortOwnerGraphicsItem, public ProcessorObserver {
Q_OBJECT
public:
    //constructor + destructor
    ProcessorGraphicsItem(Processor* processor, NetworkEditor* networkEditor);
    virtual ~ProcessorGraphicsItem();

    //---------------------------------------------------------------------------------
    //      nwebasegraphicsitem functions
    //---------------------------------------------------------------------------------
    virtual int type() const {return UserTypesProcessorGraphicsItem;}

    //style
    virtual QRectF boundingRect() const;
    virtual QPainterPath shape() const;
protected:
    virtual void initializePaintSettings();
    virtual void mainPaint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget, NWEItemSettings& setting);
public:
    virtual void layoutChildItems();
    //---------------------------------------------------------------------------------
    //      portownergraphicsitem functions
    //---------------------------------------------------------------------------------
public:
    virtual QList<Processor*> getProcessors() const;
    virtual void saveMeta();
    virtual void loadMeta();
    virtual void renameFinished();
    virtual void toggleProcessorWidget();
    virtual void togglePropertyList();

protected:
    virtual void registerPorts();

    //---------------------------------------------------------------------------------
    //      processorobserver functions
    //---------------------------------------------------------------------------------
public:
    virtual void processorWidgetCreated(const Processor* processor);
    virtual void processorWidgetDeleted(const Processor* processor);
    virtual void portsChanged(const Processor* processor);
    virtual void propertiesChanged(const PropertyOwner* owner);
    virtual void stateChanged(const Processor* processor);

    //---------------------------------------------------------------------------------
    //      processor  functions
    //---------------------------------------------------------------------------------
public:
    //getter and setter
    virtual Processor* getProcessor() const;

    //events
    virtual void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* event);
protected:
    //members
    Processor* processor_; ///core item
};

} //namespace voreen

#endif // VRN_PROCESSORGRAPHICSITEM_H
