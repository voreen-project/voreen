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

#include "voreen/qt/networkeditor/graphlayouts/nwegl_sugiyama.h"

//gi
#include "voreen/qt/networkeditor/graphicitems/core/processorgraphicsitem.h"

//core
#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/coprocessorport.h"

//qt
#include <QMap>
#include <QList>
#include <QPoint>
#include <QMessageBox>


namespace voreen{

NWEGL_Sugiyama::NWEGL_Sugiyama() : NWEGL_Base(), shiftX_(300), overlap_(false), median_(true), portflush_(true) {}

void NWEGL_Sugiyama::setSortParameter(qreal shift, bool overlap, bool median, bool portflush) {
    shiftX_ = shift;
    overlap_ = overlap;
    median_ = median;
    portflush_ = portflush;
}

void NWEGL_Sugiyama::sort(ProcessorNetwork* network, std::vector<Processor*>* processors, QMap<Processor*, ProcessorGraphicsItem*>* itemMap) {
    //topologically ordered list of all processors
    std::vector<Processor*> procOrder = getRenderingOrder(network);

    //list of all endprocessors
    std::vector<Processor*> firstLevel = getEndProcessors(network);

    //list of all source processors
    std::vector<Processor*> sourceProc = getSourceProcessors(network);

    //list of all processors
    std::vector<Processor*> allProc = network->getProcessors();

    //copy selected processors in vector
    std::vector<Processor*> selected(*processors);

    //list of list of layer (sorted)
    std::vector<std::vector<Processor*> > layer;

    //if no or just 1 processor is selected, sort whole network
    if (processors->size() <= 1) {
        layer = getGraphLayer(network->getProcessors());
    }

    //else take only selected processors
    else  {
        layer = getGraphLayer(selected);
    }

    //sort the layer
    std::vector<std::vector<Processor*> > layerSorted = sortGraphLayer(layer);
    layer = layerSorted;

    //saves position of CoProcessor
    QPointF copos;

    //list of all processors that are connected to coProcessor
    std::vector<Processor*> connect;
    connect.clear();

    //list of processors with no position (not sorted)
    std::vector<Processor*> noPos;
    noPos.clear();

    //actual Layer that needs to be drawn
    std::vector<Processor*> actLayer;

    //previous Layer that has been drawn
    std::vector<Processor*> lastLayer;

    //variables for distances between processors on x and y axis
    qreal shiftY = 100;

    //lower left corner to start sorting
    QPointF lowerLeftCorner;

    if (processors->size() <= 1) {
        qreal x;
        qreal y;
        QPointF procPos;
        QPointF oldProcPos;
        for (size_t j = 0; j < allProc.size(); j++) {
            QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
            while (i != itemMap->constEnd()) {
                if (allProc[j] == i.key()) {
                    if (j == 0) {
                        oldProcPos = i.value()->pos();
                        x = i.value()->pos().x();
                        y = i.value()->pos().y();
                    }
                    else {
                        procPos = i.value()->pos();
                        if (procPos.x() < oldProcPos.x()) {
                            x = procPos.x();
                            oldProcPos.setX(x);
                        }
                        if (procPos.y() > oldProcPos.y()) {
                            y = procPos.y();
                            oldProcPos.setY(y);
                        }
                    }
                    ++i;
                }
                else ++i;
            }
        }
        lowerLeftCorner.setX(x);
        lowerLeftCorner.setY(y);
    }
    else {
        qreal x;
        qreal y;
        QPointF procPos;
        QPointF oldProcPos;
        for (size_t j = 0; j < selected.size(); j++) {
            QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
            while (i != itemMap->constEnd()) {
                if (selected[j] == i.key()) {
                    if (j == 0) {
                        oldProcPos = i.value()->pos();
                        x = i.value()->pos().x();
                        y = i.value()->pos().y();
                    }
                    else {
                        procPos = i.value()->pos();
                        if (procPos.x() < oldProcPos.x()) {
                            x = procPos.x();
                            oldProcPos.setX(x);
                        }
                        if (procPos.y() > oldProcPos.y()) {
                            y = procPos.y();
                            oldProcPos.setY(y);
                        }
                    }
                    ++i;
                }
                else ++i;
            }
        }
        lowerLeftCorner.setX(x);
        lowerLeftCorner.setY(y);
    }

    //variable to save y position of previous layer
    qreal newY;
    qreal oldY;

    newY = lowerLeftCorner.y();
    oldY = 0;

    if(overlap_==false) {

        //calculate shiftX, so that processors will not overlap each other
        for(size_t j = 0; j < procOrder.size(); j++){
            QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
            while (i != itemMap->constEnd()) {
                if(procOrder[j] == i.key()) {
                    qreal width = i.value()->boundingRect().width();
                    if (width > shiftX_) {
                        //shiftX is summed up to next 100 of width (e.g.: 366 -> 400, 122 -> 200, ...)
                        shiftX_ = width + (100 - ((int)width)%100);
                    }

                ++i;
                }
                else ++i;
            }
        }
    }

    //setPos for processors for all layers
    for (size_t it = 0; it < layer.size(); it++) {
        actLayer.clear();
        actLayer = layer[it];
        LDEBUGC("NWEGL_Sugiyama", "# of Elements in layer " << it << ": " << actLayer.size());

        for (size_t j = 0; j < actLayer.size(); j++) {
            QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
            while (i != itemMap->constEnd()) {
                if (i.key() == actLayer[j]) {

                    //only calculate mid portflush is true
                    qreal mid = 0;
                    if(portflush_ == true) {
                        mid = (shiftX_ - i.value()->boundingRect().width()) / 2;
                    }

                    //setPos for each processor inside the actual layer

                    //if actual layer is bigger than last layer...
                    if (actLayer.size() > lastLayer.size()+1 || actLayer.size()+1 < lastLayer.size()) {

                        //in first layer start with newY position
                        if (it==0) {
                            i.value()->setPos(lowerLeftCorner.x() + mid + shiftX_*j, newY);
                            oldY = i.value()->pos().y();
                        }
                        else {
                            i.value()->setPos(lowerLeftCorner.x() + shiftX_*j + mid, (newY-(shiftY*2)));
                            oldY = i.value()->pos().y();
                        }
                    }

                    //...or if more than 3 processors are connected to outport of one processor
                    //then increase distance on y-axis (to have more space for arrows)

                    //TODO: this are just 2 example scenarios where distance on y-axis could be increased
                    else if (getConnectedToOut(i.key()).size() > 3) {

                        //in first layer start with newY position
                        if (it==0) {
                            i.value()->setPos(lowerLeftCorner.x() + shiftX_*j + mid, newY);
                            oldY = i.value()->pos().y();
                        }
                        else {
                            i.value()->setPos(lowerLeftCorner.x() + shiftX_*j + mid, newY-(shiftY*3));
                            oldY = i.value()->pos().y();
                        }
                    }
                    //else shorter distance on y-axis
                    else {

                        //in first layer start with newY position
                        if (it==0) {
                            i.value()->setPos(lowerLeftCorner.x() + shiftX_*j + mid, newY);
                            oldY = i.value()->pos().y();
                        }
                        else {
                            i.value()->setPos(lowerLeftCorner.x() + shiftX_*j + mid, (newY-shiftY));
                            oldY = i.value()->pos().y();
                        }
                    }
                    ++i;
                }
                else ++i;
            }
        }
        newY = oldY;
        lastLayer = actLayer;
    }

    // saves position of processors in layer above
    std::vector<QPointF> abovePos;
    std::vector<qreal> aboveWidth;

    // variable to save differences between x positions of two layers
    qreal xDiff;

    // to start with the top of the network reverse list first
    std::reverse(layer.begin(), layer.end());

    //saves median of x position for every layer
    if(median_ == true){
        for (size_t it = 0; it < layer.size() ; it++) {

            actLayer.clear();
            actLayer = layer[it];
            double sum = 0;

            for(size_t act = 0; act < actLayer.size(); act++) {

                QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
                while (i != itemMap->constEnd()) {
                    if(actLayer[act] ==  i.key()) {
                        sum += i.value()->pos().x();
                        ++i;
                    }
                    else ++i;
                }
            }

            //median calculation
            sum = sum/(layer[it].size());

            if ( it != 0) {

                // if lastLayer is bigger than actLayer shift all processors on the x axis dependent on the difference in size
                if (lastLayer.size() > actLayer.size()) {
                    abovePos.clear();

                    // save position of processors in layer above
                    for (size_t k = 0; k < lastLayer.size(); k++) {
                        QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
                        while (i != itemMap->constEnd()) {
                            if (lastLayer[k]==i.key()) {

                                abovePos.push_back(i.value()->pos());
                                ++i;
                            }
                            else ++i;
                        }
                    }

                    //diff saves differences of layer sizes
                    size_t diff = lastLayer.size()-actLayer.size();
                    for (size_t actIt = 0; actIt< actLayer.size(); actIt++) {
                        QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
                        while (i != itemMap->constEnd()) {
                            if (actLayer[actIt] == i.key()) {
                                QPointF actPos = i.value()->pos();

                                if ( actIt == 0) {
                                    xDiff = abovePos[0].x() - i.value()->pos().x();
                                }

                                i.value()->setX(actPos.x()+xDiff+diff*(shiftX_/2));
                                ++i;
                            }
                            else ++i;
                        }
                    }
                }

                // if lastLayer is smaller than actLayer shift all processors on the x axis dependent on the difference in size
                if (lastLayer.size() < actLayer.size()) {
                    abovePos.clear();

                    // save position of processors in layer above
                    for (size_t k = 0; k < lastLayer.size(); k++) {
                        QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
                        while (i != itemMap->constEnd()) {
                            if (lastLayer[k]==i.key()) {
                                abovePos.push_back(i.value()->pos());
                                ++i;
                            }
                            else ++i;
                        }
                    }

                    //diff saves differences of layer sizes
                    size_t diff = actLayer.size()-lastLayer.size();

                    for (size_t actIt = 0; actIt< actLayer.size(); actIt++) {
                        QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
                        while (i != itemMap->constEnd()) {
                            if (actLayer[actIt] == i.key()) {
                                QPointF actPos = i.value()->pos();

                                //xDiff for the first element of layers is enough
                                if ( actIt == 0) {
                                    xDiff = i.value()->pos().x() - abovePos[0].x();
                                }
                                i.value()->setX(actPos.x()-xDiff-diff*(shiftX_/2));
                                ++i;
                            }
                            else ++i;
                        }
                    }
                }

                // if layers have same size, put processors on same x-axis like processors above
                if (lastLayer.size() == actLayer.size()) {
                    abovePos.clear();
                    aboveWidth.clear();

                    // save position and width of processors in layer above
                    for (size_t k = 0; k < lastLayer.size(); k++) {
                        QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
                        while (i != itemMap->constEnd()) {
                            if (lastLayer[k]==i.key()) {
                                abovePos.push_back(i.value()->pos());
                                aboveWidth.push_back(i.value()->boundingRect().width());
                                ++i;
                            }
                            else ++i;
                        }
                    }

                    // set x position for actual layer
                    for (size_t actIt = 0; actIt< actLayer.size(); actIt++) {
                        QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
                        while (i != itemMap->constEnd()) {
                            if (actLayer[actIt] == i.key()) {
                                //QPointF actPos = i.value()->pos();
                                qreal actWidth = i.value()->boundingRect().width();

                                //positioning with portflush
                                if (portflush_==true) {
                                    if(actWidth < aboveWidth[actIt]) {
                                        qreal mid = (aboveWidth[actIt] - actWidth) / 2;

                                        i.value()->setX(abovePos[actIt].x()+mid);
                                        ++i;
                                    }
                                    else if(actWidth > aboveWidth[actIt]) {
                                        qreal mid = (actWidth - aboveWidth[actIt]) / 2;

                                        i.value()->setX(abovePos[actIt].x()-mid);
                                        ++i;
                                    }
                                    else {
                                        i.value()->setX(abovePos[actIt].x());
                                        ++i;
                                    }
                                }
                                //positioning without portflush
                                else {
                                    i.value()->setX(abovePos[actIt].x());
                                    ++i;
                                }
                            }
                            else ++i;
                        }
                    }
                }
            }
            lastLayer = actLayer;
        }
    }
    // reverse AGAIN to start at bottom again
    std::reverse(layer.begin(), layer.end());



    //calculate newLowerLeftCorner after positioning of all !SELECTED! processors (only selected processors)
    QPointF newLowerLeftCorner;
    if(processors->size() > 1) {
        qreal x;
        qreal y;
        QPointF procPos;
        QPointF oldProcPos;
        for (size_t j = 0; j < selected.size(); j++) {
            QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
            while (i != itemMap->constEnd()) {
                if (selected[j] == i.key()) {
                    if (j == 0) {
                        oldProcPos = i.value()->pos();
                        x = i.value()->pos().x();
                        y = i.value()->pos().y();
                    }
                    else {
                        procPos = i.value()->pos();
                        if (procPos.x() < oldProcPos.x()) {
                            x = procPos.x();
                            oldProcPos.setX(x);
                        }
                        if (procPos.y() > oldProcPos.y()) {
                            y = procPos.y();
                            oldProcPos.setY(y);
                        }
                    }
                    ++i;
                }
                else ++i;
            }
        }
        newLowerLeftCorner.setX(x);
        newLowerLeftCorner.setY(y);
    }

    //all sorted processors in one list
    std::vector<Processor*> allSorted;
    allSorted.clear();

    //fill up sorted processors
    for (size_t it = 0; it < layer.size(); it++) {
        actLayer.clear();
        actLayer = layer[it];

        for (size_t it2 = 0; it2 < actLayer.size(); it2++) {
            allSorted.push_back(actLayer[it2]);
        }
    }

    //calculate newLowerLeftCorner after positioning of all !SORTED! processors (whole network)
    if(processors->size() <= 1) {
        qreal x;
        qreal y;
        QPointF procPos;
        QPointF oldProcPos;
        for (size_t j = 0; j < allSorted.size(); j++) {
            QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
            while (i != itemMap->constEnd()) {
                if (allSorted[j] == i.key()) {
                    if (j == 0) {
                        oldProcPos = i.value()->pos();
                        x = i.value()->pos().x();
                        y = i.value()->pos().y();
                    }
                    else {
                        procPos = i.value()->pos();
                        if (procPos.x() < oldProcPos.x()) {
                            x = procPos.x();
                            oldProcPos.setX(x);
                        }
                        if (procPos.y() > oldProcPos.y()) {
                            y = procPos.y();
                            oldProcPos.setY(y);
                        }
                    }
                    ++i;
                }
                else ++i;
            }
        }
        newLowerLeftCorner.setX(x);
        newLowerLeftCorner.setY(y);
    }

    //if processor has CoProcessorInports, set ProcessorGraphicsItem of connected Coprocessor left to partner
    for (size_t it = 0; it < procOrder.size(); it++) {
         QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
         while (i != itemMap->constEnd()) {
             if(procOrder[it] == i.key()) {
                if(i.key()->getCoProcessorInports().size()!=0) {

                    //save position of ProcessorGraphicsItem from Processor with CoProcessorInport
                    copos = i.value()->pos();

                    //get connected Coprocessor
                    std::vector<CoProcessorPort*> connectedCoPort = i.key()->getCoProcessorInports();
                    for(size_t portIt=0; portIt < connectedCoPort.size(); portIt++) {
                        connect = connectedCoPort[portIt]->getConnectedProcessors();
                    }

                    //search for connected Coprocessor and setPos left to its partner
                    for(size_t connIt = 0; connIt < connect.size(); connIt++){

                        //only calculate position of coprocessor if whole network is sorted, or coprocessor is selected
                        if(std::find(selected.begin(), selected.end(), connect[connIt])!=selected.end() || processors->size() <= 1) {
                            QMap<Processor*,ProcessorGraphicsItem*>::const_iterator j = itemMap->constBegin();
                            while (j != itemMap->constEnd()) {
                                if (j.key() == connect[connIt]) {
                                    QPointF newPoint;

                                    //diff saves differenc between lowerleft corner and coprocessor
                                    qreal diff = std::fabs(copos.x() - newLowerLeftCorner.x());

                                    newPoint.setX(copos.x()- diff -shiftX_);
                                    newPoint.setY(copos.y()+shiftY*connIt);
                                    j.value()->setPos(newPoint);
                                    ++j;
                                }
                                else ++j;
                            }
                        }
                    }
                }
                ++i;
            }
            else ++i;
        }
    }

    //do this just when whole network is sorted

    if (processors->size() <= 1) {
        //if coprocessor is connected to source, then draw source above coprocessors
        for (size_t it = 0; it < procOrder.size(); it++) {
            if (procOrder[it]->getCoProcessorOutports().size()!=0 && procOrder[it]->getInports().size()!=0) {

                //find source processors that are connected to coprocessor
                std::vector<Port*> in = procOrder[it]->getInports();
                for (size_t inIt = 0; inIt < in.size(); inIt++) {
                    std::vector<const Port*> connPort = in[inIt]->getConnected();
                    for (size_t portIt = 0; portIt < connPort.size(); portIt++) {

                        // only use source processors
                        if (connPort[portIt]->getProcessor()->isSource()) {
                            Processor* connProc = connPort[portIt]->getProcessor();
                            QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
                            while (i != itemMap->constEnd()) {
                                if (i.key() == connProc) {

                                    //only shift on x-axis to set new position
                                    QPointF newPos;
                                    QPointF oldPos;
                                    oldPos = i.value()->pos();

                                    //diff saves differenc between lowerleft corner and coprocessor
                                    qreal diff = std::fabs(copos.x() - newLowerLeftCorner.x());

                                    newPos.setX(copos.x()- diff - shiftX_);
                                    newPos.setY(oldPos.y());
                                    i.value()->setPos(newPos);
                                    ++i;
                                }
                                else ++i;
                            }
                        }
                    }
                }
            }
        }


        //look for processors that are not in any layer and save it in noPos vector
        for (size_t it = 0; it < allProc.size(); it++) {
            bool found = false;
            for (size_t it2 = 0; it2 < allSorted.size(); it2++) {
                if (allProc[it] == allSorted[it2]) {
                    found = true;
                }
            }
            if (found == false) {
                //make sure that the processor without a layer is not a coprocessor
                if (allProc[it]->getCoProcessorOutports().size()==0) {
                    noPos.push_back(allProc[it]);
                }
            }
        }

        //reverse noPos vector for layout reasons
        std::reverse(noPos.begin(), noPos.end());

        // if there are processors that have no position yet, sort them in the lower left corner
        if (!noPos.empty()) {
            for (size_t it = 0; it < noPos.size(); it++) {
                QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
                while (i != itemMap->constEnd()) {
                     if(noPos[it] == i.key()) {
                         QPointF newPoint;
                         newPoint.setX(newLowerLeftCorner.x()-shiftX_);
                         newPoint.setY(newLowerLeftCorner.y()+shiftY*it);
                         i.value()->setPos(newPoint);
                         ++i;
                     }
                     else ++i;
                }
            }
        }
    }

    //if the lower left corner differs after the sorting, put processors back in right position (so sort button can be used multiple times)
    //for selected processors only
    if(processors->size() > 1) {
        qreal x;
        qreal y;
        QPointF procPos;
        QPointF oldProcPos;
        for (size_t j = 0; j < selected.size(); j++) {
            QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
            while (i != itemMap->constEnd()) {
                if (selected[j] == i.key()) {
                    if (j == 0) {
                        oldProcPos = i.value()->pos();
                        x = i.value()->pos().x();
                        y = i.value()->pos().y();
                    }
                    else {
                        procPos = i.value()->pos();
                        if (procPos.x() < oldProcPos.x()) {
                            x = procPos.x();
                            oldProcPos.setX(x);
                        }
                        if (procPos.y() > oldProcPos.y()) {
                            y = procPos.y();
                            oldProcPos.setY(y);
                        }
                    }
                    ++i;
                }
                else ++i;
            }
        }
        newLowerLeftCorner.setX(x);
        newLowerLeftCorner.setY(y);

        if (newLowerLeftCorner.x() < lowerLeftCorner.x()) {
            qreal diff = std::fabs(lowerLeftCorner.x() - newLowerLeftCorner.x());
            for (size_t j = 0; j < selected.size(); j++) {
                QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
                while (i != itemMap->constEnd()) {
                    if (selected[j] == i.key()) {
                        QPointF pos = i.value()->pos();
                        i.value()->setX(pos.x()+diff);
                        ++i;
                    }
                    else ++i;
                }
            }
        }
    }
    //if the lower left corner differs after the sorting, put processors back in right position (so sort button can be used multiple times)
    //for whole network
    else {
        qreal x;
        qreal y;
        QPointF procPos;
        QPointF oldProcPos;
        for (size_t j = 0; j < allProc.size(); j++) {
            QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
            while (i != itemMap->constEnd()) {
                if (allProc[j] == i.key()) {
                    if (j == 0) {
                        oldProcPos = i.value()->pos();
                        x = i.value()->pos().x();
                        y = i.value()->pos().y();
                    }
                    else {
                        procPos = i.value()->pos();
                        if (procPos.x() < oldProcPos.x()) {
                            x = procPos.x();
                            oldProcPos.setX(x);
                        }
                        if (procPos.y() > oldProcPos.y()) {
                            y = procPos.y();
                            oldProcPos.setY(y);
                        }
                    }
                    ++i;
                }
                else ++i;
            }
        }
        newLowerLeftCorner.setX(x);
        newLowerLeftCorner.setY(y);
        if (newLowerLeftCorner.x() < lowerLeftCorner.x()) {
            qreal diff = std::fabs(lowerLeftCorner.x() - newLowerLeftCorner.x());
            for (size_t j = 0; j < allProc.size(); j++) {
                QMap<Processor*,ProcessorGraphicsItem*>::const_iterator i = itemMap->constBegin();
                while (i != itemMap->constEnd()) {
                    if (allProc[j] == i.key()) {
                        QPointF pos = i.value()->pos();
                        i.value()->setX(pos.x()+diff);
                        ++i;
                    }
                    else ++i;
                }
            }
        }
    }
}

std::vector<std::vector<Processor*> > NWEGL_Sugiyama::sortGraphLayer(std::vector<std::vector<Processor*> > graphLayer) {

    //sortedLayers for return
    std::vector<std::vector<Processor*> > sortedLayers;

    //saves actual and last layer
    std::vector<Processor*> actLayer;
    std::vector<Processor*> lastLayer;

    //saves sorted layer
    std::vector<Processor*> sortedLayer;

    //initialise
    double bary = 0;
    lastLayer = graphLayer.front();
    sortedLayers.push_back(lastLayer);
    int count = 0;
    Processor* save = 0;

    //saves bary for each processor
    std::vector<std::pair<double,Processor*> > baryProcList;

    for (size_t layerIt = 1; layerIt < graphLayer.size(); layerIt++) {
        actLayer = graphLayer[layerIt];
        //if(actLayer.size() >= lastLayer.size()) {
            for (size_t actIt = 0; actIt < actLayer.size(); actIt++) {
                std::vector<Processor*> connect = getConnectedToOut(actLayer[actIt]);
                for (size_t connIt = 0; connIt < connect.size(); connIt++) {
                    for (size_t lastIt = 0; lastIt < lastLayer.size(); lastIt++) {

                        //if connected with lastLayer increase bary and count
                        if (connect[connIt] == lastLayer[lastIt]) {

                            //check if processor already has a connection to proccessor in last layer
                            if (save != lastLayer[lastIt]) {
                                bary = (double)(bary+(lastIt+1));
                                count++;
                                save = lastLayer[lastIt];
                            }
                        }
                    }
                    save = 0;
                }
                bary = (double) (bary / count);

                //save bary for each process and push in baryProcList
                std::pair<double,Processor*> baryProc = std::make_pair(bary,actLayer[actIt]);
                baryProcList.push_back(baryProc);
                bary = 0;
                count = 0;
            }

            //sort baryProcList based on bary
            std::stable_sort(baryProcList.begin(), baryProcList.end(), pairCompare);

            //push sorted processors
            for (size_t baryIt = 0; baryIt < baryProcList.size(); baryIt++) {
                sortedLayer.push_back(baryProcList[baryIt].second);
            }
        //}

        //TODO: sort upside down
        /*else{
            for (size_t j = 0; j < lastLayer.size(); j++) {
            std::vector<Processor*> connect = getConnectedToIn(lastLayer[j]);
                for (size_t k = 0; k < connect.size(); k++) {
                    bary = bary+(k+1);
                }
                bary = bary / connect.size();
                std::pair<int,Processor*> baryProc(bary,lastLayer[j]);
                baryProcMap.insert(baryProc);
                bary = 0;
            }
            //std::sort(baryProcMap.begin(),baryProcMap.end());
            for (std::map<int,Processor*>::iterator it=baryProcMap.begin(); it!=baryProcMap.end(); ++it)
                it->second();
                sortedLayer.push_back(it->second());
            for (size_t l = 0; l < baryProcMap.size(); l++) {
                sortedLayer.push_back(baryProcMap[l]);
            }
        }*/
        lastLayer = actLayer;
        sortedLayers.push_back(sortedLayer);

        sortedLayer.clear();
        baryProcList.clear();
        actLayer.clear();
    }
    return sortedLayers;
}

bool NWEGL_Sugiyama::pairCompare(const std::pair<double, Processor*>& firstElem, const std::pair<double, Processor*>& secondElem) {
    return firstElem.first < secondElem.first;
}

} // namespace voreen
