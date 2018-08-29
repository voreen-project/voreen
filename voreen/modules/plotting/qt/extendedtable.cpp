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

#include "extendedtable.h"

namespace voreen {


ExtendedTable::ExtendedTable(QWidget* parent)
    : QTableView(parent)
{
    setLocale(QLocale(QLocale::English));
}

ExtendedTable::~ExtendedTable() {
}


std::pair<QPoint,QPoint> ExtendedTable::selectedArea() const {
    QModelIndexList list = QTableView::selectedIndexes();
    int x1 = -1;
    int y1 = -1;
    int y2 = -1;
    int x2 = -1;
    for (int i = 0; i < list.size(); ++i) {
        if (i == 0) {
            x1 = list.at(i).column();
            x2 = list.at(i).column();
            y1 = list.at(i).row();
            y2 = list.at(i).row();
        }
        x1 = std::min<int>(x1,list.at(i).column());
        y1 = std::min<int>(y1,list.at(i).row());
        x2 = std::max<int>(x2,list.at(i).column());
        y2 = std::max<int>(y2,list.at(i).row());
    }
    std::pair<QPoint,QPoint> result;
    result.first = QPoint(x1,y1);
    result.second = QPoint(x2,y2);
    return result;
}

std::vector<int> ExtendedTable::selectedColumns() const {
    QModelIndexList list = QTableView::selectedIndexes();
    std::vector<int> result;
    bool match;
    for (int i = 0; i < list.size(); ++i) {
        match = true;
        for (size_t j = 0; j < result.size(); ++j) {
            if (result.at(j) == list.at(i).column()) {
                match = false;
                break;
            }
        }
        if (match) {
            result.push_back(list.at(i).column());
        }
    }
    return result;
}

std::vector<int> ExtendedTable::selectedRows() const {
    QModelIndexList list = QTableView::selectedIndexes();
    std::vector<int> result;
    bool match;
    for (int i = 0; i < list.size(); ++i) {
        match = true;
        for (size_t j = 0; j < result.size(); ++j) {
            if (result.at(j) == list.at(i).row()) {
                match = false;
                break;
            }
        }
        if (match) {
            result.push_back(list.at(i).row());
        }
    }
    return result;
}

QModelIndexList ExtendedTable::selectedIndexes() const {
    return QTableView::selectedIndexes();
}

} // namespace voreen
