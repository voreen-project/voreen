/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "plotbase.h"
#include "voreen/core/utils/exception.h"
#include "tgt/assert.h"
#include "tgt/logmanager.h"
#include <algorithm>


namespace voreen {

PlotBase::PlotBase(int keyColumnCount, int dataColumnCount)
    : keyColumnCount_(keyColumnCount)
    , dataColumnCount_(dataColumnCount)
{
    tgtAssert((keyColumnCount_ >= 0 && dataColumnCount_ >= 0), "PlotBase: negativ ColumnCounts are not allowed!!");
    int columnCount = keyColumnCount + dataColumnCount;
    columns_.assign(columnCount, Column("", EMPTY));
}

PlotBase::PlotBase(const PlotBase& rhs)
    : keyColumnCount_(rhs.keyColumnCount_)
    , dataColumnCount_(rhs.dataColumnCount_)
    , columns_(rhs.columns_)
{
}

PlotBase::~PlotBase() {
}

PlotBase& PlotBase::operator=(const PlotBase& rhs) {
    // prevent self assignment
    if (this == &rhs)
        return *this;

    keyColumnCount_ = rhs.keyColumnCount_;
    dataColumnCount_ = rhs.dataColumnCount_;
    columns_ = rhs.columns_;

    return *this;
}

int PlotBase::getKeyColumnCount() const {
    return keyColumnCount_;
}

int PlotBase::getDataColumnCount() const {
    return dataColumnCount_;
}

int PlotBase::getColumnCount() const {
    return keyColumnCount_ + dataColumnCount_;
}

const std::string& PlotBase::getColumnLabel(int column) const {
    tgtAssert((column >= 0 && column < keyColumnCount_ + dataColumnCount_), "PlotBase::getColumnLabel(): Column out of bounds.");
    return columns_[column].label_;
}

void PlotBase::setColumnLabel(int column, std::string label) {
    tgtAssert((column < keyColumnCount_ + dataColumnCount_), "PlotBase::setColumnLabel(): Column out of bounds.");
    columns_[column].label_ = label;
}

PlotBase::ColumnType PlotBase::getColumnType(int column) const {
    tgtAssert((column < keyColumnCount_ + dataColumnCount_), "PlotBase::getColumnType(): Column out of bounds.");
    return columns_[column].type_;
}

void PlotBase::setColumnType(int column, ColumnType type) {
    tgtAssert((column < keyColumnCount_ + dataColumnCount_), "PlotBase::setColumnType(): Column out of bounds.");
    columns_[column].type_ = type;
}

void PlotBase::setColumnColorHint(int column, tgt::Color colorHint) {
    tgtAssert((column < keyColumnCount_ + dataColumnCount_), "PlotBase::setColumnType(): Column out of bounds.");
    columns_[column].colorHint_= colorHint;
    columns_[column].hasColorHint_= true;
}

void PlotBase::removeColumnColorHint(int column) {
    tgtAssert((column < keyColumnCount_ + dataColumnCount_), "PlotBase::setColumnType(): Column out of bounds.");
    columns_[column].hasColorHint_= false;
}

bool PlotBase::hasColumnColorHint(int column) const {
    tgtAssert((column < keyColumnCount_ + dataColumnCount_), "PlotBase::setColumnType(): Column out of bounds.");
    return columns_[column].hasColorHint_;
}

tgt::Color PlotBase::getColumnColorHint(int column) const {
    tgtAssert((column < keyColumnCount_ + dataColumnCount_), "PlotBase::setColumnType(): Column out of bounds.");
    return columns_[column].colorHint_;
}

void PlotBase::reset(int keyColumnCount, int dataColumnCount) {
    columns_.clear();

    keyColumnCount_ = keyColumnCount;
    dataColumnCount_ = dataColumnCount;
    int columnCount = keyColumnCount+dataColumnCount;

    columns_.assign(columnCount, Column("", EMPTY));
}

bool PlotBase::compareStructure(const PlotBase* other) const {
    if (!other)
        return false;

    if (keyColumnCount_ != other->getKeyColumnCount() || dataColumnCount_ != other->getDataColumnCount())
        return false;

    if (columns_ != other->columns_)
        return false;

    return true;
}

}
