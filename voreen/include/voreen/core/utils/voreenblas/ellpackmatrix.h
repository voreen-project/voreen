/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_ELLPACKMATRIX_H
#define VRN_ELLPACKMATRIX_H

#include "tgt/assert.h"
#include "tgt/types.h"
#include "tgt/logmanager.h"
#include "voreen/core/utils/exception.h"

#include <string>
#include <sstream>

namespace voreen {

template<class T>
class EllpackMatrix {
public:
    EllpackMatrix();
    EllpackMatrix(size_t numRows, size_t numCols, size_t numColsPerRow);
    ~EllpackMatrix();

    EllpackMatrix(const EllpackMatrix& other) = delete;
    EllpackMatrix& operator=(const EllpackMatrix& other) = delete;

    EllpackMatrix(EllpackMatrix&& other);
    EllpackMatrix& operator=(EllpackMatrix&& other);

    /**
     * Initializes the data buffers and throws an exception on failure (bad allocation).
     */
    void initializeBuffers();
    /// Returns whether the data buffers have already been initialized.
    bool isInitialized() const;

    /**
     * Changes the matrix dimensions. Only allowed, if matrix has not been initialized.
     */
    void setDimensions(size_t numRows, size_t numCols, size_t numColsPerRow);

    inline void setValue(size_t row, size_t col, T value);
    inline T getValue(size_t row, size_t col) const;

    // Careful: an invalid index aborts the process!
    inline T& getWritableValue(size_t row, size_t col);

    inline void setValueByIndex(size_t row, size_t col, size_t colIndex, T value);
    inline T getValueByIndex(size_t row, size_t colIndex) const;

    inline size_t getColumn(size_t row, size_t colIndex) const;
    inline int getIndex(size_t row, size_t col) const;
    inline size_t getNumRowEntries(size_t row) const;

    T* getMatrix() const;
    size_t* getIndices() const;

    size_t getNumRows() const;
    size_t getNumCols() const;
    size_t getNumColsPerRow() const;

    size_t getMatrixBufferSize() const;
    size_t getIndicesBufferSize() const;

    bool isQuadratic() const;
    bool isSymmetric() const;

    std::string toString() const;

private:
    inline std::pair<size_t, size_t> rowIndexRange(size_t row) const;
    inline size_t internalIndex(size_t row, size_t col) const;
    inline size_t getNextFreeIndex(size_t row) const;

    size_t numRows_;
    size_t numCols_;
    size_t numColsPerRow_;

    T* M_;
    size_t* indices_;
};

} //namespace


// ------------------------------------------------------------------------------
// template definitions

template<class T>
voreen::EllpackMatrix<T>::EllpackMatrix() :
    numRows_(0),
    numCols_(0),
    numColsPerRow_(0),
    M_(0),
    indices_(0)
{}

template<class T>
voreen::EllpackMatrix<T>::EllpackMatrix(size_t numRows, size_t numCols, size_t numColsPerRow) :
    numRows_(numRows),
    numCols_(numCols),
    numColsPerRow_(numColsPerRow),
    M_(0),
    indices_(0)
{}

template<class T>
voreen::EllpackMatrix<T>::EllpackMatrix(voreen::EllpackMatrix<T>&& other) :
    numRows_(other.numRows_),
    numCols_(other.numCols_),
    numColsPerRow_(other.numColsPerRow_),
    M_(other.M_),
    indices_(other.indices_)
{
    M_ = nullptr;
    indices_ = nullptr;
}

template<class T>
voreen::EllpackMatrix<T>& voreen::EllpackMatrix<T>::operator=(voreen::EllpackMatrix<T>&& other) {
    if(&other != this) {
        this->~EllpackMatrix();
        new(this) EllpackMatrix(std::move(other));
    }
    return *this;
}

template<class T>
voreen::EllpackMatrix<T>::~EllpackMatrix() {
    delete[] M_;
    delete[] indices_;
}

template<class T>
void voreen::EllpackMatrix<T>::initializeBuffers() {

    tgtAssert(!M_ && !indices_, "Data buffers already initialized");
    tgtAssert(numRows_ > 0 && numCols_ > 0, "Invalid matrix dimensions");
    tgtAssert(numColsPerRow_ > 0, "Invalid column count");

    try {
        M_ = new T[numRows_*numColsPerRow_];
        indices_ = new size_t[numRows_*numColsPerRow_];
    }
    catch (std::bad_alloc& /*e*/){
        delete[] M_;
        M_ = 0;
        delete[] indices_;
        indices_ = 0;
        //LERRORC("voreen.OpenCL.VoreenBlas", "EllMatrix::initializeBuffers(): " <<
        //    "Bad allocation during initialization of data buffers");
        throw VoreenException("Bad allocation during initialization of data buffers");
    }

    size_t numElms = numRows_*numColsPerRow_;
    std::fill_n(M_, numElms, static_cast<T>(0));
    std::fill_n(indices_, numElms, -1);
}

template<class T>
bool voreen::EllpackMatrix<T>::isInitialized() const {
    return M_ && indices_;
}

template<class T>
void voreen::EllpackMatrix<T>::setDimensions(size_t numRows, size_t numCols, size_t numColsPerRow) {
    tgtAssert(!isInitialized(), "EllMatrix already initialized");

    numRows_ = numRows;
    numCols_ = numCols;
    numColsPerRow_ = numColsPerRow;
}

template<class T>
void voreen::EllpackMatrix<T>::setValue(size_t row, size_t col, T value) {

#ifdef VRN_BLAS_DEBUG
    tgtAssert(M_ && indices_, "Data buffers not initialized");
    tgtAssert(row < numRows_ && col < numCols_, "Invalid indices");
#endif

    // value already set for these coordinates?
    int index = getIndex(row, col);

    // if no value for indices assigned, get next free col index
    if (index == -1)
        index = static_cast<int>(getNextFreeIndex(row));

    if (index >= 0) {
        // next free col index available
        M_[index] = value;
        indices_[index] = col;
    }
    else {
        // row is full
        LWARNINGC("EllMatrix", "Too many entries for row " << row);
    }
}

template<class T>
T voreen::EllpackMatrix<T>::getValue(size_t row, size_t col) const {

#ifdef VRN_BLAS_DEBUG
    tgtAssert(M_ && indices_, "Data buffers not initialized");
    tgtAssert(row < numRows_ && col < numCols_, "Invalid indices");
#endif

    int index = getIndex(row, col);
    if (index >= 0) {
        return M_[index];
    }
    else {
        // No entry for matrix element
        return static_cast<T>(0);
    }
}

template<class T>
T& voreen::EllpackMatrix<T>::getWritableValue(size_t row, size_t col) {

#ifdef VRN_BLAS_DEBUG
    tgtAssert(M_ && indices_, "Data buffers not initialized");
    tgtAssert(row < numRows_ && col < numCols_, "Invalid indices");
#endif

    auto range = rowIndexRange(row);
    for (size_t index = range.first; index != range.second; ++index) {
        auto& c = indices_[index];
        // Cols are filled from 0 to .. and never freed. Thus, if we find -1, the entry does not exist, yet.
        if (c == col || c == -1) {
            c = col;
            return M_[index];
        }
    }
    //LWARNINGC("EllMatrix", "Too many entries for row " << row);
    tgtAssert(false, "Too many entries for row");
    std::abort();
}


template<class T>
size_t voreen::EllpackMatrix<T>::getColumn(size_t row, size_t colIndex) const {
#ifdef VRN_BLAS_DEBUG
    tgtAssert(M_ && indices_, "Data buffers not initialized");
    tgtAssert(row < numRows_, "Invalid row");
    tgtAssert(colIndex < numColsPerRow_, "Invalid col index");
#endif
    return indices_[internalIndex(row, colIndex)];
}

template<class T>
int voreen::EllpackMatrix<T>::getIndex(size_t row, size_t col) const {

#ifdef VRN_BLAS_DEBUG
    tgtAssert(M_ && indices_, "Data buffers not initialized");
    tgtAssert(row < numRows_ && col < numCols_, "Invalid indices");
#endif

    auto range = rowIndexRange(row);
    for (size_t index = range.first; index != range.second; ++index) {
        if (indices_[index] == col)
            return static_cast<int>(index);
    }
    return -1;
}

template<class T>
size_t voreen::EllpackMatrix<T>::getNumRowEntries(size_t row) const {

#ifdef VRN_BLAS_DEBUG
    tgtAssert(M_ && indices_, "Data buffers not initialized");
    tgtAssert(row >= 0 && row < numRows_, "Invalid row");
#endif

    size_t entries = 0;
    for (size_t i=0; i<numColsPerRow_; i++) {
        if (indices_[internalIndex(row, i)] == -1) {
            entries++;
        }
    }
    return entries;
}

template<class T>
T voreen::EllpackMatrix<T>::getValueByIndex(size_t row, size_t colIndex) const {

#ifdef VRN_BLAS_DEBUG
    tgtAssert(M_ && indices_, "Data buffers not initialized");
    tgtAssert(row < numRows_, "Invalid row");
    tgtAssert(colIndex < numColsPerRow_, "Invalid col index");
#endif

    return M_[internalIndex(row, colIndex)];
}

template<class T>
void voreen::EllpackMatrix<T>::setValueByIndex(size_t row, size_t col, size_t colIndex, T value) {

#ifdef VRN_BLAS_DEBUG
    tgtAssert(M_ && indices_, "Data buffers not initialized");
    tgtAssert(row < numRows_ && col < numCols_, "Invalid indices");
    tgtAssert(colIndex < numColsPerRow_, "Invalid col index");
#endif

    M_[internalIndex(row, colIndex)] = value;
    indices_[internalIndex(row, colIndex)] = col;
}

template<class T>
T* voreen::EllpackMatrix<T>::getMatrix() const {
    return M_;
}

template<class T>
size_t* voreen::EllpackMatrix<T>::getIndices() const {
    return indices_;
}

template<class T>
size_t voreen::EllpackMatrix<T>::getNumRows() const {
    return numRows_;
}

template<class T>
size_t voreen::EllpackMatrix<T>::getNumCols() const {
    return numCols_;
}

template<class T>
size_t voreen::EllpackMatrix<T>::getNumColsPerRow() const {
    return numColsPerRow_;
}

template<class T>
size_t voreen::EllpackMatrix<T>::getMatrixBufferSize() const {
    return (sizeof(T) * numRows_ * numColsPerRow_);
}

template<class T>
size_t voreen::EllpackMatrix<T>::getIndicesBufferSize() const {
    return (sizeof(size_t) * numRows_ * numColsPerRow_);
}

template<class T>
bool voreen::EllpackMatrix<T>::isQuadratic() const {
    return (numRows_ == numCols_);
}

template<class T>
bool voreen::EllpackMatrix<T>::isSymmetric() const {
    if (!isQuadratic())
        return false;

    for (size_t row=0; row < numRows_; row++) {
        for (size_t colIndex=0; colIndex<getNumColsPerRow(); colIndex++) {
            size_t col = getColumn(row, colIndex);
            if(col == -1) {
                break;
            }
            if (getValue(row, col) != getValue(col, row))
                return false;
        }
    }

    return true;
}

template<class T>
std::string voreen::EllpackMatrix<T>::toString() const {
    std::ostringstream stream;
    for (size_t row=0; row < numRows_; row++) {
        stream << "Row " << row << ": ";
        for (size_t colIndex=0; colIndex < numColsPerRow_; colIndex++) {
            float value = getValueByIndex(row, colIndex);
            if (value != static_cast<T>(0))
                stream << getColumn(row, colIndex) << "=" << value << " ";
        }
        stream << "\n";
    }
    return stream.str();
}

template<class T>
std::pair<size_t, size_t> voreen::EllpackMatrix<T>::rowIndexRange(size_t row) const {
    size_t start = row*numColsPerRow_;
    size_t end = start + numColsPerRow_;
    return std::make_pair(start, end);
}

template<class T>
size_t voreen::EllpackMatrix<T>::internalIndex(size_t row, size_t col) const {
#ifdef VRN_BLAS_DEBUG
    tgtAssert(row < numRows_ && col < numColsPerRow_, "Invalid indices");
#endif
    return row*numColsPerRow_ + col;
}

template<class T>
size_t voreen::EllpackMatrix<T>::getNextFreeIndex(size_t row) const {
#ifdef VRN_BLAS_DEBUG
    tgtAssert(M_ && indices_, "Data buffers not initialized");
    tgtAssert(row < numRows_ , "Invalid row");
#endif
    auto range = rowIndexRange(row);
    for (size_t index = range.first; index != range.second; ++index) {
        if (indices_[index] == -1)
            return index;
    }
    return -1;
}

#endif
