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

#ifndef VRN_SIMILARITYMATRIX_H
#define VRN_SIMILARITYMATRIX_H

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/datastructures/datainvalidationobserver.h"
#include "voreen/core/io/serialization/serializable.h"

#include <map>

namespace voreen {

class EnsembleDataset;

class VRN_CORE_API SimilarityMatrix : public Serializable {
public:

    SimilarityMatrix(); // For deserialization only.
    explicit SimilarityMatrix(size_t size);
    SimilarityMatrix(const SimilarityMatrix& other);
    SimilarityMatrix(SimilarityMatrix&& other);

    SimilarityMatrix& operator=(const SimilarityMatrix& other);
    SimilarityMatrix& operator=(SimilarityMatrix&& other);

    size_t getSize() const;
    float& operator() (size_t i, size_t j);
    float operator() (size_t i, size_t j) const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    inline size_t index(size_t i, size_t j) const {
        // return (i < j) ? ((i*i-3*i)/2+j) : ((j*j-3*j)/2+i); // TODO: implement symmetric case.
        return i*size_ + j;
    }

    std::vector<float> data_;
    size_t size_;
};


class VRN_CORE_API SimilarityMatrixList : public DataInvalidationObservable, public Serializable {
public:

    SimilarityMatrixList(); // For deserialization only.
    explicit SimilarityMatrixList(const EnsembleDataset& dataset);
    SimilarityMatrixList(const SimilarityMatrixList& other);
    SimilarityMatrixList(SimilarityMatrixList&& other);

    size_t getSize() const;
    const std::string& getHash() const;
    const std::vector<std::string> getFieldNames() const;
    SimilarityMatrix& getSimilarityMatrix(const std::string& fieldName);
    const SimilarityMatrix& getSimilarityMatrix(const std::string& fieldName) const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    std::map<std::string, SimilarityMatrix> matrices_;
    std::string ensembleHash_;
};

}   // namespace

#endif //VRN_SIMILARITYMATRIX_H
