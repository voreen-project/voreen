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

#ifndef VRN_SIMILARITYMATRIX_H
#define VRN_SIMILARITYMATRIX_H

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/datastructures/datainvalidationobserver.h"
#include "voreen/core/io/serialization/serializable.h"

#include <map>

namespace voreen {

class EnsembleDataset;

/**
 * This class provides storage for and access to a symmetric similarity matrix (aka. distance matrix matrix).
 * Such a matrix is typically used to perform pairwise comparisons or, to apply an MDS approach.
 * @see SimilarityPlot
 */
class VRN_CORE_API SimilarityMatrix : public Serializable {
public:

    SimilarityMatrix(); // For deserialization only.

    /**
     * Constructs a similarity matrix where size equals number of columns/rows
     */
    explicit SimilarityMatrix(size_t size);

    SimilarityMatrix(const SimilarityMatrix& other);
    SimilarityMatrix(SimilarityMatrix&& other);

    SimilarityMatrix& operator=(const SimilarityMatrix& other);
    SimilarityMatrix& operator=(SimilarityMatrix&& other);

    /**
     * Returns the number of rows/columns (for a symmetric matrix, number of rows = number of columns).
     */
    size_t getSize() const;

    /**
     * Returns the entry at the specified location.
     * @note i and j must be less than size
     */
    float& operator() (size_t i, size_t j);
    float operator() (size_t i, size_t j) const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    size_t index(size_t i, size_t j) const;

    std::vector<float> data_; ///< Actual data storage.
    size_t size_;
};

/**
 * This class provides storage for and access to multiple similarity matrices that are created for an ensemble.
 * Each field (as referenced by its name) thereby is represented by a separate similarity matrix.
 * All contained matrices are guarenteed to have the same size.
 * @see SimilarityMatrix
 * @see EnsembleDataset
 * @see SimilarityMatrixCreator
 */
class VRN_CORE_API SimilarityMatrixList : public DataInvalidationObservable, public Serializable {
public:

    SimilarityMatrixList(); // For deserialization only.

    /**
     * Constructs a similarity matrix list for the specified ensemble.
     * For each field of the ensemble, a separate, uninitialized similarity matrix will be allocated.
     */
    explicit SimilarityMatrixList(const EnsembleDataset& dataset);

    SimilarityMatrixList(const SimilarityMatrixList& other);
    SimilarityMatrixList(SimilarityMatrixList&& other);

    /**
     * Returns the number of rows/columns of the contained similarity matrices.
     */
    size_t getSize() const;

    /**
     * Returns the ensemble hash for which this similarity matrix list was created.
     */
    const std::string& getHash() const;

    /**
     * Returns all field names for each of which a similarity matrix is stored.
     */
    std::vector<std::string> getFieldNames() const;

    /**
     * Returns the similarity matrix stored for the specified field name.
     * @note the field name must be contained in the vector returned by getFieldNames()
     */
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
