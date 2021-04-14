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

#ifndef VRN_PARALLELCOORDINATESAXES_H
#define VRN_PARALLELCOORDINATESAXES_H

#include "tgt/bounds.h"
#include "tgt/tgt_gl.h"
#include "tgt/vector.h"

#include <vector>
#include <string>

namespace voreen {

class ParallelCoordinatesAxes {
public:

    ParallelCoordinatesAxes(std::string ensembleHash,
                            std::vector<std::string> members,
                            std::vector<std::pair<std::string, int>> fields,
                            std::vector<std::string> axesLabels,
                            std::vector<tgt::vec2> ranges,
                            std::vector<float> values,
                            tgt::Bounds bounds,
                            size_t numTimesteps,
                            size_t numSamples );
    ParallelCoordinatesAxes( const std::string& filepath );
    ~ParallelCoordinatesAxes();

    void serialize( const std::string& filepath, bool binary = false ) const;

    size_t members() const noexcept;
    size_t fields() const noexcept;
    size_t timesteps() const noexcept;
    size_t samples() const noexcept;

    const std::string& getEnsembleHash() const;

    const std::string& getMemberName( size_t i ) const;
    const std::string& getFieldName( size_t i ) const;
    int getChannel( size_t i) const;

    const std::vector<std::string>& getMemberNames() const noexcept;
    const std::vector<std::pair<std::string, int>>& getFields() const noexcept;
    const std::vector<std::string>& getAxesLabels() const noexcept;

    tgt::vec2 getRange( size_t field ) const;
    const std::vector<tgt::vec2>& getRanges() const noexcept;

    float getValue( size_t field, size_t sample, size_t timestep = 0, size_t member = 0 ) const;
    const std::vector<float>& getValues() const noexcept;

    const tgt::Bounds& getBounds() const;

    size_t getStrideMember() const noexcept;
    size_t getStrideTimestep() const noexcept;
    size_t memorySize() const noexcept;

    /**
     * Returns the vertex buffer name.
     * If none has been assigned yet, a new name will be assigned.
     * Hence, an OpenGL context is required to be active.
     */
    GLuint getVertexBuffer() const;

private:
    size_t timesteps_, samples_;
    std::string ensembleHash_;
    std::vector<std::string> members_;
    std::vector<std::pair<std::string, int>> fields_;
    std::vector<std::string> axesLabels_;
    std::vector<tgt::vec2> ranges_;
    std::vector<float> values_;
    tgt::Bounds bounds_;
    mutable GLuint vertexBuffer_;
};
}

#endif // VRN_PARALLELCOORDINATESAXES_H
