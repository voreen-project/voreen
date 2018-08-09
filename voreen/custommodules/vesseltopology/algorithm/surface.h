/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
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
#pragma once
#include <string>
#include <vector>
#include <array>
#include <fstream>

namespace voreen {
struct StoredSurface {
    StoredSurface(std::string filename, size_t numVoxels)
        : filename_(filename)
        , numVoxels_(numVoxels)
    {
    }

    std::string filename_;
    size_t numVoxels_;
};

typedef std::vector<uint64_t> SurfaceSlice;

class SurfaceBuilder {
public:
    SurfaceBuilder();

    SurfaceBuilder(SurfaceBuilder&& other);

    StoredSurface finalize() &&;

    ~SurfaceBuilder() {}

    void push(uint64_t linearVoxelPos);
    void push_all(SurfaceSlice linearVoxelPositions);

private:
    std::string filename_;
    size_t numVoxelsStored_;
    std::ofstream file_;
};

class SurfaceReader {
public:
    SurfaceReader(StoredSurface surface);

    ~SurfaceReader();

    bool /* success or not */ read(uint64_t& val);

    size_t numVoxels() const;

private:
    StoredSurface surface_;
    std::ifstream file_;
};

template<int N>
struct SurfaceSlices {
    std::array<SurfaceSlice, N> slices_;

    void advance(SurfaceBuilder& builder) {
        builder.push_all(m<N-1>());
        m<N-1>().clear();

        for(int i=N-1; i>0; --i) {
            std::swap(slices_[i], slices_[i-1]);
        }
    }

    template<int i>
    SurfaceSlice& m() {
        static_assert(0 <= i && i < N, "Invalid index");
        return std::get<i>(slices_);
    }
};
}
