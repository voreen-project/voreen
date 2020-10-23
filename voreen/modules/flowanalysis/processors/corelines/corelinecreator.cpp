/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "corelinecreator.h"
#include "parallelvectors.h"

#include "tgt/vector.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

#include <chrono>

namespace voreen {
size_t findNeighbourSolution(size_t triangleIndex, const std::vector<int32_t> &triangleSolutionIndices, const std::vector<tgt::vec3> &solutions, std::vector<bool> &discardTriangle, tgt::svec3 dim) // returns the index of the neighbouring tet's triangle containing a solution or the maximum size_t value
{
    discardTriangle[triangleIndex] = true;

    const auto triangleIndexInTet = triangleIndex % ParallelVectors::TrianglesPerTetrahedron;

    // if there is a solution in the same tet, return it

    for (int triIndexInSameTet = 0; triIndexInSameTet < ParallelVectors::TrianglesPerTetrahedron; ++triIndexInSameTet)
    {
        const auto neighbourTriangleIndex = triangleIndex - triangleIndexInTet + triIndexInSameTet;
        if (!discardTriangle[neighbourTriangleIndex] && solutions[triangleSolutionIndices[neighbourTriangleIndex]] != solutions[triangleSolutionIndices[triangleIndex]])
        { // other solution in neighbouring tet
            return neighbourTriangleIndex;
        }
        discardTriangle[neighbourTriangleIndex] = true;
    }

    // otherwise, search for a solution in the neighbouring tet

    const auto cubeIndex = triangleIndex / (ParallelVectors::TrianglesPerTetrahedron * ParallelVectors::TetrahedraPerCube);

    auto x = cubeIndex % (dim.x - 1);
    auto y = (cubeIndex / (dim.x - 1)) % (dim.y - 1);
    auto z = (cubeIndex / ((dim.x - 1) * (dim.y - 1))) % (dim.z - 1);

    size_t neighbourTetIndexInCube;

    switch ((triangleIndex / ParallelVectors::TrianglesPerTetrahedron) % ParallelVectors::TetrahedraPerCube)
    {
    case 0: // front top left tet
    {
        switch (triangleIndexInTet)
        {
        case 0: // front
        {
            if (z == 0)
                return std::numeric_limits<size_t>::max();

            neighbourTetIndexInCube = 4;
            --z;
            break;
        }
        case 1: // back
        {
            neighbourTetIndexInCube = 5;
            break;
        }
        case 2: // right
        {
            neighbourTetIndexInCube = 1;
            break;
        }
        case 3: // top
        {
            if (y == dim.y - 2)
                return std::numeric_limits<size_t>::max();

            neighbourTetIndexInCube = 2;
            ++y;
            break;
        }
        }
        break;
    }
    case 1: // front bottom right tet
    {
        switch (triangleIndexInTet)
        {
        case 0: // front
        {
            if (z == 0)
                return std::numeric_limits<size_t>::max();

            neighbourTetIndexInCube = 3;
            --z;
            break;
        }
        case 1: // back
        {
            neighbourTetIndexInCube = 2;
            break;
        }
        case 2: // left
        {
            neighbourTetIndexInCube = 0;
            break;
        }
        case 3: // right
        {
            if (x == dim.x - 2)
                return std::numeric_limits<size_t>::max();

            neighbourTetIndexInCube = 5;
            ++x;
            break;
        }
        }
        break;
    }
    case 2: // middle right bottom tet
    {
        switch (triangleIndexInTet)
        {
        case 0: // bottom
        {
            if (y == 0)
                return std::numeric_limits<size_t>::max();

            neighbourTetIndexInCube = 0;
            --y;
            break;
        }
        case 1: // front
        {
            neighbourTetIndexInCube = 1;
            break;
        }
        case 2: // left
        {
            neighbourTetIndexInCube = 3;
            break;
        }
        case 3: // right
        {
            if (x == dim.x - 2)
                return std::numeric_limits<size_t>::max();

            neighbourTetIndexInCube = 4;
            ++x;
            break;
        }
        }
        break;
    }
    case 3: // back right bottom tet
    {
        switch (triangleIndexInTet)
        {
        case 0: // bottom
        {
            if (y == 0)
                return std::numeric_limits<size_t>::max();

            neighbourTetIndexInCube = 5;
            --y;
            break;
        }
        case 1: // left
        {
            neighbourTetIndexInCube = 4;
            break;
        }
        case 2: // right
        {
            neighbourTetIndexInCube = 2;
            break;
        }
        case 3: // back
        {
            if (z == dim.z - 2)
                return std::numeric_limits<size_t>::max();

            neighbourTetIndexInCube = 1;
            ++z;
            break;
        }
        }
        break;
    }
    case 4: // back left bottom tet
    {
        switch (triangleIndexInTet)
        {
        case 0: // left
        {
            if (x == 0)
                return std::numeric_limits<size_t>::max();

            neighbourTetIndexInCube = 2;
            --x;
            break;
        }
        case 1: // right
        {
            neighbourTetIndexInCube = 3;
            break;
        }
        case 2: // front
        {
            neighbourTetIndexInCube = 5;
            break;
        }
        case 3: // back
        {
            if (z == dim.z - 2)
                return std::numeric_limits<size_t>::max();

            neighbourTetIndexInCube = 0;
            ++z;
            break;
        }
        }
        break;
    }
    case 5: // middle left top tet
    {
        switch (triangleIndexInTet)
        {
        case 0: // left
        {
            if (x == 0)
                return std::numeric_limits<size_t>::max();

            neighbourTetIndexInCube = 1;
            --x;
            break;
        }
        case 1: // front
        {
            neighbourTetIndexInCube = 0;
            break;
        }
        case 2: // back
        {
            neighbourTetIndexInCube = 4;
            break;
        }
        case 3: // top
        {
            if (y == dim.y - 2)
                return std::numeric_limits<size_t>::max();

            neighbourTetIndexInCube = 3;
            ++y;
            break;
        }
        }
        break;
    }
    }

    const auto neighbourCubeIndex = (dim.x - 1) * ((dim.y - 1) * z + y) + x;
    const auto neighbourTetIndex = ParallelVectors::TetrahedraPerCube * neighbourCubeIndex + neighbourTetIndexInCube;

    size_t neighbourSolutionTriIndex = std::numeric_limits<size_t>::max();
    float dist = -1.0f;

    for (int triIndexInNeighbourTet = 0; triIndexInNeighbourTet < ParallelVectors::TrianglesPerTetrahedron; ++triIndexInNeighbourTet)
    {
        const auto neighbourTriangleIndex = ParallelVectors::TrianglesPerTetrahedron * neighbourTetIndex + triIndexInNeighbourTet;

        if (!discardTriangle[neighbourTriangleIndex])
        { // other solution in neighbouring tet
            auto d = tgt::lengthSq(solutions[triangleSolutionIndices[triangleIndex]] - solutions[triangleSolutionIndices[neighbourTriangleIndex]]);
            if (d > dist)
            {
                dist = d;
                neighbourSolutionTriIndex = neighbourTriangleIndex;
            }
        }
        discardTriangle[neighbourTriangleIndex] = true;
    }

    return neighbourSolutionTriIndex;
}

CorelineCreator::CorelineCreator()
    : Processor()
    , _in(Port::INPORT, "inport", "Parallel vectors solutions")
    , _out(Port::OUTPORT, "outport", "List of corelines")
    , _lengthThreshold("lengthThreshold", "Min. length of coreline", 20, 2, 1000, Processor::VALID)
{
    this->addPort(_in);
    this->addPort(_out);
    this->addProperty(_lengthThreshold);
    _lengthThreshold.onChange(MemberFunctionCallback<CorelineCreator>(this, &CorelineCreator::process));
}

void CorelineCreator::Process( const ParallelVectorSolutions& solutions, int lengthThreshold, std::vector<std::vector<tgt::vec3>>& corelines )
{
    const auto dim = solutions.dimensions;

    auto discardTriangle = std::vector<bool>(solutions.triangleSolutionIndices.size(), true);

    // filter such that only 2 solutions of a tet are used
    for (size_t tetIndex = 0; tetIndex < solutions.triangleSolutionIndices.size() / ParallelVectors::TrianglesPerTetrahedron; ++tetIndex)
    {
        int numSolutions = 0;
        for (int triIndexInTet = 0; triIndexInTet < ParallelVectors::TrianglesPerTetrahedron; ++triIndexInTet)
        {
            if (solutions.triangleSolutionIndices[tetIndex * ParallelVectors::TrianglesPerTetrahedron + triIndexInTet] >= 0)
            {
                ++numSolutions;
                discardTriangle[tetIndex * ParallelVectors::TrianglesPerTetrahedron + triIndexInTet] = false;
            }
        }

        if (numSolutions != 2)
        {
            for (int triIndexInTet = 0; triIndexInTet < ParallelVectors::TrianglesPerTetrahedron; ++triIndexInTet)
                discardTriangle[tetIndex * ParallelVectors::TrianglesPerTetrahedron + triIndexInTet] = true;
        }
    }

    // connect solutions
    for (size_t triangleIndex = 0; triangleIndex < solutions.triangleSolutionIndices.size(); ++triangleIndex)
    {
        if (discardTriangle[triangleIndex])
            continue;

        auto coreline = std::vector<tgt::vec3>{solutions.solutions[solutions.triangleSolutionIndices[triangleIndex]]};

        size_t neighbouringSolutionTriangleIndex = triangleIndex;
        while ((neighbouringSolutionTriangleIndex = findNeighbourSolution(neighbouringSolutionTriangleIndex, solutions.triangleSolutionIndices, solutions.solutions, discardTriangle, dim)) != std::numeric_limits<size_t>::max())
        {
            coreline.push_back(solutions.solutions[solutions.triangleSolutionIndices[neighbouringSolutionTriangleIndex]]);
        }

        for (size_t i = 0; i < coreline.size() / 2; ++i)
        {
            std::swap(coreline[i], coreline[coreline.size() - 1 - i]);
        }

        neighbouringSolutionTriangleIndex = triangleIndex;
        while ((neighbouringSolutionTriangleIndex = findNeighbourSolution(neighbouringSolutionTriangleIndex, solutions.triangleSolutionIndices, solutions.solutions, discardTriangle, dim)) != std::numeric_limits<size_t>::max())
        {
            coreline.push_back(solutions.solutions[solutions.triangleSolutionIndices[neighbouringSolutionTriangleIndex]]);
        }

        if (coreline.size() >= lengthThreshold)
        {
            coreline.shrink_to_fit();
            corelines.push_back(std::move(coreline));
        }
    }
}

void CorelineCreator::process()
{
    if (!_in.hasData())
        return;


    auto corelines = std::vector<std::vector<tgt::vec3>>();
    CorelineCreator::Process( *_in.getData(), _lengthThreshold.get(), corelines);

    auto geometry = new PointSegmentListGeometryVec3();
    geometry->setData(std::move(corelines));
    _out.setData(geometry);
}

} // namespace voreen
