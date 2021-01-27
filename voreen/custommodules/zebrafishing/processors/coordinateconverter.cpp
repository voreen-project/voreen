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

#include "coordinateconverter.h"

#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/datastructures/plotrow.h"
#include "modules/plotting/datastructures/plotcell.h"

#include "tgt/plane.h"
#include "tgt/line.h"

namespace voreen {

const std::string CoordinateConverter::loggerCat_("voreen.zebrafishing.CoordinateConverter");

CoordinateConverter::CoordinateConverter()
    : Processor()
    , inport_(Port::INPORT, "plotdata", "PlotData Input", false)
    , volumeInport_(Port::INPORT, "volumedata", "Volume Inport", false)
    , spineOutport_(Port::OUTPORT, "spine", "Spine (line)")
    , gutOutport_(Port::OUTPORT, "gut", "Gut (line)")
    , yolkOnsetOutport_(Port::OUTPORT, "yolkOnset", "Yolk Onset on gut (point)")
    , planeNormal_("planeNormal", "Plane Normal", tgt::vec3(0, 1, 0), tgt::vec3(-1), tgt::vec3(1))
    , planePosition_("planePosition", "Plane Position", 0.0f, -FLT_MAX, FLT_MAX, Processor::INVALID_RESULT, NumericProperty<float>::STATIC)
    , transformMatrix_("transformMatrix", "Transformation Matrix", tgt::mat4::identity, tgt::mat4(-1e10), tgt::mat4(1e10))
    //, minimumLength_("minLength", "Minimum Track Length (in steps)", 1, 1, 1, Processor::INVALID_RESULT, NumericProperty<int>::DYNAMIC)
    //, columnRange_("timeInterval", "Time step interval")
{
    addPort(inport_);
    addPort(volumeInport_);

    addPort(spineOutport_);
    addPort(gutOutport_);
    addPort(yolkOnsetOutport_);

    planeNormal_.setReadOnlyFlag(true);
    planePosition_.setReadOnlyFlag(true);
    addProperty(planeNormal_);
    addProperty(planePosition_);

    transformMatrix_.setReadOnlyFlag(true);
    addProperty(transformMatrix_);

    //addProperty(minimumLength_);
    //addProperty(columnRange_);
}

Processor* CoordinateConverter::create() const {
    return new CoordinateConverter();
}

void CoordinateConverter::process() {
    spineOutport_.clear();
    gutOutport_.clear();
    yolkOnsetOutport_.clear();

    planeNormal_.set(tgt::vec3(1.f, 0.f, 0.f));
    planePosition_.set(0.f);
    transformMatrix_.set(tgt::mat4::identity);

    const PlotData* data = dynamic_cast<const PlotData*>(inport_.getData());
    if (!data) {
        LERROR("No valid plot data");
        return;
    }

    // no output for empty plot data
    if (data->rowsEmpty())
        return;

    // if this can be not be vec3 data: abort
    if (data->getColumnCount() != 3) {
        LERROR("Plot data does not contain 3 coordinate dimensions!");
        return;
    }

    // read in all vec3 data points
    std::vector<tgt::vec3> dataPoints;

    for (int i = 0; i < data->getRowsCount(); ++i) {
        const PlotRowValue& row = data->getRow(i);
        const std::vector<PlotCellValue>& values = row.getCells();
        // we have to iterate over the coordinate dimensions
        tgt::vec3 currentValue = tgt::vec3::zero;
        for (size_t dim = 0; dim < 3; ++dim) {
            currentValue.elem[dim] = (float) values[dim].getValue();
        }
        dataPoints.push_back(currentValue);
    }

    // data layout:
    // 0: point on spine line
    // 1: direction of spine line
    //
    // 2: point on gut line
    // 3: direction of gut line
    //
    // 4: point on spine plane
    // 5: normal of spine plane
    //
    // 6: point on gut plane
    // 7: normal of gut plane
    //
    // 8: point on mean plane
    // 9: normal of mean plane
    //
    // 10: point of yolk onset
    
    tgt::mat4 voxelToWorldMatrix = volumeInport_.getData()->getVoxelToWorldMatrix();
    tgt::mat3 mInv;
    tgt::mat3 m = voxelToWorldMatrix.getRotationalPartMat3();
    m.invert(mInv);
    tgt::mat3 mInvTr = transpose(mInv);
    
    //vec3 newNorm = normalize((mInvTr * vec4(n, T(0))).xyz());
    
    // create output geometry for rendering the spine line
    tgt::vec3 spinePoint = (voxelToWorldMatrix * tgt::vec4(dataPoints.at(0), 1.f)).xyz();
    tgt::vec3 spineDirection = tgt::normalize(mInvTr * dataPoints.at(1)); 
    // find intersections with y-z planes -> x = 0 and x = volumeDim.x - 1 (transformed) 
    tgt::vec3 v_p0 = (voxelToWorldMatrix * tgt::vec4(tgt::vec3::zero, 1.f)).xyz();
    tgt::vec3 v_p1 = (voxelToWorldMatrix * tgt::vec4(volumeInport_.getData()->getDimensions().x - 1, 0.f, 0.f, 1.f)).xyz();
    tgt::vec3 n_p = normalize(mInvTr * tgt::vec3(1.f, 0.f, 0.f));
    float t0_spine = tgt::dot(n_p, v_p0 - spinePoint) / tgt::dot(n_p, spineDirection);
    float t1_spine = tgt::dot(n_p, v_p1 - spinePoint) / tgt::dot(n_p, spineDirection);

    PointListGeometryVec3* spineList = new PointListGeometryVec3();
    spineList->addPoint(spinePoint + t0_spine * spineDirection);
    spineList->addPoint(spinePoint + t1_spine * spineDirection);
    spineOutport_.setData(spineList, true);

    // create output geometry for rendering the gut line
    tgt::vec3 gutPoint = (voxelToWorldMatrix * tgt::vec4(dataPoints.at(2), 1.f)).xyz();
    tgt::vec3 gutDirection = tgt::normalize(mInvTr * dataPoints.at(3));
    float t0_gut = tgt::dot(n_p, v_p0 - gutPoint) / tgt::dot(n_p, gutDirection);
    float t1_gut = tgt::dot(n_p, v_p1 - gutPoint) / tgt::dot(n_p, gutDirection);

    PointListGeometryVec3* gutList = new PointListGeometryVec3();
    gutList->addPoint(gutPoint + t0_gut * gutDirection);
    gutList->addPoint(gutPoint + t1_gut * gutDirection);
    gutOutport_.setData(gutList, true);

    // set properties for the mean plane
    planeNormal_.set(normalize(mInvTr * normalize(dataPoints.at(9))));
    tgt::vec3 worldPoint = (voxelToWorldMatrix * tgt::vec4(0.5f * (dataPoints.at(0) + dataPoints.at(2)), 1.f)).xyz();

    float d = tgt::dot(worldPoint, planeNormal_.get());
    planePosition_.set(d);

    // create output geometry for rendering the yolk onset
    tgt::vec4 transformedPoint = voxelToWorldMatrix * tgt::vec4(dataPoints.at(10), 1.f);
    PointListGeometryVec3* pointList = new PointListGeometryVec3();
    pointList->addPoint(transformedPoint.xyz());
    yolkOnsetOutport_.setData(pointList, true);

    // compute the transformation to normalize orientation and position of the data set
    // 1. rotate the mean plane into the x-z plane
    tgt::vec3 eY = tgt::vec3(0.f, 1.f, 0.f);
    float cosY = tgt::dot(planeNormal_.get(), eY);
    tgt::mat4 firstRotation = tgt::mat4::identity;
    if (std::abs(cosY) > 0.0001) {
        // compute rotation axis
        tgt::vec3 axisY = tgt::cross(planeNormal_.get(),eY);
        firstRotation = tgt::mat4::createRotation(std::acos(cosY), axisY);
    }

    // 2. now rotate the gutLine into the x-axis
    tgt::mat3 firstRotation3 = firstRotation.getRotationalPartMat3();
    tgt::vec3 gutDirectionTransformed = firstRotation3 * gutDirection;
    // project into x-z plane
    gutDirectionTransformed.y = 0;
    tgt::vec3 eX = tgt::vec3(-1.f, 0.f, 0.f);
    tgt::mat4 secondRotation = tgt::mat4::identity;
    float cosX = tgt::dot(gutDirectionTransformed, eX);
    if (std::abs(cosX) > 0.0001) {
        // compute rotation axis
        tgt::vec3 axisX = tgt::cross(gutDirectionTransformed, eX);
        secondRotation = tgt::mat4::createRotation(std::acos(cosX), axisX);
    }

    // 3. compute translation
    tgt::vec4 translationVector = -1.f * (secondRotation * (firstRotation * transformedPoint));
    tgt::mat4 translation = tgt::mat4::createTranslation(translationVector.xyz());

    // 4. compute composite transformation
    tgt::mat4 compositeMatrix = translation * (secondRotation * firstRotation);

    transformMatrix_.set(compositeMatrix);
}

} // namespace
