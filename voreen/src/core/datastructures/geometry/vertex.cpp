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

#include "voreen/core/datastructures/geometry/vertex.h"

namespace voreen {

bool VertexBase::equals(const VertexBase& other, double epsilon) const {
    if (distance(pos_, other.pos_) > epsilon)
        return false;
    else
        return true;
}

void VertexBase::setupVertexAttributePointers(size_t stride) {
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, static_cast<GLsizei>(stride), 0);
}

void VertexBase::disableVertexAttributePointers() {
    glDisableVertexAttribArray(0);
}

double VertexBase::getDistanceToPlane(const tgt::plane& plane, double epsilon) const {
    double distance = plane.distance(pos_);
    if (std::abs(distance) <= epsilon)
        return 0;
    else
        return distance;
}

VertexBase VertexBase::interpolate(const VertexBase& v1, const VertexBase& v2, float t) {
    return VertexBase((v1.pos_ * (1.0f - t)) + (v2.pos_ * t));
}

//-------------------------------------------------------------------------------------------------

VertexNormal::VertexNormal(tgt::vec3 pos, tgt::vec3 normal)
    : VertexBase(pos)
    , normal_(normal)
{}

VertexNormal::VertexNormal(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex)
    : VertexBase(pos)
    , normal_(normal)
{}

bool VertexNormal::equals(const VertexNormal& other, double epsilon) const {
    if (distance(pos_, other.pos_) > epsilon)
        return false;

    if (distance(normal_, other.normal_) > epsilon)
        return false;

    return true;
}

void VertexNormal::setupVertexAttributePointers() {
    VertexBase::setupVertexAttributePointers(sizeof(VertexNormal));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(VertexNormal), (void*)sizeof(tgt::vec3));
}

void VertexNormal::disableVertexAttributePointers() {
    VertexBase::disableVertexAttributePointers();

    glDisableVertexAttribArray(2);
}

VertexNormal VertexNormal::interpolate(const VertexNormal& v1, const VertexNormal& v2, float t) {
    return VertexNormal((v1.pos_ * (1.0f - t)) + (v2.pos_ * t), (v1.normal_* (1.0f - t)) + (v2.normal_* t));
}

//-------------------------------------------------------------------------------------------------

VertexColor::VertexColor(tgt::vec3 pos, tgt::vec4 color)
    : VertexBase(pos)
    , color_(color)
{}

VertexColor::VertexColor(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex)
    : VertexBase(pos)
    , color_(color)
{}

bool VertexColor::equals(const VertexColor& other, double epsilon) const {
    if (distance(pos_, other.pos_) > epsilon)
        return false;

    if (distance(color_, other.color_) > epsilon)
        return false;

    return true;
}

void VertexColor::setupVertexAttributePointers() {
    VertexBase::setupVertexAttributePointers(sizeof(VertexColor));

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(VertexColor), (void*)sizeof(tgt::vec3));
}

void VertexColor::disableVertexAttributePointers() {
    VertexBase::disableVertexAttributePointers();

    glDisableVertexAttribArray(1);
}

VertexColor VertexColor::interpolate(const VertexColor& v1, const VertexColor& v2, float t) {
    return VertexColor((v1.pos_ * (1.0f - t)) + (v2.pos_ * t), (v1.color_ * (1.0f - t)) + (v2.color_ * t));
}

//-------------------------------------------------------------------------------------------------

VertexTexCoord::VertexTexCoord(tgt::vec3 pos, tgt::vec2 texCoord, tgt::ivec2 texIndex)
    : VertexBase(pos)
    , texCoord_(texCoord)
    , texIndex_(texIndex)
{}

VertexTexCoord::VertexTexCoord(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex)
    : VertexBase(pos)
    , texCoord_(texCoord)
    , texIndex_(texIndex)
{}

bool VertexTexCoord::equals(const VertexTexCoord& other, double epsilon) const {
    if (distance(pos_, other.pos_) > epsilon)
        return false;

    if (distance(texCoord_, other.texCoord_) > epsilon)
        return false;

    return true;
}

void VertexTexCoord::setupVertexAttributePointers() {
    VertexBase::setupVertexAttributePointers(sizeof(VertexTexCoord));

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, sizeof(VertexTexCoord), (void*)sizeof(tgt::vec3));
    glEnableVertexAttribArray(4);
    glVertexAttribIPointer(4, 2, GL_INT, sizeof(VertexTexCoord), (void*)(sizeof(tgt::vec3) + sizeof(tgt::vec2)));
}

void VertexTexCoord::disableVertexAttributePointers() {
    VertexBase::disableVertexAttributePointers();

    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);
}

VertexTexCoord VertexTexCoord::interpolate(const VertexTexCoord& v1, const VertexTexCoord& v2, float t) {
    return VertexTexCoord((v1.pos_ * (1.0f - t)) + (v2.pos_ * t), (v1.texCoord_ * (1.0f - t)) + (v2.texCoord_ * t), v1.texIndex_);
}

//-------------------------------------------------------------------------------------------------

VertexColorNormal::VertexColorNormal(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal)
    : VertexBase(pos)
    , color_(color)
    , normal_(normal)
{}

VertexColorNormal::VertexColorNormal(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex)
    : VertexBase(pos)
    , color_(color)
    , normal_(normal)
{}

bool VertexColorNormal::equals(const VertexColorNormal& other, double epsilon) const {
    if (distance(pos_, other.pos_) > epsilon)
        return false;

    if (distance(color_, other.color_) > epsilon)
        return false;

    if (distance(normal_, other.normal_) > epsilon)
        return false;

    return true;
}

void VertexColorNormal::setupVertexAttributePointers() {
    VertexBase::setupVertexAttributePointers(sizeof(VertexColorNormal));

    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(VertexColorNormal), (void*)sizeof(tgt::vec3));
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(VertexColorNormal), (void*)(sizeof(tgt::vec3) + sizeof(tgt::vec4)));
}

void VertexColorNormal::disableVertexAttributePointers() {
    VertexBase::disableVertexAttributePointers();

    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
}

VertexColorNormal VertexColorNormal::interpolate(const VertexColorNormal& v1, const VertexColorNormal& v2, float t) {
    return VertexColorNormal((v1.pos_ * (1.0f - t)) + (v2.pos_ * t), (v1.color_* (1.0f - t)) + (v2.color_* t), (v1.normal_* (1.0f - t)) + (v2.normal_* t));
    //return VertexColorNormal((v1.pos_ * (1.0f - t)) + (v2.pos_ * t), (v1.color_* (1.0f - t)) + (v2.color_* t), tgt::vec3(0.0f));
}

//-------------------------------------------------------------------------------------------------

VertexNormalTexCoord::VertexNormalTexCoord(tgt::vec3 pos, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex)
    : VertexBase(pos)
    , normal_(normal)
    , texCoord_(texCoord)
    , texIndex_(texIndex)
{}

VertexNormalTexCoord::VertexNormalTexCoord(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex)
    : VertexBase(pos)
    , normal_(normal)
    , texCoord_(texCoord)
    , texIndex_(texIndex)
{}

bool VertexNormalTexCoord::equals(const VertexNormalTexCoord& other, double epsilon) const {
    if (distance(pos_, other.pos_) > epsilon)
        return false;

    if (distance(normal_, other.normal_) > epsilon)
        return false;

    if (distance(texCoord_, other.texCoord_) > epsilon)
        return false;

    return true;
}

void VertexNormalTexCoord::setupVertexAttributePointers() {
    VertexBase::setupVertexAttributePointers(sizeof(VertexNormalTexCoord));

    glEnableVertexAttribArray(2);
    glEnableVertexAttribArray(3);
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(VertexNormalTexCoord), (void*)sizeof(tgt::vec3));
    glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, sizeof(VertexNormalTexCoord), (void*)(sizeof(tgt::vec3) + sizeof(tgt::vec3)));
    glVertexAttribIPointer(4, 2, GL_INT, sizeof(VertexNormalTexCoord), (void*)(sizeof(tgt::vec3) + sizeof(tgt::vec3) + sizeof(tgt::vec2)));
}

void VertexNormalTexCoord::disableVertexAttributePointers() {
    VertexBase::disableVertexAttributePointers();

    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);
}

VertexNormalTexCoord VertexNormalTexCoord::interpolate(const VertexNormalTexCoord& v1, const VertexNormalTexCoord& v2, float t) {
    return VertexNormalTexCoord((v1.pos_ * (1.0f - t)) + (v2.pos_ * t), (v1.normal_* (1.0f - t)) + (v2.normal_* t), (v1.texCoord_ * (1.0f - t)) + (v2.texCoord_* t), v1.texIndex_);
    //return VertexNormalTexCoord((v1.pos_ * (1.0f - t)) + (v2.pos_ * t), (v1.color_* (1.0f - t)) + (v2.color_* t), tgt::vec3(0.0f));
}

//-------------------------------------------------------------------------------------------------

VertexColorTexCoord::VertexColorTexCoord(tgt::vec3 pos, tgt::vec4 color, tgt::vec2 texCoord, tgt::ivec2 texIndex)
    : VertexBase(pos)
    , color_(color)
    , texCoord_(texCoord)
    , texIndex_(texIndex)
{}

VertexColorTexCoord::VertexColorTexCoord(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex)
    : VertexBase(pos)
    , color_(color)
    , texCoord_(texCoord)
    , texIndex_(texIndex)
{}

bool VertexColorTexCoord::equals(const VertexColorTexCoord& other, double epsilon) const {
    if (distance(pos_, other.pos_) > epsilon)
        return false;

    if (distance(color_, other.color_) > epsilon)
        return false;

    if (distance(texCoord_, other.texCoord_) > epsilon)
        return false;

    return true;
}

void VertexColorTexCoord::setupVertexAttributePointers() {
    VertexBase::setupVertexAttributePointers(sizeof(VertexColorTexCoord));

    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(3);
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(VertexColorTexCoord), (void*)sizeof(tgt::vec3));
    glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, sizeof(VertexColorTexCoord), (void*)(sizeof(tgt::vec3) + sizeof(tgt::vec4)));
    glVertexAttribIPointer(4, 2, GL_INT, sizeof(VertexColorTexCoord), (void*)(sizeof(tgt::vec3) + sizeof(tgt::vec4) + sizeof(tgt::vec2)));
}

void VertexColorTexCoord::disableVertexAttributePointers() {
    VertexBase::disableVertexAttributePointers();

    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);
}

VertexColorTexCoord VertexColorTexCoord::interpolate(const VertexColorTexCoord& v1, const VertexColorTexCoord& v2, float t) {
    return VertexColorTexCoord((v1.pos_ * (1.0f - t)) + (v2.pos_ * t), (v1.color_ * (1.0f - t)) + (v2.color_ * t), (v1.texCoord_ * (1.0f - t)) + (v2.texCoord_* t), v1.texIndex_);
    //return VertexColorTexCoord((v1.pos_ * (1.0f - t)) + (v2.pos_ * t), (v1.color_* (1.0f - t)) + (v2.color_* t), tgt::vec3(0.0f));
}

//-------------------------------------------------------------------------------------------------

VertexColorNormalTexCoord::VertexColorNormalTexCoord(tgt::vec3 pos, tgt::vec4 color, tgt::vec3 normal, tgt::vec2 texCoord, tgt::ivec2 texIndex)
    : VertexBase(pos)
    , color_(color)
    , normal_(normal)
    , texCoord_(texCoord)
    , texIndex_(texIndex)
{}

bool VertexColorNormalTexCoord::equals(const VertexColorNormalTexCoord& other, double epsilon) const {
    if (distance(pos_, other.pos_) > epsilon)
        return false;

    if (distance(color_, other.color_) > epsilon)
        return false;

    if (distance(normal_, other.normal_) > epsilon)
        return false;

    if (distance(texCoord_, other.texCoord_) > epsilon)
        return false;

    return true;
}

void VertexColorNormalTexCoord::setupVertexAttributePointers() {
    VertexBase::setupVertexAttributePointers(sizeof(VertexColorNormalTexCoord));

    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);
    glEnableVertexAttribArray(3);
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(VertexColorNormalTexCoord), (void*)sizeof(tgt::vec3));
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(VertexColorNormalTexCoord), (void*)(sizeof(tgt::vec3) + sizeof(tgt::vec4)));
    glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, sizeof(VertexColorNormalTexCoord), (void*)(sizeof(tgt::vec3) + sizeof(tgt::vec4) + sizeof(tgt::vec3)));
    glVertexAttribIPointer(4, 2, GL_INT, sizeof(VertexColorNormalTexCoord), (void*)(sizeof(tgt::vec3) + sizeof(tgt::vec4) + sizeof(tgt::vec3) + sizeof(tgt::vec2)));
}

void VertexColorNormalTexCoord::disableVertexAttributePointers() {
    VertexBase::disableVertexAttributePointers();

    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glDisableVertexAttribArray(4);
}

VertexColorNormalTexCoord VertexColorNormalTexCoord::interpolate(const VertexColorNormalTexCoord& v1, const VertexColorNormalTexCoord& v2, float t) {
    return VertexColorNormalTexCoord((v1.pos_ * (1.0f - t)) + (v2.pos_ * t), (v1.color_ * (1.0f - t)) + (v2.color_ * t), (v1.normal_* (1.0f - t)) + (v2.normal_* t), (v1.texCoord_ * (1.0f - t)) + (v2.texCoord_* t), v1.texIndex_);
    //return VertexColorNormalTexCoord((v1.pos_ * (1.0f - t)) + (v2.pos_ * t), (v1.color_* (1.0f - t)) + (v2.color_* t), tgt::vec3(0.0f));
}

} // namespace
