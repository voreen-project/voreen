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

#include "voreen/core/properties/link/linkevaluatortransferfunctionproperties.h"
#include "voreen/core/properties/transfunc/1d/transfunc1dproperty.h"
#include "voreen/core/properties/floatproperty.h"

namespace voreen {

void LinkEvaluatorTransFunc1DPropertyBase::eval(Property* src, Property* dst) {
    if (dynamic_cast<FloatProperty*>(src)){
        FloatProperty* floatProp = dynamic_cast<FloatProperty*>(src);
        TransFunc1DProperty* tfProp = dynamic_cast<TransFunc1DProperty*>(dst);
        tgtAssert(floatProp && tfProp, "Invalid properties to link");

        TransFunc1D* tf = tfProp->get();
        tgtAssert(tf, "No value in tf property");
        setTFProperty(*tf, floatProp->get());
    } else {
        FloatProperty* floatProp = dynamic_cast<FloatProperty*>(dst);
        TransFunc1DProperty* tfProp = dynamic_cast<TransFunc1DProperty*>(src);
        tgtAssert(floatProp && tfProp, "Invalid properties to link");

        TransFunc1D* tf = tfProp->get();
        tgtAssert(tf, "No value in tf property");
        floatProp->set(getTFProperty(*tf));
    }
}

bool LinkEvaluatorTransFunc1DPropertyBase::arePropertiesLinkable(const Property* src, const Property* dst) const
{
    return dynamic_cast<const TransFunc1DProperty*>(src) && dynamic_cast<const FloatProperty*>(dst)
        || dynamic_cast<const TransFunc1DProperty*>(dst) && dynamic_cast<const FloatProperty*>(src);
}

float LinkEvaluatorTransFunc1DThresholdMin::getTFProperty(TransFunc1D& tf) {
    tgt::vec2 threshold = tf.getThreshold();
    return tf.normalizedToRealWorld(threshold.x);
}
void LinkEvaluatorTransFunc1DThresholdMin::setTFProperty(TransFunc1D& tf, float value) {
    tgt::vec2 threshold = tf.getThreshold();
    threshold.x = tf.realWorldToNormalized(value);
    tf.setThreshold(threshold);
}

float LinkEvaluatorTransFunc1DThresholdMax::getTFProperty(TransFunc1D& tf) {
    tgt::vec2 threshold = tf.getThreshold();
    return tf.normalizedToRealWorld(threshold.y);
}
void LinkEvaluatorTransFunc1DThresholdMax::setTFProperty(TransFunc1D& tf, float value) {
    tgt::vec2 threshold = tf.getThreshold();
    threshold.y = tf.realWorldToNormalized(value);
    tf.setThreshold(threshold);
}

float LinkEvaluatorTransFunc1DDomainMin::getTFProperty(TransFunc1D& tf) {
    tgt::vec2 domain = tf.getDomain();
    return domain.x;
}
void LinkEvaluatorTransFunc1DDomainMin::setTFProperty(TransFunc1D& tf, float value) {
    tgt::vec2 domain = tf.getDomain();
    domain.x = value;
    tf.setDomain(domain);
}

float LinkEvaluatorTransFunc1DDomainMax::getTFProperty(TransFunc1D& tf) {
    tgt::vec2 domain = tf.getDomain();
    return domain.y;
}
void LinkEvaluatorTransFunc1DDomainMax::setTFProperty(TransFunc1D& tf, float value) {
    tgt::vec2 domain = tf.getDomain();
    domain.y = value;
    tf.setDomain(domain);
}

float LinkEvaluatorTransFunc1DGamma::getTFProperty(TransFunc1D& tf) {
    return tf.getGammaValue();
}
void LinkEvaluatorTransFunc1DGamma::setTFProperty(TransFunc1D& tf, float value) {
    tf.setGammaValue(value);
}

}
