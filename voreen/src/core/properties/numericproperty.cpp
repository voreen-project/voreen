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

#include "voreen/core/properties/numericproperty.h"
#include <sstream>

namespace voreen {

//template<typename T>
//NumericProperty<T>::NumericProperty(const std::string& id, const std::string& guiText, const T& value,
//                                    const T& minValue, const T& maxValue, const T& stepping,
//                                    Processor::InvalidationLevel invalidationLevel)
//    : TemplateProperty<T>(id, guiText, value, invalidationLevel),
//    minValue_(minValue),
//    maxValue_(maxValue),
//    stepping_(stepping),
//    tracking_(true),
//    numDecimals_(2)
//{
//    addValidation(NumericPropertyValidation<T>(this));
//    std::string errorMsg;
//    validate(value, errorMsg);
//    if (!errorMsg.empty())
//        LWARNINGC("voreen.TemplateProperty", errorMsg);
//}
//
//template<typename T>
//void NumericProperty<T>::setMaxValue( const T& maxValue ) {
//    maxValue_ = maxValue;
//    this->updateWidgets();
//}
//
//template<typename T>
//const T& NumericProperty<T>::getMaxValue() const {
//    return maxValue_;
//}
//
//template<typename T>
//void NumericProperty<T>::setMinValue(const T& minValue) {
//    minValue_ = minValue;
//    this->updateWidgets();
//}
//
//template<typename T>
//const T& NumericProperty<T>::getMinValue() const {
//    return minValue_;
//}
//
//template<typename T>
//void NumericProperty<T>::setStepping(const T stepping) {
//    stepping_ = stepping;
//    this->updateWidgets();
//}
//
//template<typename T>
//T NumericProperty<T>::getStepping() const {
//    return stepping_;
//}
//
//template<typename T>
//int NumericProperty<T>::getNumDecimals() const {
//    return numDecimals_;
//}
//
//template<typename T>
//void NumericProperty<T>::setNumDecimals(int numDecimals) {
//    tgtAssert(numDecimals <= 64 && numDecimals >= 0, "Invalid number of decimals");
//    numDecimals_ = numDecimals;
//    this->updateWidgets();
//}
//
//template<typename T>
//bool NumericProperty<T>::hasTracking() const {
//    return tracking_;
//}
//
//template<typename T>
//void NumericProperty<T>::setTracking(bool tracking) {
//    tracking_ = tracking;
//    this->updateWidgets();
//}
//
//template<typename T>
//void NumericProperty<T>::increase() {
//    set(value_ + stepping_);
//}
//
//template<typename T>
//void NumericProperty<T>::decrease() {
//    set(value_ - stepping_);
//}
//
//template<typename T>
//void NumericProperty<T>::serialize(Serializer& s) const {
//    Property::serialize(s);
//
//    s.serialize("value", value_);
//
//    // serialize tracking mode, if it differs from default value
//    if (!tracking_)
//        s.serialize("tracking", tracking_);
//}
//
//template<typename T>
//void NumericProperty<T>::deserialize(Deserializer& s) {
//    Property::deserialize(s);
//
//    // deserialize value
//    T value;
//    s.deserialize("value", value);
//    try {
//        set(value);
//    }
//    catch (Condition::ValidationFailed& e) {
//        s.addError(e);
//    }
//
//    // deserialize tracking mode, if available
//    try {
//        s.deserialize("tracking", tracking_);
//    }
//    catch (SerializationNoSuchDataException&) {
//        s.removeLastError();
//    }
//}
//
//// explicit template instantiation to enable distribution of
//// implementation of template class methods over .h and .cpp files
//// and speed up compilation time.
////
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<int>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<float>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<double>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<tgt::vec2>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<tgt::vec3>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<tgt::vec4>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<tgt::ivec2>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<tgt::ivec3>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<tgt::ivec4>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<tgt::dvec2>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<tgt::dvec3>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<tgt::dvec4>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<tgt::mat2>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<tgt::mat3>;
//#endif
//#ifdef DLL_TEMPLATE_INST
//template class VRN_CORE_API NumericProperty<tgt::mat4>;
//#endif

} // namespace voreen
