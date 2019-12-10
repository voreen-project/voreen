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

#ifndef VRN_BOUNDINGBOXPROPERTY_H
#define VRN_BOUNDINGBOXPROPERTY_H
#include "voreen/core/properties/property.h"
#include "voreen/core/properties/templateproperty.h"
#include "tgt/bounds.h"
#include "tgt/interval.h"
#include <cfloat>
#include <limits>

namespace voreen {


template<typename T>
class TemplateBoundingBoxProperty :public TemplateProperty<tgt::TemplateBounds<T> >{
public:
    TemplateBoundingBoxProperty<T>(const std::string& id, const std::string& guiText,  const tgt::TemplateBounds<T>& value=tgt::TemplateBounds<T>(tgt::Vector3<T>(0.0f), tgt::Vector3<T>(0.0f)),
                  const typename tgt::Vector3<T>& minimum = tgt::Vector3<T>(0.0f), const typename tgt::Vector3<T>& maximum = tgt::Vector3<T>(100.0f),
                  const typename tgt::Vector3<T>& minRange = tgt::Vector3<T>(0.0f), const typename tgt::Vector3<T>& maxRange = tgt::Vector3<T>(std::numeric_limits<T>::max()),
                  int invalidationLevel=Processor::INVALID_RESULT,
                  Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    TemplateBoundingBoxProperty<T>();
    virtual ~TemplateBoundingBoxProperty<T>() {}

    void set(const tgt::TemplateBounds<T>& val);
    const tgt::TemplateBounds<T>& get() const;

    void reset();


    typename tgt::Vector3<T> getMinValue() const;
    typename tgt::Vector3<T> getMaxValue() const;
    typename tgt::Vector3<T> getMinRange() const;
    typename tgt::Vector3<T> getMaxRange() const;

    void setMinValue(typename tgt::Vector3<T> value);
    void setMaxValue(typename tgt::Vector3<T> value);
    void setMinRange(typename tgt::Vector3<T> value);
    void setMaxRange(typename tgt::Vector3<T> value);

    void setStepping(tgt::Vector3<T>);
    tgt::Vector3<T> getStepping() const;

    virtual std::string getTypeDescription() const { return "Bounds"; }

    /// @see Property::serialize
    virtual void serialize(Serializer& s) const;

    /// @see Property::deserialize
    virtual void deserialize(Deserializer& s);


protected:
    tgt::Interval<T> xInterval_;
    tgt::Interval<T> yInterval_;
    tgt::Interval<T> zInterval_;

    tgt::Interval<T> xDefaultInterval_;
    tgt::Interval<T> yDefaultInterval_;
    tgt::Interval<T> zDefaultInterval_;

    tgt::Vector3<T> stepping_;
};


class VRN_CORE_API FloatBoundingBoxProperty :public TemplateBoundingBoxProperty<float>{
public:
    FloatBoundingBoxProperty(const std::string& id, const std::string& guiText,  const tgt::Bounds& value=tgt::Bounds(tgt::vec3(0.0f), tgt::vec3(0.0f)),
                  const tgt::vec3& minimum = tgt::vec3(0.0f), const tgt::vec3& maximum = tgt::vec3(100.0f),
                  const tgt::vec3& minRange = tgt::vec3(0.0f), const tgt::vec3& maxRange = tgt::vec3(FLT_MAX),
                  int invalidationLevel=Processor::INVALID_RESULT,
                  Property::LevelOfDetail lod = Property::LOD_DEFAULT)
        :TemplateBoundingBoxProperty<float>(id, guiText, value, minimum, maximum, minRange, maxRange, invalidationLevel)
    {
        stepping_ = tgt::vec3(0.01f);
    }
    FloatBoundingBoxProperty()
        :TemplateBoundingBoxProperty<float>(){}
    virtual Property* create() const{ return new FloatBoundingBoxProperty();}
    virtual std::string getClassName() const       { return "FloatBoundingBoxProperty"; }
    virtual std::string getTypeDescription() const { return "Bounds"; }
};

class VRN_CORE_API IntBoundingBoxProperty :public TemplateBoundingBoxProperty<int>{
public:
    IntBoundingBoxProperty(const std::string& id, const std::string& guiText,  const tgt::IntBounds& value=tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(0)),
                  const tgt::ivec3& minimum = tgt::ivec3(0), const tgt::ivec3& maximum = tgt::ivec3(100),
                  const tgt::ivec3& minRange = tgt::ivec3(0), const tgt::ivec3& maxRange = tgt::ivec3(std::numeric_limits<int>::max()),
                  int invalidationLevel=Processor::INVALID_RESULT,
                  Property::LevelOfDetail lod = Property::LOD_DEFAULT)
        :TemplateBoundingBoxProperty<int>(id, guiText, value, minimum, maximum, minRange, maxRange, invalidationLevel)
    { }

    IntBoundingBoxProperty()
        :TemplateBoundingBoxProperty<int>(){}
    virtual Property* create() const{ return new IntBoundingBoxProperty();}
    virtual std::string getClassName() const       { return "IntBoundingBoxProperty"; }
    virtual std::string getTypeDescription() const { return "Bounds"; }
};



template<typename T>
TemplateBoundingBoxProperty<T>::TemplateBoundingBoxProperty(const std::string& id, const std::string& guiText,
                             const tgt::TemplateBounds<T>& value, const typename tgt::Vector3<T>& minValue,
                             const typename tgt::Vector3<T>& maxValue, const typename tgt::Vector3<T>& minRange,
                             const typename tgt::Vector3<T>& maxRange, int invalidationLevel, Property::LevelOfDetail lod)
    : TemplateProperty<tgt::TemplateBounds<T> >(id, guiText, value, invalidationLevel, lod)
    , xInterval_(tgt::Vector2<T>(value.getLLF().x, value.getURB().x), minValue.x, maxValue.x, minRange.x, maxRange.x)
    , yInterval_(tgt::Vector2<T>(value.getLLF().y, value.getURB().y), minValue.y, maxValue.y, minRange.y, maxRange.y)
    , zInterval_(tgt::Vector2<T>(value.getLLF().z, value.getURB().z), minValue.z, maxValue.z, minRange.z, maxRange.z)
    , stepping_(tgt::Vector3<T>::one)
{
    TemplateProperty<tgt::TemplateBounds<T> >::set(value);
    xDefaultInterval_ = xInterval_;
    yDefaultInterval_ = yInterval_;
    zDefaultInterval_ = zInterval_;
}

template<typename T>
TemplateBoundingBoxProperty<T>::TemplateBoundingBoxProperty()
    : TemplateProperty<tgt::TemplateBounds<T> > ()
{}


template<typename T>
void TemplateBoundingBoxProperty<T>::set(const tgt::TemplateBounds<T>& val){
    tgt::TemplateBounds<T> old = this->get();
    if (val == old) return;
    tgt::Vector2<T> x = tgt::Vector2<T>(val.getLLF().x, val.getURB().x);
    tgt::Vector2<T> y = tgt::Vector2<T>(val.getLLF().y, val.getURB().y);
    tgt::Vector2<T> z = tgt::Vector2<T>(val.getLLF().z, val.getURB().z);
    bool valid = xInterval_.isValid(x) && yInterval_.isValid(y) && zInterval_.isValid(z);
    if (valid){
        xInterval_.set(tgt::Vector2<T>(val.getLLF().x, val.getURB().x));
        yInterval_.set(tgt::Vector2<T>(val.getLLF().y, val.getURB().y));
        zInterval_.set(tgt::Vector2<T>(val.getLLF().z, val.getURB().z));

        // we want the same value in the template property
        TemplateProperty<tgt::TemplateBounds<T> >::set(val);

        // template property does this invalidate for us
        //this->invalidate();
    }else{
        LWARNINGC("voreen.TemplateBoundingBoxProperty", "Set invalid value.");
    }


}

template<typename T>
const tgt::TemplateBounds<T>& TemplateBoundingBoxProperty<T>::get() const{
    tgt::Vector3<T> llf(xInterval_.get().x, yInterval_.get().x, zInterval_.get().x);
    tgt::Vector3<T> urb(xInterval_.get().y, yInterval_.get().y, zInterval_.get().y);
    tgt::TemplateBounds<T> bounds(llf, urb);

    // HACK(apv): because of animation system we need to inherit from
    // TemplateProperty and implement it's in our case stupid interface

    // this should (in theory) do nothing, because the templateproperty should
    // always have the same value as we have
    // But to make sure, here we set it again
    const_cast<TemplateBoundingBoxProperty*>(this)->TemplateProperty<tgt::TemplateBounds<T> >::set(bounds);

    // because of the definition of templateproperty this method has to
    // return a refernce, so we use the value in the templateproperty,
    // that we don't have to use a variable on the stack (undefinied by c++std)
    const tgt::TemplateBounds<T>& ref = TemplateProperty<tgt::TemplateBounds<T> >::get();
    tgtAssert(ref == bounds, "BoundingBoxProperty/TemplateProperty incoherence");
    return ref;
}

template<typename T>
void TemplateBoundingBoxProperty<T>::reset(){
    xInterval_ = xDefaultInterval_;
    yInterval_ = yDefaultInterval_;
    zInterval_ = zDefaultInterval_;
    this->invalidate();
}

template<typename T>
typename tgt::Vector3<T> TemplateBoundingBoxProperty<T>::getMinValue() const{
    return tgt::Vector3<T>(xInterval_.getMinValue(), yInterval_.getMinValue(), zInterval_.getMinValue());
}

template<typename T>
typename tgt::Vector3<T> TemplateBoundingBoxProperty<T>::getMaxValue() const{
    return tgt::Vector3<T>(xInterval_.getMaxValue(), yInterval_.getMaxValue(), zInterval_.getMaxValue());
}

template<typename T>
typename tgt::Vector3<T> TemplateBoundingBoxProperty<T>::getMinRange() const{
    return tgt::Vector3<T>(xInterval_.getMinRange(), yInterval_.getMinRange(), zInterval_.getMinRange());
}

template<typename T>
typename tgt::Vector3<T> TemplateBoundingBoxProperty<T>::getMaxRange() const{
    return tgt::Vector3<T>(xInterval_.getMaxRange(), yInterval_.getMaxRange(), zInterval_.getMaxRange());
}

template<typename T>
void TemplateBoundingBoxProperty<T>::setMinValue(typename tgt::Vector3<T> value){
    if (value== getMinValue()) return;
    xInterval_.setMinValue(value.x);
    yInterval_.setMinValue(value.y);
    zInterval_.setMinValue(value.z);
    this->invalidate();
}

template<typename T>
void TemplateBoundingBoxProperty<T>::setMaxValue(typename tgt::Vector3<T> value){
    if (value== getMaxValue()) return;
    xInterval_.setMaxValue(value.x);
    yInterval_.setMaxValue(value.y);
    zInterval_.setMaxValue(value.z);
    this->invalidate();
}

template<typename T>
void TemplateBoundingBoxProperty<T>::setMinRange(typename tgt::Vector3<T> value){
    if (value== getMinRange()) return;
    xInterval_.setMinRange(value.x);
    yInterval_.setMinRange(value.y);
    zInterval_.setMinRange(value.z);
    this->invalidate();
}

template<typename T>
void TemplateBoundingBoxProperty<T>::setMaxRange(typename tgt::Vector3<T> value){
    if (value== getMaxRange()) return;
    xInterval_.setMaxRange(value.x);
    yInterval_.setMaxRange(value.y);
    zInterval_.setMaxRange(value.z);
    this->invalidate();
}

template<typename T>
void TemplateBoundingBoxProperty<T>::setStepping(tgt::Vector3<T> stepping){
    if (stepping_ != stepping){
        stepping_ = stepping;
        this->invalidate();
    }
}

template<typename T>
tgt::Vector3<T> TemplateBoundingBoxProperty<T>::getStepping() const{
    return stepping_;
}

template<typename T>
void TemplateBoundingBoxProperty<T>::serialize(Serializer& s) const{
    tgt::TemplateBounds<T> b = this->get();
    s.serialize("bounds", b);
    s.serialize("minValue", this->getMinValue());
    s.serialize("maxValue", this->getMaxValue());
    s.serialize("minRange", this->getMinRange());
    s.serialize("maxRange", this->getMaxRange());
    s.serialize("stepping", this->getStepping());
}

template<typename T>
void TemplateBoundingBoxProperty<T>::deserialize(Deserializer& s){
    typename tgt::Vector3<T> minValue, maxValue, minRange, maxRange, stepping;
    tgt::TemplateBounds<T> bounds;
    s.optionalDeserialize("minValue", minValue,tgt::Vector3<T>::zero);
    s.optionalDeserialize("maxValue", maxValue,tgt::Vector3<T>::one);
    s.optionalDeserialize("minRange", minRange,tgt::Vector3<T>::zero);
    s.optionalDeserialize("maxRange", maxRange,tgt::Vector3<T>(std::numeric_limits<T>::max()));
    s.optionalDeserialize("stepping", stepping,tgt::Vector3<T>::one);

    try{
        // try new serialization
        s.deserialize("bounds", bounds);
    }catch(...){
        s.removeLastError();
        // fallback to old serialization
        typename tgt::Vector3<T> llf;
        typename tgt::Vector3<T> urb;
        s.optionalDeserialize("llf", llf,minValue);
        s.optionalDeserialize("urb", urb,maxValue);
        bounds = tgt::TemplateBounds<T>(llf, urb);
    }

    this->setMinValue(minValue);
    this->setMaxValue(maxValue);
    this->setMinRange(minRange);
    this->setMaxRange(maxRange);
    this->setStepping(stepping);

    this->set(bounds);
}

} // namespace voreen

#endif // VRN_BOUNDINGBOXPROPERTY_H
