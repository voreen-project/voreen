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

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/transfunc/1d/1dgaussian/transfunc1dgaussian.h"
#include "voreen/core/datastructures/transfunc/1d/1dgaussian/utils/transfuncmappingcurve.h"
#include "voreen/core/datastructures/volume/histogramutils.h"
#include <math.h>

#ifdef VRN_MODULE_DEVIL
    #include <IL/il.h>
#endif

#include "tgt/logmanager.h"

#include <tinyxml/tinyxml.h>
#include <fstream>
#include <sstream>


namespace voreen {

// some stream helper functions
inline int readInt(std::ifstream& stream);
inline short readShort(std::ifstream& stream);
inline double readDouble(std::ifstream& stream);

const std::string TransFunc1DGaussian::loggerCat_("voreen.TransFunc1DGaussian");

TransFunc1DGaussian::TransFunc1DGaussian(int width)
    : TransFunc1D(width)
{
    setToStandardFunc();
}

TransFunc1DGaussian::TransFunc1DGaussian(const TransFunc1DGaussian& tf)
    : TransFunc1D(tf.dimensions_.x, tf.dataType_, tf.filter_)
{
    setMemberValuesFrom(&tf);
    //pre-integration table map is not copied
    preIntegrationTableMap_.setTransFunc(this);
}

TransFunc1DGaussian::~TransFunc1DGaussian() {
    for (size_t i = 0; i < curves_.size(); ++i)
        delete curves_[i];

    preIntegrationTableMap_.clear();

    if (preIntegrationProgram_)
        ShdrMgr.dispose(preIntegrationProgram_);
}


TransFunc1DGaussian* TransFunc1DGaussian::clone() const {
    TransFunc1DGaussian* func = new TransFunc1DGaussian();

    func->setMemberValuesFrom(this);

    return func;
}

void TransFunc1DGaussian::setMemberValuesFrom(const TransFuncBase* transfunc) {
    tgtAssert(transfunc, "null pointer passed");

    TransFunc1D::setMemberValuesFrom(transfunc);

    const TransFunc1DGaussian* tf1DGaussian = dynamic_cast<const TransFunc1DGaussian*>(transfunc);
    if (!tf1DGaussian) {
        LWARNING("updateFrom(): passed parameter is not of type TransFunc1DGaussian");
        return;
    }

    clearCurves();
    std::vector<TransFuncMappingCurve*>::const_iterator it;
    for (it = tf1DGaussian->curves_.begin(); it!=tf1DGaussian->curves_.end(); it++) {
        curves_.push_back((*it)->clone());
    }

    invalidateTexture();
}

bool TransFunc1DGaussian::compareTo(const TransFuncBase& tf) const {
    if(const TransFunc1DGaussian* tf1d = dynamic_cast<const TransFunc1DGaussian*>(&tf)) {
        bool res = TransFunc1D::compareTo(tf);
        if(!res || (curves_.size() != tf1d->curves_.size()))
            return false;
        for(size_t i = 0; i < curves_.size(); i++){
            if(curves_[i] != tf1d->curves_[i])
                return false;
        }
        return true;
    }
    return false;
}

    //--------------------------------------
    //  handle texture
    //--------------------------------------
tgt::Vector4<GLubyte> TransFunc1DGaussian::getMappingForValueUByte(float value) {

    // If there are no curves, any further calculation is meaningless
    if (curves_.empty())
        return tgt::Vector4<GLubyte>(0, 0, 0, 0);

    // temporary variable for the sum of every value derived from the different curves
    tgt::Vector4<float> color = tgt::Vector4<float>(0.f, 0.f, 0.f, 0.f);
    float alphaSum = 0.f;    //< the sum of the opacities of all added colors

    // calculate the color at position value for all curves
    std::vector<TransFuncMappingCurve*>::const_iterator curveIter = curves_.begin();
    while ((curveIter != curves_.end())) {

        TransFuncMappingCurve* curve = *curveIter;

        // deactivated curves will be ignored
        if (curve->isActive()) {
            tgt::Vector4<GLubyte> curveColor = curve->getColorAt(value);
            // Color values are mixed additively weighted by their opacity.
            // The highest alpha value of all curves is the result alpha.
            //
            // To avoid a completely black result color use a minimum
            // positive alpha value for the color weights.
            float alpha = curve->getOpacityAt(value);
            if (alpha < std::numeric_limits<float>::min())
                alpha = std::numeric_limits<float>::min();
            alphaSum += alpha;

            color.r += (curveColor.r * alpha);
            color.g += (curveColor.g * alpha);
            color.b += (curveColor.b * alpha);
            color.a = (curveColor.a > color.a ? curveColor.a : color.a);
        }

        curveIter++;
    }

    // normalise the color based on the sum of the opacities of all added colors
    // this counters darkening of colors with low opacities in the mixing process
    if (alphaSum > 1.f)
        alphaSum = 1.f;
    color.r = color.r / alphaSum;
    color.g = color.g / alphaSum;
    color.b = color.b / alphaSum;
    // by adding the values up for all curves, values outside of [0;255] are possible
    color = tgt::clamp<float>(color, 0.f, 255.f);

    // convert values from float to unsigned byte
    tgt::Vector4<uint8_t> ret;
    ret.r = static_cast<uint8_t>(color.r);
    ret.g = static_cast<uint8_t>(color.g);
    ret.b = static_cast<uint8_t>(color.b);
    ret.a = static_cast<uint8_t>(color.a);

    return ret;
}

tgt::Vector4<GLfloat> TransFunc1DGaussian::getMappingForValueFloat(float value) {
    // not used. 1D tf is ubyte by default
    return tgt::Vector4<GLfloat>(getMappingForValueUByte(value)) / 255.f;
}

    //--------------------------------------
    //  handle curves
    //--------------------------------------
int TransFunc1DGaussian::getNumCurves() const {
    return static_cast<int>(curves_.size());
}

const TransFuncMappingCurve* TransFunc1DGaussian::getCurve(int i) const {
    return curves_.at(i);
}

TransFuncMappingCurve* TransFunc1DGaussian::getCurve(int i) {
    return curves_.at(i);
}

const std::vector<TransFuncMappingCurve*> TransFunc1DGaussian::getCurves() const {
    return curves_;
}

void TransFunc1DGaussian::setCurves(std::vector<TransFuncMappingCurve*> curves) {
    std::vector<TransFuncMappingCurve *>::iterator curveIter = curves_.begin();
    // First delete all the referenced objects in the heap
    while (curveIter != curves_.end()) {
        delete (*curveIter);
        ++curveIter;
    }
    // then delete the entries in the vector
    curves_.clear();
    curves_ = curves;
    invalidateTexture();
}

void TransFunc1DGaussian::reset(){
    setToStandardFunc();
}

void TransFunc1DGaussian::sortCurves(TransFuncMappingCurve* selectedCurve) {

    if (!selectedCurve)
        return;

    std::sort(curves_.begin(), curves_.end(),
        [selectedCurve](TransFuncMappingCurve* a, TransFuncMappingCurve* b) -> bool
        {
            if (a->isActive() != b->isActive()) // activation state
                return b->isActive();
            if (a == selectedCurve)  // a selected
                return false;
            else if (b == selectedCurve) // b selected
                return true;
            else
                return false; // a,b are equal
        });
}

void TransFunc1DGaussian::addCurve(TransFuncMappingCurve* curve) {
    curves_.push_back(curve);
    invalidateTexture();
}

void TransFunc1DGaussian::updateCurve(TransFuncMappingCurve* /*curve*/) {
    // sort the curves by their activation state
    std::sort(curves_.begin(), curves_.end(),
        [](TransFuncMappingCurve* a, TransFuncMappingCurve* b) -> bool
        {
            if (a->isActive() != b->isActive())
                return b->isActive();
            else
                return false; // a,b are equal
        });
    invalidateTexture();
}

void TransFunc1DGaussian::removeCurve(TransFuncMappingCurve* curve) {
    std::vector<TransFuncMappingCurve *>::iterator curveIter = find(curves_.begin(), curves_.end(), curve);
    if (curveIter != curves_.end())
        curves_.erase(curveIter);
    delete curve;

    invalidateTexture();
}

void TransFunc1DGaussian::clearCurves() {

    // First delete all the referenced objects in the heap
    std::vector<TransFuncMappingCurve *>::iterator curveIter = curves_.begin();
    while (curveIter != curves_.end()) {
        delete (*curveIter);
        ++curveIter;
    }
    // then delete the entries in the vector
    curves_.clear();

    invalidateTexture();
}

bool TransFunc1DGaussian::isEmpty() const {
    return curves_.empty();
}

tgt::vec4 TransFunc1DGaussian::getMeanValue(float segStart, float segEnd) const {
    tgt::ivec4 result(0);
    if (!tex_ || textureInvalid_) {
        LWARNING("getMeanValue(): texture is invalid");
        return result;
    }

    float width = static_cast<float>(tex_->getWidth());
    for (int i = static_cast<int>(segStart*width); i < segEnd*width; ++i)
        result += tgt::ivec4(tex_->texel<tgt::col4>(i));

    return static_cast<tgt::vec4>(result)/(segEnd*width-segStart*width);
}

    //--------------------------------------
    //  function helper
    //--------------------------------------

void TransFunc1DGaussian::autoTreshold(Histogram1D* hist) {

    // obtain a fitting histogram for the current domain
    Histogram1D* subHisto = subHistogram(hist, domain_.x, domain_.y);

    // if the histogram hast <2 buckets with samples, no curves can be generated
    size_t nonEmptyBuckets = 0;
    for (size_t i = 0; i < subHisto->getNumBuckets(); i++)
        if (subHisto->getBucket(i) > 0) nonEmptyBuckets++;
    if (nonEmptyBuckets < 2){
        VoreenApplication::app()->showMessageBox(getGuiName(), "The histogram contains too few filled buckets for an Otsu classification. Maybe try a different domain?", true);
        return;
    }

    // classification using Otsu's method
    GaussianCurve* resultCurves = otsuTreshold(subHisto);

    // map the curves from the real world intensity subspace of the histogram to the normalised interval
    // (note that subHisto->getMinValue() == domain_.x and subHisto->getMaxValue() == domain_.y)
    GaussianCurve c0 = mapCurve(resultCurves[0], domain_, tgt::vec2(0.f, 1.f));
    GaussianCurve c1 = mapCurve(resultCurves[1], domain_, tgt::vec2(0.f, 1.f));

    float m0 = c0.mean;
    float m1 = c1.mean;
    // clamp variance values so that the markers are still accessible in the editor
    // (base markers are located at a distance of 3*variance away from the mean marker)
    float v0 = tgt::clamp(c0.variance, 0.0000001f, 0.33333333f);
    float v1 = tgt::clamp(c1.variance, 0.0000001f, 0.33333333f);;

    // means outside the normalised interval should not occur due to the usage of the sub-interval
    tgtAssert(m0 >= 0.f && m0 <= 1.f && m1 >= 0.f && m1 <= 1.f, "gaussian distributions outside the normalized intensity interval");

    // add two curves for the calculated gaussian distributions
    clearCurves();
    // the first curve is usually the 'background' of the dataset and will therefore be mapped to a low opacity
    TransFuncMappingCurve* curve = new TransFuncMappingCurve(m0, v0, 0.f, tgt::col4(255, 0, 0, 2));
    curves_.push_back(curve);
    // the second curve is usually the 'foreground' object
    curve = new TransFuncMappingCurve(m1, v1, 0.f, tgt::col4(0, 255, 0, 255));
    curve->setUnicolor(false);
    curve->setBaseColorL(tgt::col4(0, 0, 0, 255));
    curve->setBaseColorR(tgt::col4(0, 0, 0, 255));
    curves_.push_back(curve);

    // cleanup
    delete[] resultCurves;
    delete subHisto;

    invalidateTexture();
}

void TransFunc1DGaussian::setToStandardFunc() {
    clearCurves();
    // the standard function is a single white curve in the middle of the intensity space
    curves_.push_back(new TransFuncMappingCurve(0.5f, 0.01f, 0.f, tgt::col4(255)));
    setThreshold(0.f,1.f);
    setGammaValue(1.f);
    setAlphaMode(TF_USE_ALPHA);
    invalidateTexture();
}

    //--------------------------------------
    //  load and save
    //--------------------------------------
const std::vector<std::string> TransFunc1DGaussian::getLoadFileFormats() const {
    std::vector<std::string> res;
    res.push_back("tfg");
    return res;
}

const std::vector<std::string> TransFunc1DGaussian::getSaveFileFormats() const {
    std::vector<std::string> res;
    res.push_back("tfg");
    return res;
}

void TransFunc1DGaussian::serialize(Serializer& s) const {
    TransFunc1D::serialize(s);

    // serialize curves...
    s.serialize("Curves", curves_, "curve");
}

void TransFunc1DGaussian::deserialize(Deserializer& s) {
    TransFunc1D::deserialize(s);

    // deserialize curves...
    s.deserialize("Curves", curves_, "curve");
}

bool TransFunc1DGaussian::save(const std::string& filename) const {
    //look for fileExtension
    std::string fileExtension;
    size_t dotPosition = filename.rfind(".");
    if (dotPosition == std::string::npos)
        return false;
    else
        fileExtension = filename.substr(dotPosition+1);

    if (fileExtension == "tfg")
        return saveTfg(filename);
    else {
        LWARNING("File format '" << fileExtension << "' not allowed.");
        return false;
    }
}

bool TransFunc1DGaussian::saveTfg(const std::string& filename) const {

    // open file stream
    std::ofstream stream(filename.c_str(), std::ios_base::out);
    if (stream.fail()) {
        LWARNING("Unable to open file " << filename << " for writing.");
        return false;
    }

    // serialize to stream
    bool success = true;
    try {
        XmlSerializer s(filename);
        s.serialize("TransFuncGaussian", this);

        s.write(stream);
        if (stream.bad()) {
            LWARNING("Unable to write to file: " << filename);
            success = false;
        }
        stream.close();
    }
    catch (SerializationException &e) {
        LWARNING("SerializationException: " << e.what());
        stream.close();
        success = false;
    }

    // log result
    if (success)
        LINFO("Saved transfer function to file: " << filename);
    else
        LWARNING("Saving transfer function failed.");

    return success;
}



bool TransFunc1DGaussian::load(const std::string& filename) {
    // Extract the file extension
    std::string fileExtension;
    size_t dotPosition = filename.rfind(".");
    if (dotPosition != std::string::npos)
        // => the last (seperating) dot was found
        fileExtension = filename.substr(dotPosition+1);
    else
        return false;


    if (fileExtension == "tfg")
        return loadTfg(filename);
    else {
        LWARNING("File format '" << fileExtension << "' not allowed.");
        return false;
    }
}

bool TransFunc1DGaussian::loadTfg(const std::string& filename) {
    // open file stream
    std::ifstream stream(filename.c_str(), std::ios_base::in);
    if (stream.fail()) {
        LWARNING("Unable to open file " << filename << " for reading.");
        return false;
    }

    // deserialize from stream
    bool success = true;
    try {
        XmlDeserializer d(filename);
        d.read(stream);
        d.deserialize("TransFuncGaussian", *this);
        stream.close();
    }
    catch (SerializationException &e) {
        LWARNING("SerializationException: " << e.what());
        stream.close();
        success = false;
    }

    // log result
    if (success)
        LDEBUG("Loaded transfer function from file: " << filename);
    else
        LWARNING("Loading transfer function failed.");

    invalidateTexture();
    return success;
}


// --------------------------------------------------------------------------------
// some stream helper functions

/**
 * Used not only to extract neccessary information but also to fast-forward through the header
 * Buffersize is given static, because the procedure shouldn't be edited if 64bit compiler is used
 * (and therefor the sizes i.e. int, double might change)
 * These functions are 'inline', because it might be faster if they are inlined by the compiler
 * (okay, not quite _that_ convincing i admit).
 *
 * @param stream the already opened file stream used to extract the data
 * @return the read int value
 */
inline int readInt(std::ifstream& stream) {
    char* buffer = new char[4];

    // The bytes are inserted backwards because of a discrepancy of
    // most significant bit <-> least significant bit between Java and C++
    for (int i = 3; i >= 0; --i)
        stream >> buffer[i];

    int ret = *(reinterpret_cast<int*>(buffer));
    delete[] buffer;
    return ret;
}

/**
 * \sa TransFuncIntensity::readInt(std::ifstream& stream)
 */
inline short readShort(std::ifstream& stream) {
    char* buffer = new char[2];

    for (int i = 1; i >= 0; --i)
        stream >> buffer[i];

    short ret = *(reinterpret_cast<short*>(buffer));
    delete[] buffer;
    return ret;
}

/**
 * \sa TransFuncIntensity::readInt(std::ifstream& stream)
 */
inline double readDouble(std::ifstream& stream) {
    char* buffer = new char[8];

    for (int i = 7; i >= 0; --i)
        stream >> buffer[i];

    double ret = *(reinterpret_cast<double*>(buffer));
    delete[] buffer;
    return ret;
}

//------------------------------------------------------------------
//------------------------------------------------------------------
//  Meta Data
//------------------------------------------------------------------
//------------------------------------------------------------------
TransFunc1DGaussianMetaData::TransFunc1DGaussianMetaData()
    : TransFuncMetaDataGeneric<TransFunc1DGaussian>()
{}
TransFunc1DGaussianMetaData::TransFunc1DGaussianMetaData(TransFunc1DGaussian* transfunc)
    : TransFuncMetaDataGeneric<TransFunc1DGaussian>(transfunc)
{}
TransFunc1DGaussianMetaData::TransFunc1DGaussianMetaData(const std::vector<TransFunc1DGaussian*>& transfunc)
    : TransFuncMetaDataGeneric<TransFunc1DGaussian>(transfunc)
{}

MetaDataBase* TransFunc1DGaussianMetaData::clone() const {
    if (transFunc_.empty())
        return new TransFunc1DGaussianMetaData();
    else {
        std::vector<TransFunc1DGaussian*> newFunc;
        for(size_t i = 0; i < transFunc_.size(); i++) {
            if(transFunc_[i])
                newFunc.push_back(static_cast<TransFunc1DGaussian*>(transFunc_[i]->clone()));
            else
                newFunc.push_back(0);
        }
        return new TransFunc1DGaussianMetaData(newFunc);
    }
}

} // namespace voreen
