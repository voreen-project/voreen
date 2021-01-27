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

#include "voreen/core/datastructures/transfunc/2d/2dprimitives/transfunc2dprimitives.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/utils/transfuncprimitive.h"

#include "tgt/glmath.h"
#include "tgt/gpucapabilities.h"
#include "tgt/framebufferobject.h"
#include "tgt/matrixstack.h"
#include "tgt/shadermanager.h"

#include "tgt/glcontextmanager.h"

#ifdef VRN_MODULE_DEVIL
    #include <IL/il.h>
#endif

#include <fstream>


namespace voreen {

const std::string TransFunc2DPrimitives::loggerCat_("voreen.TransFunc2DPrimitives");

TransFunc2DPrimitives::TransFunc2DPrimitives(int width, int height)
    : TransFunc2D(width, height, TF_FLOAT, tgt::Texture::LINEAR)
    , helpTextureInvalid_(true)
    , fbo_(0)
    , helpTex_(0)
    , pickingTex_(0)
    , tfProgram_(0)
{
}

TransFunc2DPrimitives::~TransFunc2DPrimitives() {

    if (tfProgram_)
        ShdrMgr.dispose(tfProgram_);

    if (fbo_) {
        delete fbo_;
        fbo_ = 0;
        LGL_ERROR;
    }

    if (helpTex_) {
        delete helpTex_;
        helpTex_ = 0;
        LGL_ERROR;
    }

    if (pickingTex_) {
        delete pickingTex_;
        pickingTex_ = 0;
        LGL_ERROR;
    }

    std::vector<TransFuncPrimitive*>::iterator it;
    while (!primitives_.empty()) {
        it = primitives_.begin();
        delete *it;
        primitives_.erase(it);
    }
}

TransFunc2DPrimitives* TransFunc2DPrimitives::clone() const {
    TransFunc2DPrimitives* func = new TransFunc2DPrimitives();
    func->setMemberValuesFrom(this);
    return func;
}

void TransFunc2DPrimitives::setMemberValuesFrom(const TransFuncBase* transfunc) {
    tgtAssert(transfunc, "null pointer passed");

    TransFunc2D::setMemberValuesFrom(transfunc);

    const TransFunc2DPrimitives* tf2D = dynamic_cast<const TransFunc2DPrimitives*>(transfunc);
    if (!tf2D) {
        LWARNING("setMemberValuesFrom(): passed parameter is not of type TransFunc2DPrimitives");
        return;
    }

    primitives_.clear();
    std::vector<TransFuncPrimitive*>::const_iterator it;
    for(it = tf2D->primitives_.begin(); it!=tf2D->primitives_.end(); it++) {
        primitives_.push_back((*it)->clone());
    }

    helpTextureInvalid_ = true;
}

bool TransFunc2DPrimitives::compareTo(const TransFuncBase& tf) const {
    if(TransFunc2DPrimitives* tf2d = dynamic_cast<TransFunc2DPrimitives*>(const_cast<TransFuncBase*>(&tf))) {

        // textures available?
        if (!tex_ || !tf2d->tex_ || !tex_->getCpuTextureData() || !tf2d->tex_->getCpuTextureData())
            return false;

        // compare texture data

        return ((TransFunc2D::compareTo(tf)) &&
            (memcmp(tex_->getCpuTextureData(), tf2d->tex_->getCpuTextureData(), tex_->getSizeOnCPU()) == 0));
    }
    return false;
}

void TransFunc2DPrimitives::reset() {
    clear();
}

std::string TransFunc2DPrimitives::getShaderDefines(const std::string& defineName) const {
    return TransFunc2D::getShaderDefines(defineName) + "#define CLASSIFICATION_REQUIRES_GRADIENT\n";
}

void TransFunc2DPrimitives::invalidateTexture() {
    TransFunc2D::invalidateTexture();
    helpTextureInvalid_ = true;
}

bool TransFunc2DPrimitives::isTextureValid() const {
    return TransFunc2D::isTextureValid() && !helpTextureInvalid_;
}

    //--------------------------------------
    //  handle texture
    //--------------------------------------
tgt::Vector4<GLfloat> TransFunc2DPrimitives::getMappingForValueFloat(float x, float y) {
    if(helpTextureInvalid_ || helpTex_->getDimensions() != dimensions_)
        updateHelpTexture();
    //check if help texture has been updated correctly
    tgtAssert(helpTex_, "no help texture!");
    tgtAssert(!helpTextureInvalid_, "help texture still invalid!");
    //calculate correct position
    tgt::ivec2 texPos = tgt::clamp(tgt::ivec2(x*dimensions_.x,y*dimensions_.y),
                                   tgt::ivec2::zero, tgt::ivec2(dimensions_.x-1, dimensions_.y-1));
    return helpTex_->texel<tgt::Vector4<GLfloat> >(texPos); //always GL_RGBA
}

tgt::Vector4<GLubyte> TransFunc2DPrimitives::getMappingForValueUByte(float x, float y) {
    // not used. 2D tf is float by default
    return tgt::Vector4<GLubyte>(getMappingForValueFloat(x,y)*255.f);
}

void TransFunc2DPrimitives::updateHelpTexture() {

    // Context sensitive code.
    tgt::GLContextStateGuard guard;

    // load the shader program if necessary
    if (!tfProgram_) {
        tfProgram_ = ShdrMgr.loadSeparate("transfunc2dprimitives/transfunc2dprimitives.vert", "transfunc2dprimitives/transfunc2dprimitives.frag", generateShaderHeader(), false);
        tgtAssert(tfProgram_, "no transfunc2dprimitives shader program");
    }

    // create FBO
    if (!fbo_) {
        LDEBUG("Creating FBO...");
        fbo_ = new tgt::FramebufferObject();
        LGL_ERROR;
        if (!fbo_) {
            LERROR("Failed to initialize framebuffer object");
            delete fbo_;
            fbo_ = 0;
            return;
        }
    }
    fbo_->activate();
    fbo_->detachAll();

    //create help texture
    if(!helpTex_ || helpTex_->getDimensions() != dimensions_) {
        delete helpTex_;
        helpTex_ = new tgt::Texture(dimensions_, (GLint)GL_RGBA, (GLint)GL_RGBA, (GLenum)(dataType_ == TF_FLOAT ? GL_FLOAT : GL_UNSIGNED_BYTE), filter_);
        helpTex_->uploadTexture();
    }

    // create picking texture
    if (!pickingTex_ || pickingTex_->getDimensions() != dimensions_) {
        delete pickingTex_;
        pickingTex_ = new tgt::Texture(dimensions_, (GLint)GL_RED_INTEGER, (GLint) GL_R32I, (GLenum) GL_INT, tgt::Texture::NEAREST);
        pickingTex_->uploadTexture();
    }

    // attach textures to fbo

    fbo_->attachTexture(helpTex_);
    if (!fbo_->isComplete()) {
        LERROR("Invalid framebuffer object - could not attach help texture");
        fbo_->deactivate();
        return;
    }
    fbo_->attachTexture(pickingTex_, GL_COLOR_ATTACHMENT1);
    if (!fbo_->isComplete()) {
        LERROR("Invalid framebuffer object - could not attach picking texture");
        fbo_->deactivate();
        return;
    }

    GLenum* buffers = new GLenum[2];
    buffers[0] = GL_COLOR_ATTACHMENT0;
    buffers[1] = GL_COLOR_ATTACHMENT1;
    glDrawBuffers(2, buffers);
    // render primitives to fbo
    glViewport(0, 0, tex_->getWidth(), tex_->getHeight());

    // clear previous content
    //glClearColor(0.f,0.f,0.f,0.0f); //for debug
    glClear(GL_COLOR_BUFFER_BIT);

    // set correct projection and modelview matrices
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    MatStack.multMatrix(tgt::mat4::createOrtho(0.f, 1.f, 0.f, 1.f, -2.f, 1.f));


    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    // paint primitives
    paint();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();

    fbo_->deactivate();
    LGL_ERROR;

    helpTex_->downloadTexture();
    pickingTex_->downloadTexture();

    helpTextureInvalid_ = false;
}

std::string TransFunc2DPrimitives::generateShaderHeader() const {
    std::stringstream header;
    header << "#version 150 core\n";
    header << "#extension GL_ARB_explicit_attrib_location : enable\n";

    // matrix stack uniforms
    header << "uniform mat4 projectionMatrixStack_;\n";
    header << "uniform mat4 modelViewMatrixStack_;\n";
    header << "uniform mat4 modelViewProjectionMatrixStack_;\n";
    header << "uniform mat4 projectionMatrixInverseStack_;\n";
    header << "uniform mat4 modelViewMatrixInverseStack_;\n";
    header << "uniform mat4 modelViewProjectionMatrixInverseStack_;\n";

    return header.str();
}

    //--------------------------------------
    //  handle primitives
    //--------------------------------------
void TransFunc2DPrimitives::addPrimitive(TransFuncPrimitive* p) {
    primitives_.push_back(p);
    invalidateTexture();
}

void TransFunc2DPrimitives::removePrimitive(TransFuncPrimitive* p) {
    std::vector<TransFuncPrimitive*>::iterator it;
    for (it = primitives_.begin(); it != primitives_.end(); ++it) {
        if (*it == p) {
            delete *it;
            primitives_.erase(it);
            invalidateTexture();
            break;
        }
    }
}

void TransFunc2DPrimitives::clear() {
    std::vector<TransFuncPrimitive*>::iterator it;
    while (!primitives_.empty()) {
        it = primitives_.begin();
        delete *it;
        primitives_.erase(it);
    }

    invalidateTexture();
}

size_t TransFunc2DPrimitives::getNumPrimitives() const {
    return primitives_.size();
}

void TransFunc2DPrimitives::paint() {
    glEnable(GL_BLEND);
    //glBlendFuncSeparate(GL_ONE, GL_ONE, GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glBlendFunc(GL_ONE, GL_ONE);
    glBlendEquation(GL_MAX);
    glDisable(GL_DEPTH_TEST);

    tfProgram_->activate();

    for (size_t i = 0; i < primitives_.size(); ++i) {
        tfProgram_->setUniform("tfPrimitiveID_", static_cast<int>(i+1));
        primitives_[i]->paint();
    }

    tfProgram_->deactivate();

    glEnable(GL_DEPTH_TEST);
    glBlendFunc(GL_ONE,GL_ZERO);
    glBlendEquation(GL_FUNC_ADD);
    glDisable(GL_BLEND);
}

void TransFunc2DPrimitives::paintForSelection() {
    /*for (size_t i = 0; i < primitives_.size(); ++i)
        primitives_[i]->paintForSelection(static_cast<GLubyte>(i));*/
}

void TransFunc2DPrimitives::paintInEditor() {
    /*for (size_t i = 0; i < primitives_.size(); ++i)
        primitives_[i]->paintInEditor();*/
}

const TransFuncPrimitive* TransFunc2DPrimitives::getPrimitive(int i) const {
    if ((i < 0) || (i >= static_cast<int>(primitives_.size())))
        return 0;
    else
        return primitives_[i];
}

TransFuncPrimitive* TransFunc2DPrimitives::getPrimitive(int i) {
    if ((i < 0) || (i >= static_cast<int>(primitives_.size())))
        return 0;
    else
        return primitives_[i];
}

const TransFuncPrimitive* TransFunc2DPrimitives::getPrimitive(tgt::vec2 pos) const {
    return const_cast<TransFunc2DPrimitives*>(this)->getPrimitive(pos);
}

TransFuncPrimitive* TransFunc2DPrimitives::getPrimitive(tgt::vec2 pos) {
    if(helpTextureInvalid_ || pickingTex_->getDimensions() != dimensions_)
        updateHelpTexture();
    //check if picking texture has been updated correctly
    tgtAssert(pickingTex_, "no picking texture!");
    tgtAssert(!helpTextureInvalid_, "help texture still invalid!");
    //calculate correct position
    tgt::ivec2 texPos = tgt::clamp(tgt::ivec2(tgt::iround(pos.x * static_cast<float>(pickingTex_->getDimensions().x - 1)), tgt::iround(pos.y * static_cast<float>(pickingTex_->getDimensions().y - 1))),
                                   tgt::ivec2::zero, tgt::ivec2(dimensions_.x-1, dimensions_.y-1));

    int index = pickingTex_->texel<int>(texPos) - 1;

    //LERROR("Picked Index at " << pos << " " << texPos << ": " << index);
    if (index >= 0)
        return getPrimitive(index);
    else
        return 0;
}

TransFuncPrimitive* TransFunc2DPrimitives::getPrimitiveForClickedControlPoint(const tgt::vec2& pos) const {
    if (primitives_.empty())
        return 0;

    size_t min = 0;
    // A distance of 2 can never happen because the canvas is normalized to [0,1]x[0,1],
    // so this value is huge enough.
    float mindist = 2.f;
    float d;
    for (size_t i = 0; i < primitives_.size(); ++i) {
         d = primitives_[i]->getClosestControlPointDist(pos);
         if ((d < mindist) /*&& (d < primitives_[i]->getControlPointSize())*/) {
            mindist = d;
            min = i;
         }
    }
    if (mindist == 2.f)
        return 0;
    else
        return primitives_[min];
}
    //--------------------------------------
    //  load and save
    //--------------------------------------
const std::vector<std::string> TransFunc2DPrimitives::getLoadFileFormats() const {
    std::vector<std::string> res;
    res.push_back("tfig");
    return res;
}

const std::vector<std::string> TransFunc2DPrimitives::getSaveFileFormats() const {
    std::vector<std::string> res;
    res.push_back("tfig");
#ifdef VRN_MODULE_DEVIL
    res.push_back("png");
#endif
    return res;
}

void TransFunc2DPrimitives::serialize(Serializer& s) const {
    TransFunc2D::serialize(s);
    s.serialize("Primitives", primitives_, "Primitive");
}

void TransFunc2DPrimitives::deserialize(Deserializer& s) {
    TransFunc2D::deserialize(s);
    s.deserialize("Primitives", primitives_, "Primitive");
}

bool TransFunc2DPrimitives::save(const std::string& filename) const {
    //look for fileExtension
    std::string fileExtension;
    size_t dotPosition = filename.rfind(".");
    if (dotPosition == std::string::npos)
        return false;
    else
        fileExtension = filename.substr(dotPosition+1);

    if (fileExtension == "tfig")
        return saveTfig(filename);
    else
        return saveImage(filename);
}

bool TransFunc2DPrimitives::saveTfig(const std::string& filename) const {
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
        s.serialize("TransFunc2DPrimitives", this);

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

bool TransFunc2DPrimitives::saveImage(const std::string& filename) const {
#ifdef VRN_MODULE_DEVIL
    //extract file extension
    std::string fileExtension;
    size_t dotPosition = filename.rfind(".");
    fileExtension = filename.substr(dotPosition+1);

    //IL does _NOT_ overwrite files by default
    ilEnable(IL_FILE_OVERWRITE);
    //download texture and save as png:
    ILuint img;
    ilGenImages(1, &img);
    ilBindImage(img);

    GLubyte* im = new GLubyte[tex_->getWidth()*tex_->getHeight()*4];

    tex_->bind();
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, im);

    ilTexImage(tex_->getWidth(), tex_->getHeight(), 1, 4, IL_RGBA, IL_UNSIGNED_BYTE, im);

    if (fileExtension == "png")
        ilSave(IL_PNG, (ILstring)filename.c_str());
    else {
        ilDeleteImages(1, &img);
        delete[] im;
        return false;
    }

    ilDeleteImages(1, &img);
    delete[] im;
    return true;
#else
    LERROR("Saving of " << filename  << " failed: No DevIL support.");
    return false;
#endif // VRN_MODULE_DEVIL
}

bool TransFunc2DPrimitives::load(const std::string& filename) {
    // Extract the file extension
    std::string fileExtension;
    size_t dotPosition = filename.rfind(".");
    if (dotPosition != std::string::npos)
        // => the last (seperating) dot was found
        fileExtension = filename.substr(dotPosition+1);
    else
        return false;

    if (fileExtension == "tfig")
        return loadTfig(filename);
    else
        return false;
}

bool TransFunc2DPrimitives::loadTfig(const std::string& filename) {

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
        d.deserialize("TransFunc2DPrimitives", *this);
        stream.close();
    }
    catch (SerializationException &e) {
        LWARNING("SerializationException: " << e.what());
        stream.close();
        success = false;
    }

    // log result
    if (success)
        LINFO("Loaded transfer function from file: " << filename);
    else
        LWARNING("Loading transfer function failed.");

    return success;
}

//------------------------------------------------------------------
//------------------------------------------------------------------
//  Meta Data
//------------------------------------------------------------------
//------------------------------------------------------------------
TransFunc2DPrimitivesMetaData::TransFunc2DPrimitivesMetaData()
    : TransFuncMetaDataGeneric<TransFunc2DPrimitives>()
{}
TransFunc2DPrimitivesMetaData::TransFunc2DPrimitivesMetaData(TransFunc2DPrimitives* transfunc)
    : TransFuncMetaDataGeneric<TransFunc2DPrimitives>(transfunc)
{}
TransFunc2DPrimitivesMetaData::TransFunc2DPrimitivesMetaData(const std::vector<TransFunc2DPrimitives*>& transfunc)
    : TransFuncMetaDataGeneric<TransFunc2DPrimitives>(transfunc)
{}

MetaDataBase* TransFunc2DPrimitivesMetaData::clone() const {
    if (transFunc_.empty())
        return new TransFunc2DPrimitivesMetaData();
    else {
        std::vector<TransFunc2DPrimitives*> newFunc;
        for(size_t i = 0; i < transFunc_.size(); i++) {
            if(transFunc_[i])
                newFunc.push_back(static_cast<TransFunc2DPrimitives*>(transFunc_[i]->clone()));
            else
                newFunc.push_back(0);
        }
        return new TransFunc2DPrimitivesMetaData(newFunc);
    }
}

} // namespace voreen
