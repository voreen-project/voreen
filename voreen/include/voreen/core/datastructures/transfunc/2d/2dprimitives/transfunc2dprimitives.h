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

#ifndef VRN_TRANSFUNC2DPRIMITIVES_H
#define VRN_TRANSFUNC2DPRIMITIVES_H

#include "voreen/core/datastructures/transfunc/2d/transfunc2d.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/utils/transfuncprimitive.h"
#include "tgt/texture.h"


namespace tgt {
    class FramebufferObject;
}

namespace voreen {

class TransFuncPrimitive;

/**
 * This class represents a two-dimensional intensity-gradient based transfer function.
 *
 * - x-axis: intensity
 * - y-axis: gradient magnitude
 *
 * The transfer function is defined by primitives that can be added and edited with TransFuncEditorIntensityGradient.
 * Internally, the transfer function is represented by a two-dimensional RGBA texture of type GL_FLOAT,
 * which is updated through a framebuffer object.
 */
class VRN_CORE_API TransFunc2DPrimitives : public TransFunc2D {

    friend class TransFunc2DPrimitivesPainter;

public:
    /**
     * Constructor
     *
     * @param width width of the texture of the transfer function
     * @param height height of the texture of the transfer function
     */
    TransFunc2DPrimitives(int width = 256, int height = 256);

    /**
     * Destructor
     *
     * - deletes all primitives the transfer function consists of
     * - deletes the framebuffer object
     */
    virtual ~TransFunc2DPrimitives();

    virtual std::string getClassName() const { return "TransFunc2DPrimitives";     }
    virtual TransFunc2DPrimitives* create() const    { return new TransFunc2DPrimitives(); }
    virtual TransFunc2DPrimitives* clone() const;
    virtual void setMemberValuesFrom(const TransFuncBase* transfunc);
    virtual bool compareTo(const TransFuncBase& tf) const;
    virtual std::string getShaderDefines(const std::string& defineName = "TF_SAMPLER_TYPE") const;
    virtual void reset();
    virtual bool isTextureValid() const;
    virtual void invalidateTexture();
    //--------------------------------------
    //  handle texture
    //--------------------------------------
protected:
    /** @override */
    virtual tgt::Vector4<GLfloat> getMappingForValueFloat(float x, float y);
    /** @override */
    virtual tgt::Vector4<GLubyte> getMappingForValueUByte(float x, float y);
private:
    /**
     * Creates the help texture.
     */
    void updateHelpTexture();

    std::string generateShaderHeader() const;

    bool helpTextureInvalid_;              ///< used to check, if the texture is up to date
    tgt::FramebufferObject* fbo_;        ///< used for rendering the primitives to the help texture
    tgt::Texture* helpTex_;              ///< the help texture
    tgt::Texture* pickingTex_;          ///< texture for picking primitives. Is updated with the help texture

    tgt::Shader* tfProgram_; // shader to render the primitives into a color texture and a picking texture

    //--------------------------------------
    //  handle primitives
    //--------------------------------------
public:
    /**
     * Adds the given primitive to the transfer function.
     *
     * @param p the primitive that is added to the transfer function
     */
    void addPrimitive(TransFuncPrimitive* p);

    /**
     * Removes the given primitive from the transfer function. The primitive is deleted as well.
     *
     * @param p the primitve that will be removed from transfer function
     */
    void removePrimitive(TransFuncPrimitive* p);

    /**
     * Clears the transfer function, i.e. all primitives of the transfer function are
     * deleted.
     */
    void clear();

    /**
     * Returns the i.th primitive of the transfer function or 0
     * if no such primitive exists.
     */
    const TransFuncPrimitive* getPrimitive(int i) const;
    TransFuncPrimitive* getPrimitive(int i);

    /**
     * Returns the primitive of the transfer function at the given position or 0
     * if no primitive exists at this postion.
     * @param pos The position in normalized coordinates
     */
    const TransFuncPrimitive* getPrimitive(tgt::vec2 pos) const;
    TransFuncPrimitive* getPrimitive(tgt::vec2 pos);


    /**
     * Returns the number of primitives in the member vector primitives_.
     */
    size_t getNumPrimitives() const;

    /**
     * Returns the primitive that is under the mouse cursor or 0 if there is none.
     *
     * @param pos position of the mouse cursor
     * @return primitive that is under the mouse cursor or 0 if there is none
     */
    TransFuncPrimitive* getPrimitiveForClickedControlPoint(const tgt::vec2& pos) const;
protected:
    /**
     * Calls paint for all primitives
     */
    void paint();

    /**
     * Paints all primitives for selection purposes.
     */
    void paintForSelection();

    /**
     * Paints all primitives for display in an editor. Control points are added for every primitive.
     * An outline is added to the selected primitive.
     */
    void paintInEditor();

private:
    std::vector<TransFuncPrimitive*> primitives_; ///< primitives the transfer function consists of

    //--------------------------------------
    //  load and save
    //--------------------------------------
public:
    virtual const std::vector<std::string> getLoadFileFormats() const;
    virtual const std::vector<std::string> getSaveFileFormats() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    /**
     * The central entry point for loading a gradient transfer function.
     * The file extension is extracted and based on that the apropriate
     * load function is called. If there is no extension, loading will be unsuccessful.
     *
     * Currently supported extensions include:
     * tfig
     *
     * @param filename the filename, which should be opened
     * @return true, if loading was succesful, false otherwise
     */
    bool load(const std::string& filename);

    /**
     * Saves the transfer function to a file. Any data in the file will be overwritten.
     * The supported extensions include:
     * tfig, png
     *
     * @param filename the name of the file the transfer function will be saved to
     * @return true, if the operation was successful, false otherwise
     */
    bool save(const std::string& filename) const;
protected:
    /**
     * Saves the transfer function to a XML file. Returns true if the operation
     * was successful and false otherwise.
     *
     * @param filename name of the file the transfer function is saved to
     * @return true when save was successful, false otherwise
     */
    bool saveTfig(const std::string& filename) const;

    /**
     * Saves the transfer function to an image. Returns true if the operation
     * was successful and false otherwise.
     *
     * @param filename name of the file the transfer function is saved to
     * @return true when save was successful, false otherwise
     */
    bool saveImage(const std::string& filename) const;

    /**
     * Loads a gradient transfer function from a XML-File with ending tfig.
     *
     * @param filename the name of the file that should be opened
     * @return true, if loading was successful, false otherwise
     */
    bool loadTfig(const std::string& filename);

private:
    static const std::string loggerCat_; ///< the logger category
};

///------------------------------------------------------------------------------------------------

/**
 * Metadata encapsulating TransFunc1DKeys.
 */
class VRN_CORE_API TransFunc2DPrimitivesMetaData : public TransFuncMetaDataGeneric<TransFunc2DPrimitives> {
public:
    TransFunc2DPrimitivesMetaData();
    TransFunc2DPrimitivesMetaData(TransFunc2DPrimitives* transfunc);
    TransFunc2DPrimitivesMetaData(const std::vector<TransFunc2DPrimitives*>& transfunc);
    virtual MetaDataBase* clone() const;
    virtual MetaDataBase* create() const     { return new TransFunc2DPrimitivesMetaData(); }
    virtual std::string getClassName() const { return "TransFunc2DPrimitivesMetaData";     }
    virtual std::string toString() const {return "TransFunc2DPrimitivesMetaData";}
};

} // namespace voreen

#endif // VRN_TRANSFUNC2DPRIMITIVES_H
