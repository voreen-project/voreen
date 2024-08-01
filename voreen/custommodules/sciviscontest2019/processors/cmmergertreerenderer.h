/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_CMMERGERTREERENDERER_H
#define VRN_CMMERGERTREERENDERER_H


#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/renderport.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/colorproperty.h"

#include "tgt/vector.h"
#include "tgt/texture.h"

#include "../ports/cmhaloport.h"
#include "../interaction/camera2dinteractionhandler.h"
#include "../interaction/cmselectionmanager.h"
#include "../interaction/cmpropertyanimator.h"

//use namespace voreen
namespace voreen {

/**
 * A RenderProcessor that displays halos in an abstract fashion. Halos from
 * all time steps related to the currently selected halo will be displayed
 * as a merger tree (see literature).
 */
class VRN_CORE_API CMMergerTreeRenderer : public RenderProcessor {

public:
    /**
     * Constructor
     */
    CMMergerTreeRenderer();
    ~CMMergerTreeRenderer();

    //------------------------------------------
    //  Pure virtual functions of base classes
    //------------------------------------------
    virtual Processor* create() const { return new CMMergerTreeRenderer();     }
    virtual std::string getClassName() const { return "CMMergerTreeRenderer";  }
    virtual std::string getCategory() const  { return "Viscontest2019";        }

protected:

    virtual void setDescriptions() { setDescription("Render processor that abstractly displays evolution of halos using merger trees."); }
    virtual void process();

    /**
     * Overwrites the base implementation of this function.
     * It is used to load the needed shader.
     * @see Processor
     */
    virtual void initialize();

    /**
     * Overwrites the base implementation of this function.
     * It is used to free the used shader.
     * @see Processor
     */
    virtual void deinitialize();

private:
    struct VertexLayout{
        tgt::vec3 pos;
        GLuint flags;
        GLuint haloID;
        VertexLayout(tgt::vec3 p, GLuint f, GLuint id)
            : pos(p)
            , flags(f)
            , haloID(id)
        { }
    };

    /**
     * Creates and sets up opengl buffers for use.
     */
    void setupBuffers();

    /**
     * Deletes opengl buffers and frees ressources.
     */
    void deleteBuffers();

    /**
     * Renders the links between halos of the merger tree
     */
    void renderLinks(tgt::mat4& projectionMatrix, tgt::mat4& viewMatrix);

    /**
     * Renders the vertices of the merger tree, i.e. the halos themselves.
     */
    void renderVerts(tgt::mat4& projectionMatrix, tgt::mat4& viewMatrix, tgt::Shader* shader);

    /**
     * Renders the ids for all vertices of the merger tree, i.e. the halos themselves.
     */
    void renderVertIDs(tgt::mat4& projectionMatrix, tgt::mat4& viewMatrix, tgt::Shader* shader);

    /**
     * Move the camera to the specified destination.
     */
    void moveCameraTo(const tgt::vec2& dest);

    /**
     * Process the given key event and take appropriate actions.
     */
    void processKeyEvent(tgt::KeyEvent*);

    /**
     * Used to recompile shaders etc. if the property changes.
     */
    void useTexturesChanged();

    /**
     * Used to update interal datastructures if the time step property changes.
     */
    void timeStepChanged();

    /**
     * Used to update interal datastructures if the selected halo changes.
     */
    void selectedHaloChanged();

    /**
     * Get a pointer to the currently selected halo.
     */
    const CMHalo* getSelectedHalo();

    /**
     * Algorithm that recursively builds the merger tree.
     * @param currentHalo halo whose ancestors still have to be added
     * @param path iterator from leaf to root of the merger tree.
     *          *path is of the same generation as *currentHalo (or == *pathEnd)
     * @param pathEnd end of the path from leaf to root of the merger tree
     * @param vertlinks output parameter that holds the indices for vertex links
     *          Indices are referring to currentVerts_
     * @param yspace y-position the current halo will be placed at
     * @param xpos x-position the current halo will be placed at
     * @param descVertPos index of the descendant of currentHalo in currentVerts_
     * @return y-position currentHalo was placed at.
     */
    float addAncestors(const CMHalo* const currentHalo, std::deque<const CMHalo*>::const_reverse_iterator& path, const std::deque<const CMHalo*>::const_reverse_iterator& pathEnd, std::vector<GLushort>& vertlinks, float& yspace, float xpos, GLushort descVertPos=0);

    /**
     * Calculate and update the halo path if the time step or halo changes.
     */
    void calculateHaloPath();

    //-------------
    //  members
    //-------------

    /// Output of the rendered image
    RenderPort outport_;
    /// Inport for halo data
    CMHaloPort inport_;
    /// Property for the radius of the halos in rendering
    FloatProperty radiusProp_;
    /// Property for width of links between halos in rendering
    FloatProperty lineWidthProp_;
    /// Property for distance between halos in rendering
    FloatProperty spacingProp_;
    /// Shader for halo rendering
    ShaderProperty vertShaderProp_;
    /// Shader for selection rendering
    ShaderProperty vertSelectionShaderProp_;
    /// Shader for link rendering
    ShaderProperty linkShaderProp_;
    /// ID of the currently selected halo
    IntProperty selectedHaloIDProp_;
    /// ID of a hovered over halo
    IntProperty mouseOverHaloIDProp_;
    /// Current time step
    FloatProperty timeStep_;

    /// Used for handling key events
    EventProperty<CMMergerTreeRenderer>* keyPressEvent_;
    /// Handler for 2d interaction
    Camera2DInteractionHandler* cameraHandler_;
    /// Position of the camera (Center of the screen)
    FloatVec2Property center_;
    /// Current width of the window (in merger tree coordinates => used to specify "zoom level")
    FloatProperty windowWidth_;
    /// Duration of animations (in ms) on changed selectedHaloIDProp_
    IntProperty animationDurationProp_;
    /// Width of the last rendered image
    size_t lastOutportWidth_;
    /// Halos that form the selected path
    std::deque<const CMHalo*> selectedPath_;
    /// Vertex positions for the currently rendered merger tree
    std::vector<VertexLayout> currentVerts_;
    /// Texture used to render halos
    tgt::Texture* haloTexture_;
    /// Texture used to render links
    tgt::Texture* linkTexture_;
    /// Enable or disable texturing for the merger tree
    BoolProperty useTexturesProp_;
    /// Zoom level for the link texture
    FloatProperty linkTextureZoomLevel_;
    /// Color property for ordinary halos
    ColorProperty vertColorProp_;
    /// Color property for halos without parents
    ColorProperty newVertColorProp_;
    /// Color property for the selected halo
    ColorProperty selectedVertColorProp_;
    /// Color property for hovered over halos
    ColorProperty mouseOverVertColorProp_;
    /// Color property for halos on the selected path
    ColorProperty onPathVertColorProp_;
    /// Color property for ordinary links
    ColorProperty linkColorProp_;
    /// Color property for links on the selected path
    ColorProperty pathLinkColorProp_;

    /// Selection manager
    CMSelectionManager* selectionManager_;
    /// Property animator
    CMPropertyAnimator propertyAnimator_;

    /// vao for halo rendering
    GLuint vao_;
    /// vao for selection rendering
    GLuint selectionVao_;
    /// vbo that holds the halo positions
    GLuint vertVbo_;
    /// vbo that holds link indices (referring to vertVbo_)
    GLuint linkIbo_;

    /// Number of links (We do not store the links in RAM after generating them,
    ///         but need to know the size of linkIbo for rendering)
    size_t linkCount_;
    /// Lock to prevent cyclic updating of property
    bool targetSelectionLock_;
};

} // namespace

#endif // VRN_COSMOLOGYMERGERTREERENDERER_H
