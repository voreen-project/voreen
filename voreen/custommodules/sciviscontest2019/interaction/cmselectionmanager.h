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

#ifndef VRN_CMSELECTIONMANAGER_H
#define VRN_CMSELECTIONMANAGER_H

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/interaction/interactionhandler.h"
#include "voreen/core/datastructures/rendertarget/rendertarget.h"
namespace voreen {

/**
 * Interaction handler that helps picking objects by rendering their IDs into a separate frame buffer.
 * The rendering itself has to be done by the user! Rendering of IDs should be done as follows:
 *
 * selectionManager_->activate();
 * selectionManager_->resize(outportSize);
 * selectionManager_->clear();
 *
 * RENDERING
 *
 * selectionManager_->deactivate();
 *
 */
class CMSelectionManager : public InteractionHandler {
    friend class EventProperty<CMSelectionManager>;
public:

    /**
     * Default constructor needed for serialization. Do not call it directly.
     */
    CMSelectionManager();

    /**
     * Constructor.
     * @param id Identifier that must be unique across all interaction handlers
     *  of a processor. Must not be empty.
     * @param guiName The string that is to be displayed in the GUI.
     * @param size Size for the render target IDs will be rendered to.
     * @param selectionIDProp Property that will hold the ID of the last element clicked.
     * @param mouseOverIDProp Property that will hold the ID of the last element hovered over.
     * @param noSelectionID ID that will be returned if nothing was rendered at a
     *          given position. Because of that IDs returned by idAt will be the
     *          id rendered + noSelectionID.
     */
    CMSelectionManager(const std::string& id, const std::string& guiName, tgt::ivec2 size, IntProperty* selectionIDProp, IntProperty* mouseOverIDProp, const int noSelectionID);
    /**
     * Destructor.
     */
    ~CMSelectionManager();

    virtual std::string getClassName() const   { return "CMSelectionManager";     }
    virtual InteractionHandler* create() const { return new CMSelectionManager(); }

    /**
     * Resize the frame buffer used for ID rendering.
     */
    void resize(tgt::ivec2 newsize);

    /**
     * Initialize. Has to be called before using the CMSelectionManager.
     */

    void initialize();
    /**
     * Deinitialize. Has to be called before destroying the CMSelectionManager.
     * It frees up additionally allocated ressources.
     */
    void deinitialize();

    /**
     * Activate the interal render target and make it ready for being renderd into.
     */
    void activate();

    /**
     * Deactivate the interal render target after rendering to it.
     */
    void deactivate();

    /**
     * Clear the internal render target.
     */
    void clear();

    /*
     * Manually get the ID last rendered at the specified position.
     * @return rendered id at pos + noSelectionID
     * @note Currently the whole texture will be downloaded, so use carefully!
     */
    int idAt(tgt::ivec2 pos);

    virtual void onEvent(tgt::Event* e);

protected:

    /**
     * Handle mouse movement events and potentially update mouseOverIDProp_.
     */
    void mouseMoveEvent(tgt::MouseEvent* e);

    /**
     * Handle mouse press events and potentially update selectionIDProp_.
     */
    void mousePressEvent(tgt::MouseEvent* e);

private:
    /// ID that will be returned if nothing was rendered at a given position.
    /// Because of that IDs returned by idAt will be the id rendered + noSelectionID.
    const int noSelectionID_;
    /// Location of the last mouse press event. Used to ignore mouse shift events
    /// and only update on actual clicks on the same position.
    tgt::ivec2 mousePressedPos_;

    /// Render target for IDs
    RenderTarget renderTarget_;
    /// Prop Property that will hold the ID of the last element clicked.
    IntProperty* selectionIDProp_;
    /// Prop Property that will hold the ID of the last element hovered over.
    IntProperty* mouseOverIDProp_;

    /// Mouse movement event
    EventProperty<CMSelectionManager>* moveEvent_;
    /// Mouse click event
    EventProperty<CMSelectionManager>* clickEvent_;
};
}
#endif //VRN_CMSELECTIONMANAGER_H
