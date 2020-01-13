/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_TOUCHTABLEPREVIEWICON_H
#define VRN_TOUCHTABLEPREVIEWICON_H

#include "tgt/bounds.h"
#include "tgt/texture.h"

namespace voreen {

    /**
    * The TouchTablePreviewIcon represents a preset or a screenshot etc in a menu for the user to select.
    */
class TouchTablePreviewIcon{

public:
    /**
    * The TouchTablePreviewIcon represents a preset or a screenshot etc. in a menu for the user to select.
    *
    *@param ll: the lower left corner of the icon
    *@param width: the width of the icon
    *@param height: the height of the icon
    *@param symbolTex: the texture on th icon representing a preset, screenshot etc.
    *@param caption: a caption to the icon, e.g. the time and date the represented was taken
    */
    TouchTablePreviewIcon(tgt::ivec2 ll, int width, int height, tgt::Texture* symbolTex, std::string filename, std::string caption = "", bool isSelected = false);

    virtual tgt::ivec2 getLL() const;
    virtual void setLL(tgt::ivec2 posLL);

    virtual int getWidth() const;
    virtual void setWidth(int width);

    virtual int getHeight() const;
    virtual void setHeight(int height);

    virtual tgt::Texture* getTexture() const;
    virtual void setTexture(tgt::Texture* symbolTex);

    virtual std::string getCaption() const;
    virtual void setCaption(std::string caption);

    virtual std::string getFilename() const;
    virtual void setFilename(std::string filename);

    virtual bool isSelected() const;
    virtual void setIsSelected(bool isSelected);

    virtual bool contains(tgt::ivec2 tp);
    virtual bool intersects(tgt::ivec2 ll, tgt::ivec2 ur) const;

protected:

     tgt::Bounds bounds_;            ///< lower left corner of the icon

    int width_;                        ///< width of the icon
    int height_;                    ///< length of the icon

    tgt::Texture* symbolTex_;        ///< texture shown by the preview icon

    std::string caption_;            ///< caption to the icon

    std::string filename_;            ///< filename of the element that is represented by the icon

    bool isSelected_;                ///< bool to check whether the icon has been selected
};
}

#endif
