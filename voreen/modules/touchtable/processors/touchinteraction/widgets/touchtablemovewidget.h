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

#ifndef VRN_TOUCHTABLEMOVEWIDGET_H
#define VRN_TOUCHTABLEMOVEWIDGET_H

#include "touchtablewidget.h"

namespace voreen {

/**
 * TouchScreenWidget that toggles the "movement" flag of the TouchScreenOverlay,
 * ie. switches between changing the camera orientation / changing the camera position.
 */
class TouchTableMoveWidget : public TouchTableWidget {

public:

    TouchTableMoveWidget();

    virtual Processor* create() const;

    virtual std::string getClassName() const;

    virtual void initialize();

    /**
     * Called by TouchScreenOverlay if a TouchPoint associated with this widget is released.
     * Toggles the "movement" flag of the overlay processor.
     *
     * @see TouchTableWidget
     */
    virtual void pressed();

    /**
     * Updates the active status depending on the camera movement status of the overlay.
     *
     * @see TouchTableWidget
    */
    virtual void setOverlay(TouchTableOverlay* overlay);


protected:

    virtual void setDescriptions();

};

} // namespace

#endif // VRN_TOUCHTALBEMOVEWIDGET_H
