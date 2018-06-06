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

#ifndef VRN_PLOTPICKINGMANAGER_H
#define VRN_PLOTPICKINGMANAGER_H

#include "tgt/vector.h"
#include "tgt/assert.h"
#include "voreen/core/datastructures/rendertarget/rendertarget.h"

namespace voreen {

/**
 * \brief Manager class for picking plot objects.
 *
 * This class can convert tgt::col3 to a PlotData table position (row, column)
 * and vice versa.
 *
 * Table positions are interpreted as follows:
 * (-1,-1)       - no cell selected
 * (row, -1)     - entire row selected
 * (-1, column)  - entire column selected
 * (row, column) - single cell selected.
 *
 * The mapping of positions to colors is done as follows:
 * (0, 0, 0) adresses no cell.
 * A tgt::col3 has a size of (3*8)=24 bits, of which we reserve the lower
 * ceil(log_2(columnCount)) bits for the column index, the remaining bits
 * are usef for the row index.
 * The maximum number of cells which can be addressed with this encoding is
 * (2^24 - 1) = 16 777 215. In cases where columnCount is not a power of 2,
 * the maximum number of addressable cells will be less, but in any way be
 * at least (2^23 - 1) = 8 388 607.
 *
 * On exceeding these limits PlotPickingManager will raise an assertion
 * and set the color to (0, 0, 0).
 */
class VRN_CORE_API PlotPickingManager {
public:
    /// default constructor
    PlotPickingManager();

    /// Initializes current render target.
    void initializeTarget();

    /**
     * \brief Calculates how many bits must be reserved for storing the column number.
     *
     * \note    If columnCount changes all previously generated colors may become
     *          invalid. PlotPickingManager does not take care of this!
     */
    void setColumnCount(int columnCount);

    /**
     * \brief Encodes Converts (row, column) to a color regarding current columnCount.
     *
     * \param   row     row to encode, -1 means that the hole column is mapped to the color
     * \param   column  column to encode, -1 means that the hole row is mapped to the color
     */
    tgt::col3 getColorFromCell(int row, int column);

    /**
     * \brief Encodes Converts (row, column) to a color regarding current columnCount and sets it in current OpenGL context.
     *
     * \param   row     row to encode, -1 means that the hole column is mapped to the color
     * \param   column  column to encode, -1 means that the hole row is mapped to the color
     */
    void setGLColor(int row, int column);

    /**
     * \brief Decodes color to cell (row, column) regarding current columnCount.
     * This is the inverse function of getColorFromCell.
     *
     * \param   color   color to decode
     *
     * \return  ivec2 (row, column) containing decoded cell
     */
    tgt::ivec2 getCellFromColor(tgt::col3 color);

    /**
     * \brief Decodes the color at viewport position \a pos of current render target regarding current columnCount.
     *
     * Returns cell at \a pos using getColorAtPos an getCellFromColor. (-1, column) means that
     * the hole column it mapped to this position. (-1,-1) means, that no cell is mapped to this position.
     *
     * \param   pos     position in viewport coordinates
     *
     * \return  ivec2 (row, column) containing decoded cell
     */
    tgt::ivec2 getCellAtPos(tgt::ivec2 pos);

    /**
     * \brief Returns the color at viewport position \a pos of current render target.
     *
     * \param   pos     position in viewport coordinates
     */
    tgt::col3 getColorAtPos(tgt::ivec2 pos);

    /// Converts tgt::col3 to tgt::Color.
    tgt::Color convertColor(tgt::col3 col);

    /// Activates given RenderTarget and sets it as OpenGL context.
    void activateTarget(std::string debugLabel = "");

    void deactivateTarget();

    /// Clears current render target.
    void clearTarget();

    /**
     * \brief   Sets current render target to \a rt.
     *
     * \note    As there are no deep copies perfomed, \a rt must exist as long as it is
     *          set to this PlotPickingManager. Violation will result in undefined behaviour,
     *          most likely a SEGFAULT.
     *
     * \param   rt      pointer to render target which shall be used in future
     **/
    void setRenderTarget(RenderTarget* rt);

    /// returns pointer to current renderTarget
    RenderTarget* getRenderTarget();

private:

    RenderTarget* rt_;              ///< render target where information will be stored
    int bitsReservedForColumns_;    ///< number of bits needed for encoding columns
    int numbersReservedForColumns_; ///< equals (2^bitsReservedForColumns_)

};


} //namespace voreen

#endif
