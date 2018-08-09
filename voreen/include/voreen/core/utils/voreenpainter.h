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

#ifndef VRN_VOREENPAINTER_H
#define VRN_VOREENPAINTER_H

#include "tgt/exception.h"
#include "tgt/vector.h"
#include "tgt/painter.h"

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/property.h"

#include <vector>
#include <list>

namespace tgt {
class GLCanvas;
}

namespace voreen {

class NetworkEvaluator;
class CanvasRenderer;

/**
 * Specialized painter for voreen. This class takes care of rendering to a given tgt::GLCanvas,
 * using the assigned NetworkEvaluator.
 */
class VRN_CORE_API VoreenPainter : public tgt::Painter {
public:

    /**
     * Constructor. Initializes the VoreenPainter.
     *
     * @param canvas The canvas that this painter renders on
     * @param evaluator the network evaluator whose process() function if to be called on canvas repaints.
     * @param canvasRenderer the canvas renderer associated with the canvas
     */
    VoreenPainter(tgt::GLCanvas* canvas, NetworkEvaluator* evaluator, CanvasRenderer* canvasRenderer);

    virtual ~VoreenPainter();

    virtual void initialize();

    void repaint();

    virtual void sizeChanged(const tgt::ivec2& size);

    /**
     * This method gets the focus of the painter's canvas (getCanvas()) , calls the \a processor_ and actually paints something.
     */
    virtual void paint();

    /**
     * Returns the associated NetworkEvaluator.
     */
    NetworkEvaluator* getEvaluator() const;

    /**
     * Returns the associated CanvasRenderer.
     */
    CanvasRenderer* getCanvasRenderer() const;

    /**
     * Renders a snapshot to a file.
     *
     * The snapshot generation is delegated to the associated CanvasRenderer.
     *
     * @note requires Voreen to be compiled with DevIL support
     *
     * @param fileName name of the snapshot file
     * @param size size of the snapshot
     */
    void renderToSnapshot(const std::string& fileName, const tgt::ivec2& size);

protected:

    NetworkEvaluator* evaluator_;

    CanvasRenderer* canvasRenderer_;

    static const std::string loggerCat_;

};

} // namespace voreen

#endif // VRN_VOREENPAINTER_H
