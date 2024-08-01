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

#ifndef VRN_STREAMLINETOBOUNDINGBOX_H
#define VRN_STREAMLINETOBOUNDINGBOX_H

#include "voreen/core/processors/processor.h"

#include "../../ports/streamlinelistport.h"
#include "voreen/core/ports/geometryport.h"

namespace voreen {

    /**
     * This processor converts a StreamlineList into a bounding box which can be rendered via the BoundingBoxRenderer.
     * The output geometry consists of only two points (LLF and URB) of the StreamlineList bounds.
     *
     * @see StreamlineCreator
     * @see Streamline
     */
class VRN_CORE_API StreamlineToBoundingBox : public Processor {
public:
    StreamlineToBoundingBox();
    virtual ~StreamlineToBoundingBox();
    virtual Processor* create() const         {return new StreamlineToBoundingBox(); }

    virtual std::string getClassName() const  { return "StreamlineToBoundingBox";  }
    virtual std::string getCategory() const   { return "Converter";          }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE; }

protected:
    virtual void setDescriptions() {
        setDescription("This processor converts a StreamlineList into a bounding box which can be rendered via the BoundingBoxRenderer. " \
                       "The output geometry consists of only two points (LLF and URB) of the StreamlineList bounds.");
    }

    virtual void process();

    //--------------
    //  Member
    //--------------
private:
    StreamlineListPort inport_;     ///< inport containing the streamlines to get the bounding box
    GeometryPort outport_;          ///< port containing the 2 points of the bounding box

    static const std::string loggerCat_;
};

}   //namespace

#endif
