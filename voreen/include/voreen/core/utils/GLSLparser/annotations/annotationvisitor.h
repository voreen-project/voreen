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

#ifndef VRN_ANNOTATIONVISITOR_H
#define VRN_ANNOTATIONVISITOR_H

#include "voreen/core/utils/GLSLparser/parsetreenode.h"

#include "voreen/core/utils/GLSLparser/glslannotation.h"
#include "voreen/core/utils/GLSLparser/annotations/annotationnodes.h"

#include <stdexcept>
#include <vector>

namespace voreen {

namespace glslparser {

class AnnotationVisitor : public ParseTreeVisitor {
public:
    AnnotationVisitor(const bool freeAnnotations = true);
    virtual ~AnnotationVisitor();

    virtual bool visit(ParseTreeNode* const node);

    const std::vector<GLSLAnnotation*>& getAnnotations(ParseTreeNode* const root);

    void freeAnnotations();

private:
    void visitNode(AnnotationNode* const node);
    std::vector<std::string> visitNode(AnnotationTagName* const node);
    void visitNode(AnnotationTags* const node);
    std::vector<Token*> visitNode(AnnotationTagValues* const node);

private:
    mutable bool freeAnnotations_;
    std::vector<GLSLAnnotation*> annotations_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
