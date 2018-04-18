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

#include "voreen/core/utils/GLSLparser/annotations/annotationvisitor.h"

#include <list>

namespace voreen {

namespace glslparser {

AnnotationVisitor::AnnotationVisitor(const bool freeAnnotations)
    : freeAnnotations_(freeAnnotations)
{
}

AnnotationVisitor::~AnnotationVisitor() {
    if (freeAnnotations_)
        freeAnnotations();
}

bool AnnotationVisitor::visit(ParseTreeNode* const node) {
    if (node == 0)
        return false;

    bool res = true;
    try {
        switch (node->getNodeType()) {
            case AnnotationNodeTypes::NODE_ANNOTATION:
                visitNode(dynamic_cast<AnnotationNode*>(node));
                break;

            case AnnotationNodeTypes::NODE_ANNOTATION_TAGS:
                visitNode(dynamic_cast<AnnotationTags*>(node));
                break;

            case AnnotationNodeTypes::NODE_ANNOTATION_TAG_NAME:
                visitNode(dynamic_cast<AnnotationTagName*>(node));
                break;

            case AnnotationNodeTypes::NODE_ANNOTATION_TAG_VALUES:
                visitNode(dynamic_cast<AnnotationTagValues*>(node));
                break;

            default:
                res = false;
                log_ << "AnnotationVisitor::visit() called for ParseTreeNode!\n";
                break;
        }   // switch
    } catch (std::runtime_error& e) {
        log_ << "  Exception: " << e.what() << "\n";
        return false;
    }

    return res;
}

const std::vector<GLSLAnnotation*>& AnnotationVisitor::getAnnotations(ParseTreeNode* const root)
{
    freeAnnotations();
    freeAnnotations_ = false;
    visitAll(root);
    return annotations_;
}

void AnnotationVisitor::freeAnnotations() {
    for (size_t i = 0; i < annotations_.size(); ++i)
        delete annotations_[i];
    annotations_.clear();
}

// private methods
//

void AnnotationVisitor::visitNode(AnnotationNode* const node) {
    if (node == 0)
        return;

    if (AnnotationTags* const n = dynamic_cast<AnnotationTags* const>(node))
        visitNode(n);
    else if (AnnotationTagName* const n = dynamic_cast<AnnotationTagName* const>(node))
        visitNode(n);
    else if (AnnotationTagValues* const n = dynamic_cast<AnnotationTagValues* const>(node)) {
        try {
            visitNode(n);
        } catch (std::runtime_error& e) {
            log_ << "  Exception: " << e.what() << "\n";
        }
    }
}

std::vector<std::string> AnnotationVisitor::visitNode(AnnotationTagName* const node) {
    std::vector<std::string> names;
    if (node == 0)
        return names;

    std::list<AnnotationTagName*> tagNameNodes;
    tagNameNodes.push_back(node);

    while (! tagNameNodes.empty()) {
        AnnotationTagName* const front = tagNameNodes.front();
        tagNameNodes.pop_front();

        IdentifierToken* const myName = front->getName();
        if (myName != 0)
            names.push_back(myName->getValue());

        const std::vector<AnnotationTagName*>& more = front->getMoreNames();
        if (! more.empty())
            tagNameNodes.insert(tagNameNodes.end(), more.begin(), more.end());
    }

    return names;
}

void AnnotationVisitor::visitNode(AnnotationTags* const node) {
    if (node == 0)
        return;

    try {
        const std::vector<std::pair<AnnotationTagName*, AnnotationTagValues*> >& tags =
        node->getTags();
        for (size_t t = 0; t < tags.size(); ++t) {
            std::vector<std::string> tagNames = visitNode(tags[t].first);
            const std::vector<Token*>& tagValue = visitNode(tags[t].second);

            for (size_t n = 0; n < tagNames.size(); ++n)
                annotations_.push_back(new GLSLAnnotation(tagNames[n], tagValue));
        }
    } catch (std::runtime_error& e) {
        log_ << "  Exception:" << e.what() << "\n";
        return;
    }
}

std::vector<Token*> AnnotationVisitor::visitNode(AnnotationTagValues* const node)
{
    if (node != 0)
        return node->getTokens();

    return std::vector<Token*>();
}

}   // namespace glslparser

}   // namespace voreen
