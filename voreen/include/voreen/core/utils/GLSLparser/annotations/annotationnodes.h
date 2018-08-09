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

#ifndef VRN_ANNOTATIONNODES_H
#define VRN_ANNOTATIONNODES_H

#include "voreen/core/utils/GLSLparser/annotations/annotationterminals.h"
#include "voreen/core/utils/GLSLparser/parsetreenode.h"

namespace voreen {

namespace glslparser {

struct AnnotationNodeTypes {
    enum {
        NODE_ANNOTATION,
        NODE_ANNOTATION_TAGS,
        NODE_ANNOTATION_TAG_NAME,
        NODE_ANNOTATION_TAG_VALUES
    };
};

// ============================================================================

class AnnotationNode : public ParseTreeNode {
public:
    AnnotationNode(const int symbolID)
        : ParseTreeNode(symbolID)
    {
    }

    virtual int getNodeType() const { return AnnotationNodeTypes::NODE_ANNOTATION; }
};

// ============================================================================

class AnnotationTagName : public AnnotationNode {
public:
    AnnotationTagName(IdentifierToken* const name)
        : AnnotationNode(name->getTokenID()),
        name_(dynamic_cast<IdentifierToken* const>(name->getCopy()))
    {
    }

    virtual ~AnnotationTagName() {
        delete name_;
        for (size_t i = 0; i < list_.size(); ++i)
            delete list_[i];
    }

    virtual int getNodeType() const { return AnnotationNodeTypes::NODE_ANNOTATION_TAG_NAME; }

    void addTagName(AnnotationTagName* const tagName) {
        if (tagName != 0)
            list_.push_back(tagName);
    }

    IdentifierToken* getName() const { return name_; }
    const std::vector<AnnotationTagName*>& getMoreNames() const { return list_; }

protected:
    IdentifierToken* const name_;
    std::vector<AnnotationTagName*> list_;
};

// ============================================================================

class AnnotationTagValues : public AnnotationNode {
public:
    AnnotationTagValues(Token* const token)
        : AnnotationNode(token->getTokenID())
    {
        addValue(token);
    }

    virtual ~AnnotationTagValues() {
        for (size_t i = 0; i < values_.size(); ++i)
            delete values_[i];
    }

    virtual int getNodeType() const { return AnnotationNodeTypes::NODE_ANNOTATION_TAG_VALUES; }

    void addValue(Token* const token) {
        if (token != 0)
            values_.push_back(token->getCopy());
    }

    void merge(AnnotationTagValues* const other) {
        if (other != 0) {
            values_.insert(values_.end(), other->values_.begin(), other->values_.end());
            other->values_.clear();
        }
    }

    const std::vector<Token*>& getTokens() const { return values_; }

protected:
    std::vector<Token*> values_;
};

// ============================================================================

class AnnotationTags : public AnnotationNode {
public:
    AnnotationTags(AnnotationTagName* const tagName, AnnotationTagValues* const tagValues)
        : AnnotationNode(AnnotationTerminals::ID_AT)
    {
        addTag(tagName, tagValues);
    }

    virtual ~AnnotationTags() {
        for (size_t i = 0; i < tags_.size(); ++i) {
            delete tags_[i].first;
            delete tags_[i].second;
        }
    }

    virtual int getNodeType() const { return AnnotationNodeTypes::NODE_ANNOTATION_TAGS; }

    void addTag(AnnotationTagName* const tagName, AnnotationTagValues* const tagValues) {
        if ((tagName != 0) && (tagValues != 0))
            tags_.push_back(std::make_pair(tagName, tagValues));
    }

    void merge(AnnotationTags* const other) {
        if (other != 0) {
            tags_.insert(tags_.end(), other->tags_.begin(), other->tags_.end());
            other->tags_.clear();
        }
    }

    const std::vector<std::pair<AnnotationTagName*, AnnotationTagValues*> >& getTags() const {
        return tags_;
    }

protected:
    typedef std::pair<AnnotationTagName*, AnnotationTagValues*> Tag;
    std::vector<Tag> tags_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
