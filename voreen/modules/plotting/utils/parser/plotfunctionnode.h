/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_PLOTFUNCTIONNODE_H
#define VRN_PLOTFUNCTIONNODE_H

#include "voreen/core/utils/GLSLparser/parsetreenode.h"
#include "plotfunctiontoken.h"

namespace voreen {

namespace glslparser {


/// struct to identify the nodetypes
struct PlotFunctionNodeTypes {
    enum {
        NODE_PLOTFUNCTION,
            NODE_PLOTFUNCTION_TERM,
                NODE_TERM_LIST,
                NODE_CENTER_OPERATION,
                NODE_PRE_OPERATION,
                NODE_CENTER_ISOLATOR,
                NODE_FUNCTION_TERM,
                NODE_RANGE_TERM,
                NODE_PLOTFUNCTION_TAGS,
                NODE_PLOTFUNCTION_TAG_NAME,
                NODE_PLOTFUNCTION_TAG_VALUES
    };
};

// ============================================================================

class PlotFunctionNode : public ParseTreeNode {
public:
    PlotFunctionNode(const int symbolID);

    virtual ~PlotFunctionNode() {};

    // returns the Type of the Node
    virtual int getNodeType() const;
    //returns a copy of the PlotFunctionObject
    virtual PlotFunctionNode* getCopy() const;

};

// ============================================================================

class PlotFunctionTerm : public PlotFunctionNode {
public:
    PlotFunctionTerm(Token* const token, bool deleteToken = false);

    virtual ~PlotFunctionTerm();

    // returns the Type of the Node
    virtual int getNodeType() const;
    //returns a copy of the PlotFunctionObject
    virtual PlotFunctionNode* getCopy() const;
    // returns the Pointer to the saved Token
    virtual Token* getToken() const;

protected:
    Token* const token_;
};

// ============================================================================

class PlotFunctionRangeTerm : public PlotFunctionTerm {
public:
    PlotFunctionRangeTerm(BracketToken* const bracket, PlotFunctionTerm* const term);

    virtual ~PlotFunctionRangeTerm();

    // returns the Type of the Node
    virtual int getNodeType() const;
    //returns a copy of the PlotFunctionObject
    virtual PlotFunctionNode* getCopy() const;
    // returns the vector of Pointer to the saved Tokens
    const std::vector<Token*>& getTokens() const;


protected:
    BracketToken* const bracket_;
    PlotFunctionTerm* const term_;
    std::vector<Token*> tokens_;
};


// ============================================================================

class PlotFunctionPreOperation : public PlotFunctionTerm {
public:
    PlotFunctionPreOperation(OperatorToken* const op,
        PlotFunctionTerm* const term, bool deleteToken = false);

    virtual ~PlotFunctionPreOperation();

    // returns the Type of the Node
    virtual int getNodeType() const;

protected:
    OperatorToken* const op_;
    PlotFunctionTerm* const term_;
};

// ============================================================================

class PlotFunctionCenterOperation : public PlotFunctionTerm {
public:
    PlotFunctionCenterOperation(OperatorToken* const op, PlotFunctionTerm* const preterm,
        PlotFunctionTerm* const postterm, bool deleteToken = false);

    virtual ~PlotFunctionCenterOperation();

    // returns the Type of the Node
    virtual int getNodeType() const;

protected:
    OperatorToken* const op_;
    PlotFunctionTerm* const preterm_;
    PlotFunctionTerm* const postterm_;
};

// ============================================================================

class PlotFunctionIsolatorTerm : public PlotFunctionTerm {
public:
    PlotFunctionIsolatorTerm(IsolatorToken* const isolator,
        PlotFunctionTerm* const preTerm,
        PlotFunctionTerm* const postTerm);

    virtual ~PlotFunctionIsolatorTerm();

    // returns the Type of the Node
    virtual int getNodeType() const;

protected:
    IsolatorToken* const isolator_;
    PlotFunctionTerm* const preTerm_;
    PlotFunctionTerm* const postTerm_;
};
// ============================================================================

class PlotFunctionFunctionTerm : public PlotFunctionTerm {
public:
    PlotFunctionFunctionTerm(
        PlotFunctionTerm* const preTerm,
        PlotFunctionTerm* const postTerm);

    virtual ~PlotFunctionFunctionTerm();

    // returns the Type of the Node
    virtual int getNodeType() const;

protected:
//    IsolatorToken* const isolator_;
    PlotFunctionTerm* const preTerm_;
    PlotFunctionTerm* const postTerm_;
};

// ============================================================================
class PlotFunctionTermList : public PlotFunctionTerm {
public:
    PlotFunctionTermList(IsolatorToken* const isolator, PlotFunctionTerm* const term);

    virtual ~PlotFunctionTermList();

    // returns the Type of the Node
    virtual int getNodeType() const;
    /// Merge to PlotFunctionTermLists
    void merge(PlotFunctionTermList* const other);

protected:
    IsolatorToken* const isolator_;
    std::vector<PlotFunctionTerm*> termlist_;
};

// ============================================================================

class PlotFunctionTagName : public PlotFunctionTerm {
public:
    PlotFunctionTagName(IdentifierToken* const name);

    virtual ~PlotFunctionTagName();

    // returns the Type of the Node
    virtual int getNodeType() const;
    // add one more TagName
    void addTagName(PlotFunctionTagName* const tagName);
    // get the Main Token
    IdentifierToken* getName() const;
    //returns the saved PlotfunctionTagnames
    const std::vector<PlotFunctionTagName*>& getMoreNames() const;

protected:
    IdentifierToken* const name_;
    std::vector<PlotFunctionTagName*> list_;
};

// ============================================================================

class PlotFunctionTagValues : public PlotFunctionTerm {
public:
    PlotFunctionTagValues(IdentifierToken* const token);

    virtual ~PlotFunctionTagValues();

    // returns the Type of the Node
    virtual int getNodeType() const;
    // add a further Token to this TagValue
    void addValue(Token* const token);
    //Merge two PlotFunctionTagValues
    void merge(PlotFunctionTagValues* const other);
    // returns the vector of Pointers to the saved Tokens
    const std::vector<Token*>& getTokens() const;

protected:
    std::vector<Token*> values_;
};

// ============================================================================

class PlotFunctionTags : public PlotFunctionTerm {
public:
    PlotFunctionTags(PlotFunctionTagName* const tagName, PlotFunctionTagValues* const tagValues);

    virtual ~PlotFunctionTags();

    // returns the Type of the Node
    virtual int getNodeType() const;
    // add a PlotFunctionstags to this
    void addTag(PlotFunctionTagName* const tagName, PlotFunctionTagValues* const tagValues);
    // Merge two PlotFunctionstags in this
    void merge(PlotFunctionTags* const other);
    // return the saved pairs of Tagsnames and Tagvalues
    const std::vector<std::pair<PlotFunctionTagName*, PlotFunctionTagValues*> >& getTags() const;

protected:
    typedef std::pair<PlotFunctionTagName*, PlotFunctionTagValues*> Tag;
    std::vector<Tag> tags_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
