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

#ifndef VRN_PARSETREENODE_H
#define VRN_PARSETREENODE_H

#include <stdexcept>
#include <sstream>
#include <vector>

namespace voreen {

namespace glslparser {

class ParseTreeNode;

class ParseTreeVisitor {
public:
    ParseTreeVisitor()
        : log_(std::ios_base::out | std::ios_base::binary)
    {
    }

    virtual ~ParseTreeVisitor() {}

    void clearLog() { log_.str(""); }
    const std::ostringstream& getLog() const { return log_; }

    virtual bool visit(ParseTreeNode* const node) = 0;

    /**
     * Calls visit for the root and all of its children and theirs and so on by
     * performing a depth-first traversal.
     */
    void visitAll(ParseTreeNode* const root);

protected:
    std::ostringstream log_;
};

// ============================================================================

class ParseTreeNode {
public:
    ParseTreeNode(const int symbolID);
    virtual ~ParseTreeNode() = 0;

    virtual int getNodeType() const = 0;

    const std::vector<ParseTreeNode*>& getChildren() const { return children_; }
    int getSymbolID() const { return symbolID_; }

    void addChild(ParseTreeNode* const child);

private:
    const int symbolID_;

    // FIXME: vector containing children is only needed in few cases, e.g. preprocessor, plotfunction
    typedef std::vector<ParseTreeNode*> ChildVector;
    ChildVector children_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
