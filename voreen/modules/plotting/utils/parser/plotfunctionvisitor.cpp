/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "plotfunctionvisitor.h"

#include "plotfunctionlexer.h"
#include "plotfunctionnode.h"
#include "plotfunctionterminals.h"

#include <stack>

namespace voreen {

namespace glslparser {

PlotFunctionVisitor::PlotFunctionVisitor()
    : ParseTreeVisitor()
    , tokens_()
{
}

PlotFunctionVisitor::~PlotFunctionVisitor() {
    tokens_.clear();
    partialTrees_.clear();
}

bool PlotFunctionVisitor::visit(ParseTreeNode* const node) {
    if (node == 0)
        return false;
    bool res = true;
    try {
        switch (node->getNodeType()) {
            case PlotFunctionNodeTypes::NODE_PLOTFUNCTION:
                visitNode(dynamic_cast<PlotFunctionNode* const>(node));
                break;

            case PlotFunctionNodeTypes::NODE_PLOTFUNCTION_TERM:
                visitNode(dynamic_cast<PlotFunctionTerm* const>(node));
                break;

            case PlotFunctionNodeTypes::NODE_RANGE_TERM:
                visitNode(dynamic_cast<PlotFunctionRangeTerm* const>(node));
                break;

            case PlotFunctionNodeTypes::NODE_FUNCTION_TERM:
                visitNode(dynamic_cast<PlotFunctionFunctionTerm* const>(node));
                break;

            case PlotFunctionNodeTypes::NODE_PLOTFUNCTION_TAGS:
                visitNode(dynamic_cast<PlotFunctionTags* const>(node));
                break;

            case PlotFunctionNodeTypes::NODE_PLOTFUNCTION_TAG_NAME:
                visitNode(dynamic_cast<PlotFunctionTagName* const>(node));
                break;

            case PlotFunctionNodeTypes::NODE_PLOTFUNCTION_TAG_VALUES:
                visitNode(dynamic_cast<PlotFunctionTagValues* const>(node));
                break;

            case PlotFunctionNodeTypes::NODE_CENTER_ISOLATOR:
                visitNode(dynamic_cast<PlotFunctionIsolatorTerm* const>(node));
                break;

            case PlotFunctionNodeTypes::NODE_CENTER_OPERATION:
                visitNode(dynamic_cast<PlotFunctionCenterOperation* const>(node));
                break;

            case PlotFunctionNodeTypes::NODE_PRE_OPERATION:
                visitNode(dynamic_cast<PlotFunctionPreOperation* const>(node));
                break;

            case PlotFunctionNodeTypes::NODE_TERM_LIST:
                visitNode(dynamic_cast<PlotFunctionTermList* const>(node));
                break;

            default:
                res = false;
                log_ << "PlotFunctionVisitor::visit() called for ParseTreeNode!\n";
                break;
        }   // switch
    } catch (std::runtime_error& e) {
        log_ << "  Exception: " << e.what() << "\n";
        return false;
    }

    return res;
}

const std::vector<PlotFunctionVisitor::TokenVector>& PlotFunctionVisitor::getPlotFunctionToken(ParseTreeNode* const root) {
    tokens_.clear();
    temptokens_.clear();
    getPartialTrees(root);
    if (partialTrees_.size() == 0) {
        visitAll(root);
        tokens_.resize(1);
        tokens_[0].function = temptokens_;
    }
    else {
        for (size_t i = 0; i < partialTrees_.size(); ++i) {
            tokens_.resize(partialTrees_.size());
            temptokens_.clear();
            visitAll(partialTrees_[i].function);
            tokens_[i].function = temptokens_;
            temptokens_.clear();
            visitAll(partialTrees_[i].interval);
            tokens_[i].interval = temptokens_;
        }
    }
    return tokens_;
}

void PlotFunctionVisitor::devideTree(ParseTreeNode* const root) {
    if (root == 0)
        return;
    ParseTreeNode* lastroot = 0;
    std::stack<ParseTreeNode*> tree;
    tree.push(root);
    while (! tree.empty()) {
        ParseTreeNode* const top = tree.top();
        tree.pop();
        if (top == 0)   // should never happen
            continue;

        if (top->getSymbolID() != PlotFunctionTerminals::ID_COLON &&
            top->getSymbolID() != PlotFunctionTerminals::ID_COMMA &&
            top->getSymbolID() != PlotFunctionTerminals::ID_SEMICOLON)
            continue;

        // Push the children in reversed order onto the stack to keep their
        // order.
        //
        const std::vector<ParseTreeNode*>& children = top->getChildren();
        for (size_t j = 0; j < children.size(); ++j)
            tree.push(children[children.size() - 1 - j]);

        if (top->getSymbolID() == PlotFunctionTerminals::ID_COLON) {
            partialTrees_.resize(partialTrees_.size()+1);
            partialTrees_.back().function = top->getChildren().at(0);
            partialTrees_.back().interval = top->getChildren().at(1);
        }
        else if (top->getSymbolID() == PlotFunctionTerminals::ID_SEMICOLON && children.size() >= 2 &&
            children.at(1)->getSymbolID() != PlotFunctionTerminals::ID_COLON) {
            lastroot = children.at(1);
        }
    }
    if (lastroot) {
        partialTrees_.resize(partialTrees_.size()+1);
        partialTrees_.back().function = lastroot;
        partialTrees_.back().interval = 0;
    }
}

void PlotFunctionVisitor::getPartialTrees(voreen::glslparser::ParseTreeNode *const mainRoot) {
    partialTrees_.clear();
    devideTree(mainRoot);
}

/// start of the private functions

void PlotFunctionVisitor::visitNode(PlotFunctionNode* const node) {
    if (node == 0)
        return;

    if (PlotFunctionTags* const n = dynamic_cast<PlotFunctionTags* const>(node))
        visitNode(n);
    else if (PlotFunctionRangeTerm* const n = dynamic_cast<PlotFunctionRangeTerm* const>(node))
        visitNode(n);
    else if (PlotFunctionFunctionTerm* const n = dynamic_cast<PlotFunctionFunctionTerm* const>(node))
        visitNode(n);
    else if (PlotFunctionTagName* const n = dynamic_cast<PlotFunctionTagName* const>(node))
        visitNode(n);
    else if (PlotFunctionTagValues* const n = dynamic_cast<PlotFunctionTagValues* const>(node))
        visitNode(n);
    else if (PlotFunctionIsolatorTerm* const n = dynamic_cast<PlotFunctionIsolatorTerm* const>(node))
        visitNode(n);
    else if (PlotFunctionCenterOperation* const n = dynamic_cast<PlotFunctionCenterOperation* const>(node))
        visitNode(n);
    else if (PlotFunctionPreOperation* const n = dynamic_cast<PlotFunctionPreOperation* const>(node))
        visitNode(n);
    else if (PlotFunctionTermList* const n = dynamic_cast<PlotFunctionTermList* const>(node))
        visitNode(n);
}

void PlotFunctionVisitor::visitNode(PlotFunctionTerm* const node) {
    if (node == 0)
        return;
    if (PlotFunctionTags* const n = dynamic_cast<PlotFunctionTags* const>(node))
        visitNode(n);
    else if (PlotFunctionRangeTerm* const n = dynamic_cast<PlotFunctionRangeTerm* const>(node))
        visitNode(n);
    else if (PlotFunctionFunctionTerm* const n = dynamic_cast<PlotFunctionFunctionTerm* const>(node))
        visitNode(n);
    else if (PlotFunctionTagName* const n = dynamic_cast<PlotFunctionTagName* const>(node))
        visitNode(n);
    else if (PlotFunctionTagValues* const n = dynamic_cast<PlotFunctionTagValues* const>(node))
        visitNode(n);
    else if (PlotFunctionIsolatorTerm* const n = dynamic_cast<PlotFunctionIsolatorTerm* const>(node))
        visitNode(n);
    else if (PlotFunctionCenterOperation* const n = dynamic_cast<PlotFunctionCenterOperation* const>(node))
        visitNode(n);
    else if (PlotFunctionPreOperation* const n = dynamic_cast<PlotFunctionPreOperation* const>(node))
        visitNode(n);
    else if (PlotFunctionTermList* const n = dynamic_cast<PlotFunctionTermList* const>(node))
        visitNode(n);
    else
        temptokens_.push_back(node->getToken());
}

void PlotFunctionVisitor::visitNode(PlotFunctionFunctionTerm* const node) {
    if (node == 0)
        return;
    temptokens_.push_back(node->getToken());
}

void PlotFunctionVisitor::visitNode(PlotFunctionRangeTerm* const node) {
    if (node == 0)
        return;
    std::vector<Token*> token = node->getTokens();
    temptokens_.push_back(token.at(1));
}


void PlotFunctionVisitor::visitNode(PlotFunctionTags* const node) {
    if (node == 0)
        return;
    temptokens_.push_back(node->getToken());
}

void PlotFunctionVisitor::visitNode(PlotFunctionTagName* const node) {
    if (node == 0)
        return;
    temptokens_.push_back(node->getToken());
}

void PlotFunctionVisitor::visitNode(PlotFunctionTagValues* const node) {
    if (node == 0)
        return;
    const std::vector<Token*> nodes = node->getTokens();
    for (size_t i = 0; i < nodes.size(); ++i) {
        temptokens_.push_back(nodes[i]);
    }
}

void PlotFunctionVisitor::visitNode(PlotFunctionIsolatorTerm* const node) {
    if (node == 0)
        return;
    temptokens_.push_back(node->getToken());
}

void PlotFunctionVisitor::visitNode(PlotFunctionCenterOperation* const node) {
    if (node == 0)
        return;
    temptokens_.push_back(node->getToken());
}

void PlotFunctionVisitor::visitNode(PlotFunctionPreOperation* const node) {
    if (node == 0)
        return;
    temptokens_.push_back(node->getToken());
}

void PlotFunctionVisitor::visitNode(PlotFunctionTermList* const node) {
    if (node == 0)
        return;
    temptokens_.push_back(node->getToken());
}


}   // namespace glslparser

}   // namespace voreen
