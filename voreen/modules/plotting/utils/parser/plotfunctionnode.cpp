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

#include "plotfunctionnode.h"
#include "plotfunctionterminals.h"

namespace voreen {

namespace glslparser {

PlotFunctionNode::PlotFunctionNode(const int symbolID)
    :ParseTreeNode(symbolID)
{}

int PlotFunctionNode::getNodeType() const
{
    return PlotFunctionNodeTypes::NODE_PLOTFUNCTION;
}

PlotFunctionNode* PlotFunctionNode::getCopy() const {
    return new PlotFunctionNode(*this);
}

// ============================================================================

PlotFunctionTerm::PlotFunctionTerm(Token* const token, bool deleteToken)
    :PlotFunctionNode(token->getTokenID())
    , token_(deleteToken? token : dynamic_cast<Token* const>(token->getCopy()))
{
}

PlotFunctionTerm::~PlotFunctionTerm() {
    delete token_;
}

int PlotFunctionTerm::getNodeType() const
{
    return PlotFunctionNodeTypes::NODE_PLOTFUNCTION_TERM;
}

PlotFunctionNode* PlotFunctionTerm::getCopy() const {
    return new PlotFunctionTerm(*this);
}

Token* PlotFunctionTerm::getToken() const {
    return token_;
}

// ============================================================================

PlotFunctionRangeTerm::PlotFunctionRangeTerm(BracketToken* const bracket,
    PlotFunctionTerm* const term)
    : PlotFunctionTerm(new VariablesToken(PlotFunctionTerminals::ID_RANGE_TERM,"RANGE-TERM"),true)
    , bracket_(dynamic_cast<BracketToken* const>(bracket->getCopy()))
    , term_(dynamic_cast<PlotFunctionTerm* const>(term->getCopy()))
{
    addChild(term_);
    tokens_.push_back(getToken());
    tokens_.push_back(bracket_);
}


PlotFunctionRangeTerm::~PlotFunctionRangeTerm() {
    delete bracket_;
    tokens_.clear();
}

int PlotFunctionRangeTerm::getNodeType() const {
    return PlotFunctionNodeTypes::NODE_RANGE_TERM;
}

PlotFunctionNode* PlotFunctionRangeTerm::getCopy() const {
    return new PlotFunctionRangeTerm(bracket_,term_);
}

const std::vector<Token*>& PlotFunctionRangeTerm::getTokens() const {
    return tokens_;
}

// ============================================================================

PlotFunctionPreOperation::PlotFunctionPreOperation(OperatorToken* const op,
    PlotFunctionTerm* const term, bool deleteToken)
    : PlotFunctionTerm(op,deleteToken)
    , op_(dynamic_cast<OperatorToken* const>(op->getCopy()))
    , term_(dynamic_cast<PlotFunctionTerm* const>(term->getCopy()))
{
    addChild(term_);
}


PlotFunctionPreOperation::~PlotFunctionPreOperation() {
    delete op_;
//    delete term_;
}

int PlotFunctionPreOperation::getNodeType() const {
    return PlotFunctionNodeTypes::NODE_PRE_OPERATION;
}

// ============================================================================

PlotFunctionCenterOperation::PlotFunctionCenterOperation(OperatorToken* const op,
    PlotFunctionTerm* const preterm,
    PlotFunctionTerm* const postterm, bool deleteToken)
    : PlotFunctionTerm(op,deleteToken)
    , op_(dynamic_cast<OperatorToken* const>(op->getCopy()))
    , preterm_(dynamic_cast<PlotFunctionTerm* const>(preterm->getCopy()))
    , postterm_(dynamic_cast<PlotFunctionTerm* const>(postterm->getCopy()))
{
    addChild(preterm_);
    addChild(postterm_);
}

PlotFunctionCenterOperation::~PlotFunctionCenterOperation() {
    delete op_;
    //delete preterm_;
    //delete postterm_;
}

int PlotFunctionCenterOperation::getNodeType() const {
    return PlotFunctionNodeTypes::NODE_CENTER_OPERATION;
}
// ============================================================================

PlotFunctionIsolatorTerm::PlotFunctionIsolatorTerm(IsolatorToken* const isolator,
    PlotFunctionTerm* const preterm,
    PlotFunctionTerm* const postTerm)
    : PlotFunctionTerm(isolator)
    , isolator_(dynamic_cast<IsolatorToken* const>(isolator->getCopy()))
    , preTerm_(dynamic_cast<PlotFunctionTerm* const>(preterm->getCopy()))
    , postTerm_(dynamic_cast<PlotFunctionTerm* const>(postTerm->getCopy()))
{
    addChild(preTerm_);
    addChild(postTerm_);
}

PlotFunctionIsolatorTerm::~PlotFunctionIsolatorTerm() {
    delete isolator_;
    //delete term_;
    //delete intervalTerm_;
}

int PlotFunctionIsolatorTerm::getNodeType() const {
    return PlotFunctionNodeTypes::NODE_CENTER_ISOLATOR;
}
// ============================================================================

PlotFunctionFunctionTerm::PlotFunctionFunctionTerm(
    PlotFunctionTerm* const preterm,
    PlotFunctionTerm* const postTerm)
    : PlotFunctionTerm(new FunctionToken(PlotFunctionTerminals::ID_FUNCTION_TERM,"FUNCTION-TERM",-1),true)
    , preTerm_(dynamic_cast<PlotFunctionTerm* const>(preterm->getCopy()))
    , postTerm_(dynamic_cast<PlotFunctionTerm* const>(postTerm->getCopy()))
{
    addChild(preTerm_);
    addChild(postTerm_);
}

PlotFunctionFunctionTerm::~PlotFunctionFunctionTerm() {
    //delete isolator_;
    //delete term_;
    //delete intervalTerm_;
}

int PlotFunctionFunctionTerm::getNodeType() const {
    return PlotFunctionNodeTypes::NODE_FUNCTION_TERM;
}


// ============================================================================

PlotFunctionTermList::PlotFunctionTermList(IsolatorToken* const isolator,
    PlotFunctionTerm* const term)
    : PlotFunctionTerm(isolator)
    , isolator_(dynamic_cast<IsolatorToken* const>(isolator->getCopy()))
{
    termlist_.push_back(dynamic_cast<PlotFunctionTerm* const>(term->getCopy()));
}

PlotFunctionTermList::~PlotFunctionTermList() {
    delete isolator_;
    for (size_t i = 0; i < termlist_.size(); ++i)
        delete termlist_[i];

    //delete term_;
    //delete intervalTerm_;
}

int PlotFunctionTermList::getNodeType() const {
    return PlotFunctionNodeTypes::NODE_TERM_LIST;
}

void PlotFunctionTermList::merge(PlotFunctionTermList* const other) {
    if (other != 0) {
        termlist_.insert(termlist_.end(), other->termlist_.begin(), other->termlist_.end());
        other->termlist_.clear();
    }

}


// ============================================================================

PlotFunctionTagName::PlotFunctionTagName(IdentifierToken* const name)
    : PlotFunctionTerm(name),
    name_(dynamic_cast<IdentifierToken* const>(name->getCopy()))
{}

PlotFunctionTagName::~PlotFunctionTagName() {
    delete name_;
    for (size_t i = 0; i < list_.size(); ++i)
        delete list_[i];
}

int PlotFunctionTagName::getNodeType() const {
    return PlotFunctionNodeTypes::NODE_PLOTFUNCTION_TAG_NAME;
}

void PlotFunctionTagName::addTagName(PlotFunctionTagName* const tagName) {
    if (tagName != 0)
        list_.push_back(tagName);
}

IdentifierToken* PlotFunctionTagName::getName() const {
    return name_;
}

const std::vector<PlotFunctionTagName*>& PlotFunctionTagName::getMoreNames() const {
    return list_;
}

// ============================================================================

PlotFunctionTagValues::PlotFunctionTagValues(IdentifierToken* const token)
    : PlotFunctionTerm(token)
{
    addValue(token);
}

PlotFunctionTagValues::~PlotFunctionTagValues() {
    for (size_t i = 0; i < values_.size(); ++i)
        delete values_[i];
}

int PlotFunctionTagValues::getNodeType() const {
    return PlotFunctionNodeTypes::NODE_PLOTFUNCTION_TAG_VALUES;
}

void PlotFunctionTagValues::addValue(Token* const token) {
    if (token != 0)
        values_.push_back(token->getCopy());
}

void PlotFunctionTagValues::merge(PlotFunctionTagValues* const other) {
    if (other != 0) {
        values_.insert(values_.end(), other->values_.begin(), other->values_.end());
        other->values_.clear();
    }
}

const std::vector<Token*>& PlotFunctionTagValues::getTokens() const {
    return values_;
}

// ============================================================================

PlotFunctionTags::PlotFunctionTags(PlotFunctionTagName* const tagName,
    PlotFunctionTagValues* const tagValues)
    : PlotFunctionTerm(new IdentifierToken(PlotFunctionTerminals::ID_FUNCTION,0),true)
{
    addTag(tagName, tagValues);
}

PlotFunctionTags::~PlotFunctionTags() {
    for (size_t i = 0; i < tags_.size(); ++i) {
        delete tags_[i].first;
        delete tags_[i].second;
    }
}

int PlotFunctionTags::getNodeType() const {
    return PlotFunctionNodeTypes::NODE_PLOTFUNCTION_TAGS;
}

void PlotFunctionTags::addTag(PlotFunctionTagName* const tagName, PlotFunctionTagValues* const tagValues) {
    if ((tagName != 0) && (tagValues != 0))
        tags_.push_back(std::make_pair(tagName, tagValues));
}

void PlotFunctionTags::merge(PlotFunctionTags* const other) {
    if (other != 0) {
        tags_.insert(tags_.end(), other->tags_.begin(), other->tags_.end());
        other->tags_.clear();
    }
}

const std::vector<std::pair<PlotFunctionTagName*, PlotFunctionTagValues*> >& PlotFunctionTags::getTags() const {
    return tags_;
}



}   // namespace glslparser
}   // namespace voreen
