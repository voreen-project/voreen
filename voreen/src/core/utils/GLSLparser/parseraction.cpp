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

#include "voreen/core/utils/GLSLparser/parseraction.h"

#include "voreen/core/utils/GLSLparser/grammarsymbol.h"

#include <sstream>

namespace voreen {

namespace glslparser {

ProductionStub::ProductionStub(const Production& production)
    : headID_(production.getHeadID())
{
    const std::vector<GrammarSymbol*>& body = production.getBody();
    for (size_t i = 0; i < body.size(); ++i)
        bodyIDs_.push_back(body[i]->getSymbolID());
}

ProductionStub::ProductionStub(const int headID, const int* bodyIDs, const int bodySize)
    : headID_(headID)
{
    for (int i = 0; i < bodySize; ++i)
        bodyIDs_.push_back(bodyIDs[i]);
}

// ============================================================================

ParserAction::~ParserAction() {}

std::string ParserAction::toString() const {
    switch (type_) {
        case ACTION_TRANSITION:
            return "transition";
        case ACTION_REDUCE:
            return "reduce";
        case ACTION_ACCEPT:
            return "accept";
        case ACTION_ERROR:
            return "error";
    }
    return "";
}

// ============================================================================

ParserActionTransition::ParserActionTransition(const int nextState)
    : ParserAction(ACTION_TRANSITION),
    nextState_(nextState)
{
}

std::string ParserActionTransition::generateCode() const {
    std::ostringstream oss;
    oss << "ParserActionTransition(" << nextState_ << ");";
    return oss.str();
}

std::string ParserActionTransition::toString() const {
    std::ostringstream oss;
    oss << ParserAction::toString() << " -> " << nextState_;
    return oss.str();
}

// ============================================================================

ParserActionReduce::ParserActionReduce(const int productionID)
    : ParserAction(ACTION_REDUCE),
    productionID_(productionID)
{
}

std::string ParserActionReduce::generateCode() const {
    std::ostringstream oss;
    oss << "ParserActionReduce(" << productionID_ << ");";
    return oss.str();
}

std::string ParserActionReduce::toString() const {
    std::ostringstream oss;
    oss << ParserAction::toString() << ": pid = " << productionID_;
    return oss.str();
}

// ============================================================================

std::string ParserActionAccept::generateCode() const {
    return "ParserActionAccept();";
}

// ============================================================================

ParserActionError::ParserActionError()
    : ParserAction(ACTION_ERROR),
    lineNumber_(0),
    posCharInLine_(0),
    errorMsg_()
{
}

ParserActionError::ParserActionError(const std::string& errorMsg)
    : ParserAction(ACTION_ERROR),
    lineNumber_(0),
    posCharInLine_(0),
    errorMsg_(errorMsg)
{
}

ParserActionError::~ParserActionError() {
}

std::string ParserActionError::getError() const {
    std::ostringstream oss;
    if ((lineNumber_ > 0) && (posCharInLine_ > 0))
        oss << "Error in line " << lineNumber_ << " near position " << posCharInLine_ << ":\n";

    oss << errorMsg_;
    return oss.str();
}

std::string ParserActionError::generateCode() const {
    return "ParserActionError();";
}

}   // namespace glslparser

}   // namespace voreen
