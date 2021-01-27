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

#ifndef VRN_PARSERACTION_H
#define VRN_PARSERACTION_H

#include "voreen/core/utils/GLSLparser/production.h"

#include <list>
#include <string>

namespace voreen {

namespace glslparser {

class ProductionStub {
public:
    ProductionStub(const Production& production);
    ProductionStub(const int headID, const int* bodyIDs, const int bodySize);

    int getHeadID() const { return headID_; }
    const std::list<int>& getBodyIDs() const { return bodyIDs_; }

private:
    int headID_;
    std::list<int> bodyIDs_;
};

// ============================================================================

class ParserAction {
public:
    enum ActionType {
        ACTION_TRANSITION,
        ACTION_REDUCE,
        ACTION_ERROR,
        ACTION_ACCEPT
    };

public:
    virtual ~ParserAction() = 0;
    ActionType getType() const { return type_; }

    virtual std::string generateCode() const = 0;
    virtual std::string toString() const;

protected:
    ParserAction(const ActionType type) : type_(type) {
    }

protected:
    const ActionType type_;
};

// ============================================================================

class ParserActionTransition : public ParserAction {
public:
    ParserActionTransition(const int nextState);
    virtual ~ParserActionTransition() {}

    int getNextState() const { return nextState_; }

    virtual std::string generateCode() const;
    virtual std::string toString() const;

protected:
    const int nextState_;
};

// ============================================================================

class ParserActionReduce : public ParserAction {
public:
    ParserActionReduce(const int productionID);
    virtual ~ParserActionReduce() {}

    int getProductionID() const { return productionID_; }

    virtual std::string generateCode() const;
    virtual std::string toString() const;

protected:
    const int productionID_;
};

// ============================================================================

class ParserActionAccept : public ParserAction {
public:
    ParserActionAccept() : ParserAction(ACTION_ACCEPT) {}
    virtual ~ParserActionAccept() {}

    virtual std::string generateCode() const;
};

// ============================================================================

class ParserActionError : public ParserAction {
public:
    ParserActionError();
    ParserActionError(const std::string& errorMsg);
    virtual ~ParserActionError();

    virtual std::string generateCode() const;
    std::string getError() const;
    const std::string& getErrorMessage() const { return errorMsg_; }

    void setErrorMessage(const std::string& errorMsg) { errorMsg_ = errorMsg; }
    void setLineNumber(const int lineNumber) { lineNumber_ = lineNumber; }
    void setCharPosition(const int charPosition) { posCharInLine_ = charPosition; }

private:
    int lineNumber_;
    int posCharInLine_;
    std::string errorMsg_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
