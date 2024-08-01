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

#include "voreen/core/utils/GLSLparser/annotations/annotationparser.h"

#include "voreen/core/utils/GLSLparser/annotations/annotationlexer.h"
#include "voreen/core/utils/GLSLparser/annotations/annotationnodes.h"

namespace voreen {

namespace glslparser {

AnnotationParser::AnnotationParser(const std::list<Token*>& tokens)
    : Parser(0, false),
    tokens_(tokens)
{
}

AnnotationParser::~AnnotationParser() {
    while (! tokens_.empty()) {
        delete tokens_.front();
        tokens_.pop_front();
    }
}

// protected
//

Token* AnnotationParser::nextToken() {
    if (! tokens_.empty()) {
        Token* const t = tokens_.front();
        tokens_.pop_front();
        return t;
    }

    return 0;
}

void AnnotationParser::expandParseTree(const int productionID,
                                       const std::vector<Parser::ParserSymbol*>& reductionBody)
{
    switch (productionID) {
        case 0:  // [$START$] ::= [annotation]
        case 1:  // [annotation] ::= [tag-list]
        case 2:  // [tag-list] ::= [tag]
            break;

        case 3:  // [tag-list] ::= [tag-list] [tag]
        case 4:  // [tag-list] ::= [tag-list] , [tag]
            {
                AnnotationTags* const tag = popNode<AnnotationTags>();
                AnnotationTags* const list = popNode<AnnotationTags>();

                if ((list != 0) && (tag != 0)) {
                    list->merge(tag);
                    pushNode(list);
                }
                delete tag;
            }
            break;

        case 5:  // [tag] ::= [tag-assignment] [ [value-list] ]
        case 6:  // [tag] ::= [tag-assignment] [single-value]
            {
                AnnotationTagValues* const values = popNode<AnnotationTagValues>();
                AnnotationTagName* const  names = popNode<AnnotationTagName>();

                if ((values != 0) && (names != 0))
                    pushNode(new AnnotationTags(names, values));
            }
            break;

        case 7:  // [tag-assignment] ::= [tag-name] =
            break;

        case 8:  // [tag-assignment] ::= [tag-assignment] [tag-name] =
            if (reductionBody.size() == 3) {
                AnnotationTagName* const name = popNode<AnnotationTagName>();
                AnnotationTagName* const list = popNode<AnnotationTagName>();

                if ((list != 0) && (name != 0)) {
                    list->addTagName(name);
                    pushNode(list);
                }
            }
            break;

        case 9:  // [tag-name] ::= @ IDENTIFIER
            if (reductionBody.size() == 2) {
                IdentifierToken* const name =
                    dynamic_cast<IdentifierToken* const>(reductionBody[0]->getToken());

                if (name != 0)
                    pushNode(new AnnotationTagName(name));
            }
            break;

        case 10:  // [value-list] ::= [single-value]
            break;

        case 11:  // [value-list] ::= [value-list] , [single-value]
            if (reductionBody.size() == 3) {
                AnnotationTagValues* const val = popNode<AnnotationTagValues>();
                AnnotationTagValues* const list = popNode<AnnotationTagValues>();

                if ((list != 0) && (val != 0)) {
                    list->merge(val);
                    pushNode(list);
                }
                delete val;
            }
            break;

        case 12:  // [single-value] ::= true
        case 13:  // [single-value] ::= false
        case 14:  // [single-value] ::= INTCONST
        case 15:  // [single-value] ::= FLOATCONST
            if (reductionBody.size() == 1) {
                Token* const token = reductionBody[0]->getToken();

                if (token != 0)
                    pushNode(new AnnotationTagValues(token));
            }
            break;

        case 16:  // [single-value] ::= " STRING "
        case 17:  // [single-value] ::= ' CHARACTER '
            if (reductionBody.size() == 3) {
                Token* const token = reductionBody[1]->getToken();

                if (token != 0)
                    pushNode(new AnnotationTagValues(token));
            }
            break;

    } // switch (productionID
}

ParserAction* AnnotationParser::action(const int stateID, const int symbolID) const {
    switch (stateID) {
        case 0:
            switch (symbolID) {
                case 4:  // @
                    return new ParserActionTransition(1);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 1:
            switch (symbolID) {
                case 10:  // IDENTIFIER
                    return new ParserActionTransition(7);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 2:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionAccept();

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 3:
            switch (symbolID) {
                case -1:  // $END$
                case 4:  // @
                case 5:  // ,
                    return new ParserActionReduce(2);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 4:
            switch (symbolID) {
                case 4:  // @
                    return new ParserActionTransition(1);


                case 12:  // FLOATCONST
                    return new ParserActionTransition(10);

                case 11:  // INTCONST
                    return new ParserActionTransition(11);

                case 1:  // [
                    return new ParserActionTransition(12);

                case 9:  // false
                    return new ParserActionTransition(13);

                case 8:  // true
                    return new ParserActionTransition(16);

                case 6:  // "
                    return new ParserActionTransition(8);

                case 7:  // '
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 5:
            switch (symbolID) {
                case -1:  // $END$
                    return new ParserActionReduce(1);

                case 4:  // @
                    return new ParserActionTransition(1);

                case 5:  // ,
                    return new ParserActionTransition(17);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 6:
            switch (symbolID) {
                case 3:  // =
                    return new ParserActionTransition(19);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 7:
            switch (symbolID) {
                case 3:  // =
                    return new ParserActionReduce(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 8:
            switch (symbolID) {
                case 13:  // STRING
                    return new ParserActionTransition(20);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 9:
            switch (symbolID) {
                case 14:  // CHARACTER
                    return new ParserActionTransition(21);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 10:
            switch (symbolID) {
                case -1:  // $END$
                case 4:  // @
                case 5:  // ,
                case 2:  // ]
                    return new ParserActionReduce(15);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 11:
            switch (symbolID) {
                case -1:  // $END$
                case 4:  // @
                case 5:  // ,
                case 2:  // ]
                    return new ParserActionReduce(14);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 12:
            switch (symbolID) {

                case 12:  // FLOATCONST
                    return new ParserActionTransition(10);

                case 11:  // INTCONST
                    return new ParserActionTransition(11);

                case 9:  // false
                    return new ParserActionTransition(13);

                case 8:  // true
                    return new ParserActionTransition(16);

                case 6:  // "
                    return new ParserActionTransition(8);

                case 7:  // '
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 13:
            switch (symbolID) {
                case -1:  // $END$
                case 4:  // @
                case 5:  // ,
                case 2:  // ]
                    return new ParserActionReduce(13);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 14:
            switch (symbolID) {
                case -1:  // $END$
                case 4:  // @
                case 5:  // ,
                    return new ParserActionReduce(6);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 15:
            switch (symbolID) {
                case 3:  // =
                    return new ParserActionTransition(24);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 16:
            switch (symbolID) {
                case -1:  // $END$
                case 4:  // @
                case 5:  // ,
                case 2:  // ]
                    return new ParserActionReduce(12);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 17:
            switch (symbolID) {
                case 4:  // @
                    return new ParserActionTransition(1);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 18:
            switch (symbolID) {
                case -1:  // $END$
                case 4:  // @
                case 5:  // ,
                    return new ParserActionReduce(3);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 19:
            switch (symbolID) {
                case 1:  // [
                case 8:  // true

                case 6:  // "
                case 7:  // '
                case 12:  // FLOATCONST
                case 11:  // INTCONST
                case 9:  // false
                case 4:  // @
                    return new ParserActionReduce(7);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 20:
            switch (symbolID) {
                case 6:  // "
                    return new ParserActionTransition(26);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 21:
            switch (symbolID) {
                case 7:  // '
                    return new ParserActionTransition(27);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 22:
            switch (symbolID) {
                case 2:  // ]
                case 5:  // ,
                    return new ParserActionReduce(10);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 23:
            switch (symbolID) {
                case 5:  // ,
                    return new ParserActionTransition(28);

                case 2:  // ]
                    return new ParserActionTransition(29);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 24:
            switch (symbolID) {
                case 1:  // [
                case 8:  // true

                case 6:  // "
                case 7:  // '
                case 12:  // FLOATCONST
                case 11:  // INTCONST
                case 9:  // false
                case 4:  // @
                    return new ParserActionReduce(8);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 25:
            switch (symbolID) {
                case -1:  // $END$
                case 4:  // @
                case 5:  // ,
                    return new ParserActionReduce(4);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 26:
            switch (symbolID) {
                case -1:  // $END$
                case 4:  // @
                case 5:  // ,
                case 2:  // ]
                    return new ParserActionReduce(16);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 27:
            switch (symbolID) {
                case -1:  // $END$
                case 4:  // @
                case 5:  // ,
                case 2:  // ]
                    return new ParserActionReduce(17);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 28:
            switch (symbolID) {

                case 12:  // FLOATCONST
                    return new ParserActionTransition(10);

                case 11:  // INTCONST
                    return new ParserActionTransition(11);

                case 9:  // false
                    return new ParserActionTransition(13);

                case 8:  // true
                    return new ParserActionTransition(16);

                case 6:  // "
                    return new ParserActionTransition(8);

                case 7:  // '
                    return new ParserActionTransition(9);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 29:
            switch (symbolID) {
                case -1:  // $END$
                case 4:  // @
                case 5:  // ,
                    return new ParserActionReduce(5);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        case 30:
            switch (symbolID) {
                case 2:  // ]
                case 5:  // ,
                    return new ParserActionReduce(11);

                default:
                    return new ParserActionError();
            }   // switch (symbolID
            break;

        default:
            return new ParserActionError();

    } // switch (stateID

    return 0;
}

int AnnotationParser::gotoState(const int stateID, const int symbolID) const {
    switch (stateID) {
        case 0:
            switch (symbolID) {
                case 15:  // [annotation]
                    return 2;

                case 17:  // [tag]
                    return 3;

                case 18:  // [tag-assignment]
                    return 4;

                case 16:  // [tag-list]
                    return 5;

                case 21:  // [tag-name]
                    return 6;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 4:
            switch (symbolID) {
                case 20:  // [single-value]
                    return 14;

                case 21:  // [tag-name]
                    return 15;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 5:
            switch (symbolID) {
                case 17:  // [tag]
                    return 18;

                case 18:  // [tag-assignment]
                    return 4;

                case 21:  // [tag-name]
                    return 6;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 12:
            switch (symbolID) {
                case 20:  // [single-value]
                    return 22;

                case 19:  // [value-list]
                    return 23;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 17:
            switch (symbolID) {
                case 17:  // [tag]
                    return 25;

                case 18:  // [tag-assignment]
                    return 4;

                case 21:  // [tag-name]
                    return 6;

                default:
                    return -1;
            }   // switch (symbolID
            break;

        case 28:
            switch (symbolID) {
                case 20:  // [single-value]
                    return 30;

                default:
                    return -1;
            }   // switch (symbolID
            break;

    } // switch (stateID

    return -1;
}

ProductionStub* AnnotationParser::findProduction(const int productionID) const {
    switch (productionID) {
        case 0:  // [$START$] ::= [annotation]
            {
                int bodyIDs[] = {15};
                return new ProductionStub(0, bodyIDs, 1);
            }
        case 1:  // [annotation] ::= [tag-list]
            {
                int bodyIDs[] = {16};
                return new ProductionStub(15, bodyIDs, 1);
            }
        case 2:  // [tag-list] ::= [tag]
            {
                int bodyIDs[] = {17};
                return new ProductionStub(16, bodyIDs, 1);
            }
        case 3:  // [tag-list] ::= [tag-list] [tag]
            {
                int bodyIDs[] = {16, 17};
                return new ProductionStub(16, bodyIDs, 2);
            }
        case 4:  // [tag-list] ::= [tag-list] , [tag]
            {
                int bodyIDs[] = {16, 5, 17};
                return new ProductionStub(16, bodyIDs, 3);
            }
        case 5:  // [tag] ::= [tag-assignment] [ [value-list] ]
            {
                int bodyIDs[] = {18, 1, 19, 2};
                return new ProductionStub(17, bodyIDs, 4);
            }
        case 6:  // [tag] ::= [tag-assignment] [single-value]
            {
                int bodyIDs[] = {18, 20};
                return new ProductionStub(17, bodyIDs, 2);
            }
        case 7:  // [tag-assignment] ::= [tag-name] =
            {
                int bodyIDs[] = {21, 3};
                return new ProductionStub(18, bodyIDs, 2);
            }
        case 8:  // [tag-assignment] ::= [tag-assignment] [tag-name] =
            {
                int bodyIDs[] = {18, 21, 3};
                return new ProductionStub(18, bodyIDs, 3);
            }
        case 9:  // [tag-name] ::= @ IDENTIFIER
            {
                int bodyIDs[] = {4, 10};
                return new ProductionStub(21, bodyIDs, 2);
            }
        case 10:  // [value-list] ::= [single-value]
            {
                int bodyIDs[] = {20};
                return new ProductionStub(19, bodyIDs, 1);
            }
        case 11:  // [value-list] ::= [value-list] , [single-value]
            {
                int bodyIDs[] = {19, 5, 20};
                return new ProductionStub(19, bodyIDs, 3);
            }
        case 12:  // [single-value] ::= true
            {
                int bodyIDs[] = {8};
                return new ProductionStub(20, bodyIDs, 1);
            }
        case 13:  // [single-value] ::= false
            {
                int bodyIDs[] = {9};
                return new ProductionStub(20, bodyIDs, 1);
            }
        case 14:  // [single-value] ::= INTCONST
            {
                int bodyIDs[] = {11};
                return new ProductionStub(20, bodyIDs, 1);
            }
        case 15:  // [single-value] ::= FLOATCONST
            {
                int bodyIDs[] = {12};
                return new ProductionStub(20, bodyIDs, 1);
            }
        case 16:  // [single-value] ::= " STRING "
            {
                int bodyIDs[] = {6, 13, 6};
                return new ProductionStub(20, bodyIDs, 3);
            }
        case 17:  // [single-value] ::= ' CHARACTER '
            {
                int bodyIDs[] = {7, 14, 7};
                return new ProductionStub(20, bodyIDs, 3);
            }
    } // switch (productionID

    return 0;
}

std::string AnnotationParser::symbolID2String(const int symbolID) const {
    switch (symbolID) {
        case 0: return "[$START$]";
        case 1: return "[";
        case 2: return "]";
        case 3: return "=";
        case 4: return "@";
        case 5: return ",";
        case 6: return "\"";
        case 7: return "'";
        case 8: return "true";
        case 9: return "false";
        case 10: return "IDENTIFIER";
        case 11: return "INTCONST";
        case 12: return "FLOATCONST";
        case 13: return "STRING";
        case 14: return "CHARACTER";
        case 15: return "[annotation]";
        case 16: return "[tag-list]";
        case 17: return "[tag]";
        case 18: return "[tag-assignment]";
        case 19: return "[value-list]";
        case 20: return "[single-value]";
        case 21: return "[tag-name]";
        case 22: return "$END$";

    }  // switch (symbolID
    return "";
}

}   // namespace glslparser

}   // namespace voreen
