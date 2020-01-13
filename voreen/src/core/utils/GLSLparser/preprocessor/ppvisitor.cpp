/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/core/utils/GLSLparser/preprocessor/ppvisitor.h"

#include "voreen/core/utils/GLSLparser/preprocessor/ppparser.h"
#include "tgt/filesystem.h"

namespace voreen {

namespace glslparser {

PreprocessorVisitor::PreprocessorVisitor()
    : translation_(std::ios_base::out | std::ios_base::in),
    symbols_("global")
{
    // TODO: replace the value returned by this macro with the correct one when it is
    // expanded or invoked.
    //
    TokenList* fileMacro = new TokenList();
    ConstantToken * fileToken = new ConstantToken(PreprocessorTerminals::ID_CONSTANT, "0", ConstantToken::TYPE_INT);
    fileMacro->addToken(fileToken);
    symbols_.insertSymbol(new PreprocessorSymbol("__FILE__", false, fileMacro));
    delete fileToken;

    // TODO: replace the value returned by this macro with the correct one when it is
    // expanded or invoked.
    //
    TokenList* lineMacro = new TokenList();
    ConstantToken* lineToken = new ConstantToken(PreprocessorTerminals::ID_CONSTANT, "0", ConstantToken::TYPE_INT);
    lineMacro->addToken(lineToken);
    symbols_.insertSymbol(new PreprocessorSymbol("__LINE__", false, lineMacro));
    delete lineToken;

    TokenList* versionMacro = new TokenList();
    ConstantToken* versionToken = new ConstantToken(PreprocessorTerminals::ID_CONSTANT, "150",
        ConstantToken::TYPE_INT);
    versionMacro->addToken(versionToken);
    symbols_.insertSymbol(new PreprocessorSymbol("__VERSION__", false, versionMacro));
    delete versionToken;
}

PreprocessorVisitor::~PreprocessorVisitor() {
}

void PreprocessorVisitor::setIncludePath(const std::string& includePath) {
    includePath_ = tgt::FileSystem::absolutePath(includePath);
}

std::ostringstream& PreprocessorVisitor::translate(std::istream* is, const std::string& directory) {

    // First attempt to parser a present additional shader header
    //
    ParseTreeNode* headerNode = 0;
    try {
        if (! shaderHeader_.empty()) {
            std::istringstream iss(shaderHeader_);
            PreprocessorParser headerParser(&iss);
            headerNode = headerParser.parse();
            log_ << "parsing of additional shader header:\n\n";
            log_ << headerParser.getLog().str() << "\n";
        }
    } catch (std::bad_alloc&) {
        log_ << "  Exception: PreprocessorVisitor::translate(): failed to created parser for shader header!\n";
    }

    try {
        setIncludePath(directory);
        PreprocessorParser parser(is);
        ParseTreeNode* root = parser.parse();

        if (headerNode && root) {
            headerNode->addChild(root);
            root = headerNode;
        }
        else {
            if (headerNode)
                root = headerNode;
            log_ << parser.getLog().str() << "\n";
        }

        translate(root);
        delete root;

    } catch (std::bad_alloc&) {
        log_ << "  Exception: PreprocessorVisitor::translate(): failed to created preprocessor parser!\n";
    }

    return translation_;
}

std::ostringstream& PreprocessorVisitor::translate(ParseTreeNode* const root) {
    translation_.str("");

    visitAll(root);

    return translation_;
}

bool PreprocessorVisitor::visit(ParseTreeNode* const node) {
    if (node == 0)
        return false;

    bool res = true;
    try {
        switch (node->getNodeType()) {
            case PreprocessorNodeTypes::NODE_EXPRESSION:
                visitNode(dynamic_cast<Expression* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_PARENTHESIS:
                visitNode(dynamic_cast<ParenthesisExpression* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_EXPRESSIONLIST:
                visitNode(dynamic_cast<ExpressionList* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_ARITHMETIC:
                visitNode(dynamic_cast<ArithmeticExpression* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_INT_CONSTANT:
                visitNode(dynamic_cast<IntConstant* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_UNARY_EXPRESSION:
                visitNode(dynamic_cast<UnaryExpression* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_BINARY_EXPRESSION:
                visitNode(dynamic_cast<BinaryExpression* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_LOGICAL_EXPRESSION:
                visitNode(dynamic_cast<LogicalExpression* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_BINARY_LOGICAL:
                visitNode(dynamic_cast<LogicalBinaryExpression* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_DEFINED_OPERATOR:
                visitNode(dynamic_cast<DefinedOperator* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_TEXT:
                visitNode(dynamic_cast<TextNode* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_TOKEN_LIST:
                visitNode(dynamic_cast<TokenList* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_MACRO:
                visitNode(dynamic_cast<Macro* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_STATEMENT:
                visitNode(dynamic_cast<Statement* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_IDENTIFIER_LIST:
                visitNode(dynamic_cast<IdentifierList* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_DEFINE:
                visitNode(dynamic_cast<DefineDirective* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_ERROR:
                visitNode(dynamic_cast<ErrorDirective* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_EXTENSION:
                visitNode(dynamic_cast<ExtensionDirective* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_INCLUDE:
                visitNode(dynamic_cast<IncludeDirective* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_LINE:
                visitNode(dynamic_cast<LineDirective* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_PRAGMA:
                visitNode(dynamic_cast<PragmaDirective* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_VERSION:
                visitNode(dynamic_cast<VersionDirective* const>(node));
                break;

            case PreprocessorNodeTypes::NODE_CONDITIONAL:
            case PreprocessorNodeTypes::NODE_IF:
            case PreprocessorNodeTypes::NODE_IFDEF:
            case PreprocessorNodeTypes::NODE_IFNDEF:
                visitNode(dynamic_cast<ConditionalDirective* const>(node));
                break;

            default:
                res = false;
                log_ << "visit() called for ParseTreeNode...\n";
                break;
        }   // switch (symbolID)
    } catch (std::runtime_error& e) {
        log_ << "  Exception: " << e.what() << "\n";
        res = false;
    }

    return res;
}

// private methods
//

int PreprocessorVisitor::visitNode(Expression* const n) {
    if (ArithmeticExpression* const arex = dynamic_cast<ArithmeticExpression* const>(n))
        return visitNode(arex);
    else if (ParenthesisExpression* const pex = dynamic_cast<ParenthesisExpression* const>(n))
        return visitNode(pex);
    else if (LogicalExpression* const logex = dynamic_cast<LogicalExpression* const>(n))
        return visitNode(logex);
    else if (Macro* const m = dynamic_cast<Macro* const>(n))
        return visitNode(m);

    throw std::runtime_error("error: unknown expression class!");
}

int PreprocessorVisitor::visitNode(ExpressionList* const) {
    throw std::runtime_error("PreprocessorVisitor::visitNode for class ExpressionList is unimplemented!");
    return 0;
}

int PreprocessorVisitor::visitNode(ParenthesisExpression* const n) {
    if (n == 0)
        return 0;

    return visitNode(n->getInterior());
}

int PreprocessorVisitor::visitNode(ArithmeticExpression* const n) {
    if (IntConstant* const c = dynamic_cast<IntConstant* const>(n))
        return c->getValue();
    else if (BinaryExpression* const b = dynamic_cast<BinaryExpression* const>(n))
        return visitNode(b);
    else if (UnaryExpression* const u = dynamic_cast<UnaryExpression* const>(n))
        return visitNode(u);

    throw std::runtime_error("error: unknown arithmetic expression class!");
}

int PreprocessorVisitor::visitNode(IntConstant* const n) {
    return n->getValue();
}

int PreprocessorVisitor::visitNode(UnaryExpression* const n) {
    switch (n->getSymbolID()) {
        case PreprocessorTerminals::ID_PLUS:
            return visitNode(n->getExpression());

        case PreprocessorTerminals::ID_DASH:
            return (-visitNode(n->getExpression()));

        case PreprocessorTerminals::ID_OP_COMPLEMENT:
            return (~ visitNode(n->getExpression()));

        case PreprocessorTerminals::ID_OP_NOT:
            return (! visitNode(n->getExpression()));
    }

    throw std::runtime_error("unknown unary operator!");
}

int PreprocessorVisitor::visitNode(BinaryExpression* const n) {
    switch (n->getSymbolID()) {
        case PreprocessorTerminals::ID_BIT_OR:
            return (visitNode(n->getLeft()) | visitNode(n->getRight()));

        case PreprocessorTerminals::ID_BIT_XOR:
            return (visitNode(n->getLeft()) ^ visitNode(n->getRight()));

        case PreprocessorTerminals::ID_BIT_AND:
            return (visitNode(n->getLeft()) & visitNode(n->getRight()));

        case PreprocessorTerminals::ID_OP_LSHIFT:
            return (visitNode(n->getLeft()) << visitNode(n->getRight()));

        case PreprocessorTerminals::ID_OP_RSHIFT:
            return (visitNode(n->getLeft()) >> visitNode(n->getRight()));

        case PreprocessorTerminals::ID_PLUS:
            return (visitNode(n->getLeft()) + visitNode(n->getRight()));

        case PreprocessorTerminals::ID_DASH:
            return (visitNode(n->getLeft()) - visitNode(n->getRight()));

        case PreprocessorTerminals::ID_OP_MUL:
            return (visitNode(n->getLeft()) * visitNode(n->getRight()));

        case PreprocessorTerminals::ID_OP_DIV:
            return (visitNode(n->getLeft()) / visitNode(n->getRight()));

        case PreprocessorTerminals::ID_OP_MOD:
            return (visitNode(n->getLeft()) % visitNode(n->getRight()));
    }   // switch

    throw std::runtime_error("unknown binary operator!");
}

int PreprocessorVisitor::visitNode(LogicalExpression* const n) {
    if (LogicalBinaryExpression* const lbex = dynamic_cast<LogicalBinaryExpression* const>(n))
        return visitNode(lbex);
    else if (DefinedOperator* const defop = dynamic_cast<DefinedOperator* const>(n))
        return visitNode(defop);

    throw std::runtime_error("unknown logical expression class!");
}

int PreprocessorVisitor::visitNode(LogicalBinaryExpression* const n) {
    switch (n->getSymbolID()) {
        case PreprocessorTerminals::ID_LOGICAL_OR:
            return (visitNode(n->getLeft()) || visitNode(n->getRight()));

        case PreprocessorTerminals::ID_LOGICAL_AND:
            return (visitNode(n->getLeft()) && visitNode(n->getRight()));

        case PreprocessorTerminals::ID_EQ:
            return (visitNode(n->getLeft()) == visitNode(n->getRight()));

        case PreprocessorTerminals::ID_NEQ:
            return (visitNode(n->getLeft()) != visitNode(n->getRight()));

        case PreprocessorTerminals::ID_LESS:
            return (visitNode(n->getLeft()) < visitNode(n->getRight()));

        case PreprocessorTerminals::ID_LEQ:
            return (visitNode(n->getLeft()) <= visitNode(n->getRight()));

        case PreprocessorTerminals::ID_GEQ:
            return (visitNode(n->getLeft()) >= visitNode(n->getRight()));

        case PreprocessorTerminals::ID_GREATER:
            return (visitNode(n->getLeft()) > visitNode(n->getRight()));
    }

    throw std::runtime_error("error: unknown binary logical operator!");
}

int PreprocessorVisitor::visitNode(DefinedOperator* const n) {
    if (n == 0)
        return 0;

    PreprocessorSymbol* const symbol = symbols_.findSymbol(n->getIdentifier());
    return (((symbol != 0) && (symbol->isDefined())) ? 1 : 0);
}

int PreprocessorVisitor::visitNode(TextNode* const n) {
    if (n == 0)
        return 0;

    typedef std::map<std::string, PreprocessorSymbol*> SymbolsMap;
    const SymbolsMap& symbols = symbols_.getSymbolsMap();

    std::string text = n->getRawText();

    // Expand all macro invocations within the body until no further
    // expansions are performed.
    //
    // FIXME: recursive macro-definitions will cause this to hook up in
    // an infinite loop (and then crash).
    //
    for (bool expanded = true; expanded == true; ) {
        expanded = false;
        for (SymbolsMap::const_iterator it = symbols.begin(); it != symbols.end(); ++it) {
            if (it->second->isDefined()) {
                if (expandMacro(text, it->second) == true)
                    expanded = true;
            }
        }
    }   // for (it

    translation_ << text;
    return 1;
}

int PreprocessorVisitor::visitNode(TokenList* const) {
    throw std::runtime_error("PreprocessorVisitor::visitNode() for class TokenList is umplemented!");
    return 0;
}

int PreprocessorVisitor::visitNode(Macro* const n) {
    if (n == 0)
        return 0;

    PreprocessorSymbol* const symbol = symbols_.findSymbol(n->getIdentifier());
    if (symbol == 0) {
        log_ << "error: call to undefined macro '" << n->getIdentifier() << "'!\n";
        return 0;
    }

    std::list<Token*> macroBody;
    if (symbol->isFunction() && n->isFunction()) {
        if (n->getNumParameters() == symbol->getNumFormals())
            macroBody = n->expandMacro(symbol);
        else {
            log_ << "error: number of formal parameters does not match the number of actual parameters!\n";
            return 0;
        }
    } else {
        macroBody = symbol->getBody()->getCopy();
    }

    // Push this token to the front of the the macro's body / token stream
    // before creating a new temporary parser to make that parser only parse
    // expressions. This could be considered to be hack as the parser becomes
    // 'hijacked'.
    //
    macroBody.push_front(new Token(PreprocessorTerminals::ID_PIRATE));

    PreprocessorParser bodyParser(macroBody);
    ParseTreeNode* const root = bodyParser.parse();
    Expression* const exp = dynamic_cast<Expression* const>(root);
    if (exp != 0)
        return visitNode(exp);
    else {
        log_ << "failed to parse macro body:\n";
        log_ << bodyParser.getLog().str() << "\n";
    }
    delete root;

    return 0;
}

int PreprocessorVisitor::visitNode(Statement* const) {
    throw std::runtime_error("PreprocessorVisitor::visitNode() for class Statement is umplemented!");
    return 0;
}

int PreprocessorVisitor::visitNode(IdentifierList* const) {
    throw std::runtime_error("PreprocessorVisitor::visitNode() for class IdentifierList is umplemented!");
    return 0;
}

bool PreprocessorVisitor::visitNode(DefineDirective* const n) {
    if (n == 0)
        return false;

    PreprocessorSymbol* symbol = 0;
    if (n->getFormals() != 0) {
        const std::vector<std::string>& formals = n->getFormals()->getIdentifiers();
        symbol = new PreprocessorSymbol(n->getIdentifier(), n->isFunction(), n->getBody(), formals);
    } else
        symbol = new PreprocessorSymbol(n->getIdentifier(), n->isFunction(), n->getBody());

    bool res = symbols_.insertSymbol(symbol);

    // Within the translation, the directive is replaced by a newline to keep line
    // numbering consistent for the GLSL parser.
    //
    translation_ << "\n";

    return res;
}

void PreprocessorVisitor::visitNode(ErrorDirective* const) {
    // Within the translation, the directive is replaced by a newline to keep line
    // numbering consistent for the GLSL parser.
    //
    translation_ << "\n";
}

void PreprocessorVisitor::visitNode(ExtensionDirective* const) {
    // Within the translation, the directive is replaced by a newline to keep line
    // numbering consistent for the GLSL parser.
    //
    translation_ << "\n";
}

bool PreprocessorVisitor::visitNode(IncludeDirective* const n) {
    if (n == 0)
        return false;

    bool res = false;
    try {
        std::string includeFile = adjustIncludeFile(n->getFileName());
        PreprocessorParser includeParser(includeFile);

        log_ << "parsing included file '" << n->getFileName() << "'...\n";
        ParseTreeNode* root = includeParser.parse();
        if (root != 0) {
            visitAll(root);
            res = true;
        } else {
            log_ << "failed to parse include directive '" << n->getFileName() << "':\n";
            log_ << includeParser.getLog().str() << "\n";
        }

        delete root;
    } catch (std::bad_alloc&) {
        log_ << "failed to open file '" << n->getFileName() << "' for inclusion!\n";
    }

    return res;
}

int PreprocessorVisitor::visitNode(LineDirective* const n) {
    if (n == 0)
        return 0;

    // Within the translation, the directive is replaced by a newline to keep line
    // numbering consistent for the GLSL parser.
    //
    translation_ << "\n";

    return n->getLineNumber();
}

void PreprocessorVisitor::visitNode(PragmaDirective* const) {
    // Within the translation, the directive is replaced by a newline to keep line
    // numbering consistent for the GLSL parser.
    //
    translation_ << "\n";
}

bool PreprocessorVisitor::visitNode(UndefineDirective* const n) {
    if (n == 0)
        return false;

    PreprocessorSymbol* const symbol = symbols_.findSymbol(n->getIdentifier());
    if ((symbol != 0) && (symbol->isDefined())) {
        symbol->undefine();
        return true;
    }

    return false;
}

int PreprocessorVisitor::visitNode(VersionDirective* const n) {
    if (n == 0)
        return 0;

    // Within the translation, the directive is replaced by a newline to keep line
    // numbering consistent for the GLSL parser.
    //
    translation_ << "\n";

    return n->getVersion();
}

void PreprocessorVisitor::visitNode(ConditionalDirective* const n) {
    if (n == 0)
        return;

    ParseTreeNode* next = 0;
    Expression* const cond = n->getCondition();
    next = (visitNode(cond) != 0) ? n->getTrue() : n->getFalse();

    if (n->getTrue() == 0)
        log_ << "warning: empty if statement!\n";

    // add newline to translation for consistency for the GLSL parser
    //
    translation_ << "\n";

    if (next != 0) {
        const std::vector<ParseTreeNode*>& children = next->getChildren();
        for (size_t i = 0; (next != 0); next = children[i++]) {
            visit(next);

            if (i >= children.size())
                break;
        }
    }
}

std::string PreprocessorVisitor::adjustIncludeFile(const std::string& fileName) {
    std::string fullPath = includePath_;
#ifdef WIN32
    if ((! fullPath.empty()) && (fullPath[fullPath.size() - 1] != '\\'))
        fullPath += "\\";

    size_t offset = fullPath.size();
    fullPath += fileName;
    for (size_t i = offset; i < fullPath.size(); ++i) {
        if (fullPath[i] == '/')
            fullPath[i] = '\\';
    }
#else
    if ((! fullPath.empty()) && (fullPath[fullPath.size() - 1] != '/'))
        fullPath += "/";
    fullPath += fileName;
#endif
    return fullPath;
}

bool PreprocessorVisitor::expandMacro(std::string& input, PreprocessorSymbol* const symbol) const
{
    if (symbol == 0)
        return false;

    const std::string& name = symbol->getID();

    bool expanded = false;
    for (size_t pos = input.find(name, 0), lastPos = 0; pos != std::string::npos;
        pos = input.find(name, lastPos))
    {
        std::string body = symbol->getBody()->toString(terminals_);
        expanded = true;
        lastPos = pos + name.size();

        if (symbol->isFunction()) {
            const std::vector<std::string>& formals = symbol->getFormals();
            std::vector<std::string> parameters = readParameters(input, lastPos);

            if (symbol->getNumFormals() == parameters.size()) {
                for (size_t i = 0; i < parameters.size(); ++i)
                    body = replace(body, formals[i], parameters[i]);
            }

            std::string macroInvocation = input.substr(pos, lastPos - pos);
            input = replace(input, macroInvocation, body);
        } else {
            std::string macroInvocation = input.substr(pos, lastPos - pos);
            input = replace(input, macroInvocation, body);
            break;
        }
    }

    return expanded;
}

std::vector<std::string> PreprocessorVisitor::readParameters(const std::string& input,
                                                             size_t& offset) const
{
    std::vector<std::string> params;

    bool readingParameter = false;
    size_t startPos = std::string::npos;
    for (size_t i = offset; i < input.size(); ++i) {
        const char top = input[i];

        if (Lexer::isWhitespace(top)) {
            if (! readingParameter)
                ++startPos;
            continue;
        }

        if ((top == '(') && (! readingParameter)) {
            startPos = i + 1;
            continue;
        }

        if (startPos != std::string::npos) {
            if (top == ',') {
                params.push_back(input.substr(startPos, (i - startPos)));
                readingParameter = false;
                startPos = i + 1;
            } else if (top == ')') {
                if (i > startPos)
                    params.push_back(input.substr(startPos, (i - startPos)));
                offset = i + 1;
                break;
            } else
                readingParameter = true;
        }
    }   // for

    return params;
}

std::string PreprocessorVisitor::replace(const std::string& in, const std::string& search,
                                         const std::string& replacement) const
{

    std::ostringstream out;

    size_t lastPos = 0;
    for (size_t pos = in.find(search, 0); pos != std::string::npos;
        pos = in.find(search, lastPos))
    {
        out << in.substr(lastPos, pos - lastPos);
        out << replacement;

        lastPos = pos + search.size();
    }

    if (lastPos < in.size())
        out << in.substr(lastPos);

    return out.str();
}

}   // namespace glslparser

}   // namespace voreen
