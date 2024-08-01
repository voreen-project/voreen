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

// This class is based on the GLSL syntax highlighter written by Markus Kramer for
// Shader Maker - a cross-platform GLSL editor

#ifndef VRN_SYNTAXHIGHLIGHTER_H
#define VRN_SYNTAXHIGHLIGHTER_H

#include <QSyntaxHighlighter>

#include "voreen/qt/voreenqtapi.h"

class QTextDocument;

namespace voreen {

class VRN_QT_API SyntaxHighlighter : public QSyntaxHighlighter {
public:
    SyntaxHighlighter(QTextDocument* document);

    virtual ~SyntaxHighlighter();

protected:
    void highlightBlock(const QString& text);

    // rules
    typedef struct highlightRule_s {
        QRegExp pattern;
        QTextCharFormat format;
    } highlightRule_t;
    QList<highlightRule_t> m_rules;

    // comment expressions
    QRegExp m_commentStartExpression;
    QRegExp m_commentEndExpression;

    // format infos
    QTextCharFormat m_keywordFormat;
    QTextCharFormat m_reservedKeywordFormat;
    QTextCharFormat m_integerFormat;
    QTextCharFormat m_floatFormat;
    QTextCharFormat m_registerFormat;
    QTextCharFormat m_builtInFunctionFormat;
    QTextCharFormat m_preprocessorFormat;
    QTextCharFormat m_singleLineCommentFormat;
    QTextCharFormat m_multiLineCommentFormat;
};

} // namespace

#endif // VRN_SYNTAXHIGHLIGHTER_H
