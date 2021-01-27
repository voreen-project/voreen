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

// This class is based on the GLSL syntax highlighter written by Markus Kramer for
// Shader Maker - a cross-platform GLSL editor

#include "voreen/qt/widgets/syntaxhighlighter.h"

#include <QTextDocument>

namespace voreen {

SyntaxHighlighter::SyntaxHighlighter(QTextDocument* document)
    : QSyntaxHighlighter(document)
{}

SyntaxHighlighter::~SyntaxHighlighter() {}

void SyntaxHighlighter::highlightBlock(const QString& text) {
    if (text == "")
        return;

    // check out all rules
    foreach (highlightRule_t rule, m_rules) {
        QRegExp expression(rule.pattern);
        int index = text.indexOf(expression);

        while (index >= 0) {
            int length = expression.matchedLength();
            setFormat(index, length, rule.format);
            index = text.indexOf(expression, index + length);
        }
    }

    setCurrentBlockState(0);

    int startIndex = 0;
    if (previousBlockState() != 1)
        startIndex = text.indexOf(m_commentStartExpression);

    while (startIndex >= 0) {
        int endIndex = text.indexOf(m_commentEndExpression, startIndex);
        int commentLength;

        if (endIndex == startIndex)
            break;

        if (endIndex == -1) {
            setCurrentBlockState(1);
            commentLength = text.length() - startIndex;
        } else
            commentLength = endIndex - startIndex + m_commentEndExpression.matchedLength();

        setFormat(startIndex, commentLength, m_multiLineCommentFormat);
        startIndex = text.indexOf(m_commentStartExpression, startIndex + commentLength);
    }
}

} // namespace
