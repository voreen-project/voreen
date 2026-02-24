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
    foreach (const highlightRule_t& rule, m_rules) {
        QRegularExpressionMatchIterator matchIterator = rule.pattern.globalMatch(text);
        while (matchIterator.hasNext()) {
            const QRegularExpressionMatch match = matchIterator.next();
            setFormat(match.capturedStart(), match.capturedLength(), rule.format);
        }
    }

    setCurrentBlockState(0);

    int startIndex = 0;
    if (previousBlockState() != 1) {
        const QRegularExpressionMatch startMatch = m_commentStartExpression.match(text);
        startIndex = startMatch.capturedStart();
    }

    while (startIndex >= 0) {
        const QRegularExpressionMatch endMatch = m_commentEndExpression.match(text, startIndex);
        int endIndex = endMatch.capturedStart();
        int commentLength;

        if (endIndex == startIndex)
            break;

        if (endIndex == -1) {
            setCurrentBlockState(1);
            commentLength = text.length() - startIndex;
        } else
            commentLength = endIndex - startIndex + endMatch.capturedLength();

        setFormat(startIndex, commentLength, m_multiLineCommentFormat);
        const QRegularExpressionMatch nextStartMatch = m_commentStartExpression.match(text, startIndex + commentLength);
        startIndex = nextStartMatch.capturedStart();
    }
}

} // namespace
