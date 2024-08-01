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

#include "voreen/qt/widgets/codeedit.h"

#include <QPainter>
#include <QPaintEvent>
#include <QTextBlock>
#include <Qt>

//#include <regex>
#include <boost/regex.hpp>
#include <sstream>

#include "tgt/assert.h"

CodeEdit::CodeEdit(bool readonly, int fontSize, QWidget* parent)
    : QTextBrowser(parent)
    , textEditLock_(false)
{
    statusArea_ = new StatusArea(this);

    auto font = currentFont();
    font.setFamily("Courier New");
    font.setStyleHint(QFont::TypeWriter);
    font.insertSubstitution("Courier New", "monospace");
    font.setFixedPitch(true);
    font.setPointSize(fontSize);
    setFont(font);
    QFontMetrics metrics(font);
    setTabStopWidth(metrics.width(" ")*4);

    // TODO Are these even necessary for QTextBrowser?
    connect(this, SIGNAL(textChanged()), this, SLOT(updateStatusAreaWidth()));
    //connect(this, SIGNAL(updateRequest(const QRect &, int)), this, SLOT(updateStatusArea(const QRect &, int)));
    connect(this, SIGNAL(cursorPositionChanged()), this, SLOT(highlightCurrentLine()));
    connect(this, SIGNAL(textChanged()), this, SLOT(processCurrentLine()));
    connect(this, SIGNAL(textChanged()), statusArea_, SLOT(update()));

    updateStatusAreaWidth();
    highlightCurrentLine();

    setWordWrapMode(QTextOption::NoWrap);
    setOpenLinks(false);
    if(readonly) {
        setTextInteractionFlags(Qt::TextBrowserInteraction | Qt::TextSelectableByKeyboard);
    } else {
        setTextInteractionFlags(Qt::TextBrowserInteraction | Qt::TextEditable | Qt::TextSelectableByKeyboard);
    }

}

CodeEdit::~CodeEdit() {
    delete statusArea_;
}

int CodeEdit::statusAreaWidth() {
    int digits = 1;
    int max = qMax(1, document()->lineCount());
    while (max >= 10) {
        max /= 10;
        ++digits;
    }

    int space = 3 + fontMetrics().width(QLatin1Char('9')) * digits;

    return space;
}

void CodeEdit::updateStatusAreaWidth() {
    setViewportMargins(statusAreaWidth(), 0, 0, 0);
}

void CodeEdit::updateFontSize(unsigned char s) {
    QFont newFont = font();
    newFont.setPointSize(s);
    setFont(newFont);
    QFontMetrics metrics(newFont);
    setTabStopWidth(metrics.width(" ")*4);
    updateStatusAreaWidth();
}

static std::string processLine(const std::string& line) {
    boost::regex includeStatementRegex("#include((?:&nbsp;)*)\"(.*)\"");
    boost::regex spaceRegex(" ");
    boost::regex utf8nonbreakingSpaceRegex("\302\240");
    boost::regex lBracketRegex("<");
    boost::regex rBracketRegex(">");
    std::string output = line;
    output = boost::regex_replace(output, spaceRegex, "&nbsp;");
    output = boost::regex_replace(output, utf8nonbreakingSpaceRegex, "&nbsp;");
    output = boost::regex_replace(output, lBracketRegex, "&#060;");
    output = boost::regex_replace(output, rBracketRegex, "&#062;");
    output = boost::regex_replace(output, includeStatementRegex, "#include$1\"<a href=\"$2\">$2</a>\"");
    return output;
}

void CodeEdit::setShaderSource(const std::string& source) {
    std::istringstream sourcestream(source);

    std::string line;

    textEditLock_ = true;

    clear();
    // Hack: Select all text (which will consequently deleted) to
    // get rid of html-hyperlinks that not have been deleted by
    // clear(). Qt-bug?
    selectAll();
    while(std::getline(sourcestream, line)) {
        insertHtml(QString(processLine(line).c_str()));
        insertPlainText("\n");
    }
    textEditLock_ = false;
    // Reset cursor position
    QTextCursor cursor = textCursor();
    cursor.movePosition(QTextCursor::MoveOperation::Start);
    setTextCursor(cursor);
}

void CodeEdit::processCurrentLine() {
    if(textEditLock_) {
        return;
    }
    textEditLock_ = true;

    QTextCursor origCursor = textCursor();
    int origPos = origCursor.position();

    QTextCursor cursor = textCursor();
    cursor.select(QTextCursor::LineUnderCursor);
    std::string currentLine = cursor.selectedText().toStdString();
    setTextCursor(cursor);
    std::string processedLine = processLine(currentLine);
    insertHtml(QString(processedLine.c_str()));

    // Somehow origCursors position changes due to the operations above.
    // => We set the previously saved position and then set the original cursor.
    origCursor.setPosition(origPos);
    setTextCursor(origCursor);

    textEditLock_ = false;
}

void CodeEdit::resizeEvent(QResizeEvent* e) {
    QTextBrowser::resizeEvent(e);

    QRect cr = contentsRect();
    statusArea_->setGeometry(QRect(cr.left(), cr.top(), statusAreaWidth(), cr.height()));
}

void CodeEdit::paintEvent(QPaintEvent* e) {
    QTextBrowser::paintEvent(e);

    statusArea_->update();
}

void CodeEdit::wheelEvent(QWheelEvent* event) {

    // Block zoom since it would break font size.
    if (event->modifiers() & Qt::ControlModifier) {
        event->setModifiers(event->modifiers() & ~Qt::ControlModifier);
    }

    QTextBrowser::wheelEvent(event);
}

void CodeEdit::highlightCurrentLine() {
    QList<QTextEdit::ExtraSelection> extraSelections;

    if (isEnabled() && !isReadOnly()) {
        QTextEdit::ExtraSelection selection;

        QColor lineColor = QColor(Qt::yellow).lighter(170);

        selection.format.setBackground(lineColor);
        selection.format.setProperty(QTextFormat::FullWidthSelection, true);
        selection.cursor = textCursor();
        selection.cursor.clearSelection();
        extraSelections.append(selection);
    }

    setExtraSelections(extraSelections);
}

void CodeEdit::statusAreaPaintEvent(QPaintEvent* event) {
    QPainter painter(statusArea_);
    painter.fillRect(event->rect(), Qt::lightGray);
    painter.setPen(Qt::black);
    // Damn you, Qt! See:
    // http://www.qtcentre.org/threads/32552-Changing-the-font-of-a-QPainter
    QFont font;
    bool success = font.fromString(currentFont().toString());
    font.setUnderline(false);
    tgtAssert(success, "Failed to set line number font from string");
    painter.setFont(font);
    tgtAssert(painter.font().pointSize() == font.pointSize(), "Font size did not propagate.");

    QFontMetrics fontMetrics(font);

    int height = fontMetrics.height();
    int leading = fontMetrics.leading();

    int lineSpacing = height + std::max(0,leading); //fix for linux font size 24
    int currentLine = textCursor().blockNumber();
    int startPos = cursorRect().top() - currentLine * lineSpacing;
    int areaWidth = statusArea_->width();
    for(int line=0; line < document()->lineCount(); ++line) {
        QString number = QString::number(line+1);
        int pos = startPos + line*lineSpacing;
        painter.drawText(0, pos, areaWidth, lineSpacing,
                         Qt::AlignRight, number);
    }
}

void CodeEdit::moveCursorToPosition(int line, int col) {
    QTextCursor cursor = textCursor();
    cursor.setPosition(0);
    cursor.movePosition(QTextCursor::NextBlock, QTextCursor::MoveAnchor, line);
    if (col > 0)
        cursor.movePosition(QTextCursor::Right, QTextCursor::MoveAnchor, col);
    setTextCursor(cursor);
}

void CodeEdit::updateHighlight() {
    highlightCurrentLine();
}
