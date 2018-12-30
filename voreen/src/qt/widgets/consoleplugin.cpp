/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreen/qt/widgets/consoleplugin.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/commandqueue.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "tgt/logmanager.h"

#include <QApplication>
#include <QThread>
#include <QTextEdit>
#include <QVBoxLayout>
#include <QScrollBar>
#include <QMenu>
#include <QSettings>

namespace voreen {

class ConsoleLogQt : public tgt::Log {
public:
    ConsoleLogQt(ConsolePlugin* plugin, std::string debugStyle = "", std::string infoStyle = "",
        std::string warnStyle = "", std::string errorStyle = "",
        bool timeStamping = false, bool dateStamping = false,
        bool showCat = true, bool showLevel = true)
        : tgt::Log(),
          plugin_(plugin),
          debugStyle_(debugStyle),
          infoStyle_(infoStyle),
          warnStyle_(warnStyle),
          errorStyle_(errorStyle)
    {
        timeStamping_ = timeStamping;
        dateStamping_ = dateStamping;
        showCat_ = showCat;
        showLevel_ = showLevel;
    }

    ~ConsoleLogQt() {
        VoreenApplication::app()->getCommandQueue()->removeAll(this);
    }

    bool isOpen() { return true; }

protected:
    void logFiltered(const std::string &cat, tgt::LogLevel level, const std::string &msg, const std::string & /*extendedInfo*/ ="") {

        // Create convenience log function.
        auto logFunc = [this, cat, level, msg] {
            std::string output;

            std::string style;
            if (level == tgt::Debug)
                style = debugStyle_;
            else if (level == tgt::Info)
                style = infoStyle_;
            else if (level == tgt::Warning)
                style = warnStyle_;
            else if (level == tgt::Error)
                style = errorStyle_;

            if (dateStamping_)
                output += "[" + getDateString() + "] ";
            if (timeStamping_)
                output += "[" + getTimeString() + "] ";
            if (showCat_)
                output += cat + " ";
            if (showLevel_)
                output += "(" + getLevelString(level) + ")";

            output += "&nbsp;&nbsp;";

            output += "<span style=\"" + style + "\">";

            size_t lookHere = 0;
            size_t foundHere;
            std::string s = msg;
            const std::string br = "<br/>";
            while ((foundHere = s.find("\n", lookHere)) != std::string::npos) {
                s.replace(foundHere, 1, br);
                lookHere = foundHere + br.size();
            }

            output += s;
            output += "</span>";

            plugin_->log(output);
        };

        // accessing Qt widgets is only allowed in the GUI thread
        bool isGuiThread = (QThread::currentThread() == QCoreApplication::instance()->thread());
        if (isGuiThread) // log immediatly
            logFunc();
        else // enqueue to command queue being executed in GUI thread
            VoreenApplication::app()->getCommandQueue()->enqueue(this, LambdaFunctionCallback(logFunc));
    }

    ConsolePlugin* plugin_;
    std::string debugStyle_;
    std::string infoStyle_;
    std::string warnStyle_;
    std::string errorStyle_;
};

ConsolePlugin::ConsolePlugin(QWidget* parent, tgt::LogLevel logLevel, bool autoScroll)
    : QWidget(parent)
    , autoScroll_(autoScroll)
{
    setObjectName(tr("Console"));

    QSettings settings;
    settings.beginGroup("ConsolePlugin");
    bool disabled = settings.value("disabled", false).toBool();
    settings.endGroup();

    consoleText_ = new QTextEdit(this);
    consoleText_->setReadOnly(true);
    consoleText_->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    consoleText_->setContextMenuPolicy(Qt::CustomContextMenu);
    this->connect(consoleText_, SIGNAL(customContextMenuRequested(const QPoint&)), this, SLOT(showContextMenu(const QPoint &)));

    clearText_ = new QAction("Erase All", consoleText_);
    connect(clearText_, SIGNAL(triggered()), consoleText_, SLOT(clear()));
    clearText_->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_E));
    ctrlButtonDown_ = false;

    disableAction_ = new QAction("Disable", this);
    disableAction_->setCheckable(true);
    disableAction_->setChecked(disabled);
    connect(disableAction_, SIGNAL(triggered()), this, SLOT(disableToggled()));

    QVBoxLayout* vboxLayout = new QVBoxLayout();
    vboxLayout->addWidget(consoleText_);
    setLayout(vboxLayout);

    //< only add log, if console log is already present (otherwise assuming that logging is disabled entirely)
    if (tgt::LogManager::isInited() && LogMgr.getConsoleLog()) {
        log_ = new ConsoleLogQt(this, "", "", "color: brown; font-weight: bold", "color: red; font-weight: bold");
        log_->addCat("", true, logLevel);
        LogMgr.addLog(log_);
    }

    disableToggled();
}

ConsolePlugin::~ConsolePlugin() {

    QSettings settings;
    settings.beginGroup("ConsolePlugin");
    settings.setValue("disabled", disableAction_->isChecked());
    settings.endGroup();

    if (tgt::LogManager::isInited()) {
        LogMgr.removeLog(log_);
        delete log_;
    }
}

void ConsolePlugin::keyPressEvent(QKeyEvent *e)
{
    switch(e->key()){
        case Qt::Key_Control:
            ctrlButtonDown_ = true;
            break;
        case Qt::Key_E:
            if(ctrlButtonDown_)
                consoleText_->clear();
            break;
    }
}

void ConsolePlugin::keyReleaseEvent(QKeyEvent *e)
{
    if(e->key() == Qt::Key_Control)
        ctrlButtonDown_ = false;
}

void ConsolePlugin::showContextMenu(const QPoint &pt)
{
    QMenu *menu = consoleText_->createStandardContextMenu();
    menu->addSeparator();
    menu->addAction(clearText_);
    menu->addSeparator();
    menu->addAction(disableAction_);
    menu->exec(consoleText_->mapToGlobal(pt));
    delete menu;
}

void ConsolePlugin::disableToggled(){
    // would also disable the context menu
    //consoleText_->setEnabled(!disableAction_->isChecked());
}

void ConsolePlugin::log(const std::string& msg) {

    if (disableAction_->isChecked())
        return;

    // write log message to text box
    consoleText_->append(msg.c_str());

    // scroll to bottom
    if (autoScroll_ && isVisible()) {
        QScrollBar* scrollBar = consoleText_->verticalScrollBar();
        if (scrollBar)
            scrollBar->setValue(scrollBar->maximum());
    }
}

} // namespace voreen
