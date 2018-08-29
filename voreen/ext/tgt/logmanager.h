/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#ifndef TGT_LOGMANAGER_H
#define TGT_LOGMANAGER_H

#include <string>
#include <sstream>
#include <vector>

#include "tgt/assert.h"
#include <stdarg.h>
#include "tgt/singleton.h"
#include "tgt/types.h"

namespace tgt {

/**
 * Specifies the severity of the log event.
 * Debug messages are not logged in release builds!
 */
enum LogLevel {
    Debug,
    Info,
    Warning,
    Error,
    Fatal
};

/**
 * Holds the information for filtering messages.
 */
struct TGT_API LogFilter {
    std::string cat_;
    bool children_;
    LogLevel level_;
};

/**
 * Abstract basis class for logging messages.
 */
class TGT_API Log {
public:
    virtual ~Log() {}

    /// Log a message in this Log (message is filtered based on cat and level)
    virtual void log(const std::string &cat, LogLevel level, const std::string &msg, const std::string &extendedInfo="");

    /**
     * Add a category that is accepted to this log.
     * @param cat All messages with category = cat are accepted and logged in this log.
     * @param Children If true all messages of subcategories of cat are also accepted.
     * example:
     * tgt.Texture.Reader is a subcategory of tgt.Texture
     * @param level All messages below this LogLevel are discarded, even if the category matches.
     */
    virtual void addCat(const std::string &cat, bool Children = true, LogLevel level = Debug);
    virtual bool isOpen() = 0;

    /// Sets the log level of all filters registered at this Log.
    virtual void setLogLevel(LogLevel level);

    /// Returns if the messages are time-stamped.
    inline bool getTimeStamping() const { return timeStamping_; }
    inline void setTimeStamping(const bool timeStamping) { timeStamping_ = timeStamping; }

    /// Returns if the messages are date-stamped.
    inline bool getDateStamping() const { return dateStamping_; }
     inline void setDateStamping(const bool dateStamping) { dateStamping_ = dateStamping; }

    /// Returns if the category is printed along with the messages
    inline bool getShowCat() const { return showCat_; }
    inline void setShowCat(const bool showCat) { showCat_ = showCat; }

    /// Returns if the LogLevel is printed along with the messages.
    inline bool getShowLevel() const { return showLevel_; }
    inline void setShowLevel(const bool showLevel) { showLevel_ = showLevel; }

protected:
    virtual bool testFilter(const std::string &cat, LogLevel level);
    virtual void logFiltered(const std::string &cat, LogLevel level, const std::string &msg,
                             const std::string &extendedInfo = "") = 0;
    virtual std::string getTimeString();
    virtual std::string getDateString();
    virtual std::string getLevelString(LogLevel level);

    bool timeStamping_, dateStamping_, showCat_, showLevel_;
    std::vector<LogFilter> filters_;
};


/// Implements logging to a plain Textfile
class TGT_API TextLog : public Log {
public:
    TextLog(const std::string &filename, bool dateStamping = true, bool timeStamping= true,
            bool showCat = true, bool showLevel = true);
    virtual ~TextLog();
    bool isOpen();

protected:
    FILE* file_;
    void logFiltered(const std::string &cat, LogLevel level, const std::string &msg, const std::string &extendedInfo="");
    std::string getLevelColor(LogLevel level);
};

/**
 * Implements logging to stdout.
 * You can only have one ConsoleLog in the LogMgr
 */
class TGT_API ConsoleLog : public Log {
public:
    ConsoleLog(bool dateStamping = false, bool timeStamping= false, bool showCat = true, bool showLevel = true);
    ~ConsoleLog() {}
    bool isOpen() { return true; }
    void enableColors(bool enable = true);

protected:
    void logFiltered(const std::string &cat, LogLevel level, const std::string &msg, const std::string &extendedInfo="");
    std::string getLevelColor(LogLevel level);

    bool colorOutput_;
};

///Implements a colored html log.
class TGT_API HtmlLog : public Log {
public:
    HtmlLog(const std::string &filename, bool dateStamping = false, bool timeStamping= true,
            bool showCat = true, bool showLevel = true, bool append = true);
    ~HtmlLog();
    bool isOpen();

    std::string getAbsFilename() const;

protected:
    void logFiltered(const std::string &cat, LogLevel level, const std::string &msg, const std::string &extendedInfo="");
    std::string getLevelString(LogLevel level);
    std::string getLevelColor(LogLevel level);

    FILE* file_;
    std::string absFilename_;
};


class LogManager;
#ifdef DLL_TEMPLATE_INST
template class TGT_API Singleton<LogManager>;
#endif

/**
 * The Logmanager distributes logmessages to all Logs registered to the manager.
 * Logmessages consist of a message, a logging category and a loglevel.
 * See Log and LogLevel classes for further information.
 *
 * Usage:
 * Define "static const std::string loggerCat_;" in every class that logs messages and use the macros defined below:
 * LWARNING("Warning!");
 * Alternatively, LWARNINGC("Cat", "Warning!") may be used, which does not require the definition of loggerCat_.
 *
 * LDEBUG statements are removed if TGT_DEBUG is not defined!
 * @author Stefan Diepenbrock
 */
class TGT_API LogManager : public Singleton<LogManager> {
public:
    LogManager(const std::string& logDir = "");
    ~LogManager();

    /// Initialize logmanager, put all logfiles in logDir
    void reinit(const std::string& logDir);
    std::string getLogDir() const { return logDir_; }

    /// Log message
    void log(const std::string& cat, LogLevel level, const std::string& msg, const std::string& extendedInfo="");

    /// Add a log to the manager, from now all messages received by the manager are also distributed to this log.
    /// All logs are deleted upon destruction of the manager.
    /// If a ConsoleLog is added it will replace an existing one, the old one will be deleted.
    void addLog(Log* log);

    // Remove a log from the manager.
    void removeLog(Log* log);

    /// Removes all registered logs.
    void clear();

    /// Sets the log level of all registered Logs.
    void setLogLevel(LogLevel level);

    /// Return the ConsoleLog (or 0 if there is none)
    ConsoleLog* getConsoleLog() { return consoleLog_; }

    /// Returns all Logs currently registered at the LogManager.
    std::vector<Log*> getLogs() const;

protected:
    std::string logDir_;
    std::vector<Log*> logs_;
    ConsoleLog* consoleLog_;

};

} // namespace

#define LogMgr tgt::LogManager::getRef()

// Use "do { ... } while (0)" to allow "if (foo) LINFO("bar"); else ...", which would fail
// otherwise.
// Compare: http://gcc.gnu.org/onlinedocs/cpp/Swallowing-the-Semicolon.html

#ifdef TGT_DEBUG
    #ifdef __GNUC__
        #define LDEBUG(msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __PRETTY_FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(loggerCat_, tgt::Debug, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LINFO(msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __PRETTY_FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(loggerCat_, tgt::Info, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LWARNING(msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __PRETTY_FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(loggerCat_, tgt::Warning, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LERROR(msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __PRETTY_FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(loggerCat_, tgt::Error, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LFATAL(msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __PRETTY_FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(loggerCat_, tgt::Fatal, _tmp.str(), _tmp2.str()); \
        } while (0)

        //with category parameter:
        #define LDEBUGC(cat, msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __PRETTY_FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(cat, tgt::Debug, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LINFOC(cat, msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __PRETTY_FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(cat, tgt::Info, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LWARNINGC(cat, msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __PRETTY_FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(cat, tgt::Warning, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LERRORC(cat, msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __PRETTY_FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(cat, tgt::Error, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LFATALC(cat, msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __PRETTY_FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(cat, tgt::Fatal, _tmp.str(), _tmp2.str()); \
        } while (0)
    #else
        #define LDEBUG(msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(loggerCat_, tgt::Debug, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LINFO(msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(loggerCat_, tgt::Info, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LWARNING(msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(loggerCat_, tgt::Warning, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LERROR(msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(loggerCat_, tgt::Error, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LFATAL(msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(loggerCat_, tgt::Fatal, _tmp.str(), _tmp2.str()); \
        } while (0)

        //with category parameter:
        #define LDEBUGC(cat, msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(cat, tgt::Debug, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LINFOC(cat, msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(cat, tgt::Info, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LWARNINGC(cat, msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(cat, tgt::Warning, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LERRORC(cat, msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(cat, tgt::Error, _tmp.str(), _tmp2.str()); \
        } while (0)

        #define LFATALC(cat, msg) \
        do { \
            std::ostringstream _tmp, _tmp2; \
            _tmp2 << __FUNCTION__  << " File: " << __FILE__ << "@" << __LINE__;\
            _tmp << msg; \
            LogMgr.log(cat, tgt::Fatal, _tmp.str(), _tmp2.str()); \
        } while (0)
     #endif
#else

    //#define LDEBUG(msg)
    #define LDEBUG(msg) \
    do { \
        std::ostringstream _tmp; \
        _tmp << msg; \
        LogMgr.log(loggerCat_, tgt::Debug, _tmp.str()); \
    } while (0)

    #define LINFO(msg) \
    do { \
        std::ostringstream _tmp; \
        _tmp << msg; \
        LogMgr.log(loggerCat_, tgt::Info, _tmp.str()); \
    } while (0)

    #define LWARNING(msg) \
    do { \
        std::ostringstream _tmp; \
        _tmp << msg; \
        LogMgr.log(loggerCat_, tgt::Warning, _tmp.str()); \
    } while (0)

    #define LERROR(msg) \
    do { \
        std::ostringstream _tmp; \
        _tmp << msg; \
        LogMgr.log(loggerCat_, tgt::Error, _tmp.str()); \
    } while (0)

    #define LFATAL(msg) \
    do { \
        std::ostringstream _tmp; \
        _tmp << msg; \
        LogMgr.log(loggerCat_, tgt::Fatal, _tmp.str()); \
    } while (0)

    //
    // with category parameter
    //

    //#define LDEBUGC(cat, msg)
    #define LDEBUGC(cat, msg) \
    do { \
        std::ostringstream _tmp; \
        _tmp << msg; \
        LogMgr.log(cat, tgt::Debug, _tmp.str()); \
    } while (0)

    #define LINFOC(cat, msg) \
    do { \
        std::ostringstream _tmp; \
        _tmp << msg; \
        LogMgr.log(cat, tgt::Info, _tmp.str()); \
    } while (0)

    #define LWARNINGC(cat, msg) \
    do { \
        std::ostringstream _tmp; \
        _tmp << msg; \
        LogMgr.log(cat, tgt::Warning, _tmp.str()); \
    } while (0)

    #define LERRORC(cat, msg) \
    do { \
        std::ostringstream _tmp; \
        _tmp << msg; \
        LogMgr.log(cat, tgt::Error, _tmp.str()); \
    } while (0)

    #define LFATALC(cat, msg) \
    do { \
        std::ostringstream _tmp; \
        _tmp << msg; \
        LogMgr.log(cat, tgt::Fatal, _tmp.str()); \
    } while (0)

#endif //TGT_DEBUG

#endif //TGT_LOGMANAGER_H
