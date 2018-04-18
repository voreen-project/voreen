/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 Visualization and Computer Graphics Group, *
 * Department of Computer Science, University of Muenster, Germany.   *
 * <http://viscg.uni-muenster.de>                                     *
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

#ifndef TGT_QTAPPLICATION_H
#define TGT_QTAPPLICATION_H

#include <qapplication.h>

#include "tgt/guiapplication.h"
#include "tgt/qt/qttimer.h"
#include "tgt/qt/qtcanvas.h"

namespace tgt {

class QtApplication : public GUIApplication {
public:
    /// Constructor
    /// @param argc number of command line arguments
    /// @param vector of command line arguments
    QtApplication(int& argc, char** argv);

    /// Destructor
    virtual ~QtApplication();

    /// initialize tgt and Qt, create a Qt canvas for every QtCanvas and stuff ...
    virtual int init(InitFeature::Features featureset = InitFeature::ALL) ;

    /// enter qt mainloop (QApplication::exec())
    virtual int run();

    /// tells qt to quit mainloop and exit program
    virtual void quit(int status = 0);

    /// return the Qt::QApplication object of our tgt::QtApplication
    virtual QApplication* getQApplication() { return qapplication_; };

protected:

    QApplication* qapplication_;

};

} // namespace tgt

#endif // TGT_QTAPPLICATION_H
