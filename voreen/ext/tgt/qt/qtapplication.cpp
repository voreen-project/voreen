/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2021 University of Muenster, Germany,           *
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

#include <qapplication.h>
#include "tgt/qt/qtapplication.h"
#include "tgt/qt/qtcanvas.h"

namespace tgt {

QtApplication::QtApplication(int& argc, char** argv)
    : GUIApplication(argc, argv)
{
    qapplication_ = new QApplication(*argc_, argv_);
}

QtApplication::~QtApplication(){
    deinitGL();
//FIXME: Handled by Qt or not?
//    delete qapplication_;
    deinit();
}

int QtApplication::init(InitFeature::Features featureset) {
    return GUIApplication::init(featureset);
}

int QtApplication::run() {
    if (hasRun) {
        std::cout << "Failure: QtApplication::run() should be called only once.\n";
        return -1;
    }
    hasRun = true;
    return qapplication_->exec();
}

void QtApplication::quit(int status) {
    qapplication_->exit(status);
}

}
