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

#include "voreen/core/utils/observer.h"

namespace voreen {

Observer::~Observer() {
    clearObserveds();
}

bool Observer::observes(const ObservableBase* observable) const {
    return (observeds_.find(const_cast<ObservableBase*>(observable)) != observeds_.end());
}

void Observer::clearObserveds() const {
    while (!observeds_.empty())
        removeObserved(*observeds_.begin());
}

bool Observer::addObserved(const ObservableBase* observed) const {

    if (!observed)
        return false;

    std::set<const ObservableBase*>::const_iterator it = observeds_.find(observed);
    if (it != observeds_.end())
        return false;

    std::pair<std::set<const ObservableBase*>::iterator, bool> pr = observeds_.insert(observed);
    if (pr.second == true) {
        observed->addObserver(this);
    }
    return pr.second;

}

bool Observer::removeObserved(const ObservableBase* observed) const {

    if (!observed)
        return false;

    if (observeds_.find(observed) != observeds_.end()) {
        observeds_.erase(observed);
        observed->removeObserver(this);
        return true;
    }
    else {
        return false;
    }
}

void Observer::stopObservation(const ObservableBase* observable) const {
    removeObserved(observable);
}

// ---------------------------------------------------------------------------

ObservableBase::~ObservableBase() {
    clearObservers();
}

bool ObservableBase::addObserver(const Observer* observer) const {

    if (!observer)
        return false;

    std::pair<std::set<const Observer*>::iterator, bool> pr = observers_.insert(observer);
    if (pr.second)
        observer->addObserved(this);

    return pr.second;
}

bool ObservableBase::removeObserver(const Observer* observer) const {

    if (!observer)
        return false;

    if (observers_.find(observer) != observers_.end()) {
        observers_.erase(observer);
        observer->removeObserved(this);
        return true;
    }
    else {
        return false;
    }
}

void ObservableBase::clearObservers() const {
    while (!observers_.empty())
        removeObserver(*observers_.begin());
}

} // namespace
