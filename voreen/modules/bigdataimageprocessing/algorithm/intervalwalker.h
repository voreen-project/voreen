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

#ifndef VRN_INTERVALWALKER_H
#define VRN_INTERVALWALKER_H

#include "tgt/assert.h"
#include<vector>
#include<algorithm>
#include<string>
#include<iostream>

template<typename Index, typename Value>
struct Interval {
    Interval() = default;
    Interval(Index begin, Index end, Value value)
        : begin(begin)
        , end(end)
        , value(value)
    { }

    Index begin; //inclusive
    Index end; //exclusive
    Value value;
};

template<typename Index, typename Value>
struct IntervalWalkerIterator {
public:
    typedef Interval<Index, Value> IV;

    IntervalWalkerIterator(Index currentPos, std::vector<IV>& current)
        : currentPos_(currentPos)
        , current_(current)
        , nextMove_(current.begin())
        , nextCheck_(nextMove_)
    { }


    typename std::vector<IV>::iterator next() {
        auto end = current_.end();
        while(nextCheck_ != end) {
            if(nextCheck_->end > currentPos_) {
                // Interval counts -> swap it to valid pos and return iter to it
                //std::iter_swap(nextCheck, nextMove);
                *nextMove_ = *nextCheck_;
                auto ret = nextMove_;
                ++nextMove_;
                ++nextCheck_;
                return ret;
            } else {
                // Interval is ignored -> removed
                ++nextCheck_;
            }
        }
        tgtAssert(nextCheck_ >= nextMove_, "Invalid check/move pos");
        int64_t toRemove = std::distance(nextMove_, nextCheck_);
        if(toRemove > 0) {
            size_t size = current_.size();
            tgtAssert(toRemove <= static_cast<int64_t>(size), "Invalid number of elements to remove");
            size_t newLen = size - static_cast<size_t>(toRemove);
            tgtAssert(newLen <= size, "Size was not reduced");
            current_.resize(newLen);
        }
        return current_.end();
    }
    typename std::vector<IV>::iterator end() {
        return current_.end();
    }

    Index currentPos() const {
        return currentPos_;
    }

private:
    const Index currentPos_;
    std::vector<IV>& current_;
    typename std::vector<IV>::iterator nextMove_;
    typename std::vector<IV>::iterator nextCheck_;
};

template<typename Index, typename Value>
struct IntervalWalker {
public:
    typedef Interval<Index, Value> IV;

    IntervalWalker(Index begin, std::vector<IV>&& intervals)
        : currentPos_(begin-1)
        , beginning_(std::move(intervals))
        , current_()
    {
        std::sort(beginning_.begin(), beginning_.end(), [](const IV& i1, const IV& i2) {
            return i1.begin > i2.begin; // See invariant of `beginning`.
        });
    }
    IntervalWalkerIterator<Index, Value> next() {
        currentPos_ += 1;
        updateCurrent();
        return IntervalWalkerIterator<Index, Value>(
            currentPos_,
            current_
        );
    }

private:
    void updateCurrent() {
        while(!beginning_.empty() && beginning_.back().begin <= currentPos_) {
            current_.push_back(beginning_.back());
            beginning_.pop_back();
        }
    }

    Index currentPos_;
    std::vector<IV> beginning_; //sorted by begin, largest to smallest
    std::vector<IV> current_;
};

#ifdef VRN_INTERVALWALKER_TEST

template<typename I, typename V>
std::ostream& operator << (std::ostream& s, const Interval<I, V>& v) {
    return (s << "[" << +v.begin << " " << +v.end << "](" << +v.value << ")");
}

#define error(msg) std::cerr << caseNum << " | Error: " << msg << std::endl; return 0

typedef Interval<size_t, size_t> TIV;
int test(int caseNum, std::vector<TIV> vals, std::vector<std::vector<size_t>> expected) {
    std::vector<TIV> valsClone { vals };
    IntervalWalker<size_t, size_t> iw(0, std::move(valsClone));

    for(size_t i=0; i<expected.size(); ++i) {
        auto intervals = iw.next();
        auto& expectedIntervals = expected[i];
        size_t j = 0;
        while(true) {
            std::vector<TIV>::iterator it = intervals.next();
            if(it == intervals.end()) {
                break;
            }
            if(j >= expectedIntervals.size()) {
                error("Outer: " << i << " Inner: " << j << " Received more Intervals than expected. At least: " << *it);
            }

            if(it->value != expectedIntervals[j]) {
                error("Outer: " << i << " Inner: " << j << " Intervals differ: " << vals[expectedIntervals[j]] << " <-> " << *it);
            }
            ++j;
        }
        if(j < expectedIntervals.size()) {
            error("Outer: " << i << " Inner: " << j << " Expected at least one more Interval: " << vals[expectedIntervals[j]]);
        }
    }

    return 0;
}

int main() {
    test(0,
        std::vector<TIV> { },
        std::vector<std::vector<size_t>> {
            std::vector<size_t> {},
            std::vector<size_t> {},
            std::vector<size_t> {},
            std::vector<size_t> {},
        }
    );
    test(1,
        std::vector<TIV> { TIV(1, 3, 0) },
        std::vector<std::vector<size_t>> {
            std::vector<size_t> {},
            std::vector<size_t> {0},
            std::vector<size_t> {0},
            std::vector<size_t> {},
        }
    );
    test(1,
        std::vector<TIV> { TIV(0, 0, 0),  },
        std::vector<std::vector<size_t>> {
            std::vector<size_t> {},
            std::vector<size_t> {},
            std::vector<size_t> {},
            std::vector<size_t> {},
        }
    );
    test(1,
        std::vector<TIV> { TIV(0, 1, 0), TIV(1, 2, 1) },
        std::vector<std::vector<size_t>> {
            std::vector<size_t> {0},
            std::vector<size_t> {1},
            std::vector<size_t> {},
            std::vector<size_t> {},
        }
    );
    test(1,
        std::vector<TIV> { TIV(0, 1, 0), TIV(1, 3, 1), TIV(2, 4, 2) },
        std::vector<std::vector<size_t>> {
            std::vector<size_t> {0},
            std::vector<size_t> {1},
            std::vector<size_t> {1, 2},
            std::vector<size_t> {2},
        }
    );
    test(1,
        std::vector<TIV> { TIV(2, 4, 2), TIV(1, 3, 1), TIV(0, 1, 0) },
        std::vector<std::vector<size_t>> {
            std::vector<size_t> {0},
            std::vector<size_t> {1},
            std::vector<size_t> {1, 2},
            std::vector<size_t> {2},
        }
    );
    test(1,
        std::vector<TIV> { TIV(2, 4, 2), TIV(1, 3, 1), TIV(0, 1, 0) },
        std::vector<std::vector<size_t>> {
            std::vector<size_t> {},
            std::vector<size_t> {1},
            std::vector<size_t> {1, 2},
            std::vector<size_t> {2},
        }
    );
}

#endif //VRN_INTERVALWALKER_TEST
#endif //VRN_INTERVALWALKER_H
