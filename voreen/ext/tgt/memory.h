/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2024 University of Muenster, Germany,           *
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

#ifndef TGT_MEMORY_H
#define TGT_MEMORY_H

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>
#include <functional>

namespace tgt {

// See https://stackoverflow.com/questions/17902405/how-to-implement-make-unique-function-in-c11/17902439#17902439
//
// The previous version (as shown in the c++ docs for std::make_unique in 'possible implementation' does not cover
// make_unique for arrays.

template<class T> struct _Unique_if {
    typedef std::unique_ptr<T> _Single_object;
};

template<class T> struct _Unique_if<T[]> {
    typedef std::unique_ptr<T[]> _Unknown_bound;
};

template<class T, size_t N> struct _Unique_if<T[N]> {
    typedef void _Known_bound;
};

template<class T, class... Args>
typename _Unique_if<T>::_Single_object
make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<class T>
typename _Unique_if<T>::_Unknown_bound
make_unique(size_t n) {
    typedef typename std::remove_extent<T>::type U;
    return std::unique_ptr<T>(new U[n]());
}

template<class T, class... Args>
typename _Unique_if<T>::_Known_bound
make_unique(Args&&...) = delete;

class ScopeGuard {
public:
    template< class Func >
    ScopeGuard( Func const& cleanup )
        : cleanup_( cleanup )
    {}

    ScopeGuard( ScopeGuard&& other )
        : cleanup_( std::move( other.cleanup_ ) )
    { other.dismiss(); }

    ScopeGuard& operator=(ScopeGuard const&) = delete;
    ScopeGuard(ScopeGuard const&) = delete;
    ScopeGuard& operator=(ScopeGuard&&) = default;

    ~ScopeGuard() { cleanup_(); }

    void dismiss() { cleanup_ = []{}; }

private:
    std::function<void()> cleanup_;
};

} //namespace tgt

#endif // TGT_MEMORY_H
