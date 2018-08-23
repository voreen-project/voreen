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

#ifndef TGT_MATRIX_H
#define TGT_MATRIX_H

#include <cstring>

#include "tgt/vector.h"
#include "tgt/assert.h"

/*
    FIXME build an SSE-optimized mat4
*/

/*
    Performance remarks:

    - Most functions and operators are programmed with loops.
    Compilers can unroll that.
    gcc for instance will do this with -funroll-loops.

    - Modern compilers (gcc 4.x and icc) can even vectorize these loops
    and do these caclulations with SIMD-instructions.

    - Almost all functions are inline!
    This can have signifcant advantages in codegeneration

    - The return-value optimization which is supported by all
    modern compilers (gcc, icc and msvc for instance)
    will optimize unnecessary constructor-overhead away.

    - it could be usefull to use options like -funroll-loops
*/

namespace tgt {

/**
    Matrix2 class. Works basicly like the matrix stuff in GLSL.
    But remember that this is a ROW-MATRIX and NOT a COL-MATRIX
    like in OpenGL or GLSL respectively.
*/
template<class T>
struct Matrix2 {
    typedef T ElemType;

    enum {
        size = 4,
        cols = 2,
        rows = 2
    };
    typedef Vector2<T> RowType;

    union {
        struct {
            T t00, t01;
            T t10, t11;
        };
        T elemRowCol[rows][cols];
        T elem[size];
        struct {
            RowType elemRows[rows];
        };
    };

/*
    constructors
*/
    /// default constructor
    Matrix2() {}
    /// Init all elements with the same value
    explicit Matrix2(T t) {
        for (size_t i = 0; i < size; ++i)
            elem[i] = t;
    }
    /// Init with another Matrix of another type
    template<class U>
    Matrix2(const Matrix2<U>& m) {
        for (size_t i = 0; i < m.size; ++i)
            elem[i] = T(m.elem[i]);
    }
    /// Init from array with equal size
    explicit Matrix2(const T* t) {
        for (size_t i = 0; i < size; ++i)
            elem[i] = t[i];
    }
    /// Init componentwisely
    Matrix2(T _t00, T _t01,
            T _t10, T _t11) {
        t00 = _t00; t01 = _t01;
        t10 = _t10; t11 = _t11;
    }
    /// Init with two Vector2
    Matrix2(const Vector2<T>& v1, const Vector2<T>& v2) {
        for (size_t i = 0; i < v1.size; ++i)
            elem[0*rows + i] = v1.elem[i];
        for (size_t i = 0; i < v2.size; ++i)
            elem[1*rows + i] = v2.elem[i];
    }

/*
    create special matrices
*/

    static const Matrix2<T> zero;
    static const Matrix2<T> identity;

    static Matrix2<T> createZero() {
        return zero;
    }

    static Matrix2<T> createIdentity() {
        return identity;
    }

    static Matrix2<T> createScale(const Vector2<T>& v) {
        Matrix2<T> result = identity;
        result.t00 = v.elem[0];
        result.t11 = v.elem[1];
        return result;
    }

    /**
     * Creates a planar rotation with a given \p angle in \em RADIAN AND NOT \em DEGREE.
     *
     * @param angle The angle of rotation in \em radian.
     */
    static Matrix2<T> createRotation(T angle) {
        Matrix2<T> result = identity;
        result.t00 =  cos(angle);
        result.t01 = -sin(angle);
        result.t10 =  sin(angle);
        result.t11 =  cos(angle);
        return result;
    }

/*
    index operators
*/
    /*
        The C++ standard forbids classes with non trivial constructors
        as union members. Thus this hack provides functionality similar
        to the GLSL matrix index operator.
    */
    /// return arbitrary row-vector similar to GLSL. <br>
    /// BUT REMEMBER: GLSL gives you COLS NOT ROWS!
    const Vector2<T>& operator [] (size_t i) const {
        tgtAssert(i < size, "i must be less than size");
        return elemRows[i];
    }
    /// return arbitrary row-vector similar to GLSL. <br>
    /// BUT REMEMBER: GLSL gives you COLS NOT ROWS!
    Vector2<T>& operator [] (size_t i) {
        tgtAssert(i < size, "i must be less than size");
        return elemRows[i];
    }
};

/*
    init statics
*/

template<class T>
const Matrix2<T> Matrix2<T>::zero = Matrix2<T>
(
    T(0), T(0),
    T(0), T(0)
);

template<class T>
const Matrix2<T> Matrix2<T>::identity = Matrix2<T>
(
    T(1), T(0),
    T(0), T(1)
);

//------------------------------------------------------------------------------

/**
    Matrix3 class. Works basicly like the matrix stuff in GLSL.
    But remember that this is a ROW-MATRIX and NOT a COL-MATRIX
    like in OpenGL or GLSL respectively.
*/
template<class T>
struct Matrix3 {
    typedef T ElemType;

    enum {
        size = 9,
        cols = 3,
        rows = 3
    };
    typedef Vector3<T> RowType;

    union {
        struct {
            T t00, t01, t02;
            T t10, t11, t12;
            T t20, t21, t22;
        };
        T elemRowCol[rows][cols];
        T elem[size];
        struct {
            RowType elemRows[rows];
        };
    };

/*
    constructors
*/

    /// default constructor
    Matrix3() {}
    /// Init all elements with the same value
    explicit Matrix3(T t) {
        for (size_t i = 0; i < size; ++i)
            elem[i] = t;
    }
    /// Init with another Matrix of another type
    template<class U>
    Matrix3(const Matrix3<U>& m) {
        for (size_t i = 0; i < m.size; ++i)
            elem[i] = T(m.elem[i]);
    }
    /// Init from array with equal size
    explicit Matrix3(const T* t) {
        for (size_t i = 0; i < size; ++i)
            elem[i] = t[i];
    }
    /// Init componentwisely
    Matrix3(T _t00, T _t01, T _t02,
            T _t10, T _t11, T _t12,
            T _t20, T _t21, T _t22) {
        t00 = _t00; t01 = _t01; t02 = _t02;
        t10 = _t10; t11 = _t11; t12 = _t12;
        t20 = _t20; t21 = _t21; t22 = _t22;
    }
    /// Init with three Vector3
    Matrix3(const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3) {
        for (size_t i = 0; i < v1.size; ++i)
            elem[0*rows + i] = v1.elem[i];
        for (size_t i = 0; i < v2.size; ++i)
            elem[1*rows + i] = v2.elem[i];
        for (size_t i = 0; i < v3.size; ++i)
            elem[2*rows + i] = v3.elem[i];
    }

/*
    create special matrices
*/
    static const Matrix3<T> zero;
    static const Matrix3<T> identity;

    static Matrix3<T> createZero() {
        return zero;
    }

    static Matrix3<T> createIdentity() {
        return identity;
    }

    static Matrix3<T> createTranslation(const Vector2<T>& v) {
        Matrix3<T> result = identity;
        result.t02 = v.elem[0];
        result.t12 = v.elem[1];
        return result;
    }

    static Matrix3<T> createScale(const Vector3<T>& v) {
        Matrix3<T> result = identity;
        result.t00 = v.elem[0];
        result.t11 = v.elem[1];
        result.t22 = v.elem[2];
        return result;
    }

    /**
     * Rotation the x-axis as rotation axis.
     * The angle is in \em RADIAN AND NOT \em DEGREE.
     *
     * @param angle The angle of rotation in \em radian.
     */
    static Matrix3<T> createRotationX(T angle) {
        Matrix3<T> result = identity;
        result.t11 =  cos(angle);
        result.t12 = -sin(angle);
        result.t21 =  sin(angle);
        result.t22 =  cos(angle);
        return result;
    }

    /**
     * Rotation the y-axis as rotation axis.
     * The angle is in \em RADIAN AND NOT \em DEGREE.
     *
     * @param angle The angle of rotation in \em radian.
     */
    static Matrix3<T> createRotationY(T angle) {
        Matrix3<T> result = identity;
        result.t00 =  cos(angle);
        result.t02 =  sin(angle);
        result.t20 = -sin(angle);
        result.t22 =  cos(angle);
        return result;
    }

    /**
     * Rotation the z-axis as rotation axis.
     * The angle is in \em RADIAN AND NOT \em DEGREE.
     *
     * @param angle The angle of rotation in \em radian.
     */
    static Matrix3<T> createRotationZ(T angle) {
        Matrix3<T> result = identity;
        result.t00 =  cos(angle);
        result.t01 = -sin(angle);
        result.t10 =  sin(angle);
        result.t11 =  cos(angle);
        return result;
    }

    /**
     * Simliar to glRotate but \p angle is in \em RADIAN AND NOT \em DEGREE.
     *
     * @param angle The angle in \em radian.
     * @param axis The axis of rotation. Needn't be in unit length.
     */
    static Matrix3<T> createRotation(T angle, Vector3<T> axis);

    /**
     * inverts the matrix.
     * @return true - if successful -- otherwise false
    */
    bool invert(Matrix3<T>& result) const;

    T det() {
        return elem[0]*elem[4]*elem[8] + elem[1]*elem[5]*elem[6] + elem[2]*elem[3]*elem[7] - elem[0]*elem[5]*elem[7] - elem[1]*elem[3]*elem[8] - elem[2]*elem[4]*elem[6];
    }
/*
    index operator
*/
    /*
        The C++ standard forbids classes with non trivial constructors
        as union members. Thus this hack provides functionality similar
        to the GLSL matrix index operator.
    */
    /// return arbitrary row-vector similar to GLSL. <br>
    /// BUT REMEMBER: GLSL gives you COLS NOT ROWS!
    const Vector3<T>& operator [] (size_t i) const {
        tgtAssert(i < size, "i must be less than size");
        return elemRows[i];
    }
    /// return arbitrary row-vector similar to GLSL. <br>
    /// BUT REMEMBER: GLSL gives you COLS NOT ROWS!
    Vector3<T>& operator [] (size_t i) {
        tgtAssert(i < size, "i must be less than size");
        return elemRows[i];
    }
};

/*
    init statics
*/

/// init statics
template<class T>
const Matrix3<T> Matrix3<T>::zero = Matrix3<T>
(
    T(0), T(0), T(0),
    T(0), T(0), T(0),
    T(0), T(0), T(0)
);

template<class T>
const Matrix3<T> Matrix3<T>::identity = Matrix3<T>
(
    T(1), T(0), T(0),
    T(0), T(1), T(0),
    T(0), T(0), T(1)
);

template<class T>
bool Matrix3<T>::invert(Matrix3<T>& result) const {
    using std::abs; // use overloaded abs
    result = Matrix3<T>::identity;

    T z = elem[0]*(elem[4]*elem[8] - elem[5]*elem[7]) + elem[1]*(elem[5]*elem[6] - elem[3]*elem[8]) + elem[2]*(elem[3]*elem[7] - elem[4]*elem[6]);

    if( abs(z) == 0.f/*< 0.00001*/)
        return false;

    //we have a row-matrix:
    result.elem[0] = (elem[4]*elem[8] - elem[5]*elem[7]);
    result.elem[1] = (elem[2]*elem[7] - elem[1]*elem[8]);
    result.elem[2] = (elem[1]*elem[5] - elem[2]*elem[4]);

    result.elem[3] = (elem[5]*elem[6] - elem[3]*elem[8]);
    result.elem[4] = (elem[0]*elem[8] - elem[2]*elem[6]);
    result.elem[5] = (elem[2]*elem[3] - elem[0]*elem[5]);

    result.elem[6] = (elem[3]*elem[7] - elem[4]*elem[6]);
    result.elem[7] = (elem[1]*elem[6] - elem[0]*elem[7]);
    result.elem[8] = (elem[0]*elem[4] - elem[1]*elem[3]);

    result /= z;

    return true;
}
//------------------------------------------------------------------------------

/**
    Matrix4 class. Works basicly like the matrix stuff in GLSL.
    But remember that this is a ROW-MATRIX and NOT a COL-MATRIX
    like in OpenGL or GLSL respectively.
*/
template<class T>
struct Matrix4 {
    typedef T ElemType;

    enum {
        size = 16,
        cols = 4,
        rows = 4
    };
    typedef Vector4<T> RowType;

    union {
        struct {
            T t00, t01, t02, t03;
            T t10, t11, t12, t13;
            T t20, t21, t22, t23;
            T t30, t31, t32, t33;
        };
        T elemRowCol[rows][cols];
        T elem[size];
        struct {
            RowType elemRows[rows];
        };
    };

/*
    constructors
*/
    /// default constructor
    Matrix4() {}
    /// Init all elements with the same value
    explicit Matrix4(T t) {
        for (size_t i = 0; i < size; ++i)
            elem[i] = t;
    }
    /// Init with another Matrix of another type
    template<class U>
    Matrix4(const Matrix4<U>& m) {
        for (size_t i = 0; i < m.size; ++i)
            elem[i] = T(m.elem[i]);
    }
    /// Init from array with equal size
    explicit Matrix4(const T* t) {
        for (size_t i = 0; i < size; ++i)
            elem[i] = t[i];
    }
    /// Init componentwisely
    Matrix4(T _t00, T _t01, T _t02, T _t03,
            T _t10, T _t11, T _t12, T _t13,
            T _t20, T _t21, T _t22, T _t23,
            T _t30, T _t31, T _t32, T _t33) {
        t00 = _t00; t01 = _t01; t02 = _t02; t03 = _t03;
        t10 = _t10; t11 = _t11; t12 = _t12; t13 = _t13;
        t20 = _t20; t21 = _t21; t22 = _t22; t23 = _t23;
        t30 = _t30; t31 = _t31; t32 = _t32; t33 = _t33;
    }
    /// Init with four Vector4
    Matrix4(const Vector4<T>& v1, const Vector4<T>& v2,
            const Vector4<T>& v3, const Vector4<T>& v4) {
        for (size_t i = 0; i < v1.size; ++i)
            elem[0*rows + i] = v1.elem[i];
        for (size_t i = 0; i < v2.size; ++i)
            elem[1*rows + i] = v2.elem[i];
        for (size_t i = 0; i < v3.size; ++i)
            elem[2*rows + i] = v3.elem[i];
        for (size_t i = 0; i < v4.size; ++i)
            elem[3*rows + i] = v4.elem[i];
    }

    /// Return the upper left 3x3-submatrix in a 4x4 matrix
    Matrix4<T> getRotationalPart() const {
        Matrix4<T> result(T(0));
        result.t00 = t00; result.t01  = t01; result.t02 = t02;
        result.t10 = t10; result.t11  = t11; result.t12 = t12;
        result.t20 = t20; result.t21  = t21; result.t22 = t22;
        result.t33 = T(1);
        return result;
    }

    /// Return the upper left 3x3-submatrix in a 3x3 matrix
    Matrix3<T> getRotationalPartMat3() const {
        Matrix3<T> result;
        result.t00 = t00; result.t01  = t01; result.t02 = t02;
        result.t10 = t10; result.t11  = t11; result.t12 = t12;
        result.t20 = t20; result.t21  = t21; result.t22 = t22;
        return result;
    }

    /// Returns a vector representing the scaling applied to unit vectors by this matrix.
    Vector3<T> getScalingPart() const {
        Vector3<T> result(T(0));
        result.elem[0] = length(Vector3<T>(t00, t10, t20));
        result.elem[1] = length(Vector3<T>(t01, t11, t21));
        result.elem[2] = length(Vector3<T>(t02, t12, t22));
        return result;
    }

    Vector3<T> getTranslationPart() const {
        Vector3<T> result(t03,t13,t23);
        return result;
    }
/*
    create special matrices
*/
    static const Matrix4<T> zero;
    static const Matrix4<T> identity;

    static Matrix4<T> createZero() {
        return zero;
    }

    static Matrix4<T> createIdentity() {
        return identity;
    }

    static Matrix4<T> createTranslation(const Vector3<T>& v) {
        Matrix4<T> result = identity;
        result.t03 = v.elem[0];
        result.t13 = v.elem[1];
        result.t23 = v.elem[2];
        return result;
    }

    static Matrix4<T> createScale(const Vector3<T>& v) {
        Matrix4<T> result = identity;
        result.t00 = v.elem[0];
        result.t11 = v.elem[1];
        result.t22 = v.elem[2];
        return result;
    }

    /**
     * Rotation the x-axis as rotation axis.
     * The angle is in \em RADIAN AND NOT \em DEGREE.
     *
     * @param angle The angle of rotation in \em radian.
     */
    static Matrix4<T> createRotationX(T angle) {
        using namespace std; // use overloaded cos and sin
        Matrix4<T> result = identity;
        result.t11 =  cos(angle);
        result.t12 = -sin(angle);
        result.t21 =  sin(angle);
        result.t22 =  cos(angle);
        return result;
    }

    /**
     * Rotation the y-axis as rotation axis.
     * The angle is in \em RADIAN AND NOT \em DEGREE.
     *
     * @param angle The angle of rotation in \em radian.
     */
    static Matrix4<T> createRotationY(T angle) {
        using namespace std; // use overloaded cos and sin
        Matrix4<T> result = identity;
        result.t00 =  cos(angle);
        result.t02 =  sin(angle);
        result.t20 = -sin(angle);
        result.t22 =  cos(angle);
        return result;
    }

    /**
     * Rotation the z-axis as rotation axis.
     * The angle is in \em RADIAN AND NOT \em DEGREE.
     *
     * @param angle The angle of rotation in \em radian.
     */
    static Matrix4<T> createRotationZ(T angle) {
        using namespace std; // use overloaded cos and sin
        Matrix4<T> result = identity;
        result.t00 =  cos(angle);
        result.t01 = -sin(angle);
        result.t10 =  sin(angle);
        result.t11 =  cos(angle);
        return result;
    }

    static Matrix4<T> createRotationXDegree(T angle);
    static Matrix4<T> createRotationYDegree(T angle);
    static Matrix4<T> createRotationZDegree(T angle);

    /**
     * Simliar to glRotate but \p angle is in \em RADIAN AND NOT \em DEGREE.
     *
     * @param angle The angle in \em radian.
     * @param axis The axis of rotation. Needn't be in unit length.
     */
    static Matrix4<T> createRotation(T angle, Vector3<T> axis);

    static Matrix4<T> createRotationDegree(T angle, Vector3<T> axis);

    /**
    * This function returns the Matrix that gluLookAt would create and put onto the OpenGL
    * Matrix-Stack, given the position, focus-point, and up-Vector of a camera.
    */
    static Matrix4<T> createLookAt(const Vector3<T>& eye, const Vector3<T>& focus, const Vector3<T>& up);

    /**
     * This function provides a glFrustum() replacement.  In case you just need the matrix that glFrustum would
     * produce, you can use this function instead of putting the matrix on the OpenGL stack and then reading it.
     */
    static Matrix4<T> createFrustum(T left, T right, T top, T bottom, T pnear, T pfar);

    /**
     * This function provides a gluPerspective() replacement.
     * See its documentation for an explanation of the parameters. \em BUT:
     * \p fov expects \em RADIAN and not \em DEGREE.
     */
    static Matrix4<T> createPerspective(T fov, T aspect, T pnear, T pfar);

    static Matrix4<T> createOrtho(T left, T right, T bottom, T top, T pnear, T pfar);

    /**
     * We are assuming the look vector of the object is the positive z-axis
     */
    static Matrix4<T> rigidBodyTranformation(Vector3<T> finalLook, Vector3<T> finalUp, Vector3<T> finalTranslation);

/*
    index operator
*/

   /*
        The C++ standard forbids classes with non trivial constructors
        as union members. Thus this hack provides functionality similar
        to the GLSL matrix index operator.
    */
    /// return arbitrary row-vector similar to GLSL. <br>
    /// BUT REMEMBER: GLSL gives you COLS NOT ROWS!
    const Vector4<T>& operator [] (size_t i) const {
        tgtAssert(i < size, "i must be less than size");
        return elemRows[i];
    }
    /// return arbitrary row-vector similar to GLSL. <br>
    /// BUT REMEMBER: GLSL gives you COLS NOT ROWS!
    Vector4<T>& operator [] (size_t i) {
        tgtAssert(i < size, "i must be less than size");
        return elemRows[i];
    }

    /**
     * inverts the matrix.
     * @return true - if successful -- otherwise false
    */
    bool invert(Matrix4<T>& result) const;

    /**
     * Interpret this matrix as a transformation with a rotational/translational part
     * and use it to transform the given vector.
     */
    Vector3<T> transform(const Vector3<T>) const;
};

/*
    init statics
*/
template<class T>
const Matrix4<T> Matrix4<T>::zero = Matrix4<T>
(
    T(0), T(0), T(0), T(0),
    T(0), T(0), T(0), T(0),
    T(0), T(0), T(0), T(0),
    T(0), T(0), T(0), T(0)
);

template<class T>
const Matrix4<T> Matrix4<T>::identity = Matrix4<T>
(
    T(1), T(0), T(0), T(0),
    T(0), T(1), T(0), T(0),
    T(0), T(0), T(1), T(0),
    T(0), T(0), T(0), T(1)
);


/*
    typedefs for easy usage
*/

typedef Matrix2<float>  Matrix2f;
typedef Matrix2<double> Matrix2d;
typedef Matrix3<float>  Matrix3f;
typedef Matrix3<double> Matrix3d;
typedef Matrix4<float>  Matrix4f;
typedef Matrix4<double> Matrix4d;

typedef Matrix4f        Matrix;

/*
    types like in GLSL
*/
typedef Matrix2f        mat2;
typedef Matrix3f        mat3;
typedef Matrix4f        mat4;

typedef Matrix2d        dmat2;
typedef Matrix3d        dmat3;
typedef Matrix4d        dmat4;

#ifdef DLL_TEMPLATE_INST
template struct TGT_API Matrix2<float>;
template struct TGT_API Matrix3<float>;
template struct TGT_API Matrix4<float>;
template struct TGT_API Matrix2<double>;
template struct TGT_API Matrix3<double>;
template struct TGT_API Matrix4<double>;
#endif

/*
    Prepare to implement 3 * 5 * 7 = 105 operators and a couple of functions
    the lazy way with evil voodoo macro magic
*/

#define TGT_MAT_TRANSPOSE \
template<class T> inline TGT_BASE_TYPE<T> transpose(const TGT_BASE_TYPE<T>& m) { \
    TGT_BASE_TYPE<T> mRes; \
    for (size_t row = 0; row < m.rows; ++row) \
        for (size_t col = 0; col < m.cols; ++col) \
            mRes.elemRowCol[col][row] = m.elemRowCol[row][col]; \
    return mRes; \
}

/*
    Similar to GLSL make an exception with the
    Matrix-Matrix, Vector-Matrix and Matrix-Vector multiplication.
    These are "correct" operations and do not go componentwisely
*/

#define TGT_MAT_MUL_VEC \
template<class T> inline typename TGT_BASE_TYPE<T>::RowType operator * (const TGT_BASE_TYPE<T>& m, const typename TGT_BASE_TYPE<T>::RowType& v) { \
    typename TGT_BASE_TYPE<T>::RowType vRes; \
    for (size_t i = 0; i < v.size; ++i) \
        vRes[i] = dot(m[i], v); \
    return vRes; \
}

#define TGT_VEC_MUL_MAT \
template<class T> inline typename TGT_BASE_TYPE<T>::RowType operator * (const typename TGT_BASE_TYPE<T>::RowType& v, const TGT_BASE_TYPE<T>& m) { \
    typename TGT_BASE_TYPE<T>::RowType vRes; \
    TGT_BASE_TYPE<T> transposed = transpose(m); \
    for (size_t i = 0; i < v.size; ++i) \
        vRes[i] = dot(v, transposed[i]); \
    return vRes; \
}

#define TGT_MAT_MUL_MAT \
template<class T> inline TGT_BASE_TYPE<T> operator * (const TGT_BASE_TYPE<T>& m1, const TGT_BASE_TYPE<T>& m2) { \
    TGT_BASE_TYPE<T> mRes; \
    TGT_BASE_TYPE<T> transposed = transpose(m2); \
    for (size_t row = 0; row < m1.rows; ++row) \
        for (size_t col = 0; col < m1.cols; ++col) \
            mRes.elemRowCol[row][col] = \
                dot(m1[row], transposed[col]); \
    return mRes; \
}

#define TGT_MAT_MULEQ_MAT \
template<class T> inline TGT_BASE_TYPE<T>& operator *= (TGT_BASE_TYPE<T>& m1, const TGT_BASE_TYPE<T>& m2) { \
    TGT_BASE_TYPE<T> temp = m1; \
    m1 = temp * m2; \
    return m1; \
}

#define TGT_MAT_COMP_MUL \
template<class T> inline TGT_BASE_TYPE<T> matrixCompMult (const TGT_BASE_TYPE<T>& v1, const TGT_BASE_TYPE<T>& v2) { \
    TGT_BASE_TYPE<T> vRes; \
    for (size_t i = 0; i < v1.size; ++i) \
        vRes.elem[i] = v1.elem[i] * v2.elem[i]; \
    return vRes; \
}

/*
det
m[0][0] * m[1][1] * m[2][2]
+ m[0][1] * m[1][2] * m[2][0]
+ m[0][2] * m[1][0] * m[2][1]
- m[0][2] * m[1][1] * m[2][0]
- m[0][1] * m[1][0] * m[2][2]
- m[0][0] * m[1][2] * m[2][1]
*/

#define TGT_IMPLEMENT_MAT_FUNCTIONS \
    TGT_VEC_UNARY_MINUS \
    TGT_VEC_MIN \
    TGT_VEC_MAX \
    TGT_VEC_FLOOR \
    TGT_VEC_CEIL \
    TGT_VEC_MIN_SINGLE \
    TGT_VEC_MAX_SINGLE \
    TGT_VEC_MIN_SELF \
    TGT_VEC_MAX_SELF \
    TGT_VEC_CLAMP \
    TGT_VEC_CLAMP_SINGLE \
    TGT_VEC_HADD \
    TGT_VEC_HSUB \
    TGT_VEC_HMUL \
    TGT_VEC_HDIV \
    TGT_VEC_HMOD \
    TGT_VEC_HAND \
    TGT_VEC_HOR \
    TGT_VEC_HXOR \
    TGT_VEC_REL_OP_EQUAL  \
    TGT_VEC_REL_OP_NOT_EQUAL \
    TGT_VEC_LESS_THAN \
    TGT_VEC_LESS_THAN_EQUAL \
    TGT_VEC_GREATER_THAN \
    TGT_VEC_GREATER_THAN_EQUAL \
    TGT_VEC_EQUAL \
    TGT_VEC_NOT_EQUAL \
    TGT_MAT_TRANSPOSE \
    TGT_MAT_COMP_MUL \
    TGT_MAT_MUL_VEC \
    TGT_VEC_MUL_MAT \
    TGT_MAT_MUL_MAT \
    TGT_MAT_MULEQ_MAT

/*
    Implementation of Matrix2<T> operators
*/

#define TGT_BASE_TYPE Matrix2
    #define TGT_VEC_OP +
    #define TGT_VEC_OPEQ +=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP -
    #define TGT_VEC_OPEQ -=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP *
    #define TGT_VEC_OPEQ *=
        /*
            According to GLSL matrix * matrix should be exceptionally a
            "correct" algebraic multiplication.
            Instead use matrixCompMult if you need this operatrion.
        */
        TGT_VEC_OP_BASE
        TGT_BASE_OP_VEC
        TGT_VEC_OPEQ_BASE
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP /
    #define TGT_VEC_OPEQ /=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP &
    #define TGT_VEC_OPEQ &=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP |
    #define TGT_VEC_OPEQ |=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP ^
    #define TGT_VEC_OPEQ ^=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    TGT_IMPLEMENT_MAT_FUNCTIONS
#undef TGT_BASE_TYPE

/*
    Implementation of Matrix3<T> operators
*/

#define TGT_BASE_TYPE Matrix3
    #define TGT_VEC_OP +
    #define TGT_VEC_OPEQ +=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP -
    #define TGT_VEC_OPEQ -=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP *
    #define TGT_VEC_OPEQ *=
        /*
            According to GLSL matrix * matrix should be exceptionally a
            "correct" algebraic multiplication.
            Instead use matrixCompMult if you need this operatrion.
        */
        TGT_VEC_OP_BASE
        TGT_BASE_OP_VEC
        TGT_VEC_OPEQ_BASE
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP /
    #define TGT_VEC_OPEQ /=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP &
    #define TGT_VEC_OPEQ &=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP |
    #define TGT_VEC_OPEQ |=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP ^
    #define TGT_VEC_OPEQ ^=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    TGT_IMPLEMENT_MAT_FUNCTIONS
#undef TGT_BASE_TYPE

/*
    Implementation of Matrix4<T> operators
*/

#define TGT_BASE_TYPE Matrix4
    #define TGT_VEC_OP +
    #define TGT_VEC_OPEQ +=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP -
    #define TGT_VEC_OPEQ -=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP *
    #define TGT_VEC_OPEQ *=
        /*
            According to GLSL matrix * matrix should be exceptionally a
            "correct" algebraic multiplication.
            Instead use matrixCompMult if you need this operatrion.
        */
        TGT_VEC_OP_BASE
        TGT_BASE_OP_VEC
        TGT_VEC_OPEQ_BASE
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP /
    #define TGT_VEC_OPEQ /=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP &
    #define TGT_VEC_OPEQ &=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP |
    #define TGT_VEC_OPEQ |=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    #define TGT_VEC_OP ^
    #define TGT_VEC_OPEQ ^=
        TGT_IMPLEMENT_OPERATORS
    #undef  TGT_VEC_OP
    #undef  TGT_VEC_OPEQ

    TGT_IMPLEMENT_MAT_FUNCTIONS
#undef TGT_BASE_TYPE

/*
    Undefine macros
*/

#undef TGT_IMPLEMENT_MAT_FUNCTIONS
#undef TGT_MAT_COMP_MUL
#undef TGT_MAT_MULEQ_MAT
#undef TGT_MAT_MUL_MAT
#undef TGT_MAT_MUL_VEC
#undef TGT_MAT_TRANSPOSE
#undef TGT_VEC_MUL_MAT


//------------------------------------------------------------------------------
// Non inline implementation
//------------------------------------------------------------------------------

/*
    createRotation for Matrix3<T> and Matrix4<T>
*/

template<class T>
Matrix4<T> Matrix4<T>::createRotationXDegree(T angle) {
    T cosValue, sinValue;
    double doubleAngle = std::fmod(static_cast<double>(angle),360.0);

    if(doubleAngle ==   0.0) {
        cosValue = T(1); sinValue = T(0);
    } else
    if(doubleAngle ==  90.0) {
        cosValue = T(0); sinValue = T(1);
    } else
    if(doubleAngle == 180.0) {
        cosValue = T(-1); sinValue = T(0);
    } else
    if(doubleAngle == 270.0) {
        cosValue = T(0); sinValue = T(-1);
    } else {
        cosValue = T(cos(deg2rad(doubleAngle))); sinValue = T(sin(deg2rad(doubleAngle)));
    }

    return Matrix4<T>(
            T(1),       T(0),       T(0),   T(0),
            T(0),   cosValue, - sinValue,   T(0),
            T(0),   sinValue,   cosValue,   T(0),
            T(0),        T(0),      T(0),   T(1)
           );
}

template<class T>
Matrix4<T> Matrix4<T>::createRotationYDegree(T angle) {
    T cosValue, sinValue;
    double doubleAngle = std::fmod(static_cast<double>(angle),360.0);

    if(doubleAngle ==   0.0) {
        cosValue = T(1); sinValue = T(0);
    } else
    if(doubleAngle ==  90.0) {
        cosValue = T(0); sinValue = T(1);
    } else
    if(doubleAngle == 180.0) {
        cosValue = T(-1); sinValue = T(0);
    } else
    if(doubleAngle == 270.0) {
        cosValue = T(0); sinValue = T(-1);
    } else {
        cosValue = T(cos(deg2rad(doubleAngle))); sinValue = T(sin(deg2rad(doubleAngle)));
    }

    return Matrix4<T>(
              cosValue,   T(0),   sinValue,   T(0),
                  T(0),   T(1),       T(0),   T(0),
            - sinValue,   T(0),   cosValue,   T(0),
                  T(0),   T(0),       T(0),   T(1)
           );
}

template<class T>
Matrix4<T> Matrix4<T>::createRotationZDegree(T angle) {
    T cosValue, sinValue;
    double doubleAngle = std::fmod(static_cast<double>(angle),360.0);

    if(doubleAngle ==   0.0) {
        cosValue = T(1); sinValue = T(0);
    } else
    if(doubleAngle ==  90.0) {
        cosValue = T(0); sinValue = T(1);
    } else
    if(doubleAngle == 180.0) {
        cosValue = T(-1); sinValue = T(0);
    } else
    if(doubleAngle == 270.0) {
        cosValue = T(0); sinValue = T(-1);
    } else {
        cosValue = T(cos(deg2rad(doubleAngle))); sinValue = T(sin(deg2rad(doubleAngle)));
    }

    return Matrix4<T>(
            cosValue, - sinValue,   T(0),   T(0),
            sinValue,   cosValue,   T(0),   T(0),
                T(0),       T(0),   T(1),   T(0),
                T(0),       T(0),   T(0),   T(1)
           );
}

template<class T>
Matrix3<T> Matrix3<T>::createRotation(T angle, Vector3<T> axis) {
    using namespace std; // use overloaded cos and sin
    axis = normalize(axis);
    T s = T(sin(angle));
    T c = T(cos(angle));
    T x = axis.x;
    T y = axis.y;
    T z = axis.z;

    Matrix3<T> result
    (
        x*x*(1-c) + c,      x*y*(1-c) - z*s,    x*z*(1-c) + y*s,
        x*y*(1-c) + z*s,    y*y*(1-c) + c,      y*z*(1-c) - x*s,
        x*z*(1-c) - y*s,    y*z*(1-c) + x*s,    z*z*(1-c) + c
    );

    return result;
}

template<class T>
Matrix4<T> Matrix4<T>::createRotation(T angle, Vector3<T> axis) {
    using namespace std; // use overloaded cos and sin
    axis = normalize(axis);
    T s = T(sin(angle));
    T c = T(cos(angle));
    T x = axis.x;
    T y = axis.y;
    T z = axis.z;

    Matrix4<T> result
    (
        x*x*(1-c) + c,      x*y*(1-c) - z*s,    x*z*(1-c) + y*s, T(0),
        x*y*(1-c) + z*s,    y*y*(1-c) + c,      y*z*(1-c) - x*s, T(0),
        x*z*(1-c) - y*s,    y*z*(1-c) + x*s,    z*z*(1-c) + c,   T(0),
        T(0),               T(0),               T(0),            T(1)
    );

    return result;
}

template<class T>
Matrix4<T> Matrix4<T>::createRotationDegree(T angle, Vector3<T> axis) {
    using namespace std; // use overloaded cos and sin
    axis = normalize(axis);
    T s = T(sin(deg2rad(angle)));
    T c = T(cos(deg2rad(angle)));
    T x = axis.x;
    T y = axis.y;
    T z = axis.z;

    Matrix4<T> result
    (
        x*x*(1-c) + c,      x*y*(1-c) - z*s,    x*z*(1-c) + y*s, T(0),
        x*y*(1-c) + z*s,    y*y*(1-c) + c,      y*z*(1-c) - x*s, T(0),
        x*z*(1-c) - y*s,    y*z*(1-c) + x*s,    z*z*(1-c) + c,   T(0),
        T(0),               T(0),               T(0),            T(1)
    );

    return result;
}

/*
    createFrustum, createPerspective and createLookAt
*/

template<class T>
Matrix4<T> Matrix4<T>::createLookAt(const Vector3<T>& eye, const Vector3<T>& focus, const Vector3<T>& up) {

    Vector3<T> look   = normalize(focus - eye);
    Vector3<T> strafe = normalize(cross(look, normalize(up)));
    Vector3<T> up2    = cross(strafe, look);

    Matrix4<T> m(Vector4<T>( strafe, T(0) ),
                 Vector4<T>( up2, T(0) ),
                 Vector4<T>( -look, T(0) ),
                 Vector4<T>( T(0), T(0), T(0), T(1) ));

    return m * Matrix4<T>::createTranslation(-eye);
}

template<class T>
Matrix4<T> Matrix4<T>::createFrustum(T left, T right, T bottom, T top, T pnear, T pfar) {
#ifdef TGT_DEBUG
    if ((pnear == pfar) || (left == right) || (top == bottom)) {
        tgtAssert(false, "The parameters passed to createFrustum cannot be used to form a projection matrix.");
        return Matrix4<T>::identity;
    }
#endif // TGT_DEBUG

    Matrix4<T> m(
        T(2)*pnear/(right-left),        T(0),               (right+left)/(right-left),             T(0),
                T(0),           T(2)*pnear/(top-bottom),    (top+bottom)/(top-bottom),             T(0),
                T(0),                   T(0),               (pnear+pfar)/(pnear-pfar),  (T(2)*pfar*pnear)/(pnear-pfar),
                T(0),                   T(0),                          -T(1),                      T(0)
    );

    return m;
}

template<class T>
Matrix4<T> Matrix4<T>::createPerspective(T fov, T aspect, T pnear, T pfar) {
#ifdef TGT_DEBUG
    if ((pnear == pfar) || (rad2deg(fov) > T(355)) || (rad2deg(fov) < T(5)) || (aspect == T(0))) {
        tgtAssert(false, "The parameters passed to createPerspective cannot be used to form a projection matrix.");
        return Matrix4<T>::identity;
    }
#endif // TGT_DEBUG

    T f = T(1) / tan( fov/T(2) );
    Matrix4<T> m(
        f / aspect,         T(0),            T(0),                          T(0),
            T(0),             f,             T(0),                          T(0),
            T(0),           T(0),   (pnear+pfar)/(pnear-pfar),  (T(2)*pfar*pnear)/(pnear-pfar),
            T(0),           T(0),           -T(1),                          T(0)
    );

    return m;
}

template<class T>
Matrix4<T> Matrix4<T>::createOrtho(T left, T right, T bottom, T top, T pnear, T pfar) {
#ifdef TGT_DEBUG
    if ((pnear == pfar) || (left == right) || (top == bottom)) {
        tgtAssert(false, "The parameters passed to createOrthogonal cannot be used to form a projection matrix.");
        return Matrix4<T>::identity;
    }
#endif // TGT_DEBUG

    Matrix4<T> m(
        T(2)/(right-left),      T(0),              T(0),          -(right+left)/(right-left),
             T(0),         T(2)/(top-bottom),      T(0),          -(top+bottom)/(top-bottom),
             T(0),              T(0),         T(-2)/(pfar-pnear), -(pfar+pnear)/(pfar-pnear),
             T(0),              T(0),              T(0),                      T(1)
    );

    return m;
}

template<class T>
Matrix4<T> Matrix4<T>::rigidBodyTranformation(Vector3<T> finalLook, Vector3<T> finalUp, Vector3<T> finalTranslation) {
    Vector3<T> l = normalize(finalLook);
    Vector3<T> uTmp = normalize(finalUp);
    Vector3<T> s = normalize(cross(uTmp,l));
    Vector3<T> u = cross(l,s);

    return Matrix4<T>(s.x,u.x,l.x,finalTranslation.x,
                      s.y,u.y,l.y,finalTranslation.y,
                      s.z,u.z,l.z,finalTranslation.z,
                      T(0),T(0),T(0),T(1));

}


/*
    ostream operators
*/

/// ostream-operator
template<class T>
std::ostream& operator << (std::ostream& s, const Matrix2<T>& m) {
    return (s
        << "| " << m.elem[0] << " " << m.elem[1] << " |" << std::endl
        << "| " << m.elem[2] << " " << m.elem[3] << " |" << std::endl);
}

/// ostream-operator
template<class T>
std::ostream& operator << (std::ostream& s, const Matrix3<T>& m) {
    return (s
        << "| " << m.elem[0] << " " << m.elem[1] << " " << m.elem[2] << " |" << std::endl
        << "| " << m.elem[3] << " " << m.elem[4] << " " << m.elem[5] << " |" << std::endl
        << "| " << m.elem[6] << " " << m.elem[7] << " " << m.elem[8] << " |" << std::endl);
}

/// ostream-operator
template<class T>
std::ostream& operator << (std::ostream& s, const Matrix4<T>& m) {
    return (s
        << "| " << m.elem[ 0] << " " << m.elem[ 1] << " " << m.elem[ 2] << " " << m.elem[ 3] << " |" << std::endl
        << "| " << m.elem[ 4] << " " << m.elem[ 5] << " " << m.elem[ 6] << " " << m.elem[ 7] << " |" << std::endl
        << "| " << m.elem[ 8] << " " << m.elem[ 9] << " " << m.elem[10] << " " << m.elem[11] << " |" << std::endl
        << "| " << m.elem[12] << " " << m.elem[13] << " " << m.elem[14] << " " << m.elem[15] << " |" << std::endl);
}

/// This operator is not available in GLSL but very usefull: A mat4 * vec3 operator, returning a vec3
/// note: the resulting vec3 is divided by the w component of the temporary result!
template<class T>
Vector3<T> operator * (const Matrix4<T>& m, const Vector3<T>& v) {
    Vector4<T> v4(v, 1);
    v4 = m * v4;
    v4 /= v4.w;
    return v4.xyz();
}

// addtional matrix functions

template<class T>
bool Matrix4<T>::invert(Matrix4<T>& result) const {
    using std::abs; // use overloaded abs
    T t;
    Matrix4<T> tmp = *this;
    result = Matrix4<T>::identity;

    for (size_t i = 0; i < 4; ++i) {
        // Look for largest element in column
        size_t swap = i;
        for (size_t j = i + 1; j < 4; ++j) {
            if (abs(tmp.elemRowCol[j][i]) > abs(tmp.elemRowCol[i][i]))
                swap = j;
        }

        if (swap != i) {
            // Swap rows.
            for (size_t k = 0; k < 4; ++k) {
                t = tmp.elemRowCol[i][k];
                tmp.elemRowCol[i][k] = tmp.elemRowCol[swap][k];
                tmp.elemRowCol[swap][k] = t;

                t = result.elem[i*4+k];
                result.elem[i*4+k] = result.elem[swap*4+k];
                result.elem[swap*4+k] = t;
            }
        }

        if (tmp.elemRowCol[i][i] == 0) {
            // The matrix is singular
            return false;
        }

        t = tmp.elemRowCol[i][i];
        for (size_t k = 0; k < 4; k++) {
            tmp.elemRowCol[i][k] /= t;
            result.elem[i*4+k] /= t;
        }
        for (size_t j = 0; j < 4; j++) {
            if (j != i) {
                t = tmp.elemRowCol[j][i];

                for (size_t k = 0; k < 4; k++) {
                    tmp.elemRowCol[j][k] -= tmp.elemRowCol[i][k]*t;
                    result.elem[j*4+k] -= result.elem[i*4+k]*t;
                }
            }
        }
    }
    return true;
}

template<class T>
Vector3<T> Matrix4<T>::transform(const Vector3<T> v) const {
    return (*this * Vector4<T>(v, T(1))).xyz();
}

};

#endif //TGT_MATRIX_H
