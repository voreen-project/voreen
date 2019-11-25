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

#include "pyvoreenobjects.h"

namespace voreen {

//////////////////////////////////////////////////////////////////
// VolumeObject
//////////////////////////////////////////////////////////////////


PyObject* VolumeObject_new(PyTypeObject *type, PyObject */*args*/, PyObject */*kwds*/) {
    VolumeObject *self;
    self = (VolumeObject *) type->tp_alloc(type, 0);
    if (self != NULL) {
        self->format = PyUnicode_FromString("");
        if (self->format == NULL) {
            Py_DECREF(self);
            return NULL;
        }
        self->data = PyList_New(0);
        if (self->data == NULL) {
            Py_DECREF(self);
            return NULL;
        }
        self->dimX = self->dimY = self->dimZ = 0u;
        self->spacingX = self->spacingY = self->spacingZ = 1.0f;
        self->offsetX = self->offsetY = self->offsetZ = 0.0f;
        self->rwmScale = 1.0f;
        self->rwmOffset = 0.0f;
    }
    return (PyObject *) self;
}

int VolumeObject_init(VolumeObject *self, PyObject *args, PyObject *kwds) {
    // Meta data should be set using properties!
    static const char *kwlist[] = {"format", "data", "dimension", NULL};
    PyObject *format = NULL, *data = NULL, *tmp;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO(III)", (char **) kwlist,
                                     &format,
                                     &data,
                                     &self->dimX, &self->dimY, &self->dimZ))
        return -1;

    if (format) {
        tmp = self->format;
        Py_INCREF(format);
        self->format = format;
        Py_XDECREF(tmp);
    }
    if (data) {
        tmp = self->data;
        Py_INCREF(data);
        self->data = data;
        Py_XDECREF(tmp);
    }
    return 0;
}

void VolumeObject_dealloc(VolumeObject *self) {
    Py_XDECREF(self->format);
    Py_XDECREF(self->data);
    Py_TYPE(self)->tp_free((PyObject *) self);
}


//////////////////////////////////////////////////////////////////
// RenderTargetObject
//////////////////////////////////////////////////////////////////

PyObject* RenderTargetObject_new(PyTypeObject *type, PyObject */*args*/, PyObject */*kwds*/) {
    RenderTargetObject *self;
    self = (RenderTargetObject *) type->tp_alloc(type, 0);
    if (self != NULL) {
        self->internalColorFormat = -1;
        self->internalDepthFormat = -1;
        self->colorTexture = PyList_New(0);
        if (self->colorTexture == NULL) {
            Py_DECREF(self);
            return NULL;
        }
        self->depthTexture = PyList_New(0);
        if (self->depthTexture == NULL) {
            Py_DECREF(self);
            return NULL;
        }
        self->width = self->height = 0u;
    }
    return (PyObject *) self;
}

int RenderTargetObject_init(RenderTargetObject *self, PyObject *args, PyObject *kwds) {
    static const char *kwlist[] = {"internalColorFormat", "internalDepthFormat", "colorTexture", "depthTexture", "width", "height", NULL};
    PyObject *colorTexture = NULL, *depthTexture = NULL, *tmp;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "iiOO(II)", (char **) kwlist,
                                     &self->internalColorFormat,
                                     &self->internalDepthFormat,
                                     &colorTexture,
                                     &depthTexture,
                                     &self->width, &self->height
                                     ))
        return -1;

    if (colorTexture) {
        tmp = self->colorTexture;
        Py_INCREF(colorTexture);
        self->colorTexture = colorTexture;
        Py_XDECREF(tmp);
    }
    if (depthTexture) {
        tmp = self->depthTexture;
        Py_INCREF(depthTexture);
        self->depthTexture = depthTexture;
        Py_XDECREF(tmp);
    }
    return 0;
}

void RenderTargetObject_dealloc(RenderTargetObject *self) {
    Py_XDECREF(self->colorTexture);
    Py_XDECREF(self->depthTexture);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

}
