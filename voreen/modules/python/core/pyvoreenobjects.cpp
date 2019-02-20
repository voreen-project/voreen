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
        self->spacingX = self->spacingY = self->spacingZ = 0.0f;
        self->offsetX = self->offsetY = self->offsetZ = 0.0f;
        self->rwmScale = 1.0f;
        self->rwmOffset = 0.0f;
    }
    return (PyObject *) self;
}

int VolumeObject_init(VolumeObject *self, PyObject *args, PyObject *kwds) {
    static const char *kwlist[] = {"format", "data", "dimension", "spacing", "offset", "rwm scale", "rwmoffset", NULL};
    PyObject *format = NULL, *data = NULL, *tmp;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO(III)(fff)(fff)|ff", (char **) kwlist,
                                     &format,
                                     &data,
                                     &self->dimX, &self->dimY, &self->dimZ,
                                     &self->spacingX, &self->spacingY, &self->spacingZ,
                                     &self->offsetX, &self->offsetY, &self->offsetZ,
                                     &self->rwmScale, &self->rwmOffset))
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
        self->width = self->height = 0u;
    }
    return (PyObject *) self;
}

int RenderTargetObject_init(RenderTargetObject *self, PyObject *args, PyObject *kwds) {
    static const char *kwlist[] = {"format", "data", "dimension", NULL};
    PyObject *format = NULL, *data = NULL, *tmp;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO(II)", (char **) kwlist,
                                     &format,
                                     &data,
                                     &self->width, &self->height
                                     ))
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

void RenderTargetObject_dealloc(RenderTargetObject *self) {
    Py_XDECREF(self->format);
    Py_XDECREF(self->data);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

}