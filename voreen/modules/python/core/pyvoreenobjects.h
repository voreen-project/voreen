/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_PYVOREENOBJECTS_H
#define VRN_PYVOREENOBJECTS_H

#include <Python.h>
#include <structmember.h>

namespace voreen {

//////////////////////////////////////////////////////////////////
// VolumeObject
//////////////////////////////////////////////////////////////////

typedef struct {
    PyObject_HEAD
    PyObject* format;
    //unsigned int numChannels; // Encoded in format
    PyObject* data;
    unsigned int dimX, dimY, dimZ;
    float spacingX, spacingY, spacingZ;
    float offsetX, offsetY, offsetZ;
    float rwmScale;
    float rwmOffset;
} VolumeObject;

/*
 * The following struct defines the python equivalents the members of the VolumeObject struct.
 */
static PyMemberDef VolumeObject_members[] = {
    {(char*)"format",      T_OBJECT_EX,    offsetof(VolumeObject, format),     0, (char*)"Format: (Vector(2|3|4))?(float | double | u?int(8|16|32|64))" },
    //{(char*)"numChannels", T_OBJECT_EX,    offsetof(VolumeObject, numChannels),0, (char*)"Number of Channels [1, 4]"},
    {(char*)"data",        T_OBJECT_EX,    offsetof(VolumeObject, data),       0, (char*)"data"       },
    {(char*)"dimX",        T_UINT,         offsetof(VolumeObject, dimX),       0, (char*)"Dimension X"},
    {(char*)"dimY",        T_UINT,         offsetof(VolumeObject, dimY),       0, (char*)"Dimension Y"},
    {(char*)"dimZ",        T_UINT,         offsetof(VolumeObject, dimZ),       0, (char*)"Dimension Z"},
    {(char*)"spacingX",    T_FLOAT,        offsetof(VolumeObject, spacingX),   0, (char*)"Spacing X"  },
    {(char*)"spacingY",    T_FLOAT,        offsetof(VolumeObject, spacingY),   0, (char*)"Spacing Y"  },
    {(char*)"spacingZ",    T_FLOAT,        offsetof(VolumeObject, spacingZ),   0, (char*)"Spacing Z"  },
    {(char*)"offsetX",     T_FLOAT,        offsetof(VolumeObject, offsetX),    0, (char*)"Offset X"   },
    {(char*)"offsetY",     T_FLOAT,        offsetof(VolumeObject, offsetY),    0, (char*)"Offset Y"   },
    {(char*)"offsetZ",     T_FLOAT,        offsetof(VolumeObject, offsetZ),    0, (char*)"Offset Z"   },
    {(char*)"rwmScale",    T_FLOAT,        offsetof(VolumeObject, rwmScale),   0, (char*)"Real World Mapping Scale"  },
    {(char*)"rwmOffset",   T_FLOAT,        offsetof(VolumeObject, rwmOffset),  0, (char*)"Real World Mapping Offset" },
    {NULL}  /* Sentinel */
};

/*
 * Construction and destruction.
 */
PyObject* VolumeObject_new(PyTypeObject *type, PyObject */*args*/, PyObject */*kwds*/) ;
int VolumeObject_init(VolumeObject *self, PyObject *args, PyObject *kwds);
void VolumeObject_dealloc(VolumeObject *self);

/*
 * The following struct defines the python equivalent for a Voreen Volume.
 * Currently, the implementation is very rudimentary and only supports:
 *   - format
 *   - data (single channel, float, only!)
 *   - dimension
 *   - spacing
 *   - offset
 *   - real world mapping scaling
 *   - real world mapping offset
 */
static PyTypeObject VolumeObjectType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "voreen.Volume",                        /*tp_name*/
    sizeof(VolumeObject),                   /*tp_basicsize*/
    0,                                      /*tp_itemsize*/
    (destructor) VolumeObject_dealloc,      /*tp_dealloc*/
    NULL,                                   /*tp_print*/
    NULL,                                   /*tp_getattr*/
    NULL,                                   /*tp_setattr*/
    NULL,                                   /*tp_as_async*/
    NULL,                                   /*tp_repr*/
    NULL,                                   /*tp_as_number*/
    NULL,                                   /*tp_as_sequence*/
    NULL,                                   /*tp_as_mapping*/
    NULL,                                   /*tp_hash*/
    NULL,                                   /*tp_call*/
    NULL,                                   /*tp_str*/
    NULL,                                   /*tp_getattro*/
    NULL,                                   /*tp_setattro*/
    NULL,                                   /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,                     /*tp_flags*/
    "Volume, contained by a port",          /*tp_doc*/
    NULL,                                   /*tp_traverse*/
    NULL,                                   /*tp_clear*/
    NULL,                                   /*tp_richcompare*/
    0,                                      /*tp_weaklistoffset*/
    NULL,                                   /*tp_iter*/
    NULL,                                   /*tp_iternext*/
    NULL,                                   /*tp_methods*/
    VolumeObject_members,                   /*tp_members*/
    NULL,                                   /*tp_getset*/
    NULL,                                   /*tp_base*/
    NULL,                                   /*tp_dict*/
    NULL,                                   /*tp_descr_get*/
    NULL,                                   /*tp_descr_set*/
    0,                                      /*tp_dictoffset*/
    (initproc) VolumeObject_init,           /*tp_init*/
    NULL,                                   /*tp_alloc*/
    VolumeObject_new,                       /*tp_new*/
    NULL,                                   /*tp_free*/
    NULL,                                   /*tp_is_gc*/
    NULL,                                   /*tp_bases*/
    NULL,                                   /*tp_mro*/
    NULL,                                   /*tp_cache*/
    NULL,                                   /*tp_subclasses*/
    NULL,                                   /*tp_weaklist*/
    NULL,                                   /*tp_del*/
    0,                                      /*tp_version_tag*/
    NULL                                    /*tp_finalize*/
};


//////////////////////////////////////////////////////////////////
// RenderTargetObject
//////////////////////////////////////////////////////////////////

typedef struct {
    PyObject_HEAD
    int internalColorFormat;
    int internalDepthFormat;
    PyObject* colorTexture;
    PyObject* depthTexture;
    unsigned int width, height;
} RenderTargetObject;

/*
 * The following struct defines the python equivalents the members of the RenderTargetObject struct.
 */
static PyMemberDef RenderTargetObject_members[] = {
    {(char*)"internalColorFormat",  T_INT,          offsetof(RenderTargetObject, internalColorFormat),  0, (char*)"Internal color format"   },
    {(char*)"internalDepthFormat",  T_INT,          offsetof(RenderTargetObject, internalDepthFormat),  0, (char*)"Internal depth format"   },
    {(char*)"colorTexture",         T_OBJECT_EX,    offsetof(RenderTargetObject, colorTexture),         0, (char*)"Color Texture data"      },
    {(char*)"depthTexture",         T_OBJECT_EX,    offsetof(RenderTargetObject, depthTexture),         0, (char*)"Depth Texture data"      },
    {(char*)"width",                T_UINT,         offsetof(RenderTargetObject, width),                0, (char*)"Width"                   },
    {(char*)"height",               T_UINT,         offsetof(RenderTargetObject, height),               0, (char*)"Height"                  },
    {NULL}  /* Sentinel */
};

/*
 * Construction and destruction.
 */
PyObject* RenderTargetObject_new(PyTypeObject *type, PyObject* /*args*/, PyObject* /*kwds*/) ;
int RenderTargetObject_init(RenderTargetObject *self, PyObject *args, PyObject *kwds);
void RenderTargetObject_dealloc(RenderTargetObject *self);

/*
 * The following struct defines the python equivalent for a Voreen RenderTarget.
 * Currently, the implementation is very rudimentary and only supports:
 *   - format
 *   - data (single channel, float, only!)
 *   - dimension
 */
static PyTypeObject RenderTargetObjectType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "voreen.RenderTarget",                  /*tp_name*/
    sizeof(RenderTargetObject),             /*tp_basicsize*/
    0,                                      /*tp_itemsize*/
    (destructor) RenderTargetObject_dealloc,/*tp_dealloc*/
    NULL,                                   /*tp_print*/
    NULL,                                   /*tp_getattr*/
    NULL,                                   /*tp_setattr*/
    NULL,                                   /*tp_as_async*/
    NULL,                                   /*tp_repr*/
    NULL,                                   /*tp_as_number*/
    NULL,                                   /*tp_as_sequence*/
    NULL,                                   /*tp_as_mapping*/
    NULL,                                   /*tp_hash*/
    NULL,                                   /*tp_call*/
    NULL,                                   /*tp_str*/
    NULL,                                   /*tp_getattro*/
    NULL,                                   /*tp_setattro*/
    NULL,                                   /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,                     /*tp_flags*/
    "RenderTarget, contained by a port",    /*tp_doc*/
    NULL,                                   /*tp_traverse*/
    NULL,                                   /*tp_clear*/
    NULL,                                   /*tp_richcompare*/
    0,                                      /*tp_weaklistoffset*/
    NULL,                                   /*tp_iter*/
    NULL,                                   /*tp_iternext*/
    NULL,                                   /*tp_methods*/
    RenderTargetObject_members,             /*tp_members*/
    NULL,                                   /*tp_getset*/
    NULL,                                   /*tp_base*/
    NULL,                                   /*tp_dict*/
    NULL,                                   /*tp_descr_get*/
    NULL,                                   /*tp_descr_set*/
    0,                                      /*tp_dictoffset*/
    (initproc) RenderTargetObject_init,     /*tp_init*/
    NULL,                                   /*tp_alloc*/
    RenderTargetObject_new,                 /*tp_new*/
    NULL,                                   /*tp_free*/
    NULL,                                   /*tp_is_gc*/
    NULL,                                   /*tp_bases*/
    NULL,                                   /*tp_mro*/
    NULL,                                   /*tp_cache*/
    NULL,                                   /*tp_subclasses*/
    NULL,                                   /*tp_weaklist*/
    NULL,                                   /*tp_del*/
    0,                                      /*tp_version_tag*/
    NULL                                    /*tp_finalize*/
};

} // namespace

#endif // VRN_PYVOREENOBJECTS_H
