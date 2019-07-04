#! /usr/bin/env python3

# A python script to automatically generate the deletetemplates.h file.
# The resulting header file includes functions to check if a voxel matches
# a deleting template defined by Ma and Sonka (1996) and improved by
# Wang and Basu (2004).

import numpy as np
def readTemplate(template, pos):
    v = template[pos[2]+1, pos[1]+1, pos[0]+1]
    return v
def writeTemplate(template, pos, value):
    template[pos[2]+1, pos[1]+1, pos[0]+1] = value

def transformTemplate(template, transformation):
    ntemp = np.empty((3,3,3))
    for z in range(-1,2):
        for y in range(-1,2):
            for x in range(-1,2):
                pos = np.array([x,y,z])
                writeTemplate(ntemp, np.dot(transformation,pos), readTemplate(template, pos))
    return ntemp


def rotateX(s):
    return np.array([
        [ 1, 0, 0],
        [ 0, 0,-s],
        [ 0, s, 0],
        ])
def rotateY(s):
    return np.array([
        [ 0, 0, s],
        [ 0, 1, 0],
        [-s, 0, 0],
        ])
def rotateZ(s):
    return np.array([
        [ 0,-s, 0],
        [ s, 0, 0],
        [ 0, 0, 1],
        ])
def mirrorX():
    return np.array([
        [-1, 0, 0],
        [ 0, 1, 0],
        [ 0, 0, 1],
        ])
def mirrorY():
    return np.array([
        [ 1, 0, 0],
        [ 0,-1, 0],
        [ 0, 0, 1],
        ])
def mirrorZ():
    return np.array([
        [ 1, 0, 0],
        [ 0, 1, 0],
        [ 0, 0,-1],
        ])

coreA = np.array([
    [[ 0, 0, 0], [ 0, 1, 0], [ 0, 0, 0]],
    [[ 0, 0, 0], [ 0, 1, 0], [ 0, 0, 0]],
    [[-1,-1,-1], [-1,-1,-1], [-1,-1,-1]],
    ])
coreB = np.array([
    [[ 0, 0, 0], [ 0, 1, 0], [ 0, 0, 0]],
    [[ 0, 0,-1], [ 1, 1,-1], [ 0, 0,-1]],
    [[ 0,-1,-1], [ 0,-1,-1], [ 0,-1,-1]],
    ])
coreC = np.array([
    [[ 0, 0, 0], [ 0, 1, 0], [ 0, 0, 0]],
    [[ 0, 1, 0], [ 1, 1,-1], [ 0,-1,-1]],
    [[ 0, 0, 0], [ 0,-1,-1], [ 0,-1,-1]],
    ])
coreD1 = np.array([ #d1-1
    [[-1,-1,-1], [-1,-1,-1], [-1,-1,-1]],
    [[ 0,-1, 0], [-1, 1,-1], [-1,-1,-1]],
    [[ 0, 1, 0], [ 0,-1, 0], [-1,-1,-1]],
    ])
coreD2 = np.array([ #d1-2
    [[-1,-1, 0], [-1,-1, 0], [-1,-1, 0]],
    [[ 0,-1, 0], [-1, 1, 1], [-1,-1, 0]],
    [[ 0, 1, 0], [ 0,-1, 0], [-1,-1, 0]],
    ])
coreD3 = transformTemplate(coreD2, mirrorX()) #d1-3

templateNamesString = "const std::string templateNames[] = {\n"
arrayString = "const int8_t deleteTemplates[][3][3][3] = {\n"
switchString = "    switch(templateNum) {\n"
templateNum = 0
def genCheckOutsideSingle(template, x, y, z):
    pos = np.array([x, y, z])
    if readTemplate(template, pos) == 1 and (x<0 or y<0 or z<0) and x+y+z==-1:
        return """
            if(readRelative(vol,p,{0},{1},{2}) == -1) {{
                return false;
            }}
        """.format(2*x,2*y,2*z)
    return ""

def genCheckOutsideRing(template, x, y, z):
    pos = np.array([x, y, z])
    if readTemplate(template, pos) == 1:
        if x==0 and y<0 and z!= 0:
            return """
            if(   readRelative(vol,p,0,{0}*0,{1}*2) == -1
               && readRelative(vol,p,0,{0}*1,{1}*2) == -1
               && readRelative(vol,p,0,{0}*2,{1}*2) == -1
               && readRelative(vol,p,0,{0}*2,{1}*1) == -1
               && readRelative(vol,p,0,{0}*2,{1}*0) == -1)
            {{
                return false;
            }}
            """.format(y,z)
        if y==0 and x<0 and z!= 0:
            return """
            if(   readRelative(vol,p,{0}*0,0,{1}*2) == -1
               && readRelative(vol,p,{0}*1,0,{1}*2) == -1
               && readRelative(vol,p,{0}*2,0,{1}*2) == -1
               && readRelative(vol,p,{0}*2,0,{1}*1) == -1
               && readRelative(vol,p,{0}*2,0,{1}*0) == -1)
            {{
                return false;
            }}
            """.format(x,z)
        if z==0 and y<0 and x!=0:
            return """
            if(   readRelative(vol,p,{0}*0,{1}*2,0) == -1
               && readRelative(vol,p,{0}*1,{1}*2,0) == -1
               && readRelative(vol,p,{0}*2,{1}*2,0) == -1
               && readRelative(vol,p,{0}*2,{1}*1,0) == -1
               && readRelative(vol,p,{0}*2,{1}*0,0) == -1)
            {{
                return false;
            }}
            """.format(x,y)
    return ""

def add(template, identifier):
    global arrayString, templateNum, switchString, templateNamesString
    checkString = ""


    templateNamesString += '"' + identifier + '",\n'
    arrayString += "{ //"+identifier +"\n"
    for z in range(-1,2):
        arrayString += "{"
        for y in range(-1,2):
            arrayString += "{"
            for x in range(-1,2):
                pos = np.array([x,y,z])
                val = readTemplate(template, pos)
                #add concrete value to array string
                arrayString += ("%2d" % val) + ","

                #check if we have to add rules to switchstring
                if(identifier[0] == 'd'):
                    checkString += genCheckOutsideRing(template, x, y, z)
                else:
                    checkString += genCheckOutsideSingle(template, x, y, z)
            arrayString += "},"
        arrayString += "},\n"
    arrayString += "},\n\n"

    if(checkString):
        switchString += "       case " + str(templateNum) + ": // "+identifier
        switchString += checkString
        switchString += "break;\n"
    templateNum += 1


###################### A ###############################

#a1-6
add(transformTemplate(coreA, rotateY(-1)), "a1")
add(transformTemplate(coreA, rotateY( 1)), "a2")
add(transformTemplate(coreA, rotateX( 1)), "a3")
add(transformTemplate(coreA, rotateX(-1)), "a4")
add(transformTemplate(coreA, mirrorZ(  )), "a5")
add(coreA, "a6")

###################### B ###############################

#b1-b4
rot = rotateX( 1)         ; add(transformTemplate(coreB, rot), "b1 ")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreB, rot), "b2 ")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreB, rot), "b3 ")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreB, rot), "b4 ")

#b9, b6, b10, b5
rot = mirrorZ(  )         ; add(transformTemplate(coreB, rot), "b9 ")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreB, rot), "b6 ")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreB, rot), "b10")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreB, rot), "b5 ")

#b11, b8, b12, b7
add(coreB, "b11")
rot = rotateZ(-1)         ; add(transformTemplate(coreB, rot), "b8 ")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreB, rot), "b12")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreB, rot), "b7 ")

###################### C ###############################

#c4, c1, c2, c3
rot = mirrorZ(  )         ; add(transformTemplate(coreC, rot), "c4")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreC, rot), "c1")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreC, rot), "c2")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreC, rot), "c3")

#c8, c5, c6, c7
add(coreC, "c8")
rot = rotateZ(-1)         ; add(transformTemplate(coreC, rot), "c5")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreC, rot), "c6")
rot = rotateZ(-1).dot(rot); add(transformTemplate(coreC, rot), "c7")

###################### D1 ##############################

#d1, d2, d10, d9
add(coreD1, "d2 -1")
rot = rotateX( 1)         ; add(transformTemplate(coreD1, rot), "d2 -1")
rot = rotateX( 1).dot(rot); add(transformTemplate(coreD1, rot), "d10-1")
rot = rotateX( 1).dot(rot); add(transformTemplate(coreD1, rot), "d9 -1")

#d3, d4, d7, d8
rot = rotateY(-1)         ; add(transformTemplate(coreD1, rot), "d3 -1")
rot = rotateZ( 1).dot(rot); add(transformTemplate(coreD1, rot), "d4 -1")
rot = rotateZ( 1).dot(rot); add(transformTemplate(coreD1, rot), "d7 -1")
rot = rotateZ( 1).dot(rot); add(transformTemplate(coreD1, rot), "d8 -1")

#d5, d6, d11, d12
rot = rotateZ(-1)         ; add(transformTemplate(coreD1, rot), "d5 -1")
rot = rotateY(-1).dot(rot); add(transformTemplate(coreD1, rot), "d6 -1")
rot = rotateY(-1).dot(rot); add(transformTemplate(coreD1, rot), "d11-1")
rot = rotateY(-1).dot(rot); add(transformTemplate(coreD1, rot), "d12-1")

###################### D2 ##############################

#d1, d2, d10, d9
add(coreD2, "d1 -2")
rot = rotateX( 1)         ; add(transformTemplate(coreD2, rot), "d2 -2")
rot = rotateX( 1).dot(rot); add(transformTemplate(coreD2, rot), "d10-2")
rot = rotateX( 1).dot(rot); add(transformTemplate(coreD2, rot), "d9 -2")

#d3, d4, d7, d8
rot = rotateY(-1)         ; add(transformTemplate(coreD2, rot), "d3 -2")
rot = rotateZ( 1).dot(rot); add(transformTemplate(coreD2, rot), "d4 -2")
rot = rotateZ( 1).dot(rot); add(transformTemplate(coreD2, rot), "d7 -2")
rot = rotateZ( 1).dot(rot); add(transformTemplate(coreD2, rot), "d8 -2")

#d5, d6, d11, d12
rot = rotateZ(-1)         ; add(transformTemplate(coreD2, rot), "d5 -2")
rot = rotateY(-1).dot(rot); add(transformTemplate(coreD2, rot), "d6 -2")
rot = rotateY(-1).dot(rot); add(transformTemplate(coreD2, rot), "d11-2")
rot = rotateY(-1).dot(rot); add(transformTemplate(coreD2, rot), "d12-2")

###################### D3 ##############################

#d1, d2, d10, d9
add(coreD3, "d1 -3")
rot = rotateX( 1)         ; add(transformTemplate(coreD3, rot), "d2 -3")
rot = rotateX( 1).dot(rot); add(transformTemplate(coreD3, rot), "d10-3")
rot = rotateX( 1).dot(rot); add(transformTemplate(coreD3, rot), "d9 -3")

#d3, d4, d7, d8
rot = rotateY(-1)         ; add(transformTemplate(coreD3, rot), "d3 -3")
rot = rotateZ( 1).dot(rot); add(transformTemplate(coreD3, rot), "d4 -3")
rot = rotateZ( 1).dot(rot); add(transformTemplate(coreD3, rot), "d7 -3")
rot = rotateZ( 1).dot(rot); add(transformTemplate(coreD3, rot), "d8 -3")

#d5, d6, d11, d12
rot = rotateZ(-1)         ; add(transformTemplate(coreD3, rot), "d5 -3")
rot = rotateY(-1).dot(rot); add(transformTemplate(coreD3, rot), "d6 -3")
rot = rotateY(-1).dot(rot); add(transformTemplate(coreD3, rot), "d11-3")
rot = rotateY(-1).dot(rot); add(transformTemplate(coreD3, rot), "d12-3")

templateNamesString +="};"
arrayString += "};"
switchString += """
    }
    return true;
"""

header = """
#ifndef VRN_VOLUMETHINNING_TEMPLATES_H
#define VRN_VOLUMETHINNING_TEMPLATES_H

#include "tgt/vector.h"
#include "volumethinning.h"

namespace { // anonymous name space
""" + arrayString + """
""" + templateNamesString + """
int8_t readRelative(const voreen::VolumeMask* vol, const tgt::svec3& p, int dx, int dy, int dz);

bool fitsExpandedTemplate(const voreen::VolumeMask* vol, tgt::svec3 p, size_t templateNum) {
""" + switchString + """
}
} // end of anonymous name space
#endif  // VRN_VOLUMETHINNING_TEMPLATES_H
"""
with open('deletetemplates.h', 'w') as f:
    f.write(header)
