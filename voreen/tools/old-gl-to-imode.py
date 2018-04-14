#! /usr/bin/python3

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Note:
# This script is still work in progress. So far it does not convert
# gl* Matrix functions to use Matstack.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# This script is intended to ease the process of porting code using
# old compatibility gl calls to API tgt::ImmediateMode.
# However, it is not guaranteed to word in all cases. Therefore make
# sure to grep for warningMsg (defined below) after using this script
# on a number of files and review these changes manually.

# Usage: All parameters are interpreted as paths to files that are
# supposed to be processed by this script.

import sys
from re import sub, search

warningMsg = "TODO: IMODE-CONVERSION CHECK"

replacements = [
        #First replace primitives with matching or similar
        #They will be further replaced by IMode constants
        #['GL_QUADS', 'GL_TRIANGLE_FAN'],
        #['GL_QUAD_STRIP', 'GL_TRIANGLE_STRIP'],
        #['GL_POLYGON', 'GL_TRIANGLE_FAN'],

        [r"glBegin\s*\(GL_([A-Z_]*)\)", r"IMode.begin(tgt::ImmediateMode::\1)"],
        [r"glEnd\s*\(\)", r"IMode.end()"],

        #Vertices etc.
        [r"glVertex[2,3,4][i,f,d]\s*\(", r"IMode.vertex("],
        [r"glVertex([2,3,4])fv\s*\(", r"IMode.vertex(/*" + warningMsg + r"*/ *(tgt::vec\1*)"], #untested
        [r"glVertex([2,3,4])dv\s*\(", r"IMode.vertex(/*" + warningMsg + r"*/ *(tgt::dvec\1*)"], #untested
        [r"glTexCoord[2,3,4][f,d]\s*\(", r"IMode.texcoord("],
        [r"glColor[3,4][f,d]\s*\(", r"IMode.color("],
        [r"glColor([2,3,4])fv\s*\(", r"IMode.color(/*" + warningMsg + r"*/ *(tgt::vec\1*)"], #untested
        [r"glColor([2,3,4])fv\s*\(", r"IMode.color(/*" + warningMsg + r"*/ *(tgt::dvec\1*)"], #untested
        [r"glNormal[3][i,f,d]\s*\(", r"IMode.normal("],

        #Emit warnings and hints on how to port
        #TODO: Maybe this conversion could be done without regexes and counting brackets...
        [r"glVertex[2,3,4]..?\s*\(", r"/*" + warningMsg + r" and use IMode.normalize */ \0"],
        [r"glTexCoord[2,3,4]..?\s*\(", r"/*" + warningMsg + r" and use IMode.normalize */ \0"],
        [r"glColor[2,3,4]..?\s*\(", r"/*" + warningMsg + r" and use IMode.normalize */ \0"],
        [r"glNormal[2,3,4]..?\s*\(", r"/*" + warningMsg + r" and use IMode.normalize */ \0"],


        #Clipping
        [r"glEnable\s*\(\s*GL_CLIP_PLANE(\d+)", r"IMode.enableClipPlane(\1"],
        [r"glDisable\s*\(\s*GL_CLIP_PLANE(\d+)", r"IMode.disableClipPlane(\1"],
        [r"glClipPlane\s*\(\s*GL_CLIP_PLANE(\d+)\s*,\s*", r"IMode.setModelSpaceClipPlaneEquation(\1, /*" + warningMsg + r"*/ *(tgt::dplane*)"],

        #[r"GL_TEXTURE_([1,2,3])D", r"tgt::ImmediateMode::TEX\1D"],

        #TODO matstack stuff
        ]

def needsProcessing(content):
    for replacement in replacements:
        if search(replacement[0], content):
            return True

    return False

def addHeader(content):
    imodeInclude = r'#include "tgt/immediatemode/immediatemode.h"'
    imodeIncludeRegex = r'#include +"tgt/immediatemode/immediatemode\.h"'
    if search(imodeIncludeRegex, content) is None:
        content = sub(r'(#include \".*\")', r"\1\n\n{}".format(imodeInclude), content, count=1)
    #TODO matstack header

    return content

def replaceOldGL(content):
    for replacement in replacements:
        content = sub(replacement[0], replacement[1], content)
    return content


def processFile(filename):
    data = ""
    with open(filename, 'rb') as f:
        try:
            data = f.read().decode('utf8')
        except:
            print("WARNING: Cannot decode {} as utf8, trying iso-8859-1".format(filename))
            data = f.read().decode('iso-8859-1')

    if needsProcessing(data):
        print("Making {} use IMode".format(filename))
        coreProfileCompatibleData = addHeader(replaceOldGL(data))
        with open(filename, 'w') as f:
            f.write(coreProfileCompatibleData)
    else:
        print("{} does not appear to be using old GL".format(filename))

def main():
    for filename in sys.argv[1:]:
        print("Processing file {}".format(filename))
        processFile(filename)
        print("Done: {}\n".format(filename))

if __name__ == "__main__":
    main()
