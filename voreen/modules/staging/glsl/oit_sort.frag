/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "mod_oit.frag"

uniform mat4 viewMatrixInverse_;
uniform mat4 projectionMatrixInverse_;
uniform vec2 screenDimRCP_;

void main() {
    vec2 coord = gl_FragCoord.xy - vec2(0.5);

    uint head = imageLoad(headOffset, ivec2(coord)).x;

    if(head == 0U)
        discard;
    else {
        LinkedListStruct[MAX_DEPTHCOMPLEXITY] list;
        uint[MAX_DEPTHCOMPLEXITY] listIndices;
        int numFragments = 0;

        // copy linked list to local buffer:
        while(head != 0U) {
            list[numFragments] = linkedList_[head];
            listIndices[numFragments] = head;
            numFragments++;
            head = linkedList_[head].next_;
        }

#ifdef BUBBLE_SORT
        for(int j=numFragments; j>1; j--) {
            for(int i=0; i<(j-1); i++) {
                if(list[i].depth_ > list[i+1].depth_) {
                    LinkedListStruct tmp = list[i];
                    list[i] = list[i+1];
                    list[i+1] = tmp;
                }
            }
        }
#endif
#ifdef SELECTION_SORT
        for(int j=0; j<numFragments; j++) {
            int minIndex = j;
            for(int i=j+1; i<numFragments; i++) {
                if(list[i].depth_ < list[minIndex].depth_)
                    minIndex = i;
            }

            if(minIndex != j) {
                LinkedListStruct tmp = list[j];
                list[j] = list[minIndex];
                list[minIndex] = tmp;
            }
        }
#endif
#ifdef INSERTION_SORT
        for(int j=1; j<numFragments; j++) {
            LinkedListStruct toInsert = list[j];
            int hole = j;

            while((hole > 0) && (toInsert.depth_ < list[hole-1].depth_)) {
                list[hole] = list[hole-1];
                hole--;
            }
            list[hole] = toInsert;
        }
#endif

        // write back sorted fragments:
        for(int i=0; i<numFragments; i++) {
            LinkedListStruct cur = list[i];

            if((i+1) < numFragments)
                cur.next_ = listIndices[i+1];
            else
                cur.next_ = 0U;

            linkedList_[listIndices[i]] = cur;
        }

        int vol1Intersections = 0;
        uint pgId1 = 0; // a PG with volumeId = 1
        int vol2Intersections = 0;
        uint pgId2 = 0;
        int vol3Intersections = 0;
        uint pgId3 = 0;
        int vol4Intersections = 0;
        uint pgId4 = 0;

        for(int i=0; i<numFragments; i++) {
            LinkedListStruct cur = list[i];
            ProxyGeometryStruct pg = proxyGeometries_[cur.proxyGeometryId_];

            switch(pg.volumeId_) {
                case 1:
                    vol1Intersections++;
                    pgId1 = cur.proxyGeometryId_;
                    break;
                case 2:
                    vol2Intersections++;
                    pgId2 = cur.proxyGeometryId_;
                    break;
                case 3:
                    vol3Intersections++;
                    pgId3 = cur.proxyGeometryId_;
                    break;
                case 4:
                    vol4Intersections++;
                    pgId4 = cur.proxyGeometryId_;
                    break;
            }
        }

        vec2 p = gl_FragCoord.xy;
        p *= screenDimRCP_;

        // inject entry positions at near plane for camera-inside-volume:
        if((vol1Intersections % 2) != 0) {
            LinkedListStruct lls;

            setNormal(lls, vec3(1.0, 0.0, 0.0)); //TODO
            //setNormal(lls, vec3(0.0)); //TODO
            lls.depth_ = 0.0;
            lls.proxyGeometryId_ = pgId1;

            pushFront(lls);
        }

        if((vol2Intersections % 2) != 0) {
            LinkedListStruct lls;

            setNormal(lls, vec3(0.0));
            lls.depth_ = 0.0;
            lls.proxyGeometryId_ = pgId2;

            pushFront(lls);
        }

        if((vol3Intersections % 2) != 0) {
            LinkedListStruct lls;

            setNormal(lls, vec3(0.0));
            lls.depth_ = 0.0;
            lls.proxyGeometryId_ = pgId3;

            pushFront(lls);
        }

        if((vol4Intersections % 2) != 0) {
            LinkedListStruct lls;

            setNormal(lls, vec3(0.0));
            lls.depth_ = 0.0;
            lls.proxyGeometryId_ = pgId4;

            pushFront(lls);
        }
    }

}
