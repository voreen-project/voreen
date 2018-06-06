/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "pfskelwrapper.h"
#include "../pfskelmodule.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/roi/roiraster.h"
#include "modules/core/io/datvolumewriter.h"
#include "modules/core/io/datvolumereader.h"
#include "voreen/core/voreenapplication.h"

#include "tgt/filesystem.h"

namespace voreen {

using tgt::mat4;
using tgt::vec3;
using tgt::svec3;
using tgt::ivec3;

PFSkelWrapper::PFSkelWrapper()
{
}

std::vector<tgt::vec3> PFSkelWrapper::calculateSkeletonGeometry(const Volume* v, int fieldst, int highd, int branchTh) {
    if(!v)
        return std::vector<tgt::vec3>();

    const VolumeRAM_UInt8* source = v->getRepresentation<VolumeRAM_UInt8>();
    if(!source)
        return std::vector<tgt::vec3>();

    ivec3 dims = source->getDimensions();
    VolumeRAM_UInt8* target = new VolumeRAM_UInt8(dims+ivec3(2));
    target->clear();

    tgt::ivec3 i;
    for(i.z=0; i.z<dims.z; i.z++) {
        for(i.y=0; i.y<dims.y; i.y++) {
            for(i.x=0; i.x<dims.x; i.x++) {
                if(source->voxel(i))
                    target->voxel(i+ivec3(1)) = 255;
            }
        }
    }

    Volume* vh = new Volume(target, vec3(1.0f), vec3(0.0f));
    oldVolumePosition(vh);

    std::string basename = VoreenApplication::app()->getTemporaryPath("PFSkelWrapper");
    basename = FileSys.absolutePath(basename) + "/";
    FileSys.createDirectoryRecursive(basename);

    DatVolumeWriter dw;
    dw.write(basename+"orig.dat", vh);

    Graph g = calculateSkeleton(basename, "orig.raw", dims+ivec3(2), fieldst, highd, branchTh);
    tgt::mat4 vtop = v->getVoxelToWorldMatrix();

    std::vector<tgt::vec3> result;
    for(size_t i=0; i<g.getNumNodes(); i++) {
        result.push_back(vtop * g.getNode(i)->getPosition());
    }

    return result;
}

std::vector<tgt::vec3> PFSkelWrapper::calculateSkeletonVoxels(const Volume* v, int fieldst, int highd, int branchTh) {
    if(!v)
        return std::vector<tgt::vec3>();

    const VolumeRAM_UInt8* source = v->getRepresentation<VolumeRAM_UInt8>();
    if(!source)
        return std::vector<tgt::vec3>();

    ivec3 dims = source->getDimensions();
    VolumeRAM_UInt8* target = new VolumeRAM_UInt8(dims+ivec3(2));
    target->clear();

    tgt::ivec3 i;
    for(i.z=0; i.z<dims.z; i.z++) {
        for(i.y=0; i.y<dims.y; i.y++) {
            for(i.x=0; i.x<dims.x; i.x++) {
                if(source->voxel(i))
                    target->voxel(i+ivec3(1)) = 255;
            }
        }
    }

    Volume* vh = new Volume(target, vec3(1.0f), vec3(0.0f));
    oldVolumePosition(vh);

    std::string basename = VoreenApplication::app()->getTemporaryPath("PFSkelWrapper");
    basename = FileSys.absolutePath(basename) + "/";
    FileSys.createDirectoryRecursive(basename);

    DatVolumeWriter dw;
    dw.write(basename+"orig.dat", vh);

    Graph g = calculateSkeleton(basename, "orig.raw", dims+ivec3(2), fieldst, highd, branchTh);

    std::vector<tgt::vec3> result;
    for(size_t i=0; i<g.getNumNodes(); i++) {
        result.push_back(g.getNode(i)->getPosition());
    }

    return result;
}

Graph PFSkelWrapper::calculateSkeleton(ROIRaster* rr, int fieldst, int highd, int branchTh) {
    if(rr) {
        ivec3 dims = rr->getDimensions();
        VolumeRAM_UInt8* target = new VolumeRAM_UInt8(dims+ivec3(2));

        tgt::ivec3 i;
        for(i.x=0; i.x<dims.x+2; i.x++) {
            for(i.y=0; i.y<dims.y+2; i.y++) {
                for(i.z=0; i.z<dims.z+2; i.z++) {
                    target->voxel(i) = 0;
                }
            }
        }

        for(i.x=0; i.x<dims.x; i.x++) {
            for(i.y=0; i.y<dims.y; i.y++) {
                for(i.z=0; i.z<dims.z; i.z++) {
                    if(rr->getVoxel(i))
                        target->voxel(i+ivec3(1)) = 255;
                    else
                        target->voxel(i+ivec3(1)) = 0;
                }
            }
        }

        Volume* vh = new Volume(target, vec3(1.0f), vec3(0.0f));
        oldVolumePosition(vh);

        std::string basename = VoreenApplication::app()->getTemporaryPath("PFSkelWrapper");
        basename = FileSys.absolutePath(basename) + "/";
        FileSys.createDirectoryRecursive(basename);

        DatVolumeWriter dw;
        dw.write(basename+"orig.dat", vh);

        Graph g = calculateSkeleton(basename, "orig.raw", dims+ivec3(2), fieldst, highd, branchTh);
        tgt::mat4 vtop = tgt::mat4::createScale(rr->getGrid().getSpacing()) * tgt::mat4::createTranslation(rr->getLLF());
            //rr->getGrid()->get
        for(size_t i=0; i<g.getNumNodes(); i++) {
            GraphNode* gn = g.getNode(i);

            vec3 p = vtop * gn->getPosition();
            gn->setPosition(p);
        }

        return g;
    }

    return Graph();
}

Graph PFSkelWrapper::calculateSkeleton(const std::string& basename, const std::string& rawFile, tgt::ivec3 dims, int fieldst, int highd, int branchTh) {
    std::string call;

    std::stringstream fs_str;
    fs_str << fieldst;
    std::stringstream hd_str;
    hd_str << highd;

    std::string pfskelDir = dynamic_cast<const PFSkelModule*>(VoreenApplication::app()->getModule("pfskel"))->getPFSkelDir();

    std::string skelfn = basename + "skel" + fs_str.str() + "-" + hd_str.str();
    std::string critfn = basename + "crit" + fs_str.str();
    std::string highdfn = basename + "highd" + fs_str.str() + "-" + hd_str.str();
    std::string potffn = basename + "potf" + fs_str.str();

    std::stringstream size_call;
    size_call << " " << dims.x << " " << dims.y << " " << dims.z << " ";

    call = pfskelDir+"/PotField/driver " + basename + rawFile + size_call.str() + fs_str.str() + " " + potffn;
    system(call.c_str());

    call = pfskelDir+"/CritPts/driver " + basename + rawFile + size_call.str() + potffn + " " + critfn;
    system(call.c_str());

    call = pfskelDir+"/HighDiverg/driver " + basename + rawFile + size_call.str() + potffn + " " + hd_str.str() + " " + highdfn;
    system(call.c_str());

    call = pfskelDir+"/StreamLn/driver " + potffn + size_call.str() + critfn + " " + highdfn + " " + skelfn;
    system(call.c_str());

    std::vector<std::vector<vec3> > segments = readSkeleton(skelfn);
    std::vector<vec3> critical = readCritical(critfn);
    return postProcessSkeleton(segments, critical, branchTh);
}

std::vector<std::vector<vec3> > PFSkelWrapper::readSkeleton(std::string fname) {
    std::ifstream inFile;

    inFile.open(fname.c_str());
    if (!inFile) {
        std::cout << "Unable to open file";
        return std::vector<std::vector<vec3> >();
    }

    float x;
    int y;
    std::vector<SkelPoint> v;
    int maxId = -1;

    while (true) {
        SkelPoint sp;

        //-1.0 to remove padding:
        if(!(inFile >> x))
            goto end;
        sp.pos_.x = x-1.f;

        if(!(inFile >> x))
            goto end;
        sp.pos_.y = x-1.f;

        if(!(inFile >> x))
            goto end;
        sp.pos_.z = x-1.f;

        if(!(inFile >> y))
            goto end;
        sp.seg_ = y;
        if(y > maxId)
            maxId = y;

        //ignored float:
        if(!(inFile >> x))
            goto end;

        v.push_back(sp);
    }
end:inFile.close();
    std::cout << "read " << v.size() << " points in " << maxId << " segments." << std::endl;

    std::vector<std::vector<vec3> > segments;

    for(int i=0; i<=maxId; ++i) {
        segments.push_back(std::vector<vec3>());
    }

    for(size_t i=0; i<v.size(); ++i) {
        segments[v[i].seg_].push_back(v[i].pos_);
    }
    //for(size_t i=0; i<segments.size(); ++i) {
        //std::cout << "segment " << i << " " << segments[i].size() << std::endl;
    //}
    return segments;
}

std::vector<vec3> PFSkelWrapper::readCritical(std::string fname) {
    std::ifstream inFile;

    std::vector<vec3> critPoints;

    inFile.open(fname.c_str());
    if (!inFile) {
        std::cout << "Unable to open file";
        return critPoints;
    }

    float x;
    int y;

    while (true) {
        vec3 p;

        //-1.0 to remove padding:
        //Coordinate:
        if(!(inFile >> x))
            goto end;
        p.x = x-1.f;

        if(!(inFile >> x))
            goto end;
        p.y = x-1.f;

        if(!(inFile >> x))
            goto end;
        p.z = x-1.f;

        //---------------
        //Type(ignored):
        if(!(inFile >> y))
            goto end;
        //---------------
        //Jacobian matrix eigenvectors and values (ignored):
        for(int k=0; k<12; ++k) {
            if(!(inFile >> x))
                goto end;
        }

        critPoints.push_back(p);
    }
end:    inFile.close();
        std::cout << "read " << critPoints.size() << " critical points" << std::endl;

        return critPoints;
}

bool areNeigbors(GraphNode* p1, GraphNode* p2) {
    tgt::ivec3 a = p1->getPosition();
    tgt::ivec3 b = p2->getPosition();
    return ( (abs(a.x - b.x) <= 1) && (abs(a.y - b.y) <= 1) && (abs(a.z - b.z) <= 1) );
}

bool isRedundant(GraphNode* p) {
    if(p->getNumNeighbors() != 2)
        return false;

    return areNeigbors(p->getNeighbor(0), p->getNeighbor(1));
}

Graph PFSkelWrapper::postProcessSkeleton(const std::vector<std::vector<vec3> >& segments, const std::vector<tgt::vec3>& critPoints, int branchTh) {
    std::vector<std::vector<GraphNode*> > centerline;

    for(size_t i=0; i<segments.size(); ++i) {
        if(segments[i].size() > 0) {
            centerline.push_back(std::vector<GraphNode*>());
            centerline.back().push_back(new GraphNode(segments[i][0]));

            for(size_t j=1; j<segments[i].size(); ++j) {
                GraphNode* cur = new GraphNode(segments[i][j]);

                centerline.back().back()->connect(cur);
                centerline.back().push_back(cur);
            }
        }
    }

    for(size_t i=0; i<critPoints.size(); i++) {
        centerline.push_back(std::vector<GraphNode*>());
        centerline.back().push_back(new GraphNode(critPoints[i]));
    }

    //construct minimal spanning tree using prims algorithm:
    std::vector<GraphNode*> leafs;

    //start with last segment:
    std::vector<GraphNode*> centerTree;
    centerTree.insert(centerTree.end(), centerline.back().begin(), centerline.back().end());
    if(centerline.back().size() == 1) {
        leafs.push_back(centerline.back().back());
    }
    else {
        leafs.push_back(centerline.back().front());
        leafs.push_back(centerline.back().back());
    }
    centerline.pop_back();

    while(!centerline.empty()) {
        float mindist = FLT_MAX;
        int mini = -1;
        int minj = -1;
        bool start = true;
        for(size_t i=0; i<leafs.size(); ++i) {
            vec3 p = leafs[i]->getPosition();
            for(size_t j=0; j<centerline.size(); ++j) {
                if(centerline[j].size() == 1) {
                    if(distance(p, centerline[j][0]->getPosition()) < mindist) {
                        mindist = distance(p, centerline[j][0]->getPosition());
                        mini = static_cast<int>(i);
                        minj = static_cast<int>(j);
                        start = true;
                    }
                }
                else {
                    if(distance(p, centerline[j][0]->getPosition()) < mindist) {
                        mindist = distance(p, centerline[j][0]->getPosition());
                        mini = static_cast<int>(i);
                        minj = static_cast<int>(j);
                        start = true;
                    }
                    if(distance(p, centerline[j].back()->getPosition()) < mindist) {
                        mindist = distance(p, centerline[j].back()->getPosition());
                        mini = static_cast<int>(i);
                        minj = static_cast<int>(j);
                        start = false;
                    }
                }
            }
        }
        //found minimum distance:
        centerTree.insert(centerTree.end(), centerline[minj].begin(), centerline[minj].end());
        if(centerline[minj].size() == 1) {
            leafs.push_back(centerline[minj].back());
            centerline[minj].back()->connect(leafs[mini]);
        }
        else {
            leafs.push_back(centerline[minj].front());
            leafs.push_back(centerline[minj].back());
            if(start)
                centerline[minj].front()->connect(leafs[mini]);
            else
                centerline[minj].back()->connect(leafs[mini]);
        }
        centerline.erase(centerline.begin()+minj);
    }

    //remove dead critical points:
    bool clean = false;
    while(!clean) {
        clean = true;
        int shortest = static_cast<int>(centerTree.size())+1;
        int s = -1;

        for(size_t i=0; i<centerTree.size(); ++i) {
            if(centerTree[i]->getNumNeighbors() == 1) {
                GraphNode* cur = centerTree[i];
                GraphNode* next = cur->getNeighbor(0);
                int steps = 1;

                while(next->getNumNeighbors() == 2) {
                    steps++;
                    GraphNode* old = cur;
                    cur = next;
                    if(cur->getNeighbor(0) == old)
                        next = cur->getNeighbor(1);
                    else
                        next = cur->getNeighbor(0);

                }
                //std::cout << "leaf " << steps << std::endl;
                if(steps < shortest) {
                    shortest = steps;
                    s = static_cast<int>(i);
                }
            }
        }
        //std::cout << "shortest: " << shortest << std::endl;
        if((s > -1) && (shortest < branchTh)) {
            clean = false;

            GraphNode* cur = centerTree[s];

            cur->removeLinks();
            delete cur;
            centerTree.erase(centerTree.begin()+s);
        }
    }

    //clean = false;
    //while(!clean) {
        //clean = true;

        //for(size_t i=0; i<centerTree.size(); ++i) {
            //if(isRedundant(centerTree[i])) {
                //centerTree[i]->removeNode();
                //delete centerTree[i];
                //centerTree.erase(centerTree.begin()+i);
                //clean = false;
                //break;
            //}
        //}

        //for(size_t i=0; i<centerTree.size(); ++i) {
            //if(centerTree[i]->getNumNeighbors() == 1) {
                //if(!areNeigbors(centerTree[i], centerTree[i]->getNeighbor(0))) {
                    //centerTree[i]->removeLinks();
                    //delete centerTree[i];
                    //centerTree.erase(centerTree.begin()+i);
                    //clean = false;
                    //break;
                //}
            //}
        //}
    //}

    //for(size_t i=0; i<centerTree.size(); ++i) {
        //for(size_t j=0; j<centerTree[i]->getNumNeighbors(); ++j) {
            //if(!areNeigbors(centerTree[i], centerTree[i]->getNeighbor(j))) {
                //std::cout << "unhandeld gap in CL!\n";
            //}
        //}
    //}

    Graph g;
    for(size_t i=0; i<centerTree.size(); ++i) {
        g.addNode(centerTree[i]);
    }
    return g;
}

} // namespace voreen
