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

struct BVHNode{
    int Children;
    float min, max;
    union{
        struct{
            vec3 LLF, URB;
        };
        struct{
            vec3 aabb[2];
        };
    };
};

struct Ray{
    vec3 origin;
    vec3 direction;
    vec3 inv_direction;
    int sign[3];
};

vec2 intersection_distances_no_if(Ray ray, vec3 aabb[2])
{
    float tmin;
    float tmax;

    float tymin, tymax, tzmin, tzmax;
    tmin = (aabb[ray.sign[0]].x - ray.origin.x) * ray.inv_direction.x;
    tmax = (aabb[1-ray.sign[0]].x - ray.origin.x) * ray.inv_direction.x;
    tymin = (aabb[ray.sign[1]].y - ray.origin.y) * ray.inv_direction.y;
    tymax = (aabb[1-ray.sign[1]].y - ray.origin.y) * ray.inv_direction.y;
    tzmin = (aabb[ray.sign[2]].z - ray.origin.z) * ray.inv_direction.z;
    tzmax = (aabb[1-ray.sign[2]].z - ray.origin.z) * ray.inv_direction.z;
    tmin = std::max(std::max(tmin, tymin), tzmin);
    tmax = std::min(std::min(tmax, tymax), tzmax);


    return vec2(tmin, tmax);
    // post condition:
    // if tmin > tmax (in the code above this is represented by a return value of INFINITY)
    //     no intersection
    // else
    //     front intersection point = ray.origin + ray.direction * tmin (normally only this point matters)
    //     back intersection point  = ray.origin + ray.direction * tmax
}


Ray makeRay(vec3 origin, vec3 direction) {
    vec3 inv_direction = vec3(1.0) / direction;
    Ray r;
    r.origin = origin;
    r.direction = direction;
    r.inv_direction = inv_direction;
    r.sign[0] = (inv_direction.x < 0) ? 1 : 0;
    r.sign[1] = (inv_direction.y < 0) ? 1 : 0;
    r.sign[2] = (inv_direction.z < 0) ? 1 : 0;
    return r;
}

float minVec3(vec3 v){
    return std::min(v.x, std::min(v.y, v.z));

}

vec3 computeNodeExitPoint(const vec3 nodeLLF, const vec3 nodeURB, const vec3 nodeEntry, const vec3 rayDir) {
    //tgtAssert(inRange(nodeEntry, nodeLLF, nodeURB), "node entry point outside node");

    vec3 exitPlanes = vec3(rayDir.x >= 0.f ? nodeURB.x : nodeLLF.x,
                                    rayDir.y >= 0.f ? nodeURB.y : nodeLLF.y,
                                    rayDir.z >= 0.f ? nodeURB.z : nodeLLF.z);
    vec3 tNodeExit;
    tNodeExit.x = rayDir.x != 0.f ? ((exitPlanes.x - nodeEntry.x) / rayDir.x) : 1e6f;
    tNodeExit.y = rayDir.y != 0.f ? ((exitPlanes.y - nodeEntry.y) / rayDir.y) : 1e6f;
    tNodeExit.z = rayDir.z != 0.f ? ((exitPlanes.z - nodeEntry.z) / rayDir.z) : 1e6f;
    //tgtAssert(tgt::hand(tgt::greaterThanEqual(tNodeExit, vec3::zero)), "at least one negative node exit parameter");

    float tNodeExitMin = minVec3(tNodeExit);
    //tgtAssert(inRange(tNodeExitMin, 0.f, 1.f), "minimum node exit parameter outside range [0.0;1.0]");

    vec3 nodeExit = nodeEntry + (tNodeExitMin - 1e-6f)*rayDir;
    //tgtAssert(inRange(nodeExit, nodeLLF, nodeURB), "node exit point outside node");
    return nodeExit;
}

bool pointInBox(const vec3 nodeLLF, const vec3 nodeURB, const vec3 node){
    return node.x >= nodeLLF.x && node.x <= nodeURB.x
        && node.y >= nodeLLF.y && node.y <= nodeURB.y
        && node.z >= nodeLLF.z && node.z <= nodeURB.z;
}

int depth = 0;

void traverseBVH(std::vector<BVHNode> &nodes, int pos, vec3 origin, vec3 dir){
    depth++;

    vec3 mid = (nodes[pos].LLF+nodes[pos].URB)*0.5f;
    for(int i = 0; i != depth;i++) printf(" ");
    int children = nodes[pos].Children;
    if (!children){
        // traverse
        printf(" (%f %f %f)[leaf]\n", mid.x, mid.y, mid.z);
        depth--;
        return;
    }else{
        printf(" (%f %f %f)\n", mid.x, mid.y, mid.z);
    }





    int first = (pointInBox(nodes[children].LLF, nodes[children].URB, origin))?0:1;
    assert(pointInBox(nodes[children+first].LLF, nodes[children+first].URB, origin));
    assert(!pointInBox(nodes[children+1-first].LLF, nodes[children+1-first].URB, origin));
    BVHNode node = nodes[children+first];
    vec3 out = computeNodeExitPoint(node.LLF, node.URB, origin, dir)+5e-6f*dir;


    traverseBVH(nodes, children+first, origin, dir);

    if (pointInBox(nodes[children+1-first].LLF, nodes[children+1-first].URB, out)){
        traverseBVH(nodes, children+1-first, out, dir);
    }
    depth--;
}

vec2 volumeRangeMinMax(const VolumeRAM* vol, svec3 LLF, svec3 URB){
    float min = 1;
    float max = -1;
    for(size_t x = LLF.x; x != URB.x; x++)
        for(size_t y = LLF.y; y != URB.y; y++)
            for(size_t z = LLF.z; z != URB.z; z++){
                float voxel = vol->getVoxelNormalized(x, y, z);
                min = std::min(voxel, min);
                max = std::max(voxel, max);
            }
    return vec2(min, max);
}


//==============================================================
//                                                            //
//               Top-Down BVH Creation                        //
//                                                            //
//==============================================================


void createBVHEx(const VolumeRAM* vol, svec3 LLF, svec3 URB, size_t leafsize, size_t pos, std::vector<BVHNode>& nodes){
    depth++;
    vec3 boxsize = URB-LLF;

    BVHNode node;
    node.URB = URB;
    node.LLF = LLF;

    // handle creation of leaf
    if (boxsize.x <= leafsize && boxsize.y <= leafsize && boxsize.z <= leafsize){
        vec2 minmax = volumeRangeMinMax(vol, LLF, URB);
        node.Children = 0;

        node.min = minmax.x;
        node.max = minmax.y;
        nodes[pos] = node;
        depth--;
        return;
    }

    // volume to big to create a leaf
    // we create a node

    vec3 midLLF = LLF;
    vec3 midURB = URB;
    int maxdim;
    if (boxsize.x > boxsize.y && boxsize.x > boxsize.z){
        maxdim = 0;
        midLLF.x = (float)((LLF.x+URB.x)/2);
        midURB.x = (float)((LLF.x+URB.x)/2);
    }else if(boxsize.y > boxsize.z){
        maxdim = 1;
        midLLF.y = (float)((LLF.y+URB.y)/2);
        midURB.y = (float)((LLF.y+URB.y)/2);
    }else{
        maxdim = 2;
        midLLF.z = (float)((LLF.z+URB.z)/2);
        midURB.z = (float)((LLF.z+URB.z)/2);
    }



    size_t newpos = nodes.size();

    node.Children = (int)newpos;

    nodes.resize(nodes.size()+2);
    createBVHEx(vol, LLF, midURB, leafsize, newpos, nodes);
    createBVHEx(vol, midLLF, URB, leafsize, newpos+1, nodes);

    node.min = std::min(nodes[newpos].min, nodes[newpos+1].min);
    node.max = std::max(nodes[newpos].max, nodes[newpos+1].max);

    nodes[pos] = node;

    for(int i = 0; i != depth; i++){
        printf("==");
    }
    printf("%i (%f %f)\n", maxdim, node.min, node.max);

    depth--;
}

void createBVHEx(const VolumeRAM* vol, int leafsize, std::vector<BVHNode>& nodes){
    svec3 LLF = svec3(0, 0, 0);
    svec3 URB = vol->getDimensions();
    nodes.resize(1);
    createBVHEx(vol, LLF, URB, leafsize, 0, nodes);
}


void enumerateLeafs(const VolumeRAM* vol, size_t leafsize, std::vector<BVHNode> &nodes){
    svec3 size = vol->getDimensions();
    svec3 node_count = size/svec3(leafsize);
    for(int x = 0; x != node_count.x; x++)
        for(int y = 0; y != node_count.y; y++)
            for(int z = 0; z != node_count.z; z++){
                svec3 LLF = leafsize*svec3(x, y, z);
                svec3 URB = tgt::max((leafsize+1)*svec3(x, y, z), size);
                vec2 minmax = volumeRangeMinMax(vol, LLF, URB);

                BVHNode node;
                node.URB = URB;
                node.LLF = LLF;

                node.Children = 0;
                node.min = minmax.x;
                node.max = minmax.y;

                nodes.push_back(node);
            }
}


//==============================================================
//                                                            //
//              Bottom-up BVH Creation                        //
//                                                            //
//==============================================================
svec3 min3(svec3 a, svec3 b, svec3 c){
    return tgt::min(a, tgt::min(b, c));
}

svec3 max3(svec3 a, svec3 b, svec3 c){
    return tgt::max(a, tgt::max(b, c));
}

int getMaxDim(std::vector<BVHNode> &leafs, size_t begin, size_t end){
    assert(end > begin);
    svec3  min = svec3(-100000000);
    svec3  max = svec3( 100000000);

    for(size_t i = begin; i != end; i++){
        svec3 llf = svec3(leafs[i].LLF);
        svec3 urb = svec3(leafs[i].URB);
        min = min3(llf, urb , min);
        max = max3(llf, urb , max);
    }

    svec3 size = max-min;
    if (size.x > size.y && size.x > size.z)
        return 1;
    else if (size.y > size.z)
        return 2;
    else
        return 3;
}



bool x_less(BVHNode a, BVHNode b){
    return a.LLF.x < b.LLF.x;
}

bool y_less(BVHNode a, BVHNode b){
    return a.LLF.y < b.LLF.y;
}

bool z_less(BVHNode a, BVHNode b){
    return a.LLF.z < b.LLF.z;
}

void createBVHBottomUpEx(std::vector<BVHNode> &leafs, size_t begin, size_t end, std::vector<BVHNode> &tree, size_t pos){
    assert(end > begin);
    bool create_leaf = (end-begin) == 1;

    if (create_leaf){
        tree[pos] = leafs[0];
        return;
    }

    // we need to create a inner node

    // sort in biggest dimension
    int maxdim = getMaxDim(leafs, begin, end);
    assert(maxdim >= 1 && maxdim <= 3);
    if (maxdim == 1){
        std::sort(leafs.begin()+begin, leafs.begin()+end, x_less);
    }else if (maxdim == 2){
        std::sort(leafs.begin()+begin, leafs.begin()+end, y_less);
    }else{
        std::sort(leafs.begin()+begin, leafs.begin()+end, z_less);
    }

    size_t mid = (begin+end)/2;
    size_t new_pos = tree.size();
    tree.resize(new_pos+2);

    createBVHBottomUpEx(leafs, begin, mid, tree, new_pos);
    createBVHBottomUpEx(leafs, mid, end, tree, new_pos+1);

    BVHNode node;
    node.LLF = min(leafs[new_pos].LLF, leafs[new_pos+1].LLF);
    node.URB = max(leafs[new_pos].URB, leafs[new_pos+1].URB);
    node.Children = (int)pos;
    node.min = std::min(leafs[new_pos].min, leafs[new_pos+1].min);
    node.max = std::max(leafs[new_pos].max, leafs[new_pos+1].max);

    tree[pos] = node;

}

void createBVHBottomUp(const VolumeRAM* vol, int leafsize, std::vector<BVHNode>& nodes){
    std::vector<BVHNode> leafs;
    enumerateLeafs(vol, leafsize, leafs);
    nodes.resize(1);
    createBVHBottomUpEx(leafs, 0, leafs.size(), nodes, 0);
}
