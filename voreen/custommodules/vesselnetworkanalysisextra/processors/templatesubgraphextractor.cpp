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

#include "templatesubgraphextractor.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "modules/vesselnetworkanalysis/algorithm/vesselgraphrefinement.h"

#include "tgt/tgt_math.h"
#include "tgt/quaternion.h"
#include "tgt/memory.h"
#include "tgt/immediatemode/immediatemode.h"

#include "yaml-cpp/yaml.h"

#include <numeric>
#include <queue>
#include <unordered_set>
#include <unordered_map>

namespace voreen {

struct TemplateNode;

struct Direction {
    tgt::vec3 direction_; //Always: length == 1

    Direction(const tgt::vec3& direction)
        : direction_(tgt::normalize(direction))
    {
        if(direction == tgt::vec3::zero) {
            throw tgt::Exception("Trying to initialize Direction vector zero length");
        }
    }

    // from-to
    Direction(const tgt::vec3& start, const tgt::vec3& end)
        : direction_(tgt::normalize(end-start))
    {
    }

    float angleTo(const Direction& other) const {
        return std::acos(tgt::dot(other.direction_, direction_));
    }

    Direction transformed_by(const tgt::mat3& transformation) const {
        return Direction(transformation*direction_);
    }
};


struct DirectionRange {
    Direction direction_;
    float max_deviation_;

    DirectionRange(const tgt::vec3& direction, float max_deviation)
        : DirectionRange(Direction(direction), max_deviation)
    {
    }

    DirectionRange(const Direction& direction, float max_deviation)
        : direction_(direction)
        , max_deviation_(max_deviation)
    {
    }

    DirectionRange(const DirectionRange& other)
        : direction_(other.direction_)
        , max_deviation_(other.max_deviation_)
    {
    }

    bool contains(const Direction& dir) const {
        return direction_.angleTo(dir) < max_deviation_;
    }

    DirectionRange transformed_by(const tgt::mat3& transformation) const {
        return DirectionRange(direction_.transformed_by(transformation), max_deviation_);
    }
};

struct TemplateBranch {
    DirectionRange branch_direction_;
    float minLength_;
    float maxLength_;
    float minStraightness;
    mutable float maxBenefitCache_;
    mutable float missingCostCache_;

    // May be null
    std::unique_ptr<TemplateNode> edge_;

    TemplateBranch(float minLength, float maxLength, float minStraightness, DirectionRange branch_direction, TemplateNode* edge = nullptr)
        : minLength_(minLength)
        , maxLength_(maxLength)
        , minStraightness(minStraightness)
        , branch_direction_(branch_direction)
        , edge_(edge)
        , maxBenefitCache_(std::numeric_limits<float>::quiet_NaN())
        , missingCostCache_(std::numeric_limits<float>::quiet_NaN())
    {
        tgtAssert(maxLength_ >= minLength_, "max < min");
    }

    float missing_cost() const;
    float max_benefit() const;
};


struct TemplateNode {
    TemplateBranch main_branch_;
    TemplateBranch off_branch_;

    TemplateNode(TemplateBranch&& main_branch, TemplateBranch&& off_branch)
        : main_branch_(std::move(main_branch))
        , off_branch_(std::move(off_branch))
    {
    }
};

float TemplateBranch::missing_cost() const {
    //TODO: maybe we want to incorporate the max benefit?
    if(std::isnan(missingCostCache_)) {
        missingCostCache_ = minLength_;
        if(edge_) {
            missingCostCache_ += edge_->main_branch_.missing_cost();
            missingCostCache_ += edge_->off_branch_.missing_cost();
        }
    }
    return missingCostCache_;
}

float TemplateBranch::max_benefit() const {
    if(std::isnan(maxBenefitCache_)) {
        maxBenefitCache_ = maxLength_;
        if(edge_) {
            maxBenefitCache_ += edge_->main_branch_.max_benefit();
            maxBenefitCache_ += edge_->off_branch_.max_benefit();
        }
    }
    return maxBenefitCache_;
}

static float CONST_DEGREE_TO_RADIAN = 2*tgt::PI/360;

static tgt::vec3 tryReadVec(const YAML::Node& node) {
    tgt::vec3 out;
    try {
        out.x = node["x"].as<float>();
        out.y = node["y"].as<float>();
        out.z = node["z"].as<float>();
    } catch(...) {
        float roll = node["roll"].as<float>() * CONST_DEGREE_TO_RADIAN;
        float pitch = node["pitch"].as<float>() * CONST_DEGREE_TO_RADIAN;

        tgt::mat3 pitch_mat = tgt::mat3::createRotationX(-pitch);
        tgt::mat3 roll_mat = tgt::mat3::createRotationZ(roll);

        out = roll_mat * pitch_mat * tgt::vec3(0,0,1);
    }
    return out;
}

static DirectionRange tryReadDirectionRange(const YAML::Node& node) {
    return DirectionRange(tryReadVec(node["dir"]), node["angle"].as<float>() * CONST_DEGREE_TO_RADIAN);
}

static TemplateBranch tryReadBranch(const YAML::Node& node) {
    tgtAssert(node.size() > 0, "Empty node");

    const auto& element = node[node.size()-1];
    DirectionRange dir = tryReadDirectionRange(element);
    float min_length = 0;
    float max_length = std::numeric_limits<float>::infinity();
    float min_straightness = 0;
    try {
        min_length = element["min_length"].as<float>();
    } catch(...) { }
    try {
        max_length = element["max_length"].as<float>();
    } catch(...) { }
    try {
        min_straightness = element["min_straightness"].as<float>();
    } catch(...) { }
    TemplateBranch currentBranch(min_length, max_length, min_straightness, dir);

    for(size_t i=node.size()-1; i != 0; --i) {
        const auto& element = node[i-1];

        DirectionRange dir = tryReadDirectionRange(element);
        float min_length = 0;
        float max_length = std::numeric_limits<float>::infinity();
        float min_straightness = 0;
        try {
            min_length = element["min_length"].as<float>();
        } catch(...) { }
        try {
            max_length = element["max_length"].as<float>();
        } catch(...) { }
        try {
            min_straightness = element["min_straightness"].as<float>();
        } catch(...) { }

        if(element["branch"]) {
            TemplateBranch main_branch = std::move(currentBranch);
            TemplateBranch off_branch = tryReadBranch(element["branch"]);
            currentBranch = TemplateBranch(min_length, max_length, min_straightness, dir, new TemplateNode(std::move(main_branch), std::move(off_branch)));
        } else {
            tgtAssert(false, "Invalid branch element");
        }
    }
    return currentBranch;
}

static TemplateBranch parse_template(const std::string& template_file_path) {
    YAML::Node tpl = YAML::LoadFile(template_file_path);
    return tryReadBranch(tpl);
}

struct BranchSearchCone {
    DirectionRange dir;
    float min_length;
    float max_length;
    tgt::vec3 node_pos;

    BranchSearchCone(const TemplateBranch& template_branch, const VesselGraphNode& node, const tgt::mat3& orientation_transformation)
        : dir(template_branch.branch_direction_.transformed_by(orientation_transformation))
        , min_length(template_branch.minLength_)
        , max_length(template_branch.maxLength_)
        , node_pos(node.pos_)
    {
    }
};

const std::string TemplateSubgraphExtractor::loggerCat_("voreen.vesselnetworkanalysisextra.templatesubgraphextractor");

TemplateSubgraphExtractor::TemplateSubgraphExtractor()
    : GeometryRendererBase()
    , inport_(Port::INPORT, "graph.input", "Graph Input", false, Processor::INVALID_RESULT)
    , startingPoint_(Port::INPORT, "pointlist.seedpoints", "Starting Point", false, Processor::INVALID_RESULT)
    , outport_(Port::OUTPORT, "graph.output", "Normalized Graph Output", false, Processor::VALID)
    , enabled_("enabled", "Enabled", true)
    , templateFile_("templateFile", "Template File", "Path", "", "YAML (*.yml)", FileDialogProperty::OPEN_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT, VoreenFileWatchListener::OPTIONAL_ON)
    , reloadTemplate_("reloadTemplate", "Reload Template")
    , keepBounds_("keepBounds", "Keep Bounds of Input Graph", true)
    , solvingStrategy_("solvingStrategy", "Solving Strategy")
    , coneRenderMode_("ConeRenderMode", "Cone Render Mode")
    , missingBranchConeColor_("missingBranchConeColor_", "Missing Branch Cone Color", tgt::Color(1.0f, 0.2f, 0.1f, 1.0f))
    , locatedBranchConeColor_("locatedBranchConeColor_", "Located Branch Cone Color", tgt::Color(0.5f, 0.8f, 0.4f, 1.0f))
    , coneLength_("coneLength", "Cone Length (mm)", 1, 0, 1000)
    , lightPosition_("lightPosition", "Light Source Position", tgt::vec4(2.3f, 1.5f, 1.5f, 1.f), tgt::vec4(-10000), tgt::vec4(10000))
    , lightAmbient_("lightAmbient", "Ambient Light", tgt::Color(0.4f, 0.4f, 0.4f, 1.f))
    , lightDiffuse_("lightDiffuse", "Diffuse Light", tgt::Color(0.8f, 0.8f, 0.8f, 1.f))
    , lightSpecular_("lightSpecular", "Specular Light", tgt::Color(0.6f, 0.6f, 0.6f, 1.f))
    , materialShininess_("materialShininess", "Shininess", 60.f, 0.1f, 128.f)
{

    addPort(inport_);
    addPort(startingPoint_);
    addPort(outport_);

    addProperty(enabled_);
    addProperty(templateFile_);
    addProperty(reloadTemplate_);
    addProperty(keepBounds_);
    addProperty(solvingStrategy_);
        solvingStrategy_.addOption("greedy", "Greedy", STRATEGY_GREEDY);
        solvingStrategy_.addOption("global", "Global", STRATEGY_GLOBAL);
        solvingStrategy_.selectByValue(STRATEGY_GLOBAL);
    addProperty(coneRenderMode_);
        coneRenderMode_.addOption("none", "None", CONES_NONE);
        coneRenderMode_.addOption("missing", "Missing", CONES_MISSING);
        coneRenderMode_.addOption("all", "All", CONES_ALL);
        coneRenderMode_.selectByValue(CONES_MISSING);
        ON_CHANGE(coneRenderMode_, TemplateSubgraphExtractor, adjustColorPropertyVisibility);
    addProperty(missingBranchConeColor_);
    addProperty(locatedBranchConeColor_);
    addProperty(coneLength_);

        addProperty(lightPosition_);
            lightPosition_.setGroupID("lighting");
        addProperty(lightAmbient_);
            lightAmbient_.setGroupID("lighting");
        addProperty(lightDiffuse_);
            lightDiffuse_.setGroupID("lighting");
        addProperty(lightSpecular_);
            lightSpecular_.setGroupID("lighting");
        addProperty(materialShininess_);
            materialShininess_.setGroupID("lighting");
    setPropertyGroupGuiName("lighting", "Lighting Parameters");

    adjustColorPropertyVisibility();
}

TemplateSubgraphExtractor::~TemplateSubgraphExtractor() {
}

static bool tryExtractPoints(const Geometry* geometry, std::vector<tgt::vec3>& output) {
    auto transformation = geometry->getTransformationMatrix();
    auto transform = [&transformation] (const tgt::vec3& v) {
        return (transformation*tgt::vec4(v, 1)).xyz();
    };
    if(const PointSegmentListGeometryVec3* seedList = dynamic_cast<const PointSegmentListGeometryVec3*>(geometry)) {
        for(int i = 0; i < seedList->getNumSegments(); i++) {
            auto segment = seedList->getSegment(i);
            for(const tgt::vec3& v: segment) {
                output.push_back(transform(v));
            }
        }
    } else if(const PointListGeometryVec3* seeds = dynamic_cast<const PointListGeometryVec3*>(geometry)) {
        for(const tgt::vec3& v: *seeds) {
            output.push_back(transform(v));
        }
    } else {
        return false;
    }
    return true;
}

static bool tryExtractCenterPoint(const GeometryPort& g, tgt::vec3& output) {
    if(!g.hasData()) {
        return false;
    }
    std::vector<tgt::vec3> points;
    if(!tryExtractPoints(g.getData(), points)) {
        return false;
    }
    if(points.empty()) {
        return false;
    }
    output = std::accumulate(points.begin(), points.end(), tgt::vec3::zero) / static_cast<float>(points.size());
    return true;
}

struct NodePath {
    std::vector<const VesselGraphEdge*> edges_;
    const VesselGraphNode* head_;
    const VesselGraphNode* tail_;
    float length_;

    const VesselGraphNode* head() const {
        return head_;
    }

    const VesselGraphNode* tail() const {
        return tail_;
    }

    const NodePath extended_by(const VesselGraphEdge* new_edge) {
        tgtAssert(&new_edge->getNode1() == head_ || &new_edge->getNode2() == head_, "Invalid path extension: not connected");
        tgtAssert(!new_edge->isLoop(), "Invalid path extension: loop");

        std::vector<const VesselGraphEdge*> longer_path(edges_.begin(), edges_.end());
        longer_path.push_back(new_edge);
        return NodePath(std::move(longer_path), &new_edge->getOtherNode(*head_), tail_, length_ + new_edge->getLength());
    }

    NodePath(const VesselGraphEdge* initial, const VesselGraphNode* head)
        : NodePath(std::vector<const VesselGraphEdge*>(), head, &initial->getOtherNode(*head), initial->getLength())
    {
        edges_.push_back(initial);
    }

    float getLength() const {
        return length_;
    }

    float getDistance() const {
        return tgt::distance(head_->pos_, tail_->pos_);
    }

    float getStraightness() const {
        float length = getLength();
        tgtAssert(length > 0, "invalid (nonpositive) length");
        return getDistance()/length;
    }

    Direction direction() const {
        tgt::vec3 rw_from_to = head_->pos_ - tail_->pos_;
        return Direction(rw_from_to);
    }

    std::vector<VesselSkeletonVoxel> joinEdges() const {
        std::vector<VesselSkeletonVoxel> output;
        const VesselGraphNode* current_tail = tail_;
        for(const VesselGraphEdge* edge: edges_) {
            const VesselGraphNode* current_head;
            const auto& new_voxels = edge->getVoxels();
            if(&edge->getNode1() == current_tail) {
                current_head = &edge->getNode2();
                output.insert(output.end(), new_voxels.begin(), new_voxels.end());
            } else {
                tgtAssert(&edge->getNode2() == current_tail, "invalid node path");
                current_head = &edge->getNode1();
                for(auto it = edge->getVoxels().rbegin(); it != edge->getVoxels().rend(); ++it) {
                    output.push_back(*it);
                }
            }

            current_tail = current_head;
        }
        return output;
    }

private:
    NodePath(std::vector<const VesselGraphEdge*> edges, const VesselGraphNode* head, const VesselGraphNode* tail, float length)
        : edges_(edges)
        , head_(head)
        , tail_(tail)
        , length_(length)
    {
    }

};

struct PathLongerThan {
    constexpr bool operator()(const NodePath &lhs, const NodePath &rhs) const
    {
        return lhs.length_ > rhs.length_;
    }
};

static Direction branchDirection(const VesselGraphNode& node, const VesselGraphEdge& branch) {
    tgtAssert(!branch.isLoop(), "direction for loop");
    tgtAssert(!branch.getVoxels().empty(), "direction for loop");
    tgt::vec3 rw_from_to = branch.getOtherNode(node).pos_ - node.pos_;
    return Direction(rw_from_to);
}

static Direction leavingDirection(const VesselGraphNode& node, const VesselGraphEdge& branch) {
    tgtAssert(!branch.isLoop(), "direction for loop");
    tgtAssert(!branch.getVoxels().empty(), "direction for loop");

    const int max_voxels = 5;

    const auto& voxels = branch.getVoxels();
    tgt::vec3 dir_sum = tgt::vec3::zero;
    int i = 0;

    //TODO: maybe do some exponentional smoothing along the path?
    if(&branch.getNode1() == &node) {
        //=> we are interested in the _first_ voxels of the path
        for(auto it = voxels.begin(); it != voxels.end(); ++it, ++i) {
            dir_sum += (it->pos_ - node.pos_);//(1 << i);
            if(i == max_voxels) {
                break;
            }
        }
    } else {
        //=> we are interested in the _last_ voxels of the path
        for(auto it = voxels.rbegin(); it != voxels.rend(); ++it, ++i) {
            dir_sum += (it->pos_ - node.pos_);//(1 << i);
            if(i == max_voxels) {
                break;
            }
        }
    }
    return Direction(dir_sum);
}

struct PathTreeNode {
    NodePath path_;
    PathTreeNode* parent_;
    std::unordered_set<PathTreeNode*> children_;

    const static std::string loggerCat_;

    PathTreeNode(NodePath path, PathTreeNode* parent)
        : path_(path)
        , parent_(parent)
        , children_()
    {
        if(parent_) {
            tgtAssert(parent_->children_.count(this) == 0, "Node already added to parent!");
            parent_->children_.insert(this);
        }
    }

    ~PathTreeNode()
    {
        for(PathTreeNode* child : children_) {
            child->parent_ = nullptr;
            delete child;
        }
        if(parent_) {
            parent_->children_.erase(this);
        }
    }

    //Should probably only be called on finished trees
    //void collect_tree_components(std::unordered_set<const VesselGraphNode*>& nodes, std::vector<const NodePath*>& paths, std::vector<BranchSearchCone>& missingBranchCones, std::vector<BranchSearchCone>& allBranchCones) {
    void collect_tree_components(std::unordered_set<const VesselGraphNode*>& nodes, std::vector<const NodePath*>& paths) {
        //TODO: we may be able to collect cones here as well

        paths.push_back(&path_);

        // Will result in superfluous inserts, but thats not a problem in terms of performance or correctness
        nodes.insert(path_.head());
        nodes.insert(path_.tail());

        tgtAssert(children_.size() <= 2, "Invalid number of children (more than 2)");
        for(PathTreeNode* child : children_) {
            child->collect_tree_components(nodes, paths);
        }
    }

    tgt::mat3 createOrientationMatrixLocalToGlobal() {
        /*
        if(parent_ && parent_->parent_) {
            const PathTreeNode* grandfather = parent_->parent_;

            const NodePath& grandfatherPath = grandfather->path_;

            const tgt::vec3 ahead = path_.head_->pos_ - path_.tail_->pos_;
            tgtAssert(ahead != tgt::vec3::zero, "ahead vector is zero");
            const tgt::vec3 top = grandfatherPath.head_->pos_ - grandfatherPath.tail_->pos_;
            tgtAssert(top != tgt::vec3::zero, "top vector is zero");

            return tgt::mat4::rigidBodyTranformation(ahead, top, tgt::vec3::zero).getRotationalPartMat3();

        } else {
            //We may want to lift this restriction later, but will have to be more careful to choose the grandfather then.
            if(parent_ && !parent_->parent_) {
                LWARNING("No grandfather edge, falling back to normal estimation from curvature");
            }

            const tgt::vec3 ahead = path_.head_->pos_ - path_.tail_->pos_;
            const tgt::vec3 center = (path_.head_->pos_ + path_.tail_->pos_)*0.5f;

            tgt::vec3 top = tgt::vec3::zero;
            size_t num_voxels = 0;
            for(const VesselGraphEdge* edge: path_.edges_) {
                for(const auto& voxel: edge->getVoxels()) {
                    top += voxel.pos_ - center;
                    ++num_voxels;
                }
            }
            tgtAssert(num_voxels > 0, "No voxels in path");

            //top (and ahead) will be normalized in rigidBodyTransformation
            return tgt::mat4::rigidBodyTranformation(ahead, top, tgt::vec3::zero).getRotationalPartMat3();
        }
        */

        // Relative orientation seems to be to unstable. we try globally consistent orientations instead.
        return tgt::mat3::identity;
    }

private:
    PathTreeNode(const PathTreeNode&);
    void operator=(const PathTreeNode&);
};

const std::string PathTreeNode::loggerCat_("voreen.templatesubgraphextractor.pathtreenode");

struct UnprocessedBranch {
    const VesselGraphNode* node;
    const TemplateBranch& t_branch;
    PathTreeNode* parent_path;

    UnprocessedBranch(const VesselGraphNode* node, const TemplateBranch& branch)
        : node(node)
        , t_branch(branch)
        , parent_path(nullptr)
    {
    }
    UnprocessedBranch(PathTreeNode* parentPath, const TemplateBranch& branch)
        : node(parentPath->path_.head())
        , t_branch(branch)
        , parent_path(parentPath)
    {
    }

    std::pair<NodePath, BranchSearchCone> find_path_from_template();

    std::pair<std::vector<NodePath>, BranchSearchCone> find_all_fitting_paths(const VesselGraphEdge* parent_edge_to_ignore = nullptr, const VesselGraphEdge* sibling_edge_to_ignore = nullptr);
};

std::pair<NodePath, BranchSearchCone> UnprocessedBranch::find_path_from_template() {
    const VesselGraphNode* current = node;
    const TemplateBranch& template_branch = t_branch;

    //TODO: we could do try to get initial the orientation from the volume or metadata
    const tgt::mat3 orientation_local_to_global = parent_path ? parent_path->createOrientationMatrixLocalToGlobal() : tgt::mat3::identity;

    std::unordered_set<const VesselGraphNode*> visited_nodes;
    visited_nodes.insert(current);

    //TODO: somehow consider branches (and directions) at END of a path!
    std::priority_queue<NodePath, std::vector<NodePath>, PathLongerThan> paths;
    auto global_branch_dir_range = template_branch.branch_direction_.transformed_by(orientation_local_to_global);
    for(const auto& branch: current->getEdges()) {
        // TODO: this restriction might not make sense, we may have to split the directions into "leaving" and "final" or something like that
        //auto dir = branchDirection(*current, branch);
        //if(global_branch_dir_range.contains(dir) && !branch.get().isLoop()) {
        paths.emplace(&branch.get(), &branch.get().getOtherNode(*current));
        //}
    }

    BranchSearchCone cone(t_branch, *current, orientation_local_to_global);

    while(true) {
        if(paths.empty()) {
            throw cone;
        }
        NodePath p = paths.top();
        paths.pop();
        if(template_branch.minLength_ <= p.getDistance() && p.getDistance() <= template_branch.maxLength_ && global_branch_dir_range.contains(p.direction())) {
            return std::make_pair(p, cone);
        }
        for(const auto& edge : p.head()->getEdges()) {
            if(edge.get().isLoop()) {
                continue;
            }
            const auto& new_node = edge.get().getOtherNode(*p.head());
            if(visited_nodes.count(&new_node) == 0) {
                visited_nodes.insert(&new_node);
                paths.push(p.extended_by(&edge.get()));
            }
        }
    }
}

std::pair<std::vector<NodePath>, BranchSearchCone> UnprocessedBranch::find_all_fitting_paths(const VesselGraphEdge* parent_edge_to_ignore, const VesselGraphEdge* sibling_edge_to_ignore) {
    const VesselGraphNode* current = node;
    const TemplateBranch& template_branch = t_branch;

    //TODO: we could do try to get initial the orientation from the volume or metadata
    const tgt::mat3 orientation_local_to_global = parent_path ? parent_path->createOrientationMatrixLocalToGlobal() : tgt::mat3::identity;

    //TODO: maybe just ignoring all nodes that are already in a path is not a good idea. (Imaging an erroneous segmentation creates a shortcut to the desired node!)
    // On the other hand, we definitely do not want to walk back into nodes of previous node paths!
    std::unordered_set<const VesselGraphNode*> visited_nodes;
    visited_nodes.insert(current);

    std::priority_queue<NodePath, std::vector<NodePath>, PathLongerThan> paths;
    auto global_branch_dir_range = template_branch.branch_direction_.transformed_by(orientation_local_to_global);
    for(const auto& branch: current->getEdges()) {
        if(branch.get().isLoop() || &branch.get() == parent_edge_to_ignore || &branch.get() == sibling_edge_to_ignore) {
            continue;
        }
        auto dir = leavingDirection(*current, branch);
        if(global_branch_dir_range.contains(dir)) {
            paths.emplace(&branch.get(), &branch.get().getOtherNode(*current));
        }
    }

    std::vector<NodePath> results;

    BranchSearchCone cone(t_branch, *current, orientation_local_to_global);

    while(true) {
        if(paths.empty()) {
            return std::make_pair(results, cone);
        }
        NodePath p = paths.top();
        paths.pop();
        float length = p.getLength();
        if(length <= template_branch.maxLength_) {
            if(length >= template_branch.minLength_ && p.getStraightness() >= template_branch.minStraightness) {
                results.push_back(p);
            }
            for(const auto& edge : p.head()->getEdges()) {
                if(edge.get().isLoop()) {
                    continue;
                }
                const auto& new_node = edge.get().getOtherNode(*p.head());
                if(visited_nodes.count(&new_node) == 0) {
                    visited_nodes.insert(&new_node);
                    paths.push(p.extended_by(&edge.get()));
                }
            }
        }
    }
}

static const VesselGraphNode* find_starting_node(const VesselGraph& graph, const tgt::vec3& starting_point) {
    auto nodes = graph.getNodes();
    float min_dist_sq = std::numeric_limits<float>::infinity();
    const VesselGraphNode* starting_node = nullptr;
    for(auto& node : nodes) {
        float current_dist_sq = tgt::distanceSq(starting_point, node.pos_);
        if(current_dist_sq < min_dist_sq) {
            min_dist_sq = current_dist_sq;
            starting_node = &node;
        }
    }
    return starting_node;
}

static void fillGraph(PathTreeNode* root, VesselGraphBuilder& output) {

    std::unordered_set<const VesselGraphNode*> nodes_to_keep;
    std::vector<const NodePath*> paths;
    root->collect_tree_components(nodes_to_keep, paths);

    std::unordered_map<const VesselGraphNode*, VGNodeID> nodes_to_new_index;
    for(auto node : nodes_to_keep) {
        VGNodeID index = output.insertNode(*node);
        nodes_to_new_index.insert(std::make_pair(node, index));
    }

    for(const auto& path : paths) {
        auto voxels = path->joinEdges();
        output.insertEdge(nodes_to_new_index[path->tail()], nodes_to_new_index[path->head()], std::move(voxels));
    }
}

struct SubgraphExtractionResult {
    float cost;
    std::unique_ptr<PathTreeNode> node;
    std::vector<BranchSearchCone> locatedSearchCones;
    std::vector<BranchSearchCone> missingSearchCones;

    SubgraphExtractionResult()
        : cost(std::numeric_limits<float>::infinity())
        , node(nullptr)
        , locatedSearchCones()
        , missingSearchCones()
    {
    }

    SubgraphExtractionResult(const TemplateBranch& missing_template, BranchSearchCone cone)
        : cost(missing_template.missing_cost())
        , node(nullptr)
        , locatedSearchCones()
        , missingSearchCones()
    {
        missingSearchCones.push_back(cone);
    }

    void set_search_cones_from(SubgraphExtractionResult& main_branch_result, SubgraphExtractionResult& off_branch_result) {
        locatedSearchCones.clear();
        missingSearchCones.clear();

        locatedSearchCones.insert(locatedSearchCones.end(), main_branch_result.locatedSearchCones.begin(), main_branch_result.locatedSearchCones.end());
        missingSearchCones.insert(missingSearchCones.end(), main_branch_result.missingSearchCones.begin(), main_branch_result.missingSearchCones.end());

        locatedSearchCones.insert(locatedSearchCones.end(), off_branch_result.locatedSearchCones.begin(), off_branch_result.locatedSearchCones.end());
        missingSearchCones.insert(missingSearchCones.end(), off_branch_result.missingSearchCones.begin(), off_branch_result.missingSearchCones.end());
    }
};

static SubgraphExtractionResult extractBranch(UnprocessedBranch b, const VesselGraphEdge* parent_edge_to_ignore = nullptr, const VesselGraphEdge* sibling_edge_to_ignore = nullptr) {

    //TODO: improve runtime by cutting branches that have larger cost than the current optimum?
    // This would require passing the best cost so far into recursive function calls

    auto plausible_paths_and_cone = b.find_all_fitting_paths(parent_edge_to_ignore, sibling_edge_to_ignore);
    auto& plausible_paths = plausible_paths_and_cone.first;
    auto& cone = plausible_paths_and_cone.second;

    if(plausible_paths.empty()) {
        return SubgraphExtractionResult(b.t_branch, cone);
    }

    SubgraphExtractionResult best_result;
    for(auto path: plausible_paths) {
        float straightness = path.getStraightness();
        float cost = 1.0f - straightness;

        if(cost >= best_result.cost) {
            continue;
        }

        std::unique_ptr<PathTreeNode> candidate = tgt::make_unique<PathTreeNode>(path, b.parent_path);

        tgtAssert(!path.edges_.empty(), "No edges in parent branch");
        const VesselGraphEdge* parent_branch_end_edge = path.edges_.back();

        if(b.t_branch.edge_) {
            SubgraphExtractionResult main_branch_result = extractBranch(UnprocessedBranch(candidate.get(), b.t_branch.edge_->main_branch_), parent_branch_end_edge);
            cost += main_branch_result.cost;

            const VesselGraphEdge* main_branch_begin_edge = nullptr;

            // TODO this may be problematic because now the main branch has priority in choosing a branch.
            // Maybe (hopefully) this is not a problem in practice.
            if(main_branch_result.node) {
                const auto& edges_of_main_branch = main_branch_result.node->path_.edges_;
                tgtAssert(!edges_of_main_branch.empty(), "No edges in main branch");
                main_branch_begin_edge = edges_of_main_branch.front();
            }

            main_branch_result.node.release(); //Avoid calling destructor on branch, will be called later if candidate is destroyed
            if(cost >= best_result.cost) {
                continue;
            }


            SubgraphExtractionResult off_branch_result = extractBranch(UnprocessedBranch(candidate.get(), b.t_branch.edge_->off_branch_), parent_branch_end_edge, main_branch_begin_edge);
            cost += off_branch_result.cost;
            off_branch_result.node.release(); //Avoid calling destructor on branch, will be called later if candidate is destroyed
            if(cost >= best_result.cost) {
                continue;
            }

            best_result.set_search_cones_from(main_branch_result, off_branch_result);
        }

        best_result.cost = cost;
        best_result.node.reset(candidate.release());
    }
    best_result.locatedSearchCones.push_back(cone);
    return best_result;
}

std::unique_ptr<VesselGraph> TemplateSubgraphExtractor::extractSubgraphGlobal(const VesselGraph& input, const tgt::vec3& starting_point, const TemplateBranch& tpl, bool keepBounds) {
    VesselGraphBuilder builder(keepBounds ? VesselGraphBuilder(input.getBounds()) : VesselGraphBuilder());
    const VesselGraphNode* starting_node = find_starting_node(input, starting_point);
    if(!starting_node) {
        // No nodes in graph
        return std::move(builder).finalize();
    }

    SubgraphExtractionResult tree = extractBranch(UnprocessedBranch(starting_node, tpl));

    missingBranches_.clear();
    locatedBranches_.clear();

    missingBranches_ = std::move(tree.missingSearchCones);
    locatedBranches_ = std::move(tree.locatedSearchCones);

    if(!missingBranches_.empty()) {
        LWARNING("" << missingBranches_.size() << " branche(s) from template could not be matched");
    }

    if(tree.node) {
        // No node => did not match any branch
        fillGraph(tree.node.get(), builder);
    }
    return std::move(builder).finalize();
}

std::unique_ptr<VesselGraph> TemplateSubgraphExtractor::extractSubgraph(const VesselGraph& input, const tgt::vec3& starting_point, const TemplateBranch& tpl, bool keepBounds) {
    VesselGraphBuilder builder(keepBounds ? VesselGraphBuilder(input.getBounds()) : VesselGraphBuilder());
    const VesselGraphNode* starting_node = find_starting_node(input, starting_point);
    if(!starting_node) {
        // No nodes in graph
        return std::move(builder).finalize();
    }

    std::queue<UnprocessedBranch> branches_to_process;
    branches_to_process.emplace(starting_node, tpl);

    std::unique_ptr<PathTreeNode> root(nullptr);

    missingBranches_.clear();
    locatedBranches_.clear();

    while(!branches_to_process.empty()) {
        auto to_process = branches_to_process.front();
        branches_to_process.pop();

        auto paths_and_cone = to_process.find_all_fitting_paths();
        auto& paths = paths_and_cone.first;
        auto& cone = paths_and_cone.second;
        if(paths.empty()) {
            missingBranches_.push_back(cone);
            continue;
        }
        locatedBranches_.push_back(cone);

        NodePath* best_path = nullptr;
        float best_length = std::numeric_limits<float>::infinity();
        for(auto& path: paths) {
            float length = path.getLength();
            if(length < best_length) {
                best_path = &path;
                best_length = length;
            }
        }
        tgtAssert(best_path, "No best path");

        PathTreeNode* path = new PathTreeNode(*best_path, to_process.parent_path);
        if(!root) {
            root.reset(path);
        }

        if(to_process.t_branch.edge_) {
            branches_to_process.push(UnprocessedBranch(path, to_process.t_branch.edge_->main_branch_));
            branches_to_process.push(UnprocessedBranch(path, to_process.t_branch.edge_->off_branch_));
        }
    }

    if(!missingBranches_.empty()) {
        LWARNING("" << missingBranches_.size() << " branche(s) from template could not be matched");
    }

    if(root) {
        // No node => did not match any branch
        fillGraph(root.get(), builder);
    }
    return std::move(builder).finalize();
}

void TemplateSubgraphExtractor::initialize() {
    GeometryRendererBase::initialize();
}

void TemplateSubgraphExtractor::deinitialize() {
    GeometryRendererBase::deinitialize();
}

void TemplateSubgraphExtractor::process() {
    const VesselGraph* input = inport_.getData();
    if(!input) {
        outport_.setData(nullptr);
        return;
    }
    if(!enabled_.get()) {
        // Pass input through
        outport_.setData(input, false);
        return;
    }
    if(templateFile_.get().empty()) {
        LERROR("No template file");
        outport_.setData(nullptr);
        return;
    }

    TemplateBranch t = parse_template(templateFile_.get());

    tgt::vec3 startingPoint;
    if(!tryExtractCenterPoint(startingPoint_, startingPoint)) {
        LWARNING("No starting point!");
        outport_.setData(nullptr);
        return;
    }

    std::unique_ptr<VesselGraph> output;
    switch(solvingStrategy_.getValue()) {
        case STRATEGY_GREEDY:
            output = extractSubgraph(*input, startingPoint, t, keepBounds_.get());
            break;
        case STRATEGY_GLOBAL:
            output = extractSubgraphGlobal(*input, startingPoint, t, keepBounds_.get());
            break;
    }

    outport_.setData(output.release());
}

void TemplateSubgraphExtractor::renderTransparent() {
    renderCones(true);
}

void TemplateSubgraphExtractor::render() {
    renderCones(false);
}

void TemplateSubgraphExtractor::renderCones(bool useTransparency) {
    auto mode = coneRenderMode_.getValue();
    if(mode == CONES_MISSING || mode == CONES_ALL) {
        for(const auto& branch: missingBranches_) {
            renderBranchSearchCone(branch, missingBranchConeColor_.get(), useTransparency);
        }
    }
    if(mode == CONES_ALL) {
        for(const auto& branch: locatedBranches_) {
            renderBranchSearchCone(branch, locatedBranchConeColor_.get(), useTransparency);
        }
    }
}

void TemplateSubgraphExtractor::renderBranchSearchCone(const BranchSearchCone& branch, const tgt::vec4& color, bool useTransparency) {
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();

    tgt::ImmediateMode::LightSource l;
    l.position = lightPosition_.get();
    l.ambientColor = lightAmbient_.get().xyz();
    l.diffuseColor = lightDiffuse_.get().xyz();
    l.specularColor = lightSpecular_.get().xyz();
    IMode.setLightSource(l);

    tgt::ImmediateMode::Material m;
    m.shininess = materialShininess_.get();
    IMode.setMaterial(m);

    //float min_l = branch.min_length;
    //float max_l = std::min(branch.max_length, 1000.0f);
    float min_l = 0;
    float max_l = coneLength_.get()*cos(branch.dir.max_deviation_);
    float tan_a = tan(branch.dir.max_deviation_);
    float lower_radius = tan_a * min_l;
    float upper_radius = tan_a * max_l;

    /// triangle mesh for rendering a cone for a missing branch
    /// top: (0,0,1), bottom (0,0,0)
    GlMeshGeometryUInt16Normal coneTriangleMesh;
    coneTriangleMesh.setCylinderGeometry(tgt::vec4::one, lower_radius, upper_radius, 1, 32, 4, false, false);

    MatStack.translate(branch.node_pos);
    MatStack.multMatrix(generateMatrixFromQuat(tgt::generateQuaternionFromTo(branch.dir.direction_.direction_, tgt::vec3(0,0,1))));
    MatStack.translate(tgt::vec3(0,0,min_l));
    MatStack.scale(tgt::vec3(1,1,max_l - min_l));

    IMode.color(color);

    GlMeshGeometryUInt16Normal::renderDefault(coneTriangleMesh, useTransparency);

    MatStack.popMatrix();
}
void TemplateSubgraphExtractor::adjustColorPropertyVisibility() {
    auto m = coneRenderMode_.getValue();
    missingBranchConeColor_.setVisibleFlag(m == CONES_ALL || m == CONES_MISSING);
    locatedBranchConeColor_.setVisibleFlag(m == CONES_ALL);
}

} // namespace voreen
