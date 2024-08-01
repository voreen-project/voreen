/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include <modules/ensembleanalysis/utils/ensemblehash.h>
#include "modules/ensembleanalysis/utils/utils.h"

#include "volumeoverlaptracker.h"

#include "modules/plotting/datastructures/plotcell.h"
#include "modules/plotting/datastructures/plotdata.h"
#include "modules/plotting/utils/plotlibrary/plotlibrary.h"
#include "modules/plotting/utils/plotlibrary/plotlibraryopengl.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

namespace voreen {

    const std::string VolumeOverlapTracker::loggerCat_("voreen.FeatureExtractor");

    static const tgt::ivec2 MARGINS(75, 50);

    //filter single Feature
    bool filterFeature = false, resetted = true;
    tgt::vec2 filtlerMouseClickPos;
    int startRow = -1;
    VolumeOverlapTracker::Feature* selectedFeature = NULL;
    std::map<int, std::vector<VolumeOverlapTracker::Feature*>> featuresInRowsSpecific;
    std::map<int, std::vector<VolumeOverlapTracker::Edge*>> filteredEdgesSpecific;

    //filter necessary Features and Edges
    std::map<int, std::vector<VolumeOverlapTracker::Edge*>> filteredEdges;
    std::map<int, std::vector<VolumeOverlapTracker::Feature*>> featuresInRows;

    //sort Features in Rows
    std::vector<int> rowOffsets;
    int columnOffset = 0;

    //scroll
    float scrollOffset = 0;
    float scrollOffsetTick = 10;

    VolumeOverlapTracker::VolumeOverlapTracker()
        : RenderProcessor()
        , inport_(Port::INPORT, "volumehandle.inport", "Volumelist Inport", false)
        , originalInport_(Port::INPORT, "ensembledatastructurein", "Ensemble Datastructure Input", false)
        , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", true, voreen::Processor::INVALID_RESULT, voreen::RenderPort::RENDERSIZE_RECEIVER)
        , track_("apply", "Track Overlap", Processor::INVALID_RESULT)
        , threshold_("threshold", "Threshold", 0.5f, 0.0f, 1.0f)
        , rowsOnCanvas_("rowsOnCanvas", "Rows on Canvas", 10, 1, 20)
        , columnsOnCanvas_("columnSelector", "Columns On Canvas", tgt::vec2(0, 0), 0, 0)
        , enableIdSelector_("enableIdSelector", "enable ID selection")
        , minTimestepSelector_("minTimestepSelector", "select Min Timestep", 0, 0, 250)
        , idSelector_("idSelector", "selected IDs")
        , colorReference_("colorReference", "Color Reference")
        , colorReferenceMode_("colorReferenceMode", "Color Reference Mode")
        , plotLib_(new PlotLibraryOpenGl())
    {
        addPort(outport_);
        addPort(inport_);
        addPort(originalInport_);

        addProperty(track_);
        ON_CHANGE(track_, VolumeOverlapTracker, trackOverlap);
        addProperty(threshold_);

        addProperty(rowsOnCanvas_);
        addProperty(columnsOnCanvas_);
        rowsOnCanvas_.setGroupID("cr");
        columnsOnCanvas_.setGroupID("cr");
        setPropertyGroupGuiName("cr", "Column and Rows Properties");

        addProperty(enableIdSelector_);
        addProperty(minTimestepSelector_);
        addProperty(idSelector_);
        idSelector_.setVisibleFlag(false);
        enableIdSelector_.setGroupID("idSelection");
        minTimestepSelector_.setGroupID("idSelection");
        idSelector_.setGroupID("idSelection");
        setPropertyGroupGuiName("idSelection", "ID Selection");

        addProperty(colorReference_);
        addProperty(colorReferenceMode_);
        colorReference_.setGroupID("color");
        colorReferenceMode_.setGroupID("color");
        setPropertyGroupGuiName("color", "Color");

        colorReference_.addOption("angle", "angle");
        colorReference_.addOption("magnitude", "magnitude");
        colorReference_.addOption("radius", "radius");
        colorReference_.addOption("spin transition-induced density anomaly", "spin transition-induced density anomaly");
        colorReference_.addOption("temperature", "temperature");
        colorReference_.addOption("temperature anomaly", "temperature anomaly");
        colorReference_.addOption("thermal conductivity", "thermal conductivity");
        colorReference_.addOption("thermal expansivity", "thermal expansivity");
        colorReference_.addOption("vx", "vx");
        colorReference_.addOption("vy", "vy");
        colorReference_.addOption("vz", "vz");

        colorReferenceMode_.addOption("median", "median");
        colorReferenceMode_.addOption("mean", "mean");
    }

    VolumeOverlapTracker::~VolumeOverlapTracker() {
        for (int i = 0; i < edges_.size(); i++) {
            for (int edge = 0; edge < edges_[i].size(); edge++) {
                delete edges_[i][edge];
            }
        }
        for (int i = 0; i < features_.size(); i++) {
            for (int key : featuresKeys_[i]) {
                delete features_[i][key];
            }
        }
    }

    Processor* VolumeOverlapTracker::create() const {
        return new VolumeOverlapTracker();
    }

    void VolumeOverlapTracker::process() {
        if (enableIdSelector_.get()) {
            idSelector_.setReadOnlyFlag(false);
            minTimestepSelector_.setReadOnlyFlag(false);
        }
        else {
            idSelector_.setReadOnlyFlag(true);
            minTimestepSelector_.setReadOnlyFlag(true);
        }

        minTimestepSelector_.setMaxValue(250);

        // Change depth func to apply rendering order properly.
        glDepthFunc(GL_LEQUAL);

        // Render picking pass.
        glLineWidth(7.0f);

        outport_.activateTarget();
        renderAxes();
        if (edges_.size() > 0) {
            renderTreePlot();
        }
        outport_.deactivateTarget();

        // Restore state.
        glDepthFunc(GL_LESS);
        glLineWidth(1.0f);
        IMode.color(tgt::col4::one);
    }

    void VolumeOverlapTracker::trackOverlap() {
        //reset
        for (int i = 0; i < edges_.size(); i++) {
            for (int edge = 0; edge < edges_[i].size(); edge++) {
                delete edges_[i][edge];
            }
        }
        for (int i = 0; i < features_.size(); i++) {
            for (int key : featuresKeys_[i]) {
                delete features_[i][key];
            }
        }
        features_.clear();
        featuresKeys_.clear();
        edges_.clear();
        filterFeature = false;
        resetted = true;

        //prep Data
        const VolumeList* dataList = inport_.getData();

        tgt::vec2 minmax = tgt::vec2(1000000, -1000000);
        for (int timestep = 0; timestep < dataList->size(); timestep++) {
            VolumeRAMRepresentationLock data(dataList->at(timestep));

            RealWorldMapping rwm = dataList->at(timestep)->getRealWorldMapping();

            tgt::vec3 dimensions = data->getDimensions();

            std::map<int, Feature*> features;
            std::set<int> keys;
//#pragma omp parallel for
        //for (long z = 0; z < static_cast<long>(dimensions.z); z++) {
            //extract features
            for (size_t z = 0; z < dimensions.z; z++) {
                for (size_t y = 0; y < dimensions.y; y++) {
                    for (size_t x = 0; x < dimensions.x; x++) {
                        float value = rwm.normalizedToRealWorld(data->getVoxelNormalized(x, y, z, 0));
                        if (value != 0) {
                            if (features[value] == NULL) {
                                features[value] = new Feature(value);
                                keys.insert(value);
                            }
                            features[value]->voxels.insert(tgt::vec3(x, y, z));
                            features[value]->voxelCount++;
                        }
                    }
                }
            }

            //set color
            //daten auslesen und bearbeitbar machen
            const EnsembleDataset& ensemble = *originalInport_.getData();

            for (const EnsembleMember& member : ensemble.getMembers()) {
                if (dataList->size() == member.getTimeSteps().size()) {
                    int fieldNameSize = member.getTimeSteps()[timestep].getFieldNames().size();
                    int i = timestep * fieldNameSize;
                    int fieldnameIndex = 0;
                    for (int f = 0; f < member.getTimeSteps()[timestep].getFieldNames().size(); f++) {
                        if (member.getTimeSteps()[timestep].getFieldNames()[f] == colorReference_.get()) {
                            fieldnameIndex = f;
                        }
                    }
                    VolumeRAMRepresentationLock dataColor(ensemble.getVolumes()[i + fieldnameIndex]);

                    RealWorldMapping rwmColor = ensemble.getVolumes()[i + fieldnameIndex]->getRealWorldMapping();

                    //mean
                    if (colorReferenceMode_.get() == "mean") {
                        for (auto pair : features) {
                            float sum = 0;
                            for (vec3Comparable voxel : pair.second->voxels) {
                                float value = rwmColor.normalizedToRealWorld(dataColor->getVoxelNormalized(voxel.voxel, 0));
                                sum += value;
                            }
                            sum /= pair.second->voxelCount;
                            if (sum < minmax.x) {
                                minmax.x = sum;
                            }
                            else if (sum > minmax.y) {
                                minmax.y = sum;
                            }
                            pair.second->color = tgt::vec4(sum, 0, 0, 1);
                        }
                    }//median
                    else if (colorReferenceMode_.get() == "median") {
                        for (auto pair : features) {
                            float oneBack = 0, current = 0;
                            bool first = true;
                            for (vec3Comparable voxel : pair.second->voxels) {
                                float value = rwmColor.normalizedToRealWorld(dataColor->getVoxelNormalized(voxel.voxel, 0));
                                if (first) {
                                    current = value;
                                    first = false;
                                }
                                if (oneBack < current && value < current) {
                                    current = (oneBack < value) ? value : oneBack;
                                }
                                else if (oneBack > current && value > current) {
                                    current = (oneBack < value) ? oneBack : value;
                                }
                                oneBack = value;
                            }
                            if (current < minmax.x) {
                                minmax.x = current;
                            }
                            else if (current > minmax.y) {
                                minmax.y = current;
                            }
                            pair.second->color = tgt::vec4(current, 0, 0, 1);
                        }
                    }
                }
                else {
                    std::cout << "<FeatureTracker>: Number of Timesteps do not match! Ensemble: " << member.getTimeSteps().size() << " connected: " << dataList->size() << std::endl;
                }
            }

            for (auto pair : features) {
                float value = 0;
                if (minmax.x > 0) {
                    value = (pair.second->color.x - minmax.x) / (minmax.y - minmax.x);
                }
                else {
                    value = (pair.second->color.x + (-1 * minmax.x)) / (minmax.y + (-1 * minmax.x));
                }
                pair.second->color = tgt::vec4(1 - value, 1 - value, 1 - value, 1);
            }

            //push into global variables
            features_.push_back(features);
            featuresKeys_.push_back(keys);

            //look up intersection and generate corresponding edges
            if (features_.size() > 1) {
                for (int i = features_.size() - 1; i < features_.size(); i++) {
                    std::map<int, Feature*> features0 = features_[i - 1];
                    std::map<int, Feature*> features1 = features_[i];

                    std::vector<Edge*> edges;

                    for (auto pair0 : features0) {
                        for (auto pair1 : features1) {
                            std::vector<vec3Comparable> backinsert;

                            std::set_intersection(pair0.second->voxels.begin(), pair0.second->voxels.end(),
                                pair1.second->voxels.begin(), pair1.second->voxels.end(),
                                std::back_inserter(backinsert));

                            float is0 = (float)backinsert.size() / (float)pair0.second->voxels.size();
                            float is1 = (float)backinsert.size() / (float)pair1.second->voxels.size();

                            Edge* edge = new Edge(pair0.second, pair1.second, is0, is1);

                            pair0.second->up.push_back(edge);
                            pair1.second->down.push_back(edge);

                            edges.push_back(edge);
                        }
                        //clear voxels list to not overflow RAM
                        pair0.second->voxels.clear();
                    }
                    edges_.push_back(edges);                    
                }
            }
        } 

        //clear voxels from last Features
        for (auto pair : features_[features_.size() - 1]) {
            pair.second->voxels.clear();
        }
    }

    void VolumeOverlapTracker::renderAxes() {
        // Set Plot status.
        plotLib_->setWindowSize(outport_.getSize());
        plotLib_->setAxesWidth(1.0f);
        plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
        plotLib_->setLineWidth(1.0f);
        plotLib_->setMaxGlyphSize(1.0f);
        plotLib_->setMarginBottom(MARGINS.y);
        plotLib_->setMarginTop(MARGINS.y);
        plotLib_->setMarginLeft(MARGINS.x);
        plotLib_->setMarginRight(MARGINS.x);
        plotLib_->setMinimumScaleStep(32, PlotLibrary::X_AXIS);
        plotLib_->setMinimumScaleStep(32, PlotLibrary::Y_AXIS);
        plotLib_->setDomain(Interval<plot_t>(-1.0, 1.0), PlotLibrary::X_AXIS);
        plotLib_->setDomain(Interval<plot_t>(-1.0, 1.0), PlotLibrary::Y_AXIS);
        plotLib_->setDimension(PlotLibrary::TWO);

        if (plotLib_->setRenderStatus()) {
            plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
            plotLib_->renderAxes();
            plotLib_->setDrawingColor(tgt::Color(0, 0, 0, .5f));
            plotLib_->setFontSize(12);
            plotLib_->renderAxisLabel(PlotLibrary::X_AXIS, "Features");
            plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, "Timesteps");
        }
        plotLib_->resetRenderStatus();

        // Reset state.
        IMode.color(1.0f, 1.0f, 1.0f, 1.0f);
        LGL_ERROR;
    }

    void VolumeOverlapTracker::renderTreePlot() {
        //reset
        filteredEdges.clear();
        filteredEdgesSpecific.clear();
        featuresInRows.clear();
        featuresInRowsSpecific.clear();
        for (int i = 0; i < features_.size(); i++) {
            for (int key : featuresKeys_[i]) {
                features_[i][key]->column = -1;
            }
        }

        tgt::vec2 canvasSize = outport_.getSize();
        tgt::vec2 plotSize = tgt::vec2(canvasSize.x - (MARGINS.x * 2), canvasSize.y - (MARGINS.y * 2));

        //filter Edges with Threshold
        for (int i = 0; i < edges_.size(); i++) {
            std::vector<Edge*> edges;
            for (Edge* edge : edges_[i])
                if (edge->is0 >= threshold_.get() || edge->is1 >= threshold_.get()) {
                    edges.push_back(edge);
                }
            filteredEdges[i] = edges;
        }

        //save which Features go in what Row
        for (int i = 0; i < filteredEdges.size(); i++) {
            for (Edge* edge : filteredEdges[i]) {
                if (!containsFeature(featuresInRows[i], edge->f0)) {
                    featuresInRows[i].push_back(edge->f0);
                }
                if (!containsFeature(featuresInRows[i + 1], edge->f1)) {
                    featuresInRows[i + 1].push_back(edge->f1);
                }
            }
        }
        /*
        //set color to largest connected Feature
        for (int i = 1; i < featuresInRows.size(); i++) {
            for (Feature* feature : featuresInRows[i]) {
                int maxSize = 0;
                for (Edge* down : feature->down) {
                    if (containsEdge(filteredEdges[i - 1], down) && down->f0->voxelCount > maxSize) {
                        maxSize = down->f0->voxelCount;
                        feature->color = down->f0->color;
                    }
                }
            }
        }*/

        //sort Features in Rows
        auto compareFeatures = [](Feature* a, Feature* b) -> bool {
            return (a->voxelCount < b->voxelCount);
        };
        for (int row = 0; row < featuresInRows.size(); row++) {
            std::sort(featuresInRows[row].begin(), featuresInRows[row].end(), compareFeatures);
            std::reverse(featuresInRows[row].begin(), featuresInRows[row].end());
        }

        //sort rows before filtering
        columnOffset = 0;
        for (int row = 0; row < featuresInRows.size(); row++) {
            for (Feature* feature : featuresInRows[row]) {
                //reset OffsetVector
                rowOffsets.clear();
                for (int i = 0; i < featuresInRows.size(); i++) {
                    rowOffsets.push_back(0);
                }
                traversEdgesSort(feature, row);
                //set columnOffset
                int max = 0;
                for (int offset : rowOffsets) {
                    if (max < offset) {
                        max = offset;
                    }
                }
                columnOffset += max;
            }
        }

        //filter single Feature
        if (filterFeature) {
            if (resetted) {
                tgt::vec2 mouseClickPos = filtlerMouseClickPos;
                mouseClickPos.y += scrollOffset;
                int featuresOnCanvas = columnsOnCanvas_.get().y - columnsOnCanvas_.get().x;
                tgt::vec2 rowColumn = mapCanvasPosToRowCol(mouseClickPos, plotSize.x / featuresOnCanvas, plotSize.y / featuresInRows.size(), featuresInRows.size());
                if (rowColumn.x != -1 && rowColumn.y != -1) {
                    for (Feature* feature : featuresInRows[rowColumn.x]) {
                        if (feature->column - columnsOnCanvas_.get().x == rowColumn.y) {
                            selectedFeature = feature;
                        }
                    }
                }
                startRow = rowColumn.x;
                resetted = false;
            }
            if (selectedFeature != NULL) {
                traversEdges(selectedFeature, startRow);

                featuresInRows = featuresInRowsSpecific;
                filteredEdges = filteredEdgesSpecific;

                //set Properties for Linking
                if (enableIdSelector_.get()) {
                    idSelector_.reset();

                    if (startRow >= 0 && startRow < featuresInRows.size()) {
                        for (Feature* feature : featuresInRows[startRow]) {
                            idSelector_.addRow(std::to_string(feature->id));
                        }
                    }
                }
            }
            else {
                featuresInRows.clear();
                filteredEdges.clear();
            }
            //sort rows again
            columnOffset = 0;
            for (int row = 0; row < featuresInRows.size(); row++) {
                for (Feature* feature : featuresInRows[row]) {
                    //reset OffsetVector
                    rowOffsets.clear();
                    for (int i = 0; i < featuresInRows.size(); i++) {
                        rowOffsets.push_back(0);
                    }
                    traversEdgesSort(feature, row);
                    //set columnOffset
                    int max = 0;
                    for (int offset : rowOffsets) {
                        if (max < offset) {
                            max = offset;
                        }
                    }
                    columnOffset += max;
                }
            }
        }
        else {
            resetted = true;
            selectedFeature = NULL;
            columnsOnCanvas_.setMaxValue(columnOffset);
        }
        
        //prep canvas distance variables
        int rowCount = features_.size();
        if (rowCount > rowsOnCanvas_.get()) {
            rowCount = rowsOnCanvas_.get();
        }
        int maxFeatureSize = 0;
        for (int i = 0; i < featuresInRows.size(); i++) {
            for (Feature* feature : featuresInRows[i]) {
                if (feature->voxelCount > maxFeatureSize) {
                    maxFeatureSize = feature->voxelCount;
                }
            }
        }
        int columnCount = (filterFeature) ? columnOffset : columnsOnCanvas_.get().y - columnsOnCanvas_.get().x;
        float columnSize = plotSize.x / columnCount;
        float rowSize = plotSize.y / rowCount;
        scrollOffsetTick = rowSize;

        //place Features
        for (int row = 0; row < featuresInRows.size(); row++) {
            for (int i = 0; i < featuresInRows[row].size(); i++) {
                Feature* current = featuresInRows[row][i];
                if (filterFeature) {
                    current->posOnCanvas = positionFeature(current->column, row, columnSize, rowSize, canvasSize);
                }
                else {
                    if (current->column >= columnsOnCanvas_.get().x && current->column <= columnsOnCanvas_.get().y) {
                        current->column -= columnsOnCanvas_.get().x;
                        current->posOnCanvas = positionFeature(current->column, row, columnSize, rowSize, canvasSize);
                    }
                    else {
                        current->column = -1;
                    }
                }
            }
        }

        //prep scrolloffset
        float scrollOffsetGL = (scrollOffset / canvasSize.y) * 2;

        //render Lines between intersecting Features
        for (int i = 0; i < filteredEdges.size(); i++) {
            for (Edge* edge : filteredEdges[i]) {
                if (edge->f0->column != -1 && edge->f1->column != -1) {
                    float size0 = (edge->is0 * 6) + 1;
                    float size1 = (edge->is1 * 6) + 1;
                    tgt::vec2 f0Pos = edge->f0->posOnCanvas;
                    f0Pos.y += scrollOffsetGL;
                    tgt::vec2 f1Pos = edge->f1->posOnCanvas;
                    f1Pos.y += scrollOffsetGL;
                    renderLine(f0Pos, f1Pos, size0, tgt::vec4(edge->is0, 0, 0, 0.5));
                    renderLine(f0Pos, f1Pos, size1, tgt::vec4(0, 0, edge->is1, 0.5));
                }
            }
        }

        //render Features
        for (int row = 0; row < featuresInRows.size(); row++) {
            for (Feature* feature : featuresInRows[row]) {
                if (feature->column != -1) {
                    float size = ((((float)feature->voxelCount / (float)maxFeatureSize) * 0.5) + 0.5) * columnSize / (float)canvasSize.x;
                    float height = ((((float)feature->voxelCount / (float)maxFeatureSize) * 0.5) + 0.5) * 100 / (float)canvasSize.y;
                    if (size > height)
                        size = height;
                    tgt::vec2 fPos = feature->posOnCanvas;
                    fPos.y += scrollOffsetGL;
                    renderCircle(fPos, size, feature->color, canvasSize);
                }
            }
        }

        //render Timestep Info
        for (int row = 0; row < featuresInRows.size(); row++) {
            tgt::Font font(VoreenApplication::app()->getFontPath("Vera.ttf"));
            font.setFontSize(18);
            font.setFontColor(tgt::vec4(0, 0, 0, 1));

            tgt::vec3 pos = tgt::vec3(MARGINS.x - 33, MARGINS.y + ((row * rowSize) + (rowSize / 2)), 0);
            pos.y += scrollOffset;

            font.render(pos, std::to_string(row + minTimestepSelector_.get()), canvasSize);
        }
    }

    void VolumeOverlapTracker::traversEdges(Feature* feature, int i) {
        if (containsFeature(featuresInRows[i], feature) && !containsFeature(featuresInRowsSpecific[i], feature)) {
            featuresInRowsSpecific[i].push_back(feature);
            feature->column = -1;
            for (Edge* up : feature->up) {
                if (containsEdge(filteredEdges[i], up) && !containsEdge(filteredEdgesSpecific[i], up)) {
                    filteredEdgesSpecific[i].push_back(up);
                    traversEdges(up->f1, i + 1);
                }
            }
            for (Edge* down : feature->down) {
                if (containsEdge(filteredEdges[i - 1], down))
                    traversEdges(down->f0, i - 1);
            }
        }
    }

    void VolumeOverlapTracker::traversEdgesSort(Feature* feature, int i) {
        if (containsFeature(featuresInRows[i], feature) && feature->column == -1) {
            feature->column = columnOffset + rowOffsets[i];
            rowOffsets[i]++;
            for (Edge* up : feature->up) {
                if (containsEdge(filteredEdges[i], up)) {
                    traversEdgesSort(up->f1, i + 1);
                }
            }
            for (Edge* down : feature->down) {
                if (containsEdge(filteredEdges[i - 1], down))
                    traversEdgesSort(down->f0, i - 1);
            }
        }
    }

    tgt::vec2 VolumeOverlapTracker::positionFeature(int column, int row, float columnSize, float rowSize, tgt::vec2 canvasSize) {
        float x = MARGINS.x + ((column * columnSize) + (columnSize / 2));
        float y = MARGINS.y + ((row * rowSize) + (rowSize / 2));
        tgt::vec2 position = mapCanvasPosToGLPos(tgt::vec2(x, y), canvasSize);

        return position;
    }

    void VolumeOverlapTracker::renderQuad(tgt::vec2 position, float size, tgt::vec4 color) {
        IMode.begin(tgt::ImmediateMode::QUADS);

        IMode.color(color);
        IMode.vertex(tgt::vec2(position.x - size, position.y - size));
        IMode.vertex(tgt::vec2(position.x + size, position.y - size));
        IMode.vertex(tgt::vec2(position.x + size, position.y + size));
        IMode.vertex(tgt::vec2(position.x - size, position.y + size));

        IMode.end();
    }

    void VolumeOverlapTracker::renderCircle(tgt::vec2 position, float size, tgt::vec4 color, tgt::vec2 canvasSize) {
        tgt::vec2 aspect = (canvasSize.x > canvasSize.y) ? tgt::vec2(canvasSize.y / canvasSize.x, 1) : tgt::vec2(1, canvasSize.x / canvasSize.y);

        //render Cirlce
        size_t circleResolution = 64;
        IMode.begin(tgt::ImmediateMode::POLYGON);
        IMode.color(color);
        for (size_t i = 0; i < circleResolution; ++i) {
            float angle = (static_cast<float>(i) / static_cast<float>(circleResolution)) * (2 * tgt::PIf);
            tgt::vec2 direction(cosf(angle), sinf(angle));
            direction *= size;

            tgt::vec2 pos = direction * aspect + position;
            IMode.vertex(tgt::vec3(pos.x, pos.y, 0));
        }
        IMode.end();

        //render black Line around
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        glLineWidth(1.f);
        IMode.color(tgt::vec4(0, 0, 0, 1));
        for (size_t i = 0; i < circleResolution; ++i) {
            float angle = (static_cast<float>(i) / static_cast<float>(circleResolution)) * (2 * tgt::PIf);
            tgt::vec2 direction(cosf(angle), sinf(angle));
            direction *= size;

            tgt::vec2 pos = direction * aspect + position;
            IMode.vertex(tgt::vec3(pos.x, pos.y, 0));
        }
        IMode.end();
    }

    void VolumeOverlapTracker::renderLine(tgt::vec2 position0, tgt::vec2 position1, tgt::vec4 color) {
        IMode.begin(tgt::ImmediateMode::LINES);

        IMode.color(color);
        IMode.vertex(position0);
        IMode.vertex(position1);

        IMode.end();
    }

    void VolumeOverlapTracker::renderLine(tgt::vec2 position0, tgt::vec2 position1, float width, tgt::vec4 color) {
        IMode.begin(tgt::ImmediateMode::LINES);

        IMode.color(color);
        glLineWidth(width);
        IMode.vertex(position0);
        IMode.vertex(position1);

        IMode.end();
    }

    void VolumeOverlapTracker::onEvent(tgt::Event* e) {
        // Events will be triggered for RenderProcessors, even if the processor is not ready it seems.
        if (isReady()) {
            tgt::MouseEvent* event = dynamic_cast<tgt::MouseEvent*>(e);
            if (event) {
                if (event->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT && event->action() == tgt::MouseEvent::PRESSED
                    && event->modifiers() == tgt::MouseEvent::MODIFIER_NONE) {
                    filterFeature = !filterFeature;
                    if (filterFeature) {
                        filtlerMouseClickPos = tgt::vec2(event->x(), event->y());
                    }
                    process();
                }
                if (event->button() == tgt::MouseEvent::MOUSE_WHEEL_DOWN) {
                    scrollOffset -= scrollOffsetTick;
                    process();
                }
                if (event->button() == tgt::MouseEvent::MOUSE_WHEEL_UP) {
                    scrollOffset += scrollOffsetTick;
                    process();
                }
            }
        }

        RenderProcessor::onEvent(e);
    }

    tgt::vec2 VolumeOverlapTracker::mapCanvasPosToGLPos(tgt::vec2 position, tgt::vec2 canvas) {
        tgt::vec2 mapped = map(position, canvas);

        mapped.x = (mapped.x * 2) - 1;
        mapped.y = (mapped.y * 2) - 1;

        return mapped;
    }

    tgt::vec2 VolumeOverlapTracker::mapCanvasPosToRowCol(tgt::vec2 position, int columnSize, int rowSize, int rowCount) {
        int row = (position.y - MARGINS.y) / rowSize;
        int column = (position.x - MARGINS.x) / columnSize;
        if (row < 0 || column < 0) {
            return tgt::vec2(-1, -1);
        }
        return tgt::vec2(((row - rowCount) * -1) - 1, column);
    }

    tgt::vec2 VolumeOverlapTracker::map(tgt::vec2 value, tgt::vec2 max) {
        return tgt::vec2(value.x / max.x, value.y / max.y);
    }

    tgt::vec4 VolumeOverlapTracker::generateRandomColor() {
        return tgt::vec4((float)rand() / RAND_MAX, (float)rand() / RAND_MAX, (float)rand() / RAND_MAX, 1);
    }

    bool VolumeOverlapTracker::containsFeature(std::vector<Feature*> list, Feature* feature) {
        for (Feature* f : list) {
            if (f->id == feature->id) {
                return true;
            }
        }
        return false;
    }

    bool VolumeOverlapTracker::containsEdge(std::vector<Edge*> list, Edge* edge) {
        for (Edge* e : list) {
            if (e->f0->id == edge->f0->id && e->f1->id == edge->f1->id) {
                return true;
            }
        }
        return false;
    }
}   // namespace