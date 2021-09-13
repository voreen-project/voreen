/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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
#include "voreen/core/processors/renderprocessor.h"

#include <voreen/core/datastructures/volume/volumeatomic.h>
#include <include/voreen/core/ports/volumeport.h>
#include "voreen/core/properties/string/stringlistproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include <modules/ensembleanalysis/ports/ensembledatasetport.h>

#ifndef VRN_VOLUMEOVERLAPTRACKER_H
#define VRN_VOLUMEOVERLAPTRACKER_H

namespace voreen {

    class PlotLibrary;

    class VRN_CORE_API VolumeOverlapTracker : public RenderProcessor {
    public:
        VolumeOverlapTracker();
        ~VolumeOverlapTracker();
        virtual Processor* create() const;

        virtual std::string getCategory() const { return "Volume Processing"; }
        virtual std::string getClassName() const { return "VolumeOverlapTracker"; }
        virtual CodeState getCodeState() const { return CODE_STATE_STABLE; }

        struct vec3Comparable {
            tgt::vec3 voxel;

            vec3Comparable(int x, int y, int z) {
                voxel = tgt::vec3(x, y, z);
            }

            vec3Comparable(tgt::vec3 voxel_) {
                voxel = voxel_;
            }

            bool operator==(const vec3Comparable& voxel_) const {
                return (voxel_.voxel.x == voxel.x && voxel_.voxel.y == voxel.y && voxel_.voxel.z == voxel.z);
            }

            bool operator<(const vec3Comparable& voxel_) const {
                return (voxel_.voxel.x < voxel.x
                    || voxel_.voxel.x == voxel.x && voxel_.voxel.y < voxel.y
                    || voxel_.voxel.x == voxel.x && voxel_.voxel.y == voxel.y && voxel_.voxel.z < voxel.z);
            }
        };

        struct Edge;

        struct Feature {
            int id, voxelCount, column;

            std::vector<Edge*> up, down;

            std::set<vec3Comparable> voxels;

            tgt::vec2 posOnCanvas;

            tgt::vec4 color;

            Feature(int id_) {
                id = id_;
                color = tgt::vec4(0, 0, 0, 1);
                voxelCount = 0;
                column = -1;
            }

            bool operator<(const Feature& feature) const {
                return (voxelCount < feature.voxelCount);
            }
        };

        struct Edge {
            Feature* f0;
            Feature* f1;

            float is0, is1;

            Edge(Feature* f0_, Feature* f1_, float is0_, float is1_) {
                f0 = f0_;
                f1 = f1_;
                is0 = is0_;
                is1 = is1_;
            }
        };

    protected:

        virtual void setDescriptions() {
            setDescription("Tracks Features over multiple Timesteps via volume overlapping");
        }

        virtual void process();
        virtual void onEvent(tgt::Event* e);

    private:
        void trackOverlap();

        void traversEdges(Feature* feature, int i);
        void traversEdgesSort(Feature* feature, int i);

        void renderAxes();
        void renderTreePlot();
        void renderQuad(tgt::vec2 position, float size, tgt::vec4 color);
        void renderCircle(tgt::vec2 position, float size, tgt::vec4 color, tgt::vec2 canvasSize);
        void renderLine(tgt::vec2 position0, tgt::vec2 position1, tgt::vec4 color);
        void renderLine(tgt::vec2 position0, tgt::vec2 position1, float width, tgt::vec4 color);

        tgt::vec2 positionFeature(int column, int row, float columnSize, float rowSize, tgt::vec2 canvasSize);
        tgt::vec2 mapCanvasPosToGLPos(tgt::vec2 position, tgt::vec2 canvas);
        tgt::vec2 mapCanvasPosToRowCol(tgt::vec2 position, int columnSize, int rowSize, int rowCount);
        tgt::vec2 map(tgt::vec2 value, tgt::vec2 max);
        tgt::vec4 generateRandomColor();

        bool containsFeature(std::vector<Feature*> list, Feature* feature);
        bool containsEdge(std::vector<Edge*> list, Edge* edge);

        std::vector<std::map<int, Feature*>> features_;
        std::vector<std::set<int>> featuresKeys_;
        std::vector<std::vector<Edge*>> edges_;

        VolumeListPort inport_;

        EnsembleDatasetPort originalInport_;

        RenderPort outport_;

        ButtonProperty track_;

        FloatProperty threshold_;

        IntProperty rowsOnCanvas_, minTimestepSelector_;

        IntIntervalProperty columnsOnCanvas_;

        BoolProperty enableIdSelector_;

        StringListProperty idSelector_;

        StringOptionProperty colorReference_, colorReferenceMode_;

        std::unique_ptr<PlotLibrary> plotLib_;

        static const std::string loggerCat_; ///< category used in logging
    };

}   //namespace

#endif // VRN_VOLUMEOVERLAPTRACKER_H
