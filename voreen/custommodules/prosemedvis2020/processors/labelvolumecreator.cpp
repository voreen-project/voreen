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

#include "labelvolumecreator.h"
#include "custommodules/flowsimulation/utils/utils.h"

namespace voreen {

	LabelVolumeCreator::LabelVolumeCreator() : Processor(),
		maskInport_(Port::INPORT, "mask.inport", "Mask Input"),
		concreteVesselGraphInport_(Port::INPORT, "ConcreteVesselGraph", "Concrete Vessel Graph Input"),
		labelVolumeOutport_(Port::OUTPORT, "mask.outport", "Label Volume Outport")
	{
		addPort(maskInport_);
		addPort(concreteVesselGraphInport_);
		addPort(labelVolumeOutport_);
	}

	Processor* LabelVolumeCreator::create() const {
		return new LabelVolumeCreator();
	}

	std::vector<tgt::svec3> getVoxelInSphere(tgt::svec3 pos, tgt::vec3 radius) { 
		std::vector<tgt::svec3> voxelinSphere;
		tgt::vec3 fpos = pos;

		tgt::svec3 lfl_corner =  tgt::svec3(floor(fpos.x - radius.x), floor(fpos.y - radius.y), floor(fpos.z - radius.z));
		tgt::svec3 ubr_corner = tgt::svec3(ceil(fpos.x + radius.x), ceil(fpos.y + radius.y), ceil(fpos.z + radius.z) );

		// Check for points in Sphere, only check Points in cube around pos
		tgt::svec3 vox = { 0, 0, 0 };
		for (int x = lfl_corner.x; x <= ubr_corner.x; ++x) {
			for (int y = lfl_corner.y; y <= ubr_corner.y; ++y) {
				for (int z = lfl_corner.z; z <= ubr_corner.z; ++z) {
					vox.x = x;
					vox.y = y;
					vox.z = z;
					
					tgt::svec3 diff = (vox - pos); 

					float dist = tgt::length(diff);
					float d = tgt::hadd(tgt::vec3(diff * diff)/(radius*radius));

					if (d <= 1) {
						voxelinSphere.push_back(vox);
					}
				}
			}
		}
		return voxelinSphere;
	}

	std::vector<tgt::svec3> getVoxelInHalfSphere(tgt::svec3 pos, tgt::vec3 radius, tgt::vec3 normal) {
		tgt::vec3 posToVoxel;
		std::vector<tgt::svec3> voxelinHalfSphere;
		std::vector<tgt::svec3> voxelInSphere = getVoxelInSphere(pos, radius);

		//Half the sphere with imaginary disk with given normal, to prevent endpoint parts from overextending
		for (tgt::svec3 voxel : voxelInSphere) {
			posToVoxel = { (float)voxel.x - pos.x, (float)voxel.y - pos.y, (float)voxel.z - pos.z };
			float dot = tgt::dot(posToVoxel,normal);
			if (dot < 0) {
				voxelinHalfSphere.push_back(voxel);
			}
		}
		return voxelinHalfSphere;
	}


	void LabelVolumeCreator::process() {
		tgt::vec3 normal;
		std::vector<tgt::svec3> samples;
		std::vector<ConcreteVesselGraphEdge> renderingsequence;
		const VolumeBase* mask = maskInport_.getData();
		tgt::svec3 dimensions = mask->getDimensions();
		const ConcreteVesselGraph* cVesselGraph = concreteVesselGraphInport_.getData();
		std::vector<ConcreteVesselGraphEdge> edges = cVesselGraph->getEdges();
		// data containing the segmentation labels
		VolumeRAM_UInt8* representation = new VolumeRAM_UInt8(dimensions);
		VolumeRAMRepresentationLock maskVolume(mask);

		// copy for getNodes
		ConcreteVesselGraph cg = *cVesselGraph;

		uint8_t edgeID;
		representation->clear();

		// Renderingsequence
		for (ConcreteVesselGraphEdge edge : edges) {
			if (cg.getNode(edge.getStart()).getType() == BOUNDARY_NODE) {
				renderingsequence.push_back(edge);
			}
		}
		for (ConcreteVesselGraphEdge edge : edges) {
			if (cg.getNode(edge.getStart()).getType() == BIFURCATION_NODE) {
				renderingsequence.push_back(edge);
			}
		}
		for (ConcreteVesselGraphEdge edge : edges) {
			if (cg.getNode(edge.getStart()).getType() == END_NODE) {
				renderingsequence.push_back(edge);
			}
		}

		for (ConcreteVesselGraphEdge edge : renderingsequence) {
			std::vector<ConcreteVesselGraphEdgePoint> edgePoints = edge.getEdgePoints();
			edgeID = edge.getID();

			// Go backwards through each edge (better for bifurcation)
			size_t i = edgePoints.size()-1;
			while (i > 0) {
				
				// Aprroximate disk normal with current trajectory
				if (i == 0) {
					normal = mask->getWorldToVoxelMatrix() * (edgePoints[i].getPos() - edgePoints[i + 1].getPos());
				}
				else {
					normal = mask->getWorldToVoxelMatrix() * (edgePoints[i - 1].getPos() - edgePoints[i].getPos());
				}

				samples = getVoxelInHalfSphere(mask->getWorldToVoxelMatrix() * edgePoints[i].getPos(), edge.getRadius() / mask->getSpacing(), normal);

				// Color Volume, dont overwrite voxels
				for (tgt::svec3 sample : samples) {
					if (sample != tgt::min(sample, representation->getDimensions() - tgt::svec3::one) || maskVolume->getVoxelNormalized(sample) == 0 || representation->voxel(sample) != 0) {
						continue;
					}
					representation->voxel(sample) = edgeID;
				}
				i--;
			}		
		}
		Volume* segmentationVolume = new Volume(representation, mask);
		segmentationVolume->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(representation->getBaseType()));
		labelVolumeOutport_.setData(segmentationVolume);
	}
}