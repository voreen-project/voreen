/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCatmullClarkFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2001/11/13 14:13:55 $
  Version:   $Revision: 1.11 $
  Thanks:    This work was supported bt PHS Research Grant No. 1 P41 RR13218-01
             from the National Center for Research Resources

Copyright (c) 1993-2001 Leopold K�hschelm, Daniel Wagner, Sebastian Zambal
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "vtkCatmullClarkFilter.h"
#include "vtkEdgeTable.h"
#include "vtkObjectFactory.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"

vtkCatmullClarkFilter* vtkCatmullClarkFilter::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkCatmullClarkFilter");
  if(ret)
    {
    return (vtkCatmullClarkFilter*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkCatmullClarkFilter;
}

void vtkCatmullClarkFilter::insertEdgePoint(CellData &cell, EdgePoint *newElement) {
    newElement->next = cell.edgePoints;
    cell.edgePoints = newElement;
}

EdgePoint *vtkCatmullClarkFilter::findEdgePoint(CellData &cell, vtkIdType id1, vtkIdType id2) {
    EdgePoint *result = cell.edgePoints;
    bool found = false;
    while((result != NULL) && (found == false)) {
        if (((result->p1 == id1) && (result->p2 == id2)) ||
            ((result->p1 == id2) && (result->p2 == id1))) {     // found edge
            found = true;
        } else {
            result = result->next;
        }
    }
    return result;
}

void vtkCatmullClarkFilter::insertVertexPoint(CellData &cell, VertexPoint *newElement) {
    newElement->next = cell.vertexPoints;
    cell.vertexPoints = newElement;
}

VertexPoint *vtkCatmullClarkFilter::findVertexPoint(CellData &cell, vtkIdType id) {
    VertexPoint *result = cell.vertexPoints;
    bool found = false;
    while((result != NULL) && (found == false)) {
        if (result->v_old == id) {     // found vertex
            found = true;
        } else {
            result = result->next;
        }
    }
    return result;
}

void vtkCatmullClarkFilter::computeFacePoint(vtkCell *cell, float *facePoint) {
    facePoint[0] = 0.0f;
    facePoint[1] = 0.0f;
    facePoint[2] = 0.0f;
    vtkPoints *points = cell->GetPoints();
    int cpoints = cell->GetNumberOfPoints();
    for (int i = 0; i < cpoints; ++i) {
        auto current = points->GetPoint(i);
        facePoint[0] += (current[0] / cpoints);
        facePoint[1] += (current[1] / cpoints);
        facePoint[2] += (current[2] / cpoints);
    }
}

vtkIdList *vtkCatmullClarkFilter::getNeighborIds(vtkPolyData *input, vtkIdType i) {
    vtkIdList *neighbors = vtkIdList::New();

    vtkIdList *neighborCells = vtkIdList::New();
    input->GetPointCells(i, neighborCells);
    for (int j = 0; j < neighborCells->GetNumberOfIds(); ++j) {

        vtkIdList *cellPoints = vtkIdList::New();
        input->GetCellPoints(neighborCells->GetId(j), cellPoints);
        int cpoints = cellPoints->GetNumberOfIds();
        for (int k = 0; k < cpoints; ++k) {

            vtkIdType ptId1 = cellPoints->GetId(k);
            vtkIdType ptId2 = cellPoints->GetId((k+1) % cpoints);

            if (i == ptId1) {
                if (neighbors->IsId(ptId2) < 0) {
                    neighbors->InsertNextId(ptId2);
                }
            }
            if (i == ptId2) {
                if (neighbors->IsId(ptId1) < 0) {
                    neighbors->InsertNextId(ptId1);
                }
            }
        }
    }
    return neighbors;
}

void vtkCatmullClarkFilter::Subdivision(vtkPolyData *input, vtkPolyData *output) {
    // BuildLinks(), generate topological information
    input->BuildLinks();

    // create the building blocks of the output-polydata.
    vtkPoints    *outputPoints = vtkPoints::New();
    vtkCellArray *outputPolys  = vtkCellArray::New();

    // create the array holding the information about face-, edge- and
    // vertex-points
    CellData *cells = new CellData[input->GetNumberOfCells()];
    for (int i = 0; i < input->GetNumberOfCells(); ++i) {
        cells[i].edgePoints   = NULL;
        cells[i].vertexPoints = NULL;
    }
    int countInputCells = input->GetNumberOfCells();

    // for each cell...
    for (int i = 0; i < countInputCells;  i++) {

        // get the current cell
        vtkGenericCell *cell = vtkGenericCell::New();
        input->GetCell(i, cell);

        // get the number of points in the current cell
        int cpoints = cell->GetNumberOfPoints();

        //=============================================
        // Compute the face point and add it to the
        // output-points
        //=============================================
        float *facePoint = new float[3];
        computeFacePoint(cell, facePoint);
        cells[i].facePoint = outputPoints->InsertNextPoint(facePoint);

        //=============================================
        // Compute the edge points
        //=============================================
        for (int j = 0; j < cpoints; ++j) {

            // determine ids of the two vertices that form the edge
            vtkIdType id1 = cell->GetPointId(j);
            vtkIdType id2 = cell->GetPointId((j+1) % cpoints);

            // check: has the neighbouring cell already computed the edge-point?
            vtkIdList *neighborList = vtkIdList::New();
            input->GetCellEdgeNeighbors(i, id1, id2, neighborList);
            EdgePoint *edgePoint = NULL;
            vtkIdType neighbor = 0;
            if (neighborList->GetNumberOfIds() > 0) {
                neighbor = neighborList->GetId(0);      // assume, there is only one neighbor
                edgePoint = findEdgePoint(cells[neighbor], id1, id2);
            }

            if (edgePoint == NULL) {         // edge not found at neighbor; has to be computed

                // compute the face-point of the neighbor
                float *facePoint_n = new float[3];
                computeFacePoint(input->GetCell(neighbor), facePoint_n);

                float *edgePoint = new float[3];
                edgePoint[0] = ((input->GetPoint(id1))[0] + 
                            (input->GetPoint(id2))[0] +
                            (outputPoints->GetPoint(cells[i].facePoint))[0] +
                            facePoint_n[0]) / 4.0f;
                edgePoint[1] = ((input->GetPoint(id1))[1] +
                            (input->GetPoint(id2))[1] +
                            (outputPoints->GetPoint(cells[i].facePoint))[1] +
                            facePoint_n[1]) / 4.0f;
                edgePoint[2] = ((input->GetPoint(id1))[2] +
                            (input->GetPoint(id2))[2] +
                            (outputPoints->GetPoint(cells[i].facePoint))[2] +
                            facePoint_n[2]) / 4.0f;

                EdgePoint *newEdgePoint = new EdgePoint;
                newEdgePoint->p1 = id1;
                newEdgePoint->p2 = id2;
                newEdgePoint->e = outputPoints->InsertNextPoint(edgePoint);
                insertEdgePoint(cells[i], newEdgePoint);

            } else {                         // edge found at neighbor; insert given information

                EdgePoint *newEdgePoint = new EdgePoint;
                newEdgePoint->p1 = id1;
                newEdgePoint->p2 = id2;
                newEdgePoint->e = edgePoint->e;
                insertEdgePoint(cells[i], newEdgePoint);
            }
        }

        //=============================================
        // Compute the vertex points
        //=============================================
        for (int j = 0; j < cpoints; ++j) {

            // get the global id of the point
            vtkIdType id = cell->GetPointId(j);

            // check: has the neighboring cell already computed the vertex-point?
            vtkIdList *neighborCells = vtkIdList::New();
            input->GetPointCells(id, neighborCells);

            VertexPoint *vertexPoint = NULL;
            int k = 0;

            while ((vertexPoint == NULL) && 
                   (k < neighborCells->GetNumberOfIds())) {
                vertexPoint = findVertexPoint(cells[neighborCells->GetId(k)], id);
                ++k;
            }
            if (vertexPoint == NULL) {          // vertex point must be computed

                float *coord = new float[3];
                coord[0] = 0.0f;
                coord[1] = 0.0f;
                coord[2] = 0.0f;

                // find out the neighbors
                vtkIdList *neighbors = getNeighborIds(input, id);

                int n = neighbors->GetNumberOfIds();

                vtkIdList *neighborCells = vtkIdList::New();
                input->GetPointCells(id, neighborCells);
                int a = neighborCells->GetNumberOfIds();

                float weight_v = ((float)n-2.0f)/(float)n;
                float weight_e = 1.0f/(float)(n*n);
                //float weight_f = 1.0f/(float)(n*n);

                float *sum = new float[3];
                sum[0] = 0.0f;
                sum[1] = 0.0f;
                sum[2] = 0.0f;

                // add weighted sum of involved face-points
                float *fp = new float[3];
                for (int l = 0; l < a; ++l) {
                    computeFacePoint(input->GetCell(neighborCells->GetId(l)), fp);
                    sum[0] += fp[0];
                    sum[1] += fp[1];
                    sum[2] += fp[2];
                }
                delete fp;

                // add weighted sum of involved edge-points
                for (int l = 0; l < n; ++l) {
                    sum[0] += input->GetPoint(neighbors->GetId(l))[0];
                    sum[1] += input->GetPoint(neighbors->GetId(l))[1];
                    sum[2] += input->GetPoint(neighbors->GetId(l))[2];
                }

                // add weighted vertex
                coord[0] = (input->GetPoint(id))[0] * weight_v + sum[0]*weight_e;
                coord[1] = (input->GetPoint(id))[1] * weight_v + sum[1]*weight_e;
                coord[2] = (input->GetPoint(id))[2] * weight_v + sum[2]*weight_e;

/*
                coord[0] = (input->GetPoint(id))[0];
                coord[1] = (input->GetPoint(id))[1];
                coord[2] = (input->GetPoint(id))[2];
*/

                VertexPoint *newVertexPoint = new VertexPoint;
                newVertexPoint->v_old = id;
                newVertexPoint->v = outputPoints->InsertNextPoint(coord);
                insertVertexPoint(cells[i], newVertexPoint);
            } else {
                VertexPoint *newVertexPoint = new VertexPoint;
                newVertexPoint->v_old = id;
                newVertexPoint->v = vertexPoint->v;
                insertVertexPoint(cells[i], newVertexPoint);
            }
        }

        //=============================================
        // Create new Cell
        //=============================================
        for (int j = 0; j < cpoints; ++j) {
            vtkIdList *newPointIds = vtkIdList::New();
            vtkIdType id1 = cell->GetPointId((vtkIdType) j);
            vtkIdType id2 = cell->GetPointId((vtkIdType) ((j+1) % cpoints));
            vtkIdType id3 = cell->GetPointId((vtkIdType) ((j+2) % cpoints));
            vtkIdType e1  = (findEdgePoint(cells[i], id1, id2))->e;
            vtkIdType e2  = (findEdgePoint(cells[i], id2, id3))->e;
            vtkIdType v   = (findVertexPoint(cells[i], id2))->v;
            vtkIdType f   = cells[i].facePoint;
            newPointIds->InsertNextId(v);
            newPointIds->InsertNextId(e1);
            newPointIds->InsertNextId(f);
            newPointIds->InsertNextId(e2);
            outputPolys->InsertNextCell(newPointIds);
        }

        EdgePoint *ep = cells[i].edgePoints;
        while (ep != NULL) {
            ep = ep->next;
        }

        VertexPoint *vp = cells[i].vertexPoints;
        while (vp != NULL) {
            vp = vp->next;
        }
    }

    // We now assign the points and polys to the output-vtkPolyData.
    output->SetPoints(outputPoints);
    output->SetPolys(outputPolys);

    // free memory 
    // ToDo: remove linked lists and other stuff
    delete [] cells;
    outputPoints->Delete();
    outputPolys->Delete();
}

void vtkCatmullClarkFilter::Execute(void) {

    vtkPolyData *input;
    vtkPolyData *output; 
    vtkPolyData *intermediate;

    int subdivisions = GetNumberOfSubdivisions();

    if (subdivisions == 0) {
        // No subdivision
    } else if (subdivisions == 1) {
        input = GetPolyDataInput(0);//GetInput();
        output = GetOutput();
        Subdivision(input, output);
    } else for(int i = 0; i < subdivisions; ++i) {
        if (i == 0) {
            input = GetPolyDataInput(0); //GetInput();
            intermediate = vtkPolyData::New();
            Subdivision(input, intermediate);
        } else if (i < subdivisions-1) {
            input->Delete();
            input = intermediate;
            intermediate = vtkPolyData::New();
            Subdivision(input, intermediate);
        } else {
            input->Delete();
            input = intermediate;
            output = GetOutput();
            Subdivision(intermediate, output);
        }
    }
}


int vtkCatmullClarkFilter::GenerateSubdivisionPoints (vtkPolyData *inputDS,vtkIntArray *edgeData, vtkPoints *outputPts, vtkPointData *outputPD)
{
    return 1;// Not implemented
}
