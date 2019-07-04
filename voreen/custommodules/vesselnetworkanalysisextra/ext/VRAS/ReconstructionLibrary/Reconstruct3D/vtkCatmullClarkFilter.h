/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCatmullClarkFilter.h,v $
  Language:  C++
  Date:      $Date: 2001/10/11 13:37:08 $
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
// .NAME vtkCatmullClarkFilter - generate a subdivision surface using the Catmull-Clark Scheme
// .SECTION Description
// TODO

#ifndef __vtkCatmullClarkFilter_h
#define __vtkCatmullClarkFilter_h

#include <vtkApproximatingSubdivisionFilter.h>
#include <vtkIntArray.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

struct EdgePoint {
    vtkIdType e;
    vtkIdType p1;
    vtkIdType p2;
    EdgePoint *next;
};

struct VertexPoint {
    vtkIdType v;
    vtkIdType v_old;
    VertexPoint *next;
};

struct CellData {
    vtkIdType    facePoint;
    EdgePoint   *edgePoints;
    VertexPoint *vertexPoints;
};

class  vtkCatmullClarkFilter : public vtkApproximatingSubdivisionFilter
{
public:
  // Description:
  // Construct object with NumberOfSubdivisions set to 1.
  static vtkCatmullClarkFilter *New();
  vtkTypeMacro(vtkCatmullClarkFilter,vtkApproximatingSubdivisionFilter);

protected:
  vtkCatmullClarkFilter () {};
  ~vtkCatmullClarkFilter () {};

  virtual void Execute(void);

  void insertEdgePoint(CellData &cell, EdgePoint *newElement);
  EdgePoint *findEdgePoint(CellData &cell, vtkIdType id1, vtkIdType id2);
  void insertVertexPoint(CellData &cell, VertexPoint *newElement);
  VertexPoint *findVertexPoint(CellData &cell, vtkIdType id);

  void computeFacePoint(vtkCell *cell, float *facePoint);
  vtkIdList *getNeighborIds(vtkPolyData *input, vtkIdType i);

  void Subdivision(vtkPolyData* input, vtkPolyData *output);

  //0 Fail, otherwise 1 // Not implement
  int GenerateSubdivisionPoints (vtkPolyData *inputDS, vtkIntArray *edgeData,
                                  vtkPoints *outputPts,
                                  vtkPointData *outputPD);

private:
  vtkCatmullClarkFilter(const vtkCatmullClarkFilter&);  // Not implemented.
  void operator=(const vtkCatmullClarkFilter&);  // Not implemented.
};

#endif



