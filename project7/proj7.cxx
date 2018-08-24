/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include "tricase.cxx"
#include "TriangleList.h"
// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    return dims[0]*dims[1]*dims[2];
    // 2D
   // return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
     idx[0] = pointId%dims[0];
     idx[1] = (pointId/dims[0])%dims[1];
     idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    //idx[0] = pointId%dims[0];
    //idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
     idx[0] = cellId%(dims[0]-1);
     idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
     idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
}


class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

float interpolate(const float a, const float b, const float F_A, const float F_B, const float pt){
        float t = (pt-a)/(b-a);
        return F_A + (t * (F_B - F_A));
}


float
EvaluateFieldAtLocation(const float *pt, const int *dims,
                        const float *X, const float *Y, const float *F)
{
    float bbox[4];
    int pidx[2];
    int ptr[2], j;
    float v1, v2, v3, v4, a, b;
    for(int i = 0; i < GetNumberOfCells(dims); ++i){
        GetLogicalCellIndex(ptr, i, dims);
        bbox[0] = X[ptr[0]];
        bbox[1] = X[ptr[0]+1];
        bbox[2] = Y[ptr[1]] ;
        bbox[3] = Y[ptr[1]+1];
        if(pt[0] >= bbox[0] && pt[0] <= bbox[1]){
            if(pt[1] >= bbox[2] && pt[1] <= bbox[3]){
                pidx[0] = ptr[0];
                pidx[1] = ptr[1];
                v1 = F[GetPointIndex(pidx, dims)];
                pidx[0] = ptr[0]+1;
                pidx[1] = ptr[1];
                v2 = F[GetPointIndex(pidx, dims)];
                pidx[0] = ptr[0];
                pidx[1] = ptr[1]+1;
                v3 = F[GetPointIndex(pidx, dims)];
                pidx[0] = ptr[0]+1;
                pidx[1] = ptr[1]+1;
                v4 = F[GetPointIndex(pidx, dims)];
                a = interpolate(bbox[0], bbox[1], v1, v2, pt[0]);
                b = interpolate(bbox[0], bbox[1], v3, v4, pt[0]);
                return interpolate(bbox[2], bbox[3], a , b, pt[1]);
                }
        }
    }
    return 0;
}


vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

void setEdgeIndex(int edge[4][2], int ptr[4][2], const int *dims){
	edge[0][0] = GetPointIndex(ptr[0], dims);
	edge[0][1] = GetPointIndex(ptr[1], dims);
	edge[1][0] = GetPointIndex(ptr[1], dims);
	edge[1][1] = GetPointIndex(ptr[3], dims);
	edge[2][0] = GetPointIndex(ptr[2], dims);
	edge[2][1] = GetPointIndex(ptr[3], dims);
	edge[3][0] = GetPointIndex(ptr[0], dims);
	edge[3][1] = GetPointIndex(ptr[2], dims);

}

void getEdges(int *ptr, int edge){
	int count;
	if(edge == 0){
	    ptr[0] = 0;
	    ptr[1] = 1;
	}
	if(edge == 1){
	    ptr[0] = 1;
	    ptr[1] = 3;
	}
	if(edge == 2){
	    ptr[0] = 2;
	    ptr[1] = 3;
	}
	if(edge == 3){
	    ptr[0] = 0;
	    ptr[1] = 2;
	}
}	

void createTable(int table[][4])
{
	//case 0
	table[0][0] = table[0][1] = table[0][2] = table[0][3] = -1;
	//case 1
	table[1][0] = 0; table[1][1] = 3; table[1][2] = table[1][3] = -1; 
	//case 2
	table[2][0] = 0; table[2][1] = 1; table[2][2] = table[2][3] = -1;
	//case 3
	table[3][0] = 1; table[3][1] = 3; table[3][2] = table[3][3] = -1;
	//case 4
	table[4][0] = 2; table[4][1] = 3; table[4][2] = table[4][3] = -1;
	//case 5
	table[5][0] = 0; table[5][1] = 2; table[5][2] = table[5][3] = -1;
	//case 6
	table[6][0] = 0; table[6][1] = 1; table[6][2] = 2; table[6][3] = 3;
	//case 7
	table[7][0] = 1; table[7][1] = 2; table[7][2] = table[7][3] = -1;
	//case 8
	table[8][0] = 1; table[8][1] = 2; table[8][2] = table[8][3] = -1;
	//case 9
	table[9][0] = 0; table[9][1] = 3; table[9][2] = 1; table[9][3] = 2;
	//case 10
	table[10][0] = 0; table[10][1] = 2; table[10][2] = table[10][3] = -1;
	//case 11
	table[11][0] = 2; table[11][1] = 3; table[11][2] = table[11][3] = -1;
	//case 12
	table[12][0] = 1; table[12][1] = 3; table[12][2] = table[12][3] = -1;
	//case 13
	table[13][0] = 0; table[13][1] = 1; table[13][2] = table[13][3] = -1;
	//case 14
	table[14][0] = 0; table[14][1] = 3; table[14][2] = table[14][3] = -1;
	//case 15
	table[15][0] = table[15][1] =  table[15][2] = table[2][3] = -1;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj6B.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
     
    // Add 4 segments that put a frame around your isolines.  This also
    // documents how to use "AddSegment".
    TriangleList sl;

// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!
    int icase;
    int vertMap[8] = {2,3,0,1,6,7,4,5}; 
    vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
    static int edges[12][2] = { {0,1}, {1,3}, {2,3}, {0,2}, {4,5}, {5,7},
				{6,7}, {4,6}, {0,4}, {1,5}, {2,6}, {3,7} };
    int tri[3];
    for(i = 0; i < GetNumberOfCells(dims); i++){
	rgrid->GetCellPoints(i, ptids);
	icase = 0;
	for(j = 0; j < 8; j++){
	    if(F[ptids->GetId(vertMap[j])] > 3.2){
		icase += pow(2,j);	
	    }
	    
	}
	float t;
	float pt1[3];
	float pt2[3];
	float pt3[3];
	int idx[3], idx2[3];
	int edge1[2], edge2[2], edge3[2];	
	float low_pt;
	float high_pt;
	for(j = 0; j < 16; j+=3){
	    if(triCase[icase][j] != -1){
		edge1[0] = edges[triCase[icase][j]][0];
		edge1[1] = edges[triCase[icase][j]][1];
		GetLogicalPointIndex(idx, ptids->GetId(vertMap[edge1[0]]), dims);
		GetLogicalPointIndex(idx2, ptids->GetId(vertMap[edge1[1]]), dims);
		low_pt = F[ptids->GetId(vertMap[edges[triCase[icase][j]][0]])];
		high_pt = F[ptids->GetId(vertMap[edges[triCase[icase][j]][1]])];
		t = (3.2 - low_pt) / (high_pt - low_pt);

		pt1[0] = X[idx[0]] + (t * (X[idx2[0]] - X[idx[0]]));	
		pt1[1] = Y[idx[1]] + (t * (Y[idx2[1]] - Y[idx[1]]));	
		pt1[2] = Z[idx[2]] + (t * (Z[idx2[2]] - Z[idx[2]]));	
			
			
		edge2[0] = edges[triCase[icase][j+1]][0];
		edge2[1] = edges[triCase[icase][j+1]][1];
		GetLogicalPointIndex(idx, ptids->GetId(vertMap[edge2[0]]), dims);
		GetLogicalPointIndex(idx2, ptids->GetId(vertMap[edge2[1]]), dims);
		low_pt = F[ptids->GetId(vertMap[edges[triCase[icase][j+1]][0]])];
		high_pt = F[ptids->GetId(vertMap[edges[triCase[icase][j+1]][1]])];
		t = (3.2 - low_pt) / (high_pt - low_pt);
		pt2[0] = X[idx[0]] + (t * (X[idx2[0]] - X[idx[0]]));	
		pt2[1] = Y[idx[1]] + (t * (Y[idx2[1]] - Y[idx[1]]));	
		pt2[2] = Z[idx[2]] + (t * (Z[idx2[2]] - Z[idx[2]]));	
		
		edge3[0] = edges[triCase[icase][j+2]][0];
		edge3[1] = edges[triCase[icase][j+2]][1];
		GetLogicalPointIndex(idx, ptids->GetId(vertMap[edge3[0]]), dims);
		GetLogicalPointIndex(idx2, ptids->GetId(vertMap[edge3[1]]), dims);
		low_pt = F[ptids->GetId(vertMap[edges[triCase[icase][j+2]][0]])];
		high_pt = F[ptids->GetId(vertMap[edges[triCase[icase][j+2]][1]])];
		t = (3.2 - low_pt) / (high_pt - low_pt);
		pt3[0] = X[idx[0]] + (t * (X[idx2[0]] - X[idx[0]]));	
		pt3[1] = Y[idx[1]] + (t * (Y[idx2[1]] - Y[idx[1]]));	
		pt3[2] = Z[idx[2]] + (t * (Z[idx2[2]] - Z[idx[2]]));	
		sl.AddTriangle(pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2], pt3[0], pt3[1], pt3[2]);
	}
}
}

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */
    vtkPolyData *pd = sl.MakePolyData();
    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
