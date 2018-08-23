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
#include <vtkCubeSource.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkProperty.h>
#include <vtkHedgeHog.h>
#include <vtkRungeKutta4.h>
#include <vtkStreamTracer.h>
#include <vtkLineSource.h>


int main()
{
  vtkDataSetReader *rdr = vtkDataSetReader::New();
  rdr->SetFileName("proj8.vtk");
  rdr->Update();
 
  vtkContourFilter *cf = vtkContourFilter::New();
  cf->SetValue(0,2.5);
  cf->SetInputConnection(rdr->GetOutputPort());
  cf->Update();

  vtkContourFilter *cf2 = vtkContourFilter::New();
  cf2->SetValue(0, 5.0);
  cf2->SetInputConnection(rdr->GetOutputPort());
  cf2->Update();

  // Define viewport ranges
  double xmins[4] = {0,0,0.5,0.5};
  double xmaxs[4] = {0.5,0.5,1,1};
  double ymins[4] = {0,0.5,0,0.5};
  double ymaxs[4]= {0.5,1,0.5,1};
 
    vtkSmartPointer<vtkRenderWindow> renderWindow = 
     vtkSmartPointer<vtkRenderWindow>::New();
    //render 1 
    vtkSmartPointer<vtkRenderer> renderer1 = 
      vtkSmartPointer<vtkRenderer>::New();
 
    renderWindow->AddRenderer(renderer1);
    renderer1->SetViewport(xmins[0],ymins[0],xmaxs[0],ymaxs[0]);
    
    //render 2
    vtkSmartPointer<vtkRenderer> renderer2 =
      vtkSmartPointer<vtkRenderer>::New();

    renderWindow->AddRenderer(renderer2);
    renderer2->SetViewport(xmins[1],ymins[1],xmaxs[1],ymaxs[1]);
    
// renderer 3
    vtkSmartPointer<vtkRenderer> renderer3 =
      vtkSmartPointer<vtkRenderer>::New();

    renderWindow->AddRenderer(renderer3);
    renderer3->SetViewport(xmins[2],ymins[2],xmaxs[2],ymaxs[2]);

// renderer 4
    vtkSmartPointer<vtkRenderer> renderer4 =
      vtkSmartPointer<vtkRenderer>::New();

    renderWindow->AddRenderer(renderer4);
    renderer4->SetViewport(xmins[3],ymins[3],xmaxs[3],ymaxs[3]);

    // Create a mapper and actor

    vtkSmartPointer<vtkPolyDataMapper> win1Mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    win1Mapper->SetInputConnection(cf->GetOutputPort());
    
    vtkSmartPointer<vtkActor> win1Actor = vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);
    win1Actor->GetProperty()->SetColor(0.0,1.0,0.0);
    win1Actor->GetMapper()->ScalarVisibilityOff();
    
    vtkSmartPointer<vtkPolyDataMapper> win1Mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    win1Mapper2->SetInputConnection(cf2->GetOutputPort());
 
    vtkSmartPointer<vtkActor> win1Actor2 = vtkSmartPointer<vtkActor>::New(); 
    win1Actor2->SetMapper(win1Mapper2);
    win1Actor2->GetProperty()->SetColor(1.0,0.0,0.0);
    win1Actor2->GetMapper()->ScalarVisibilityOff();
       

    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renderWindow);
    
    renderer1->AddActor(win1Actor);
    renderer1->AddActor(win1Actor2);
    renderer1->SetBackground(0.0, 0.0, 0.0);
    renderWindow->SetSize(600, 600);


    renderer1->GetActiveCamera()->SetFocalPoint(0,0,0);
    renderer1->GetActiveCamera()->SetPosition(0,0,70);
    renderer1->GetActiveCamera()->SetViewUp(0,1,0);
    renderer1->GetActiveCamera()->SetClippingRange(20, 120);
    renderer1->GetActiveCamera()->SetDistance(70);

    vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
    reader->SetFileName("proj8.vtk");
    reader->Update();
    vtkDataSet *hardyglobal = reader->GetOutput();
    hardyglobal->GetPointData()->SetActiveScalars("hardyglobal");
    vtkDataSet *grad = reader->GetOutput();
    grad->GetPointData()->SetActiveAttribute("grad", vtkDataSetAttributes::VECTORS);
//
	
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin(0,0,0);
    plane->SetNormal(0,1,0);


//create cutter
    vtkSmartPointer<vtkCutter> cutter =  vtkSmartPointer<vtkCutter>::New();
    cutter->SetCutFunction(plane);
    cutter->SetInputData(hardyglobal);
    cutter->Update();
    vtkSmartPointer<vtkPolyDataMapper> cutterMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cutterMapper->SetInputConnection(cutter->GetOutputPort());
    cutterMapper->SetScalarRange(hardyglobal->GetScalarRange());

// plane actor
    vtkSmartPointer<vtkActor> planeActor = vtkSmartPointer<vtkActor>::New();
    planeActor->GetProperty()->SetColor(0.5, 1, 0.5);
    planeActor->GetProperty()->SetLineWidth(2);
    planeActor->SetMapper(cutterMapper);    
// plane 2
    vtkSmartPointer<vtkPlane> plane2 = vtkSmartPointer<vtkPlane>::New();
    plane2->SetOrigin(0,0,0);
    plane2->SetNormal(1,0,0);


//create cutter 2
    vtkSmartPointer<vtkCutter> cutter2 =  vtkSmartPointer<vtkCutter>::New();
    cutter2->SetCutFunction(plane2);
    cutter2->SetInputData(hardyglobal);
    cutter2->Update();
    vtkSmartPointer<vtkPolyDataMapper> cutterMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    cutterMapper2->SetInputConnection(cutter2->GetOutputPort());
    cutterMapper2->SetScalarRange(hardyglobal->GetScalarRange());
// plane actor 2
    vtkSmartPointer<vtkActor> planeActor2 = vtkSmartPointer<vtkActor>::New();
    planeActor2->GetProperty()->SetColor(0.5, 1, 0.5);
    planeActor2->GetProperty()->SetLineWidth(2);
    planeActor2->SetMapper(cutterMapper2);
// plane 3
    vtkSmartPointer<vtkPlane> plane3 = vtkSmartPointer<vtkPlane>::New();
    plane3->SetOrigin(0,0,0);
    plane3->SetNormal(0,0,1);


//create cutter 3
    vtkSmartPointer<vtkCutter> cutter3 =  vtkSmartPointer<vtkCutter>::New();
    cutter3->SetCutFunction(plane3);
    cutter3->SetInputData(hardyglobal);
    cutter3->Update();
    vtkSmartPointer<vtkPolyDataMapper> cutterMapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
    cutterMapper3->SetInputConnection(cutter3->GetOutputPort());
    cutterMapper3->SetScalarRange(hardyglobal->GetScalarRange());
// plane actor 3
    vtkSmartPointer<vtkActor> planeActor3 = vtkSmartPointer<vtkActor>::New();
    planeActor3->GetProperty()->SetColor(0.5, 1, 0.5);
    planeActor3->GetProperty()->SetLineWidth(2);
    planeActor3->SetMapper(cutterMapper3);

    
    renderer2->AddActor(planeActor);
    renderer2->AddActor(planeActor2);
    renderer2->AddActor(planeActor3);

    vtkSmartPointer<vtkHedgeHog> hedgehog = vtkSmartPointer<vtkHedgeHog>::New();
    hedgehog->SetInputData(grad);
    //hedgehog->SetScalarFactor(.001);
    
    vtkSmartPointer<vtkPolyDataMapper> hmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    hmapper->SetInputConnection(hedgehog->GetOutputPort());
    hmapper->SetScalarRange(grad->GetScalarRange());
    
    vtkSmartPointer<vtkActor> hactor = vtkSmartPointer<vtkActor>::New();
    hactor->SetMapper(hmapper);
    
    renderer3->AddActor(hactor);

//line source instead of a point source
   
    vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
    lineSource->SetPoint1(-9,0,0);
    lineSource->SetPoint2(9,0,0);
    lineSource->SetResolution(19);

    vtkSmartPointer<vtkRungeKutta4> integ = vtkSmartPointer<vtkRungeKutta4>::New();
    
    vtkSmartPointer<vtkStreamTracer> streamer = vtkSmartPointer<vtkStreamTracer>::New();
    streamer->SetInputConnection(reader->GetOutputPort());
    streamer->SetSourceConnection(lineSource->GetOutputPort());
    streamer->SetMaximumPropagation(100);
    streamer->SetIntegrator(integ);

    vtkSmartPointer<vtkPolyDataMapper> mapStreamLines = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapStreamLines->SetInputConnection(streamer->GetOutputPort());
    mapStreamLines->SetScalarRange(grad->GetScalarRange());

    vtkSmartPointer<vtkActor> streamActor = vtkSmartPointer<vtkActor>::New();
    streamActor->SetMapper(mapStreamLines);

    renderer4->AddActor(streamActor);     



  // This starts the event loop and invokes an initial render.
  //
  ((vtkInteractorStyle *)iren->GetInteractorStyle())->SetAutoAdjustCameraClippingRange(0);
  iren->Initialize();
  iren->Start();
  
}	
