#ifndef SURFACERECONSTRUCTION_H
#define SURFACERECONSTRUCTION_H

/*
    Surface Reconstruction is the overall wrapper for the
    3D Reconstruction Methods & Clipping and Navigation
    of the 3D Vascular Navigation.
*/

#include "../Core/tubNavContainer.h"
#include "../MathCore/CoordinateFrame.h"

#include "SurfaceSubdivision.h"
#include <vtkPolyDataAlgorithm.h>

struct Clipping {
   tubNav::Point<double> location;
   tubNav::Point<double> normal;
   bool clipDescendants;
};


class SurfaceReconstruction
{
   public:
      enum class ReconstructionMethod { Subdivision , CatmullClark , Custom};
      enum class SmoothingFilter {LinearSubdivision, LoopSubdivision,  ButterflySubdivision, LaplacianSubdivision, CatmullClark };

      SurfaceReconstruction();
      SurfaceReconstruction(std::shared_ptr<tubNav::Container<double> > tree);

      void AddNewClipping(tubNav::Point<double> location, tubNav::Point<double> n, bool clipDescendants = false);

      void RemoveClipping(int idx);
      Clipping GetClipping(int idx);

      void SetVascularTree(std::shared_ptr<tubNav::Container<double> > tree);
      void ChangeReconstructionMethod(ReconstructionMethod newMethod);
      void SetReconstruction(vtkSmartPointer<vtkPolyData> customReconstruction);


      void SetSmoothingFilter(SmoothingFilter filter, int numberOfSubdivisions);

      void DisableDecimation(){ enableDecimation = false;}
      void EnableDecimation(float Percentage =0.5){ enableDecimation = true; decimationPercentage = Percentage;}


      void ResetReconstruction();

      std::shared_ptr<tubNav::CoordinateFrame<double>> GetCurrentFrame();

      vtkSmartPointer<vtkPolyData> GetReconstruction();

      void MoveToNextPoint();
      void MoveToPreviousPoint();
      void MoveToStart();
      bool IsAtEnd();
   private:

      void ApplyClipping(int idx);
      void ApplyAllClippings();
      tubNav::ContainerIterator<double> it;

      ReconstructionMethod method;

      std::shared_ptr<tubNav::Container<double>> vascTree;
      std::shared_ptr<tubNav::CoordinateFrame<double>> frame;

      vtkSmartPointer<vtkPolyData> currentReconstruction;
      vtkSmartPointer<vtkPolyData> _customReconstruction;

      std::vector<Clipping> clippings;

      float decimationPercentage;
      bool enableDecimation;
      vtkSmartPointer<vtkPolyDataAlgorithm> smoothingFilterToApply;

};



#endif // SURFACERECONSTRUCTION_H
