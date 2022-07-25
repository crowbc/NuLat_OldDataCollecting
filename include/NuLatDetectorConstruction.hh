//DETECTOR CONSTRUCTION HEADER FILE

#ifndef NuLatDetectorConstruction_H
#define NuLatDetectorConstruction_H 1

#include "globals.hh"
#include <vector>

#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
//Header file to create many copied volumes
#include "G4PVParameterised.hh"
//Geant4 header files for making shapes
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4IntersectionSolid.hh"
//Geant4 header file for cutting chunks out of a shape
#include "G4SubtractionSolid.hh"
//Header file for creating region that collects data
#include "G4VSensitiveDetector.hh"
//Header files for setting attributes in visualization
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
//New Material Properties (added 2/7/22)
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "Materials.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class G4GenericMessenger;

/// Detector construction

class NuLatDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    NuLatDetectorConstruction(G4int xVoxels, G4int yVoxels, G4int zVoxels,
                              G4double vXDimension, G4double vYDimension, G4double vZDimension, 
                              G4double vSXDimension, G4double vSYDimension, G4double vSZDimension)
    {
        nOfVoxelsInX = xVoxels;
        nOfVoxelsInY = yVoxels;
        nOfVoxelsInZ = zVoxels;
        
        voxelXDimension = vXDimension;
        voxelYDimension = vYDimension;
        voxelZDimension = vZDimension;
        
        voxelSpacingXDimension = vSXDimension;
        voxelSpacingYDimension = vSYDimension;
        voxelSpacingZDimension = vSZDimension;
 // TODO pass these through parameter list
        PMTLength = 20*cm;
        LightGuideTaperLength = 3.5*cm;
        LightGuideSquareLength = 0.5*cm;
        lightGuideWithPMTLength = PMTLength+LightGuideTaperLength+LightGuideSquareLength;
    }
    ~NuLatDetectorConstruction() {}
    G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

  private:
    void buildExperimentalHall();

    /********************************/
    /* Logical Volume Declerations  */
    /********************************/
    G4LogicalVolume* WorldVolume();
    G4LogicalVolume* experimentalHallLog;

    G4LogicalVolume* HamamatsuR10533();
    G4LogicalVolume* LightGuideAndPMT(G4double dx1, G4double dx2, G4double dy1, G4double dy2, G4double dz, G4bool instrumented);
    G4LogicalVolume* LightGuideAndPMTLog;
    G4LogicalVolume* LightGuideAndPMTLog_noPMTs;
    G4LogicalVolume* NuLatLightGuideZBankPlusLogical;
    G4LogicalVolume* NuLatLightGuideZBankMinusLogical;
    G4LogicalVolume* NuLatLightGuideYBankPlusLogical;
    G4LogicalVolume* NuLatLightGuideYBankMinusLogical;
    G4LogicalVolume* NuLatLightGuideXBankPlusLogical;
    G4LogicalVolume* NuLatLightGuideXBankMinusLogical;
    G4LogicalVolume* muMetal_Log;
    G4LogicalVolume* pmtLens_Log;
    G4LogicalVolume* pmtPhotoCathode_Log;
    G4LogicalVolume* lightGuideTrdIntersLightGuideConeLog;
    G4LogicalVolume* lightGuideSquareLog;
    G4LogicalVolume* PMTLog;
    G4LogicalVolume* dividers_Log;
    G4LogicalVolume* dividers_inner_Log;
    G4LogicalVolume* VoxelLogical;
    
    // Parameter variable declarations
    G4int nOfVoxelsInX, nOfVoxelsInY, nOfVoxelsInZ, photoCathodeIDNum;
    
    // Dimension variable declarations
    G4double voxelXDimension, voxelYDimension, voxelZDimension;    
    G4double voxelSpacingXDimension, voxelSpacingYDimension, voxelSpacingZDimension;
    G4double  PMTLength, LightGuideTaperLength, LightGuideSquareLength, lightGuideWithPMTLength;

    /********************************/
    /* Physical Volume Declerations */
    /********************************/
    G4VPhysicalVolume* physWorld;


    /********************************/
    /* Create an instance of the    */
    /* Materials class              */
    /********************************/
    Materials* NuLatMaterials;

    /********************************/
    /*        Pointer to            */
    /*  Visualization Attributes    */
    /********************************/
    std::vector<G4VisAttributes*> NuLatVisAttributes;
};

#endif







