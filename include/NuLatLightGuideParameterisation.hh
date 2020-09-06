//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: NuLatLightGuideParameterisation.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file NuLatLightGuideParameterisation.hh
/// \brief Definition of the NuLatLightGuideParameterisation class

#ifndef NuLatLightGuideParameterisation_H
#define NuLatLightGuideParameterisation_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;

/// EM Calorimeter cell parameterisation

class NuLatLightGuideParameterisation : public G4VPVParameterisation
{
public:
    NuLatLightGuideParameterisation(
          G4int nOfVoxelsInX, G4int nOfVoxelsInY, G4int nOfVoxelsInZ,
          G4double voxelXDimension, G4double voxelYDimension, G4double voxelZDimension,
          G4double voxelSpacingXDimension,
          G4double voxelSpacingYDimension,
          G4double voxelSpacingZDimension,
          G4int xselector,
          G4int yselector,
          G4int zselector
          );
    virtual ~NuLatLightGuideParameterisation();
    
    virtual void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;
    
private:
    G4double xLightGuide[15*15];
    G4double yLightGuide[15*15];
    G4double zLightGuide[15*15];
    G4RotationMatrix* rotm[15*15];  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
