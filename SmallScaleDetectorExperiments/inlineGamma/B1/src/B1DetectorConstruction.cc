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
// $Id$
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume1(0),
  fScoringVolume2(0),
  fScoringVolume3(0),
  fScoringVolume4(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{

  //std::cout << "test\n";

  G4double a, z, density;
  G4int nelements;

  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  G4Element* Gd = new G4Element("Gadolinium", "Gd", z = 64, a = 157.2*g/mole);
  G4Element* Ga = new G4Element("Gallium", "Ga", z = 31, a = 69.72*g/mole);
  G4Element* Al = new G4Element("Aluminum", "Al", z = 13, a = 26.98*g/mole);
  G4Material* GAGGCe = new G4Material("GAGGCe",density=6.63*g/cm3,nelements=4);
  GAGGCe -> AddElement(Gd, 3);
  GAGGCe -> AddElement(Al, 2);
  GAGGCe -> AddElement(Ga, 3);
  GAGGCe -> AddElement(O, 12);
  //std::cout << "test\n";
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  //G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  //G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 40.;
  G4double world_sizeZ  = 100.;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //
  /*
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  */
  /*
  //     
  // Shape 1
  //  
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);
        
  // Conical section shape       
  G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  G4double shape1_hz = 3.*cm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  G4Cons* solidShape1 =    
    new G4Cons("Shape1", 
    shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
    shape1_phimin, shape1_phimax);
                      
  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                        shape1_mat,          //its material
                        "Shape1");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //     
  // Shape 2
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4ThreeVector pos2 = G4ThreeVector(0, -1*cm, 7*cm);

  // Trapezoid shape       
  G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  G4double shape2_dz  = 6*cm;      
  G4Trd* solidShape2 =    
    new G4Trd("Shape2",                      //its name
              0.5*shape2_dxa, 0.5*shape2_dxb, 
              0.5*shape2_dya, 0.5*shape2_dyb, 0.5*shape2_dz); //its size
                
  G4LogicalVolume* logicShape2 =                         
    new G4LogicalVolume(solidShape2,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos2,                    //at position
                    logicShape2,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  */

  G4double voxelDim = 7.; //mm
  //G4double voxelThicc = 1.; //mm
  double pos = 0, space = 2.;
    
  //1st voxel
  G4Box* voxel1Box = new G4Box("voxel1Box", voxelDim/2.,
				  voxelDim/2., voxelDim/2.);
  //std::cout << "test\n";
  G4LogicalVolume* voxel1Log = new G4LogicalVolume(voxel1Box,
						      GAGGCe, "voxel1Log");
  //std::cout << "test\n";
  new G4PVPlacement(0, G4ThreeVector(0,0,pos), voxel1Log, "voxel1Place",
		    logicWorld, false, 0);
  pos += voxelDim + space;

  //2nd voxel
  G4Box* voxel2Box = new G4Box("voxel2Box", voxelDim/2.,
				  voxelDim/2., voxelDim/2.);
  //std::cout << "test\n";
  G4LogicalVolume* voxel2Log = new G4LogicalVolume(voxel2Box,
						      GAGGCe, "voxel2Log");
  //std::cout << "test\n";
  new G4PVPlacement(0, G4ThreeVector(0,0,pos), voxel2Log, "voxel2Place",
		    logicWorld, false, 0);
  pos += voxelDim + space;

  //3rd voxel
  G4Box* voxel3Box = new G4Box("voxel3Box", voxelDim/2.,
				  voxelDim/2., voxelDim/2.);
  //std::cout << "test\n";
  G4LogicalVolume* voxel3Log = new G4LogicalVolume(voxel3Box,
						      GAGGCe, "voxel3Log");
  //std::cout << "test\n";
  new G4PVPlacement(0, G4ThreeVector(0,0,pos), voxel3Log, "voxel3Place",
		    logicWorld, false, 0);
  pos += voxelDim + space;

  //4th voxel
  G4Box* voxel4Box = new G4Box("voxel4Box", voxelDim/2.,
				  voxelDim/2., voxelDim/2.);
  //std::cout << "test\n";
  G4LogicalVolume* voxel4Log = new G4LogicalVolume(voxel4Box,
						      GAGGCe, "voxel4Log");
  //std::cout << "test\n";
  new G4PVPlacement(0, G4ThreeVector(0,0,pos), voxel4Log, "voxel4Place",
		    logicWorld, false, 0);
  

  G4Tubs* plasticDiscTub = new G4Tubs("plasticDiscTub", 0., 25.4/2.,
					       3./2., 0., M_PI*2.);

  G4Material* plasticDiscMat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4LogicalVolume* plasticDiscLog = new G4LogicalVolume(plasticDiscTub,
						       plasticDiscMat,
						       "plasticDiscLog");
  new G4PVPlacement(0, G4ThreeVector(0,0,-8.5), plasticDiscLog,
		    "plasticDiscPlace", logicWorld, false, 0);
  
  

  /*
  G4Material* brass = nist->FindOrBuildMaterial("G4_BRASS");
  G4LogicalVolume* bigBrassLog = new G4LogicalVolume(bigBrassTubs, brass,
						     "bigBrassLog");
  new G4PVPlacement(0, G4ThreeVector(0,0,-(18. + 25.4/2.)), bigBrassLog,
		    "bigBrassPlace", logicWorld, false, 0);

  G4Tubs* smallBrassTubs = new G4Tubs("smallBrassTubs", 4./2., 19./2.,
				    1./2., 0., 2*M_PI);
  
  G4Material* cupper = nist->FindOrBuildMaterial("G4_Cu");
  G4LogicalVolume* smallBrassLog = new G4LogicalVolume(smallBrassTubs, cupper,
						     "smallBrassLog");
  new G4PVPlacement(0, G4ThreeVector(0,0,-(18. - .5)), smallBrassLog,
		    "smallBrassPlace", logicWorld, false, 0);
  */

  
                
  // Set Shape2 as scoring volume
  //
  fScoringVolume1 = voxel1Log;
  fScoringVolume2 = voxel2Log;
  fScoringVolume3 = voxel3Log;
  fScoringVolume4 = voxel4Log;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
