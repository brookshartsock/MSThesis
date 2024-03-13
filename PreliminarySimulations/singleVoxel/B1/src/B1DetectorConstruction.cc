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
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  G4double a, z, density;
  G4int nelements;
  G4double ncomponents;

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //here is where GAGG is defined
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  G4Element* Gd = new G4Element("Gadolinium", "Gd", z = 64, a = 157.2*g/mole);
  G4Element* Ga = new G4Element("Gallium", "Ga", z = 31, a = 69.72*g/mole);
  G4Element* Al = new G4Element("Aluminum", "Al", z = 13, a = 26.98*g/mole);
  G4Material* GAGGCe = new G4Material("GAGGCe",density=6.63*g/cm3,nelements=4);
  GAGGCe -> AddElement(Gd, 3);
  GAGGCe -> AddElement(Al, 2);
  GAGGCe -> AddElement(Ga, 3);
  GAGGCe -> AddElement(O, 12);

  //here is where Indium (material) is defined
  G4Element* In = new G4Element("Indium"  , "In", z=49 , a=114.818*g/mole);
  G4Material* InFoil = new G4Material("InFoil", density=7.31*g/cm3, nelements=1);
  InFoil -> AddElement(In, 1);

  //here is the Novak Sheild with 30% of tungsten/epoxy
  G4Material* tungstenMat = nist -> FindOrBuildMaterial("G4_W"); //19.3 g/cm3
  G4Material* epoxy = nist -> FindOrBuildMaterial("G4_POLYSTYRENE"); //1.06 g/cm3
  double tungstenProportion = .3; //ratio of tungsten to total mass
  G4double NovakShieldDensity = 1./((tungstenProportion/19.3)
				  + ((1.-tungstenProportion)/1.06)) * g/cm3;
  //std::cout << "Novak shield density: " << NovakShieldDensity << "\n";
  G4Material* NovakShield = new G4Material("NovakShield", NovakShieldDensity,
					   ncomponents = 2);
  NovakShield -> AddMaterial(tungstenMat, tungstenProportion);
  NovakShield -> AddMaterial(epoxy, (1-tungstenProportion)); 

  //creating an instance of the class (?)  to access the voxel specifics struct
  B1DetectorConstruction::voxelSpecifics DetectorConstruction;

  G4double voxel_sizeX = DetectorConstruction.voxel_x;
  G4double voxel_sizeY = DetectorConstruction.voxel_y;
  G4double voxel_sizeZ = DetectorConstruction.voxel_z;
  G4String voxel_materialInt = DetectorConstruction.voxel_mat;
  
  // Envelope parameters
  //
  //G4double env_sizeXY = 10., env_sizeZ = 10.;
  
  //setting the env size .5 larger than the size of the voxel...
  G4double env_sizeX = voxel_sizeX + .5; //mm
  G4double env_sizeY = voxel_sizeY + .5; //mm
  G4double env_sizeZ = voxel_sizeZ + .5; //mm
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  //G4double world_sizeXY = 1.2*env_sizeXY;
  //G4double world_sizeZ  = 1.2*env_sizeZ;

  //setting the world size 1. lager than the size of the voxel...
  G4double world_sizeX = voxel_sizeX + 1.; //mm
  G4double world_sizeY = voxel_sizeY + 1.; //mm
  G4double world_sizeZ = voxel_sizeZ + 1.; //mm
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);     //its size
      
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
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeX, 0.5*env_sizeY, 0.5*env_sizeZ); //its size
      
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
 
  //     
  // Shape 1
  //  
  //G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  //G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);
        
  // Conical section shape       
  //G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  //G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  //G4double shape1_hz = 3.*cm;
  //G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  //G4double voxelSize = 7.;
  //G4Cons* solidShape1 =    
  //new G4Cons("Shape1", 
  //shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
  //shape1_phimin, shape1_phimax);
  /*                 
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
                    checkOverlaps);          //overlaps checking*/

  //     
  // Shape 2
  //
  //G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  //G4ThreeVector pos2 = G4ThreeVector(0, -1*cm, 7*cm);

  // Trapezoid shape       
  //G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  //G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  //G4double shape2_dz  = 6*cm;  
  //B1DetectorConstruction::voxelSpecifics DetectorConstruction;

  //G4double voxelSize = 7.;
  //G4double voxelThicc = DetectorConstruction.voxelThicc;

  G4Material* voxel_material;
  if(voxel_materialInt == 1){
    voxel_material = GAGGCe;
  }
  else if(voxel_materialInt == 2){
    voxel_material = InFoil;
  }
  else if(voxel_materialInt == 3){
    voxel_material = NovakShield;
  }
  else if(voxel_materialInt == 4){
    voxel_material = tungstenMat;
  }
  else if(voxel_materialInt == 5){
    voxel_material = epoxy;
  }
  else{
    voxel_material = NULL;
    std::cout << "Hey man this will cause a core dump btw >:)\n";
  }
  
  G4Box* voxelBox =    
    new G4Box("voxelBox",                      //its name
              .5*voxel_sizeX, .5*voxel_sizeY, .5*voxel_sizeZ); //its size
                
  G4LogicalVolume* voxelLog =                         
    new G4LogicalVolume(voxelBox,         //its solid
                        voxel_material,          //its material
                        "voxelLog");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,0),                    //at position
                    voxelLog,             //its logical volume
                    "voxel",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                
  // Set voxel log as the scoring volume
  //
  fScoringVolume = voxelLog;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
