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
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "B1DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0)
{

  //lol do NOT comment out if you want to fire particles using gun
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  /*
  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(6.*MeV);*/

  //srand(unsigned int(time(NULL))); //here is a random seed based on the time
  int RNGSeed = 42069; //or you can set it manually here
  srand(RNGSeed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  
  //G4double envSizeXY = 0;
  //G4double envSizeZ = 0;
  /*
  if (!fEnvelopeBox)
  {
    G4LogicalVolume* envLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if ( envLV ) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  }

  if ( fEnvelopeBox ) {
    //envSizeXY = fEnvelopeBox->GetXHalfLength()*2.;
    //envSizeZ = fEnvelopeBox->GetZHalfLength()*2.;
  }  
  else  {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n"; 
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("B1PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
     }*/

  //the goal here is isoptropic generation randomly inside the voxel...

  //creating an instance of B1DC in PGA so we can access voxelSpecifics struct...
  B1DetectorConstruction::voxelSpecifics PGA;

  G4double positionRange_x = PGA.voxel_x;
  G4double positionRange_y = PGA.voxel_y;
  G4double positionRange_z = PGA.voxel_z;

  //this will generate a value from -range/2. to range/2.
  G4double xpos = positionRange_x * (((double)std::rand()/RAND_MAX) - .5);
  G4double ypos = positionRange_y * (((double)std::rand()/RAND_MAX) - .5);
  G4double zpos = positionRange_z * (((double)std::rand()/RAND_MAX) - .5);

  //setting up isotropic momentum
  double theta = 2*M_PI*((double)std::rand()/RAND_MAX); //(0,2pi)
  double phi = 2*(((double)std::rand()/RAND_MAX) - .5); //(-1,1)
  phi = acos(phi); //lol this is math notation 

  G4double imo = cos(theta)*sin(phi);
  G4double jmo = sin(theta)*sin(phi);
  G4double kmo = cos(phi);

  //G4double energy = 10*((double)std::rand()/RAND_MAX);
  //energy = .174;
  
  fParticleGun->SetParticlePosition(G4ThreeVector(xpos,ypos,zpos));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(imo,jmo,kmo));
  //fParticleGun->SetParticleEnergy(energy);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

