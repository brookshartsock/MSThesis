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
/// \file optical/LXe/include/LXePrimaryGeneratorAction.hh
/// \brief Definition of the LXePrimaryGeneratorAction class
//
//
#ifndef LXePrimaryGeneratorAction_h
#define LXePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4String.hh"

class G4ParticleGun;
class G4Event;
//class LXeEventAction;

struct eData{
    G4double imo;
    G4double jmo;
    G4double kmo;
    G4double priEnergy;
    };

class LXePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  LXePrimaryGeneratorAction();
  virtual ~LXePrimaryGeneratorAction();

  /*
  struct eData{
    G4double imo;
    G4double jmo;
    G4double kmo;
    G4double priEnergy;
    };*/
  
  const eData& GetPGAData() const {return fPGAData;}

  
 
public:

  virtual void GeneratePrimaries(G4Event* anEvent);

private:

  G4ParticleGun* fParticleGun;
  //G4String PGAData;
  eData fPGAData;
};

#endif
