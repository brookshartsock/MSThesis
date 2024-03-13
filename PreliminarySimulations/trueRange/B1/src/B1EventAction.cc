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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
  : G4UserEventAction(),
    fRunAction(runAction),
    fEdep(0.)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
  ResetTotalDistances();
}

void B1EventAction::AddDistance(G4int trackID, G4double distance) {
  // Accumulate distance for the given track ID
  fTotalDistances[trackID] += distance;
}

G4double B1EventAction::GetTotalDistance(G4int trackID) const {
  // Retrieve the total distance traveled by the particle with the given track ID
  auto it = fTotalDistances.find(trackID);
  if (it != fTotalDistances.end()) {
    return it->second;
  } else {
    return 0.0;
  }
}

void B1EventAction::ResetTotalDistances() {
  // Clear the distances map
  fTotalDistances.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event*)
{   
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);
  /*
    if(fEdep > -1){
    std::fstream outfile;
    outfile.open("data.csv", std::fstream::in | std::fstream::out 
    | std::fstream::app);
    outfile << fEdep << "\n";
    outfile.close();
    }*/

  G4double totalDistance = GetTotalDistance(1);
  //G4cout << "Total distance traveled by the primary particle: " << totalDistance << " mm" << G4endl;

  std::fstream outfile;
  outfile.open("data.csv", std::fstream::in | std::fstream::out 
	       | std::fstream::app);
  outfile << totalDistance << "\n";
  outfile.close();


  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
