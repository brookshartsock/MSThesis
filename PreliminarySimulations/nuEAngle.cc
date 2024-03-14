//uses root to compare angle between nu and the  e following the Ga interaction
//nu + Ga -> e- + Ge+
//  lhs   ->   rhs  
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

//root includes
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TCanvas.h>

//globals
double amu = 931.4941024; // MeV per amu

double mGa = 70.924701 * amu; //mass Ga-71
double me = .510999; //mass e-
double mGe = (70.924954 * amu) - me; //mass Ge-71+
double energyExit = .1749; //energy of excited state


//numerically solves (guesses) for and returns pe
double momentumGuesser(double uNu){


  bool status = 1; //used for controlling while loop
  double dpe = 1.e5; //initial change and val for pe, 1.e5 for high uNu
  double pe = dpe;
  double diff;

  TLorentzVector nu(0,0,uNu,uNu); //neutrino Lorentz Vec.
  TLorentzVector Ga(0,0,0,mGa); //gallium
  TLorentzVector total = nu+Ga; //total for lhs

  double inM = total.M(); //inv. mass for lhs
  double inMguess;

  while(status){

    TLorentzVector e(0,0,pe,std::sqrt((me*me) + (pe*pe))); //e-
    TLorentzVector Ge(0,0,-pe,std::sqrt((mGe*mGe) + (pe*pe))); //Ge+

    total = e + Ge;
    inMguess = total.M(); //inv. mass for rhs

    diff = inM - inMguess;
    
    if(diff > 0){ //if rhs > lhs
      pe += dpe; //increase pe
    }
    else{
      if(std::abs(diff) < 1e-5){ //check to see if rhs ~= lhs
	status = 0; //terminate loop
      }
      else{ // rhs !~= lhs and lhs > rhs
	pe -= dpe; //decrease pe and dpe
	dpe *= .1;
	pe += dpe;
      }
    }

  }

  //std::cout << "Momentum guess: " << pe << " MeV\n";


  return pe;
}


//written with the help of Jonathan Folkerts
//computes Lorentz boost and returns dot product of nu/e vecs
double boost(double energyNu, double pe){


  TLorentzVector nu(0,0,energyNu,energyNu); //neutrino 
  TLorentzVector Ga(0,0,0,mGa); //gallium (again)
  TLorentzVector total = nu+Ga;

  double theta = 2*M_PI*((double)std::rand()/RAND_MAX); //(0,2pi)
  double phi = 2*(((double)std::rand()/RAND_MAX) - .5); //(-1,1)
  phi = acos(phi); //lol this is math notation 

  double px = cos(theta)*sin(phi); //this is
  double py = sin(theta)*sin(phi); //already
  double pz = cos(phi); //normalized >:)

  TLorentzVector e(px,py,pz,std::sqrt((me*me) + (pe*pe))); //hey look an e-

  TVector3 boost = total.BoostVector();
  e.Boost(boost); //boost

  px = e.Px(); //reset pn as vec. values
  py = e.Py();
  pz = e.Pz();

  double norm = std::sqrt((px*px) + (py*py) + (pz*pz)); //calc norm. const.

  px /= norm; //normalize
  py /= norm;
  pz /= norm;

  double val = std::acos(pz); //dot product, (0,0,1) \cdot (px,py,pz)
  

  return val;
}


//calls functions and stores data into csvs
void takeData(double uNu){


  double mean = 0.;
  double variance = 0.;
  
  std::string outFilename;
  outFilename = std::to_string(uNu) + ".pdf"; //double to string 

  double pe = momentumGuesser(uNu); //get momemtum estimate

  std::string title = std::to_string(uNu) //title string
    + " MeV Neutrino;\\acos(\\theta);Counts";
  TH1D *hist = new TH1D("hist", title.c_str(), 1e3, 0, M_PI); //create hist

  int iMax = 1e6; //no. of total events per histData csv
  for(int i = 0; i < iMax; ++i){
  
    double val = boost(uNu, pe); //call boost function

    mean += val;
    variance += val * val;

    hist -> Fill(val); //Fill hist

  }

  TCanvas *c = new TCanvas("c", "c", 1920, 1080); //name name size size
  hist -> Draw("hist");
  c -> Print(outFilename.c_str());

  delete hist;
  delete c;

  std::fstream outData; //for mean/std of each energy
  outData.open("data.csv", std::fstream::in | std::fstream::out 
	       | std::fstream::app);

  mean /= (double)iMax;
  variance /= (double)iMax;
  variance -= mean * mean;
  double stdev = std::sqrt(variance);

  outData << uNu << "," << mean << "," << stdev << "," << pe << "\n";

  outData.close();

  
}


//makes data.csv file, calculates energyNu and calls data function
int main(){


  std::fstream outData; //makes file for mean/std
  outData.open("data.csv", std::fstream::in | std::fstream::out 
	       | std::fstream::trunc);
  outData << "energy,mean,std,pe\n";
  outData.close();

  //increases energyNu on log scale
  for(double energyNu = .5; energyNu <= 1.e7; energyNu *= 1.5){
    
    takeData(energyNu);

  }
  

  return 0;
}
