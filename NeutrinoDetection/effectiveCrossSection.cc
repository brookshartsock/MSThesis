#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

const int arraySize = 75;
double energy[arraySize];
double flux[arraySize];


//written with help of GPT4
bool readStoreData(const std::string& filePath){

  
  std::ifstream file(filePath);
  if(!file.is_open()){
    return false;
  }

  std::string line;
  int index = 0;

  // Skip the header line
  std::getline(file, line);

  // Read each line from the file
  while(std::getline(file, line) && index < arraySize){
      
    std::stringstream ss(line);
    std::string value;

    // Get the energy value
    std::getline(ss, value, ',');
    energy[index] = std::stod(value);

    // Get the flux value
    std::getline(ss, value, ',');
    flux[index] = std::stod(value);

    index++;
	
  }

  file.close();

    
  return 1;
}


void coutData(){
  

  for(int i = 0; i < arraySize; i++){
    std::cout << "Energy: " << energy[i] << ", Flux: ";
    std::cout << flux[i] << std::endl;
  }

}


void normalizeFlux(){


  double fluxSum = 0;
  
  for(int i = 0; i < (arraySize-1); ++i){

    fluxSum += .5 * (flux[i] + flux[i+1]) * (energy[i+1] - energy[i]);

  }

  for(int i = 0; i < arraySize; ++i){

    flux[i] /= fluxSum;

  }

  /*
  //test normalize
  fluxSum = 0;
  for(int i = 0; i < (arraySize-1); ++i){

  fluxSum += .5 * (flux[i] + flux[i+1]) * (energy[i+1] - energy[i]);

  }
  std::cout << fluxSum << "\n";
  */
  

}


//for any value of energy, returning the linearly interpolated flux
double interpolateFlux(double energyVal){


  double estFlux;
  
  int i = 0;
  while(i < (arraySize - 1)){
    
    if(energyVal >= energy[i] && energyVal <= energy[i + 1]){
      double slope = (flux[i + 1] - flux[i]) / (energy[i + 1] - energy[i]);
      estFlux =  flux[i] + slope * (energyVal - energy[i]);

      i = arraySize;
    }

    ++i;
    
  }
  

  return estFlux;
}


double calculateFluxProportion(double threshold){


  double fluxProportion;
  bool condition = 1;
  
  for(int i = 0; i < (arraySize - 1); ++i){

    if(energy[i] > threshold){
      
      if(condition){ //first sum needs to consider interpolated point
	fluxProportion = .5 * (interpolateFlux(threshold) + flux[i+1]);
	fluxProportion *= (energy[i+1] - threshold);

	//std::cout << "check\n";
	//std::cout << .5 * (interpolateFlux(threshold) + flux[i+1]) << "\n";
	//std::cout << (energy[i+1] - threshold) << "\n";

	condition = 0;
      }
      else{
	fluxProportion += .5 * (flux[i] + flux[i+1]) * (energy[i+1] - energy[i]);
      }
      
    }
    
  }
  

  return fluxProportion;
}


double calculateCrossSection(double energy, double threshold, double A){
  return A * .01547 * std::pow((energy - threshold),1.8001);
}


double findXSecRatio(double threshold){


  double A = 1.; //cross section ratio A*\sigma

  double HConvolution = 0; //convolution for IBD on Hydrogen
  for(int i = 0; i < (arraySize-1); ++i){

    if(energy[i] > 1.806){
      HConvolution += .25 * (flux[i] + flux[i+1])
	* (calculateCrossSection(energy[i],1.806,A)
	   + calculateCrossSection(energy[i+1],1.806,A))
	* (energy[i+1]-energy[i]);
    }

  }

  //std::cout << HConvolution << "\n";

  double sum = 0;
  while(sum < HConvolution){

    sum = 0;
    for(int i = 0; i < (arraySize-1); ++i){

      if(energy[i] > threshold){
	sum += .25 * (flux[i] + flux[i+1])
	  * (calculateCrossSection(energy[i],threshold,A)
	     + calculateCrossSection(energy[i+1],threshold,A))
	  * (energy[i+1]-energy[i]);
      }

    }

    A += .00001;
  }

  A -= .00001;

  
  return A;
}


int main(){

  
  if(readStoreData("ANuSpectrum.csv")){

    //coutData();
    normalizeFlux();
    //coutData();

    //std::cout << interpolateFlux(.1923) << "\n";

    /*
      double threshold = .1923; //MeV
      double ANuProportion = calculateFluxProportion(threshold); //not MeV

      std::cout << "With a threshold of " << threshold << " MeV, you can detect ";
      std::cout << ANuProportion << " of the anti-neutrino spectrum.\n";
    */

    std::cout << findXSecRatio(2.167) << "\n";
	
  }

    
  return 0;
}
