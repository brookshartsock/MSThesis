#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <map>


// Global variables
const int nElements = 110;
const int nIsotopes = 41; //40 is max, need 41 to check for 0's ;)
double isotopeData[nElements][nIsotopes][2] = {0.};

// Map to store element names and their corresponding index in massData
std::map<std::string, int> elementIndex;


void readFillIsotopeData(){

  
  std::ifstream file("massData.csv");

  if (!file.is_open()) {
    std::cerr << "Error opening file." << std::endl;
    return;
  }

  int currentElementIndex = -1;
  int isotopeIndex = 0;

  std::string line;
  while(std::getline(file, line)){
    
    std::stringstream ss(line);
    std::string token;
    std::getline(ss, token, ',');

    if(isalpha(token[0])){ // Check if the token is an element name
      if (elementIndex.find(token) == elementIndex.end()) { // New element
	currentElementIndex++;
	elementIndex[token] = currentElementIndex;
      }
      else{ // Existing element
	currentElementIndex = elementIndex[token];
      }
      isotopeIndex = 0;
    }
    
    else if(!token.empty()){ // Token is a mass number
      double mass = std::stod(token);
      isotopeData[currentElementIndex+1][isotopeIndex][0] = mass;
            
      // Read the second column value if present
      if(std::getline(ss, token, ',')){
	double value = std::stod(token);
	isotopeData[currentElementIndex+1][isotopeIndex][1] = value;
      }
      else{
	isotopeData[currentElementIndex+1][isotopeIndex][1] = 0.0;
      }

      isotopeIndex++;
    }
    
  }

  file.close();

  
}


//test that mass data / abundance is working
void coutIsotopeData(){
  
  for(int Z = 0; Z < nElements; Z++){
    for(int isoNum = 0; isoNum < nIsotopes; ++isoNum){ 
      std::cout << Z << "," << isoNum << "," << isotopeData[Z][isoNum][0];
      std::cout << "," << isotopeData[Z][isoNum][1] << "\n"; 
    }
  }

}


//run to find max number of isotopes per atom
//lol it's 40, don't set nIsotopes = 40 or it will break
//because it won't ever find a 0 >:)
int findMaxIsotopeNumber(){


  int maxIsotopeNumber = 0;
  int isotopeCount = 0;
  int currentZ = 0;
  
  for(int Z = 0; Z < nElements; Z++){
    
    for(int isoNum = 0; isoNum < nIsotopes; ++isoNum){

      if(isotopeData[Z][isoNum][0] != 0.){
	isotopeCount += 1;
	
	if(isotopeCount > maxIsotopeNumber){
	  maxIsotopeNumber = isotopeCount;
	}
      }

      else{
	isotopeCount = 0;
      }
      
    }
    
  }


  return maxIsotopeNumber;
}


//gives index number corresponding to mass number
//will return -1 if mass doesn't exist...
//check to see if the mass data exists for an element
//based on proton number - Z and mass number - A
int findAIndex(int Z, int A){


  int AIndex = 0;

  bool condition = 1;
  while(condition){
    
    if(std::round(isotopeData[Z][AIndex][0]) == A){
      condition = 0;
    }
    else if(isotopeData[Z][AIndex][0] == 0){
      AIndex = -1;
      condition = 0;
    }
    else{
      ++AIndex;
    }
    
  }
  

  return AIndex;
}


//testing massExists function
void coutAIndex(){

  int maxElements = nElements/2.; //nElements
  for(int Z = 0; Z < maxElements; Z++){
    for(int isoNum = Z; isoNum <= (Z*3); ++isoNum){
      std::cout << Z << "," << isoNum << "," << findAIndex(Z, isoNum) << "\n";
    }
  }

}


//check to see if naturally abundant (NA)
bool isNA(int Z, int A){


  int AIndex = findAIndex(Z,A);
  
  if(AIndex == -1){ 
    return 0;
  }
  else if(isotopeData[Z][AIndex][1] == 0.){
    return 0;
  }
  else{
    return 1;
  }

  
}


//testing isNA function
void coutIsNA(){

  int maxElements = nElements/2.; //nElements
  for(int Z = 0; Z < maxElements; Z++){
    for(int isoNum = Z; isoNum <= (Z*3); ++isoNum){
      std::cout << Z << "," << isoNum << "," << isNA(Z, isoNum) << "\n";
    }
  }

}


//calculated IBD threshold energy given target X(Z,A)
double calculateIBDThreshold(int Z, int A){


  double mElectron = 0.5109989461; //MeV
  double mNeutron = 939.565346; //MeV
  double mAMU = 931.494102; //MeV

  int AIndex = findAIndex(Z,A);
  //std::cout << "Mass of target is " << isotopeData[Z][AIndex][0] << " AMU\n";
  double initialEnergy = isotopeData[Z][AIndex][0] * mAMU;

  AIndex = findAIndex(Z-1,A-1);
  //std::cout << "Mass of resultant is " << isotopeData[Z-1][AIndex][0];
  //std::cout << " AMU\n";
  double finalEnergy = isotopeData[Z-1][AIndex][0] * mAMU;
  finalEnergy += (2* mElectron) + mNeutron;

  double threshold = (finalEnergy*finalEnergy) - (initialEnergy*initialEnergy);
  threshold /= (2. * initialEnergy);


  return threshold;
}


//calculated neutrino CC threshold energy given target X(Z,A)
double calculateNuCCThreshold(int Z, int A){


  double mAMU = 931.494102; //MeV

  int AIndex = findAIndex(Z,A);
  //std::cout << "Mass of target is " << isotopeData[Z][AIndex][0] << " AMU\n";
  double initialEnergy = isotopeData[Z][AIndex][0] * mAMU;

  AIndex = findAIndex(Z+1,A);
  //std::cout << "Mass of resultant is " << isotopeData[Z-1][AIndex][0];
  //std::cout << " AMU\n";
  double finalEnergy = isotopeData[Z+1][AIndex][0] * mAMU;

  double threshold = (finalEnergy*finalEnergy) - (initialEnergy*initialEnergy);
  threshold /= (2. * initialEnergy);


  return threshold;
}


//calculated anti nu CC threshold energy given target X(Z,A)
double calculateANuCCThreshold(int Z, int A){


  double mElectron = 0.5109989461; //MeV
  double mAMU = 931.494102; //MeV

  int AIndex = findAIndex(Z,A);
  //std::cout << "Mass of target is " << isotopeData[Z][AIndex][0] << " AMU\n";
  double initialEnergy = isotopeData[Z][AIndex][0] * mAMU;

  AIndex = findAIndex(Z-1,A);
  //std::cout << "Mass of resultant is " << isotopeData[Z-1][AIndex][0];
  //std::cout << " AMU\n";
  double finalEnergy = isotopeData[Z-1][AIndex][0] * mAMU;
  finalEnergy += (2* mElectron);

  double threshold = (finalEnergy*finalEnergy) - (initialEnergy*initialEnergy);
  threshold /= (2. * initialEnergy);


  return threshold;
}


int main() {

  
  readFillIsotopeData();
  //coutIsotopeData();

  //std::cout << "Max isotopes per element: " << findMaxIsotopeNumber();
  //std::cout << "\n";

  //coutAIndex();

  //coutIsNA();

  //for IBD if(isNA(Z,A) && (findAIndex(Z-1,A-1) != -1))
  //for nuCC if(isNA(Z,A) && (findAIndex(Z+1,A) != -1))
  //for AnuCC if(isNA(Z,A) && (findAIndex(Z-1,A) != -1))

  std::fstream outData; //for mean/std of each energy
  outData.open("outThresholdsData.csv", std::fstream::in | std::fstream::out 
	       | std::fstream::trunc);

  outData << "Z,A,NuCCThreshold_MeV,IBDThreshold,ANuCCThreshold\n";

  for(int Z = 1; Z < nElements; ++Z){
    
    for(int A = Z; A <= (Z*6); ++A){

      if(isNA(Z,A)){
	outData << Z << "," << A << ",";

	if(findAIndex(Z+1,A) != -1){
	  outData << calculateNuCCThreshold(Z,A) << ",";
	}
	else{
	  outData << "NA" << ",";
	}

	if(findAIndex(Z-1,A-1) != -1){
	  outData << calculateIBDThreshold(Z,A) << ",";
	}
	else{
	  outData << "NA" << ",";
	}

	if(findAIndex(Z-1,A) != -1){
	  outData << calculateANuCCThreshold(Z,A) << "\n";
	}
	else{
	  outData << "NA\n";
	}
      }
     
    }
    
  }

  outData.close();
      

  return 0;
}
