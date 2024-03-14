#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>

const int cubeNum = 6;
const int ROWS = 100; //1e4 or 881117
const int COLS = cubeNum*2; //15
int data[ROWS][COLS];
float momentum[ROWS][3];
float energy[ROWS];


//written with the help of GPT4
bool readStoreData(std::string filename) {
  std::ifstream infile(filename.c_str());
  std::string line;
  int row = 0;
  int col = 0;

  if (!infile.is_open()) {
    std::cout << "No such file exists." << std::endl;
    return false;
  }

  const int rows = ROWS;
  const int momentumCols = 3;
  const int energyCol = 1;
  const int dataCols = COLS - energyCol - momentumCols;

  while (std::getline(infile, line) && row < ROWS) {
    std::istringstream iss(line);
    std::string val;
    col = 0;
    while (std::getline(iss, val, ',')) {
      std::istringstream convert(val);
      if (col < energyCol) {
	if (!(convert >> energy[row])) {
	  energy[row] = 0.0f;
	}
      } else if (col < energyCol + momentumCols) {
	if (!(convert >> momentum[row][col - energyCol])) {
	  momentum[row][col - energyCol] = 0.0f;
	}
      } else {
	if (!(convert >> data[row][col - energyCol - momentumCols])) {
	  data[row][col - energyCol - momentumCols] 
	    = cubeNum * cubeNum * cubeNum;
	}
      }
      col++;
    }
    while (col < COLS + energyCol + momentumCols) {
      if (col >= energyCol + momentumCols) {
	data[row][col - energyCol - momentumCols] 
	  = cubeNum * cubeNum * cubeNum;
      }
      col++;
    }
    row++;
  }

  return true;
}


//so... energy is imo an imo is now energy (oops); what a fucking dumbass
void fixDataStructureLol(){


  for(int i = 0; i < ROWS; ++i){

    double tempEnergy = momentum[i][2];
    momentum[i][2] = momentum[i][1];
    momentum[i][1] = momentum[i][0];
    momentum[i][0] = energy[i];
    energy[i] = tempEnergy;
    
  }
  

}


void writeTermData(){

  for(int i = 0; i < ROWS; ++i){
    std::cout << energy[i] << " ";
    for(int j = 0; j < 3; ++j){
      std::cout << momentum[i][j] << " ";
    }
    for(int j = 0; j < COLS-3; ++j){
      std::cout << data[i][j] << " ";
    }
    std::cout << "\n";
  }

}


std::string E2IDNo(int run, const int hits){

  //general
  int n, IDdata[cubeNum][2], IDi[cubeNum];
  std::string IDNo;
  bool iRepeats = 0;

  //G4 variables
  int CPNo, i, j, k;

  //sorts data from lowest to highest CPNo;
  for(int iSort = 0; iSort < hits - 1; ++iSort){

    for(int jSort = 0; jSort < hits - iSort - 1; ++jSort){

      if(data[run][jSort] > data[run][jSort + 1]){
	int temp = data[run][jSort];
	data[run][jSort] = data[run][jSort + 1];
	data[run][jSort + 1] = temp;
      }

    }

  }

  for(n = 0; n < cubeNum; ++n){ //store as empty to start. 23 = x
    IDdata[n][0] = 23;
    IDdata[n][1] = 23;
  }
  
  for(n = 0; n < hits; ++n){ //store data in IDdata[][]
        
    CPNo = data[run][n];

    i = 0;
    j = 0;
    k = 0;
    while(CPNo >= cubeNum*cubeNum){
      CPNo -= cubeNum*cubeNum;
      i += 1;
    }
    while(CPNo >= cubeNum){
      CPNo -= cubeNum;
      j += 1;
    }
    while(CPNo >= 1){
      CPNo -= 1;
      k += 1;
    }

    //std::cout << "Voxel: " << i << "," << j << "," << k << "\n";

    if(n == 0){ //set initial voxel
      IDi[n] = i;
      IDdata[n][0] = j;
      IDdata[n][1] = k;
    }
    else{
      IDi[n] = i - IDi[0];

      if(IDi[n] == 0){ //checking for i value equal to vertex
	return "";
      }

      if(n > 1){
	for(int m = 1; m < n; ++m){ //checking for repeats in i after n = 1
	  if(IDi[m] == IDi[n]){ 
	    return "";
	  }
	}
      }
      
      if((j - IDdata[0][0]) < 0){
	IDdata[IDi[n]][0] = -(j - IDdata[0][0]) + 10;
      }
      else{
	IDdata[IDi[n]][0] = j - IDdata[0][0];
      }
      
      if((k - IDdata[0][1]) < 0){
	IDdata[IDi[n]][1] = -(k - IDdata[0][1]) + 10;
      }
      else{
	IDdata[IDi[n]][1] = k - IDdata[0][1];
      }
   
    }
  }    
    
  //output for testing
  /*
    for(n = 0; n < cubeNum; ++n){
    std::cout << IDdata[n][0] << "," << IDdata[n][1] << ",";
    }
    std::cout << "\n";*/

  for(n = 1; n < cubeNum; ++n){
    IDNo.insert(0, 1, (IDdata[n][0] + 97));
    IDNo.insert(0, 1, (IDdata[n][1] + 97));
  }

  return IDNo;
}


bool linReg(std::string filename, int run, const int hits,
            double slopeY, double slopeZ, std::string ID) {

  std::string lookupFile = "lookupList.csv";
  std::ifstream file(lookupFile.c_str());
  if (!file.is_open()) {
    std::cerr << "Failed to open the file." << std::endl;
    return false;
  }

  std::string line;
  std::getline(file, line); // Skip the header line

  while (std::getline(file, line)) {
    std::istringstream ss(line);
    std::string id;
    if (std::getline(ss, id, ',')) {
      if (id == ID) {
        // Found the ID, print the other columns
        std::string eventsStr, slopeYStr, STDYStr, slopeZStr, STDZStr;
        if (std::getline(ss, eventsStr, ',') &&
            std::getline(ss, slopeYStr, ',') &&
            std::getline(ss, STDYStr, ',') &&
            std::getline(ss, slopeZStr, ',') &&
            std::getline(ss, STDZStr)) {
          int events = atoi(eventsStr.c_str());
          double mTheta = atof(slopeYStr.c_str());
          double STDY = atof(STDYStr.c_str());
          double mPhi = atof(slopeZStr.c_str());
          double STDZ = atof(STDZStr.c_str());

          std::fstream outData;
          outData.open(filename.c_str(), std::fstream::in | std::fstream::out 
                       | std::fstream::app);
          outData << run + 1 << "," << hits << "," << energy[run]
                  << "," << slopeY << "," << slopeZ << "," << mTheta
                  << "," << mPhi << "," << STDY << ","
                  << STDZ << "\n";
          outData.close();

          return true;
        }
      }
    }
  }

  // std::cout << "ID not found." << std::endl;
  return false; // ID not found
}


int main(){


  int hits;
  std::string filename = "data.csv";

  std::cout << "\n";
  std::cout << " --- Lookup List Analysis for VLXE --- \n";
  std::cout << "\n";

  if(readStoreData(filename)){
    fixDataStructureLol();
    //writeTermData();

    ///*
    std::fstream outData;
    filename.insert(0, "out");
    outData.open(filename.c_str(), std::fstream::in | std::fstream::out 
		 | std::fstream::trunc);
    outData << "run,hits,energy,slopeY,slopeZ,estSlopeY,estSlopeZ,STDY,STDZ\n";
    outData.close();

    for(int x = 0; x < ROWS; ++x){

      hits = 0;
      while(data[x][hits] != cubeNum*cubeNum*cubeNum){

	hits += 1;

      }

      if(hits > 1 && (hits < (cubeNum+1))){
	double slopeY = momentum[x][1]/momentum[x][0];
	double slopeZ = momentum[x][2]/momentum[x][0];

	std::string ID = E2IDNo(x, hits);

	if(linReg(filename, x, hits, slopeY, slopeZ, ID));
      }

    }//*/

  }

  //delete[] data;
  //delete[] momentum;


  return 0;
}
