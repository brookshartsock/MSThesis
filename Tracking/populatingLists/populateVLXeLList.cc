#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>

//"Globals"
int ROWS; //const int ROWS = 10000;
int COLS = 10; //const int COLS = 10;
int **data;
float **momentum;
const int cubeNum = 6;


//written with the help of GPT4
bool readStoreData(std::string filename){


  std::ifstream infile(filename.c_str());
  std::string line;
  int row = 0;
  int col = 0;

  if(!infile.is_open()){
    std::cout << "No such file exists." << std::endl;
    return false;
  }

  const int rows = ROWS;
  const int momentumCols = 3;
  const int dataCols = COLS - momentumCols;

  data = new int*[rows];
  momentum = new float*[rows];

  for(int i = 0; i < rows; ++i){

    data[i] = new int[dataCols];
    momentum[i] = new float[momentumCols];

  }

  while(std::getline(infile, line) && row < ROWS){

    std::istringstream iss(line);
    std::string val;

    col = 0;
    while(std::getline(iss, val, ',')){

      std::istringstream convert(val);
      if(col < momentumCols){
	if(!(convert >> momentum[row][col])){
	  momentum[row][col] = 0.0f;
	}
      } 
      else{
	if(!(convert >> data[row][col - momentumCols])){
	  data[row][col - momentumCols] = cubeNum*cubeNum*cubeNum;
	}
      }

      col++;
    }

    while(col < COLS + momentumCols){

      if(col >= momentumCols){
	data[row][col - momentumCols] = cubeNum*cubeNum*cubeNum;
      }
      col++;
    }
    row++;
  }


  return true;
}


void writeTermData(){


  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < 3; ++j){
      std::cout << momentum[i][j] << " ";
    }
    for(int j = 0; j < COLS-3; ++j){
      std::cout << data[i][j] << " ";
    }
    std::cout << "\n";
  }


}

//x dir can be negative UPDATE :(
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

    if(n == 0){ //set initial voxel
      IDi[n] = i;
      IDdata[n][0] = j;
      IDdata[n][1] = k;
    }
    else{
      IDi[n] = std::abs(i - IDi[0]);

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


//written with the help of GPT4
void lookupList(std::string ID, int run) {


  std::ifstream inputFile("lookupList.csv");
  std::vector<std::string> rows;

  // Read the existing rows from the file
  if(inputFile.is_open()){
    std::string line;
    while(std::getline(inputFile, line)){
      rows.push_back(line);
    }
    inputFile.close();
  }

  double slopeY, slopeZ;
  double yi = momentum[run][1] / momentum[run][0];
  double zi = momentum[run][2] / momentum[run][0];
  double STDEVY, STDEVZ;

  unsigned int events = 0;
  double oldSlopeY = 0, oldSlopeZ = 0, oldSTDEVY = 0, oldSTDEVZ = 0;

  bool found = false;
  size_t foundRowIndex = 0;

  // Search through lookupList.csv
  for(size_t i = 0; i < rows.size(); ++i){

    std::string& row = rows[i];
    std::istringstream iss(row);
    std::string cell;
    std::vector<std::string> cells;

    while(std::getline(iss, cell, ',')){
      cells.push_back(cell);
    }

    if(cells.size() >= 6 && cells[0] == ID){
      found = true;
      foundRowIndex = i;
      events = std::atoi(cells[1].c_str());
      oldSlopeY = std::atof(cells[2].c_str());
      oldSTDEVY = std::atof(cells[3].c_str());
      oldSlopeZ = std::atof(cells[4].c_str());
      oldSTDEVZ = std::atof(cells[5].c_str());
      break;
    }
  }

  if(found){
    slopeY = ((oldSlopeY * events) + yi) / (events + 1.0);
    slopeZ = ((oldSlopeZ * events) + zi) / (events + 1.0);

    double SSy = ((oldSTDEVY * oldSTDEVY) + (oldSlopeY * oldSlopeY)) 
      * events + (yi * yi);
    double SSz = ((oldSTDEVZ * oldSTDEVZ) + (oldSlopeZ * oldSlopeZ)) 
      * events + (zi * zi);

    STDEVY = std::sqrt((SSy / (events + 1.0)) - (slopeY * slopeY));
    STDEVZ = std::sqrt((SSz / (events + 1.0)) - (slopeZ * slopeZ));

    // Replace that row in .csv
    std::ostringstream oss;
    oss << ID << "," << (events + 1) << "," << slopeY << "," << STDEVY 
	<< "," << slopeZ << "," << STDEVZ;
    rows[foundRowIndex] = oss.str();
  } 
  else {
    events = 1;
    slopeY = yi;
    slopeZ = zi;
    STDEVY = 0.0;
    STDEVZ = 0.0;

    // Create a new row in the .csv
    std::ostringstream oss;
    oss << ID << "," << events << "," << slopeY << "," << STDEVY 
	<< "," << slopeZ << "," << STDEVZ;
    rows.push_back(oss.str());
  }

  // Write the modified/added rows back to the file
  std::ofstream outputFile("lookupList.csv");
  if (outputFile.is_open()) {
    for (size_t i = 0; i < rows.size(); ++i) {
      outputFile << rows[i] << "\n";
    }
    outputFile.close();
  } else {
    std::cerr << "Error opening lookupList.csv for writing." << std::endl;
  }
}


int main(){


  std::cout << "\n";
  std::cout << " --- Look-Up List v-1.0 --- \n";
  std::cout << "\n";

  int hits;
  std::string ID;

  std::string filename;
  //std::cout << "Enter the filename: ";
  //std::cin >> filename;
  filename = "data.csv"; // FIX FIX FIX FIX FIX FIX FIX FIX

  // std::cout << "Event number (rows): ";
  //std::cin >> ROWS;
  ROWS = 1000000;  // FIX FIX FIX FIX FIX FIX FIX FIX

  if(readStoreData(filename)){
    //writeTermData();
    
    for(int x = 0; x < ROWS; ++x){

      hits = 0;
      while(data[x][hits] != (cubeNum*cubeNum*cubeNum)){
	hits += 1;
      }

      if(hits > 1 && hits < (cubeNum+1)){

	ID = E2IDNo(x, hits);

	if(ID != ""){
	  lookupList(ID,x);
	}
      }
 
    }
  }

  delete[] data;


  return 0;
}
