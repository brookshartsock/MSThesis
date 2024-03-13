#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>

//globals
const int cubeNum = 8;
const int ROWS = 1e6; //1e6
bool doubleSiPM = 1;
bool inlin3 = 1;
bool fixDS = 1e-6; //reorders columns because I'm mentally deficient 

const int COLS = cubeNum*2; //15
double data[ROWS][COLS]; //should be double for double SiPM (ha)
float momentum[ROWS][3]; //0,1,2 -> i,j,k
float energy[ROWS];
int nanStatus = 0;
double RSQO = 1e-6; //r-squared offset (to help prevent NANs?)


bool readStoreData(std::string filename) {

  
  std::ifstream infile(filename.c_str());
  std::string line;
  int row = 0;
  int col = 0;

  if (!infile.is_open()) {
    std::cout << "No such file exists." << std::endl;
    return 0;
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
    

  return 1;
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
      std::cout << data[i][j] << " "; //can std::round() here to check >:)
    }
    std::cout << "\n";
  }

  
}


bool linReg(std::string filename, int run, const int hits,
	    double slopeY, double slopeZ){


  std::fstream outData;
  outData.open(filename.c_str(), std::fstream::in | std::fstream::out 
	       | std::fstream::app);

  //general
  int n;

  //G4 variables
  int CPNo, i, j, k;
  double voxelDim = 7., voxelThicc = 1., voxelSpace = 2.;
  double p0 = .4238; //SiPMDR = tanh(p0 * jpos)

  double fExpHall_x;
  double fExpHall_y;
  double fExpHall_z;
 
  if(inlin3){//for inline orientation
    fExpHall_x = (cubeNum * voxelThicc) + ((cubeNum + 1) * voxelSpace);
    fExpHall_y = (cubeNum * voxelDim) + ((cubeNum + 1) * voxelSpace);
    fExpHall_z = (cubeNum * voxelDim) + ((cubeNum + 1) * voxelSpace);
  }

  else{//offset
    fExpHall_x = (cubeNum * voxelThicc) + ((cubeNum + 1) * voxelSpace);
    fExpHall_y = (cubeNum * voxelDim) + ((cubeNum + 1) * voxelSpace)
      + (((cubeNum-1)/2.) * (voxelDim+voxelSpace));
    fExpHall_z = (cubeNum * voxelDim) + ((cubeNum + 1) * voxelSpace)
      + (((cubeNum-1)/2.) * (voxelDim+voxelSpace));
  }

  /*
    std::cout << "Hall dims: " << fExpHall_x << "," << fExpHall_y
    << "," << fExpHall_z << "\n";*/

  //stuff for lin reg and R^2
  double xSum = 0, xSqSum = 0;
  double ySum = 0, xySum = 0;
  double zSum = 0, xzSum = 0;

  double deltaTheta, mTheta, bTheta; //for y
  double yAve = 0, rSquaredTheta;
  double numerTheta = 0, denomTheta = 0;

  double deltaPhi, mPhi, bPhi; //for z
  double zAve = 0, rSquaredPhi;
  double numerPhi = 0, denomPhi = 0;
    
  double ipos[hits], jpos[hits], kpos[hits];

  for(n = 0; n < hits; ++n){
        
    CPNo = std::round(data[run][n]);
    double SiPMDifferenceRatio = (data[run][n] - (double)CPNo)*4.;
    //std::cout << SiPMDifferenceRatio << "\n";

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

    if(inlin3){//inline
      ipos[n] = -(fExpHall_x/2.) + voxelSpace + (voxelThicc/2.) 
	+ (i * (voxelThicc+voxelSpace));
      jpos[n] = -(fExpHall_y/2.) + voxelSpace + (voxelDim/2.) 
	+ (j * (voxelDim+voxelSpace));
      kpos[n] = -(fExpHall_z/2.) + voxelSpace + (voxelDim/2.) 
	+ (k * (voxelDim+voxelSpace));
    }

    else{//offset
      ipos[n] = -(fExpHall_x/2.) + voxelSpace + (voxelThicc/2.) 
	+ (i * (voxelThicc+voxelSpace));
      jpos[n] = -(fExpHall_y/2.) + voxelSpace + (voxelDim/2.) 
	+ (j * (voxelDim+voxelSpace)) 
	+ (i * ((voxelDim+voxelSpace)/2.));
      kpos[n] = -(fExpHall_z/2.) + voxelSpace + (voxelDim/2.) 
	+ (k * (voxelDim+voxelSpace))
	+ (i * ((voxelDim+voxelSpace)/2.));
    }

    if(doubleSiPM){
      if(std::abs(SiPMDifferenceRatio) == 1){
	jpos[n] += SiPMDifferenceRatio * 3.5;
      }
      else{
	jpos[n] += std::atanh(SiPMDifferenceRatio) / p0;
      }
    }
        
    xSum += ipos[n];
    xSqSum += (ipos[n]*ipos[n]);
        
    ySum += jpos[n];
    xySum += (ipos[n]*jpos[n]);
    yAve += jpos[n];
        
    zSum += kpos[n];
    xzSum += (ipos[n]*kpos[n]);
    zAve += kpos[n];
        
  }
    
  deltaTheta = (hits*xSqSum) - (xSum*xSum);
  bTheta = (1./deltaTheta) * ((xSqSum*ySum) - (xSum*xySum));
  mTheta = (1./deltaTheta) * ((hits*xySum) - (xSum*ySum));
  yAve /= hits;
    
  deltaPhi = (hits*xSqSum) - (xSum*xSum);
  bPhi = (1./deltaPhi) * ((xSqSum*zSum) - (xSum*xzSum));
  mPhi = (1./deltaPhi) * ((hits*xzSum) - (xSum*zSum));
  zAve /= hits;

  //std::cout << hits << "\n";

  for(n = 0; n < hits; ++n){
        
    numerTheta += pow((((mTheta*ipos[n])+bTheta) - jpos[n]), 2.);
    denomTheta += pow((jpos[n] - yAve), 2.);
        
    numerPhi += pow((((mPhi*ipos[n])+bPhi) - kpos[n]), 2.);
    denomPhi += pow((kpos[n] - zAve), 2.);
        
  }
    
  rSquaredTheta = 1 - (numerTheta/(denomTheta + RSQO));
  rSquaredPhi = 1 - (numerPhi/(denomPhi + RSQO));

  /*
    if(std::isnan(std::abs(rSquaredPhi))
    || std::isnan(std::abs(rSquaredTheta))){
    ++nanStatus;
    }*/

  if(!std::isnan(std::abs(mTheta))&&!std::isnan(std::abs(rSquaredTheta))
     &&!std::isnan(std::abs(mPhi))&&!std::isnan(std::abs(rSquaredPhi))){
    outData << run+1 << "," << hits << "," << energy[run] << "," << slopeY 
	    << "," << slopeZ << "," << mTheta << "," << mPhi << "," 
	    << rSquaredTheta << "," << rSquaredPhi << "\n";
    outData.close();

    
    return 1;

    
  }
  else{
    ++nanStatus;
    outData.close();

    return 0;

    
  }

  
}


int main(){


  int hits;
  std::string filename = "data.csv";

  std::cout << "\n";
  std::cout << " --- Linear Regression Analysis for VLXE --- \n";
  std::cout << "\n";

  if(readStoreData(filename)){
    if(fixDS){
      fixDataStructureLol();
    }
    //writeTermData();
    
    std::fstream outData;
    filename.insert(0, "out");
    outData.open(filename.c_str(), std::fstream::in | std::fstream::out 
		 | std::fstream::trunc);
    outData << "run,hits,energy,slopeY,slopeZ,estSlopeY,";
    outData << "estSlopeZ,RSqY,RSqZ\n";
    outData.close();

    for(int x = 0; x < ROWS; ++x){

      hits = 0;
      while(data[x][hits] != cubeNum*cubeNum*cubeNum){

	hits += 1;

      }

      if(hits > 1){
	double slopeY = momentum[x][1]/momentum[x][0];
	double slopeZ = momentum[x][2]/momentum[x][0];

	if(linReg(filename, x, hits, slopeY, slopeZ));
      }

    }

  }

  std::cout << "Total NANs: " << nanStatus << "\n";


  return 0;
}
