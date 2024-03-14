#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TAxis.h"

// Globals
const int dataSize = 10000;
double data[dataSize];
double timeData[dataSize];
bool savePicture = 0; //1 for saving pngs, 0 for no pngs saved
bool pTube = 0; // 0 for SiPM (positive landau), 1 for Ptube (negative)

struct peakData{
  double peak;
  double timing;
};


//written with the help of GPT4
bool readStoreData(const std::string& filename) {
  std::ifstream file(filename.c_str());
  std::string line;

  if (!file.is_open()) {
    std::cout << "Can't open file (bozo)" << std::endl;
    return false;
  }

  // Skip first 21 lines
  for (int i = 0; i < 12; ++i) {
    std::getline(file, line);
  }

  // Read and store time and data values
  for (int i = 0; i < dataSize && std::getline(file, line); ++i) {
    size_t pos = line.find(',');
    if (pos != std::string::npos) {
      std::string timeValue = line.substr(0, pos);
      std::string dataValue = line.substr(pos + 1);

      timeData[i] = std::atof(timeValue.c_str());
      data[i] = std::atof(dataValue.c_str());
      if(pTube){
	data[i] *= -1;
      }
    }
  }

  file.close();
  return 1;


}


void zeroData(){


  double zeroVal = 0.;
  int nMax = 250;
  for(int n = 0; n < nMax; ++n){ //find offset
    zeroVal += data[n];
  }
  zeroVal /= (double) nMax;

  for(int n = 0; n < dataSize; ++n){ //account for offset
    data[n] -= zeroVal;
  }


}


double timingEst(){


  double peak = data[0];
  double timing;
  int peaki;
  int i;

  int offset = 40; //set 0 for 1, 1 for 3, 2 for 5 points...
  int nPoints = (2*offset) + 1;
  double rollingAv;
  double peakAv;
  
  for(i = offset; i < dataSize-offset; ++i){

    rollingAv = 0;
    for(int j = i-offset; j <= i+offset; ++j){
      rollingAv += data[j];
    }
    rollingAv /= nPoints;

    if(rollingAv > peak){
      peak = data[i];
      peakAv = rollingAv;
      peaki = i;
    }

  }

  bool check = 1;
  i = peaki+1;
  while(check){

    rollingAv = 0;
    for(int j = i-offset; j <= i+offset; ++j){
      rollingAv += data[j];
    }
    rollingAv /= nPoints;

    i += 1;
    if((rollingAv < (peakAv/2.718281)) || (i == dataSize-offset-1)){
      check = 0;
    }

  }

  timing = timeData[i];
  
  return timing;


}

double landauFit(int fileNo, double timEst) {


  std::string outputFilename = std::to_string(fileNo) + ".png";

  std::vector<double> timeVec(timeData, timeData + sizeof(timeData) / sizeof(timeData[0]));
  std::vector<double> dataVec(data, data + sizeof(data) / sizeof(data[0]));

  int j = 0; //find j for initial time of peak
  while((timeVec[j] < -.15e-6) && j < dataSize){
    j += 1;
  }

  int k = 0; //find k corresponding to time esitmate, "won't" break 
  while((timeVec[k] < (timEst-.15e-6)) && (k < dataSize)){
    k +=1;
  }

  //double timeStep = timeVec[1] - timeVec[0];

  // Create a ROOT TCanvas to display the fit result
  TCanvas *c1 = new TCanvas("c1", "Landau Fit", 800, 600);

  // Create a TGraph to hold the data points
  TGraph *graph = new TGraph(timeVec.size(), &timeVec[0], &dataVec[0]);
 
  // Create a Landau distribution function with initial parameter values
  TF1 *landauFunc = new TF1("landauFunc","landau",
			    timeVec[j - 200], timeVec[k + 400]); //timeVec[j - 200]
 
  
  // Set some initial guesses for fit parameters
  //landauFunc->SetParameter(0, 1.0); // Amplitude
  //landauFunc->SetParameter(1, 2e-8); // Most probable value (MPV)
  //landauFunc->SetParameter(2, 1e-8); // Width
  //landauFunc->SetParameter(3, 0); //vertical offset
  //landauFunc->SetRange(timeVec[j - 10], timeVec[k + 20]); //range
  
  
  // Fit the Landau function to the data
  graph->Fit("landauFunc", "Q", "", timeVec[j - 200], timeVec[k + 400]); 
  //change to "R" to see things 0_0

  // Get the fit results
  double amp = landauFunc->GetParameter(0);
  double mpv = landauFunc->GetParameter(1);
  double width = landauFunc->GetParameter(2);

  // Calculate the maximum value of the fit
  double maxFitValue = landauFunc->GetMaximum();
  //double timeOfMaxFitValue = landauFunc->GetX(maxFitValue*.95,timeVec.front(),
  //timeVec.back());

  //double peakOverE = maxFitValue / TMath::E();
  //double timeOfPeakOverE = landauFunc->GetX(peakOverE, timeOfMaxFitValue,
  //timeVec.back());
  //double timeDiff = timeOfPeakOverE - timeOfMaxFitValue;

  // Print the fit results
  //std::cout << "Amplitude (A): " << amp << std::endl;
  //std::cout << "Most Probable Value (MPV): " << mpv << std::endl;
  //std::cout << "Width (sigma): " << width << std::endl;
  //std::cout << "Maximum Fit Value: " << maxFitValue << std::endl;
  //std::cout << "Time of peak: " << timeOfMaxFitValue << "\n";
  //std::cout << "Time at 1/e of Peak: " << timeOfPeakOverE << std::endl;
  //std::cout << "Time Difference (Peak to 1/e of Peak): " << timeDiff 
  //<< std::endl;

  if(savePicture){
  
    // Draw the data points and the fit function
    graph->SetTitle("Landau Fit");
    graph->GetXaxis()->SetTitle("Time");
    graph->GetYaxis()->SetTitle("Amplitude");
    graph->Draw("AP");
    landauFunc->Draw("SAME");

    // Save the canvas as a PDF
    c1->SaveAs(outputFilename.c_str());
  }
  
  // Clean up
  delete c1;
  delete graph;
  delete landauFunc;

  //peakData theReturn;
  //theReturn.peak = maxFitValue;
  //theReturn.timing = timeDiff;

  return maxFitValue;
}


int main(){
  std::ifstream inputFile("peakIn.txt");
  if(!inputFile){
    std::cout << "Error: Unable to open peakIn.txt file.\n";
    return 1;
  }

  std::string inputDirectory;
  std::string outputFilename;
  std::string line;
  while(std::getline(inputFile, line)){
    std::istringstream iss(line);
    std::string parameter;
    std::string value;
    if(std::getline(iss, parameter, ':') && std::getline(iss, value)){
      if(parameter == "Input (directory)"){
        inputDirectory = value;
      } 
      else if(parameter == "Output (file)"){
        outputFilename = value;
      }
    }
  }

  inputFile.close();

  if(inputDirectory.empty() || outputFilename.empty()){
    std::cout << "Error: Invalid input in peakIn.txt file.\n";
    return 1;
  }

  //int files = 10000;
  std::string filename;
  peakData thispeak;

  std::fstream outfile;
  outfile.open(outputFilename.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

  outfile << "event,peak,timing\n";

  int i = 0;
  int iMax = 1e4; //1e6
  bool check = 1, isOldData = 1; //for the loop and to check for data structure
  while(check && i < iMax){ // i = 0; i < files

    for(int channelIndex = 1; channelIndex < 5; ++channelIndex){
      std::ostringstream oss;
      oss << inputDirectory << "/C" << channelIndex << "autosave";
      if(isOldData){
	oss << std::setfill('0') << std::setw(5) << i;
      }
      else{
	oss << i;
      }
      oss << ".csv";
      filename = oss.str();

      check = readStoreData(filename);

      if(check){
	zeroData();
	double timEst = timingEst();
	double peak = landauFit(i, timEst);
      
	if(channelIndex == 1){
	  outfile << i+1 << ",";
	}
      
	outfile << peak <<  ",";

	if(channelIndex == 4){
	  outfile << "\n";
	}
      }
      else{
	if(!isOldData){
	  isOldData = 1;
	  std::cout << "Checking for old data type\n";
	  check = 1;
	}
      }
    }
    	i += 1;
  }

  outfile.close();

  return 0;
}
