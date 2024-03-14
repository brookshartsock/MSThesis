#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TAxis.h"

//globals
const int dataSize = 1000;
double timeData[dataSize];
double ch1Data[dataSize];
double ch2Data[dataSize];
bool savePicture = 0; //1 for saving pngs, 0 for no pngs saved
bool energyCalibration = 0;

std::string dataString;
int peakIndexEst; //index no. of where the peak of ch2 is
int peakRange = 75;


//written with the help of GPT4
bool readStoreData(const std::string& filename) {

  
  std::ifstream file(filename.c_str());
  std::string line;

  if (!file.is_open()) {
    std::cout << "Can't open file (bozo)" << std::endl;
    return false;
  }

  // Skip first 21 lines
  for (int i = 0; i < 21; ++i) {
    std::getline(file, line);
  }

  // Read and store time and data values
  for (int i = 0; i < dataSize && std::getline(file, line); ++i) {
    size_t pos1 = line.find(',');
    size_t pos2 = line.find(',', pos1 + 1);

    if (pos1 != std::string::npos && pos2 != std::string::npos) {
      std::string timeValue = line.substr(0, pos1);
      std::string ch1Value = line.substr(pos1 + 1, pos2 - pos1 - 1);
      std::string ch2Value = line.substr(pos2 + 1);

      timeData[i] = std::atof(timeValue.c_str());
      ch1Data[i] = std::atof(ch1Value.c_str());
      ch2Data[i] = std::atof(ch2Value.c_str());
    }
  }

  file.close();

  
  return 1;
}


void coutData(){


  for(int i = 0; i < dataSize; ++i){

    std::cout << timeData[i] << "," << ch1Data[i] << "," << ch2Data[i]
	      << "\n";
    
  }
  

}


void zeroData(){


  double zeroVal = 0.;
  int nMax = 25;

  //zero ch1 data
  for(int n = 0; n < nMax; ++n){ //find offset
    zeroVal += ch1Data[n];
  }
  zeroVal /= (double) nMax;

  for(int n = 0; n < dataSize; ++n){ //account for offset
    ch1Data[n] -= zeroVal;
  }

  //zero ch2 data
  zeroVal = 0.;
  for(int n = 0; n < nMax; ++n){ //find offset
    zeroVal += ch2Data[n];
  }
  zeroVal /= (double) nMax;

  for(int n = 0; n < dataSize; ++n){ //account for offset
    ch2Data[n] -= zeroVal;
  }

  
}


double checkCh2Data(){


  double peak = ch2Data[0];
  int peaki;
  int i;
  peakIndexEst = 0;

  int offset = 4; //set 0 for 1, 1 for 3, 2 for 5 points...
  int nPoints = (2*offset) + 1;
  double rollingAv;
  double peakAv;
  
  for(i = offset; i < dataSize-offset; ++i){

    rollingAv = 0;
    for(int j = i-offset; j <= i+offset; ++j){
      rollingAv += ch2Data[j];
    }
    rollingAv /= nPoints;

    if(rollingAv > peak){
      peak = ch2Data[i];
      peakAv = rollingAv;
      peaki = i;
    }

  }

  peakIndexEst = peaki;
  //std::cout << peakIndexEst << "\n";

  
  return peak;
}


double lFitCh2Data(int fileNo, std::string outfileName){


  std::string outputFilename = std::to_string(fileNo) + "-2x.png";

  std::vector<double> timeVec(timeData, timeData + sizeof(timeData)
			      / sizeof(timeData[0]));
  std::vector<double> dataVec(ch2Data, ch2Data + sizeof(ch2Data)
			      / sizeof(ch2Data[0]));

  TCanvas *c1 = new TCanvas("c1", "Landau Fit", 800, 600);

  // Create a TGraph to hold the data points
  TGraph *graph = new TGraph(timeVec.size(), &timeVec[0], &dataVec[0]);
 
  // Create a Landau distribution function with initial parameter values
  TF1 *landauFunc = new TF1("landauFunc","landau",
			    timeVec[peakIndexEst - peakRange],
			    timeVec[peakIndexEst + peakRange]);

  //landauFunc -> SetParLimits(1, timeVec[peakIndexEst - peakRange],
  //timeData[peakIndexEst + peakRange]);

  graph->Fit("landauFunc", "Q", "", timeVec[peakIndexEst - peakRange],
	     timeData[peakIndexEst + peakRange]);

  //std::cout << timeData[peakIndexEst + peakRange]
  //-timeData[peakIndexEst - peakRange] << "\n";

  //std::cout << peakIndexEst << "," << peakRange << "\n";
  //change to "R" to see things 0_0

  // Get the fit results
  double amp = landauFunc->GetParameter(0);
  double mpv = landauFunc->GetParameter(1);
  double width = landauFunc->GetParameter(2);

  // Calculate the maximum value of the fit
  double maxFitValue = landauFunc->GetMaximum();
  double timeOfMaxFitValue = landauFunc->GetX(maxFitValue*.95,timeVec.front(),
					      timeVec.back());


  // Print the fit results
  //std::cout << "Amplitude (A): " << amp << std::endl;
  //std::cout << "Most Probable Value (MPV): " << mpv << std::endl;
  //std::cout << "Width (sigma): " << width << std::endl;
  //std::cout << "Maximum Fit Value: " << maxFitValue << std::endl;
  //std::cout << "Time of peak: " << timeOfMaxFitValue << "\n";
  //std::cout << "Time at 1/e of Peak: " << timeOfPeakOverE << std::endl;
  //std::cout << "Time Difference (Peak to 1/e of Peak): " << timeDiff 
  //<< std::endl;

  //std::cout << savePicture << "\n";

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
  //else{
  //}

  if(energyCalibration){
    int stringLength = outfileName.size();
    outfileName.insert(stringLength-4, "2");//add a 2 before the .csv

    std::fstream outfile;
    outfile.open(outfileName.c_str(), std::fstream::in | std::fstream::out
		 | std::fstream::app);

    outfile << fileNo << "," << maxFitValue << "," <<  1 << "\n";

    outfile.close();  
  }

  else{
    //so there is "maxFitValue" for the peak value at timeOfMaxFitValue seconds
    double ch2Slope = .1509; //these are for 2x
    double ch2SlopeErr = 8.608e-5;
    double ch2Intercept = -.00001736;
    double ch2InterceptErr = .002509;

    double partial_peak = 1.0 / ch2Slope;
    double partial_intercept = -1.0 / ch2Slope;
    double partial_slope = -(maxFitValue - ch2Intercept) / (ch2Slope * ch2Slope);

    double ch2Edep = (maxFitValue - ch2Intercept) / ch2Slope;
    double ch2EdepErr =
      std::sqrt( std::pow(partial_peak * ch2SlopeErr, 2)
		 + std::pow(partial_intercept * ch2InterceptErr, 2)
		 + std::pow(partial_slope * ch2SlopeErr, 2) );

    dataString += std::to_string(ch2Edep) + ","
      + std::to_string(ch2EdepErr) + ",";
  }

  delete c1;
  delete graph;
  delete landauFunc;
  

  return timeOfMaxFitValue;
}


double lFitCh1Data(int fileNo, std::string outfileName){


  std::string outputFilename = std::to_string(fileNo) + "-1x.png";

  std::vector<double> timeVec(timeData, timeData + sizeof(timeData)
			      / sizeof(timeData[0]));
  std::vector<double> dataVec(ch1Data, ch1Data + sizeof(ch1Data)
			      / sizeof(ch1Data[0]));

  TCanvas *c1 = new TCanvas("c1", "Landau Fit", 800, 600);

  // Create a TGraph to hold the data points
  TGraph *graph = new TGraph(timeVec.size(), &timeVec[0], &dataVec[0]);
 
  // Create a Landau distribution function with initial parameter values
  TF1 *landauFunc = new TF1("landauFunc","landau",
			    timeVec[peakIndexEst - peakRange],
			    timeVec[peakIndexEst + peakRange]);

  graph->Fit("landauFunc", "Q", "", timeVec[peakIndexEst - peakRange],
	     timeVec[peakIndexEst + peakRange]); 
  //change to "R" to see things 0_0

  // Get the fit results
  double amp = landauFunc->GetParameter(0);
  double mpv = landauFunc->GetParameter(1);
  double width = landauFunc->GetParameter(2);

  // Calculate the maximum value of the fit
  double maxFitValue = landauFunc->GetMaximum();
  double timeOfMaxFitValue = landauFunc->GetX(maxFitValue*.95,timeVec.front(),
					      timeVec.back());

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
  //else{
  //}

  if(energyCalibration){
    int stringLength = outfileName.size();
    outfileName.insert(stringLength-4, "1");//add a 1 before the .csv

    std::fstream outfile;
    outfile.open(outfileName.c_str(), std::fstream::in | std::fstream::out
		 | std::fstream::app);

    outfile << fileNo << "," << maxFitValue << "," <<  1 << "\n";

    outfile.close();  
  }

  else{
    //so there is "peak" for the peak value at array index "peaki"
    double ch1Slope = .1471; //these are for 1x
    double ch1SlopeErr = 2.99e-6;
    double ch1Intercept = .0001227;
    double ch1InterceptErr = .002399;

    double partial_peak = 1.0 / ch1Slope;
    double partial_intercept = -1.0 / ch1Slope;
    double partial_slope = -(maxFitValue - ch1Intercept) / (ch1Slope * ch1Slope);

    double ch1Edep = (maxFitValue - ch1Intercept) / ch1Slope;
    double ch1EdepErr =
      std::sqrt( std::pow(partial_peak * ch1SlopeErr, 2)
		 + std::pow(partial_intercept * ch1InterceptErr, 2)
		 + std::pow(partial_slope * ch1SlopeErr, 2) );

    dataString.insert(0, std::to_string(ch1EdepErr) + ",");
    dataString.insert(0, std::to_string(ch1Edep) + ",");
  
    //dataString += std::to_string(ch1Edep) + ","
    //+ std::to_string(ch1EdepErr) + ",";
  }
  
  delete c1;
  delete graph;
  delete landauFunc;
  

  return timeOfMaxFitValue;
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

  //if(!energyCalibration){
  //int files = 10000;
  std::string filename;
  std::fstream outfile;
  if(!energyCalibration){
    outfile.open(outputFilename.c_str(), std::fstream::in | std::fstream::out
		 | std::fstream::trunc);

    outfile << "tek,ch1Edep_MeV,ch1EdepErr_MeV,ch2Edep_MeV,ch2EdepErr_MeV,";
    outfile << "time_CH1-CH2_ns\n";
  }

  int i = 0;
  int iMax = 1e6; //1e6, safety measure
  //double ch2Threshold = .0476;//.0476
  bool check = 1; //for the loop
  while(check && (i < iMax)){ // i = 0; i < files
    
    std::ostringstream oss;
    oss << inputDirectory << "/tek";
    oss << i;    
    oss << "ALL.csv";
    filename = oss.str();

    check = readStoreData(filename);
    dataString = std::to_string(i) + ',';

    if(check){ //if data can be read...
      
      //coutData();

      //savePicture = 0;

      checkCh2Data();
      //if(checkCh2Data() > ch2Threshold){
      zeroData();
      double ch2PeakTime = lFitCh2Data(i,outputFilename)*1e9;
      //std::cout << ch2PeakTime << "\n";
      double ch1PeakTime = lFitCh1Data(i,outputFilename)*1e9;
      //std::cout << ch1PeakTime << "\n";

      dataString += std::to_string(ch1PeakTime-ch2PeakTime) + "\n";

      if(!energyCalibration){
	outfile << dataString;
      }
      //}

    }

    else{
      std::cout << "Can't open file " << filename << "\n";
    }

    /*
      if((i % 100) == 0){
      std::cout << "Finished file index " << i << "\n";
      }*/

    ++i;
    
  }

  if(!energyCalibration){
    outfile.close();
  }

  
  return 0;
}
