#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>

// ROOT includes
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TStyle.h>

//globals
std::vector<double> ipos, energy, SiPM0, SiPM1;
int rows;

const int bins = 20;
double meanValue[bins] = {0};
double sdValue[bins] = {0};
int eventsPerBin[bins] = {0};


bool readStoreData(std::string filename){

  
  std::ifstream inputFile(filename.c_str());

  if(!inputFile.is_open()){
    std::cerr << "Error opening file: " << filename << std::endl;
    return false;
  }

  std::string line;
  std::getline(inputFile, line); // skip the header line

  while(std::getline(inputFile, line)){
      
    std::istringstream iss(line);
    std::string iposStr, energyStr, SiPM0Str, SiPM1Str;

    if(!(std::getline(iss, iposStr, ',') &&
	 std::getline(iss, energyStr, ',') &&
	 std::getline(iss, SiPM0Str, ',') &&
	 std::getline(iss, SiPM1Str))){
      std::cerr << "Error reading line: " << line << std::endl;
      return false;
    }

    // Convert string values to double
    double iposVal = std::stod(iposStr);
    double energyVal = std::stod(energyStr);
    double SiPM0Val = std::stod(SiPM0Str);
    double SiPM1Val = std::stod(SiPM1Str);

    // Store values in vectors
    ipos.push_back(iposVal);
    energy.push_back(energyVal);
    SiPM0.push_back(SiPM0Val);
    SiPM1.push_back(SiPM1Val);

    rows++;
  }
    

  return true;
}

void displayData(){

  
  for(int i = 0; i < rows; ++i){
      
    std::cout << ipos[i] << "," << energy[i] << "," << SiPM0[i] << ","
	      << SiPM1[i] << "\n";
	
  }

    
}


void trimData(){


  for(int i = 0; i < rows; ++i){

    if(ipos[i] == 10){
      ipos[i] -= .1; //bitch
    }
    
    int bin = std::floor(ipos[i]) + 10;
    //std::cout << bin << "\n";
    if((SiPM1[i] + SiPM0[i]) != 0){
    double value = (SiPM1[i] - SiPM0[i])/(SiPM1[i] + SiPM0[i]);
    
    meanValue[bin] += value;
    sdValue[bin] += value * value;
    ++eventsPerBin[bin];
    //std::cout << eventsPerBin[bin] << "\n";
    }
    
  }

  for(int i = 0; i < bins; ++i){

    //std::cout << meanValue[i] << "\n";
    meanValue[i] /= (double)eventsPerBin[i];
    sdValue[i] = std::sqrt((sdValue[i] / (double)eventsPerBin[i])
			   - (meanValue[i] * meanValue[i]));
    
  }


}


void outTrimData(){

  for(int i = 0; i < bins; ++i){
    std::cout << i - 9.5 << "," << meanValue[i] << "," << sdValue[i] << "\n";
  }

}


void makeScatterPlot() {
    TCanvas *canvas = new TCanvas("canvas", "Scatter Plot", 800, 600);

    // Create a TGraphErrors for the scatter plot
    TGraphErrors *scatterPlot = new TGraphErrors(bins);

    for (int i = 0; i < bins; ++i) {
        double xPos = i - 9.5; // Adjusting the x-position
        scatterPlot->SetPoint(i, xPos, meanValue[i]);
        scatterPlot->SetPointError(i, 0, sdValue[i]);
    }

    // Set axis labels
    scatterPlot->GetXaxis()->SetTitle("Position (mm)");
    scatterPlot->GetYaxis()->SetTitle("Mean (SiPM1-SiPM0)/(SiPM1 + SiPM0)");
    scatterPlot->SetTitle("Double SiPM NR - 20x6.65x7.05 mm Crystal");
    scatterPlot->GetYaxis()->SetTitleOffset(1.2);

    // Set marker style and color for the original points
    scatterPlot->SetMarkerStyle(20); // Adjust as needed
    scatterPlot->SetMarkerColor(kBlue); // Adjust as needed
    scatterPlot->SetLineColor(kBlue); // Line color for clear error bars

    // Fit a tanh function to the scatter plot
    TF1 *tanhFit = new TF1("tanhFit", "TMath::TanH([0] * x )", -8.5, 8.5);
    tanhFit->SetParameters(.16); // Set initial parameter values

    // Set fit parameter name
    tanhFit->SetParName(0, "p0");

    // Fit the tanh function to the data
    scatterPlot->Fit("tanhFit", "R");

    // Draw the scatter plot with the original points and clear error bars
    canvas->cd();
    
    gStyle->SetOptStat(1); // Disable default stats box
    gStyle->SetOptFit(1); // Disable default fit box

    /*
    // Add a new stats box manually
    TPaveStats *stats = new TPaveStats(0.65, 0.15, 0.85, 0.35, "NDC");
    stats->SetName("stats");
    stats->SetOptStat(111); // Show mean and RMS
    stats->SetOptFit(11);   // Show fit parameters
    stats->SetFitFormat("6.4g"); // Set fit parameter format
    scatterPlot->SetStats(1); // Enable stats box
    canvas->GetListOfPrimitives()->Add(stats);*/

    gStyle->SetStatX(0.9); // X position of stats box
    gStyle->SetStatY(0.25); // Y position of stats box
    gStyle->SetStatW(0.2);  // Width of stats box
    gStyle->SetStatH(0.2);  // Height of stats box
    
    scatterPlot->Draw("AP"); // "AP" option draws the points with axes and uncertainties
    tanhFit->Draw("same");

    canvas->SaveAs("doubleSiPMNR20mmPlot.pdf");
}


int main(){

  
  std::string filename = "doubleSiPMData.csv";//doubleSiPMData.csv

  if(readStoreData(filename)){
    //rows = 10;
    //displayData();
    trimData();
    outTrimData();
    makeScatterPlot();
  }
    

  return 0;
}
