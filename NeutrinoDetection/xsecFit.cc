#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include <TStyle.h>

double fitf(double *x,double *par) {

  
  double pol3 = par[4]*x[0]*x[0]*x[0] + par[3]*x[0]*x[0] + par[2]*x[0] + par[1];
  
  if(x[0] > par[5]){
    return pol3;
  }
  else{
    double A = pol3/pow(par[5],par[0]);
    return A*pow(x[0],par[0]);
  }
  
  
}


int main() {
  // File path
  std::string filePath = "IBDXSec.csv";

  // Vectors to store data
  std::vector<double> energy;
  std::vector<double> xsec;

  // Open file
  std::ifstream file(filePath);

  // Check if file is open
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filePath << std::endl;
    return 1;
  }

  // Read data from file
  std::string line;
  getline(file, line); // Skip header
  while (getline(file, line)) {
    std::stringstream ss(line);
    std::string value;
    std::vector<std::string> values;

    while (getline(ss, value, ',')) {
      values.push_back(value);
    }

    energy.push_back(std::stod(values[0]));
    xsec.push_back(std::stod(values[1]));
  }

  // Close file
  file.close();

  // Create a TGraph
  TGraph *graph = new TGraph(energy.size(), &energy[0], &xsec[0]);

  // Set graph title and axis labels
  graph->SetTitle("IBD Cross Section;Energy (MeV);Cross Section (10^{-41} cm^{2})");

  // Create a canvas
  TCanvas *canvas = new TCanvas("canvas", "IBD Cross Section", 800, 600);

  gStyle -> SetOptFit();

  // Fit a second degree polynomial
  //TF1 *fitFunc = new TF1("fitFunc", "expo", 0, 2.);

  //fitFunc->SetParLimits(0,.004, 1);
  //fitFunc->SetParLimits(1,5,10);
  
  //fitFunc->SetParLimits(3,0,100);
  //fitFunc->SetParLimits(0,0,100);
  //fitFunc->SetParLimits(4,0,10);

  //fitFunc->SetParameter(0,1);
  //fitFunc->SetParameter(1,1);
  //fitFunc->SetParameter(2,1);

  /*
  // Fit the first range with a third-degree polynomial
  TF1 *func = new TF1("fit",fitf,0,12,6);

  func -> SetParLimits(5, 1, 12);
  func -> SetParameter(5, 2);*/

  TF1 *fit = new TF1("fit","[0]*pow(x-1.806,[1])", 1.9, 12);

  graph->Fit("fit");
  graph->Draw("AP");

  fit->Draw("SAME");
  graph->Draw("P SAME");

  graph->SetMarkerStyle(104);

  graph->SetMinimum(1e-4); //min y
  graph->SetMaximum(1.1); //max y

  //func -> SetParLimits(5, 1.9, 12);
  //func -> SetParameter(5, 6);


  //graph->SetStats(1);

  canvas->SetLogy();

  //graph->Draw("SAME");
  //canvas->Update();

  // Save the canvas as a PNG file
  canvas->SaveAs("IBDXSec_Fit.png");

  // Clean up
  delete graph;
  delete canvas;

  return 0;
}
