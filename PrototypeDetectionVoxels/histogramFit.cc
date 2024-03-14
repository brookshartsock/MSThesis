#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>


//written with the help of GPT4
void createHistogram(const std::string& csvFilename, double binCorr, bool fitting, double amp, double mean, double std, double exp1, double exp2, double fitMin, double fitMax) {
  // Load the data from the CSV file into a TTree
  TFile *file = new TFile((csvFilename + ".root").c_str(), "RECREATE");
  TTree *tree = new TTree("tree", "Tree with peak data");

  // Variables to hold the data from the TTree
  Int_t event;
  Double_t peak, timing;

  // Set the branch address to link the variables with the TTree data
  tree->Branch("event", &event, "event/I");
  tree->Branch("peak", &peak, "peak/D");
  tree->Branch("timing", &timing, "timing/D");

  // Open the CSV file for reading
  std::ifstream inputFile((csvFilename + "_peakSiPM.csv").c_str());
  if (!inputFile.is_open()) {
    std::cerr << "Error: Could not open the input file." << std::endl;
    return;
  }

  // Read and parse the data from the CSV file line by line
  std::string line;
  bool firstLine = true;
  while (std::getline(inputFile, line)) {
    // Skip the header line (if present)
    if (firstLine) {
      firstLine = false;
      continue;
    }

    // Parse the comma-separated values
    std::istringstream ss(line);
    std::string eventStr, peakStr, timingStr;
    std::getline(ss, eventStr, ',');
    std::getline(ss, peakStr, ',');
    std::getline(ss, timingStr, ',');

    // Convert the parsed strings to appropriate data types
    event = std::stoi(eventStr);
    peak = std::stod(peakStr);
    timing = std::stod(timingStr);

    // Fill the TTree with the parsed data
    tree->Fill();
  }

  // Close the input file
  inputFile.close();

  // Get the number of entries in the TTree
  Long64_t nEntries = tree->GetEntries();

  // Create a vector to hold the "peak" values
  std::vector<double> peakValues;
  peakValues.reserve(nEntries);

  // Loop over the entries in the TTree and store the "peak" values
  for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
    tree->GetEntry(iEntry);

    // Ignore "inf" peak values
    if (!std::isinf(peak)) {
      peakValues.push_back(peak);
    }

  }

  // Sort the "peak" values to find the minimum and maximum
  std::sort(peakValues.begin(), peakValues.end());
  Double_t minPeak = peakValues.front();
  Double_t maxPeak = peakValues.back();

  // Calculate the bin size based on differences between consecutive peak values
  Double_t minDiff = std::numeric_limits<Double_t>::max();
  for (size_t i = 1; i < peakValues.size(); ++i) {
    Double_t diff = std::abs(peakValues[i] - peakValues[i - 1]);
    if (diff > 0 && diff < minDiff) {
      minDiff = diff;
    }
  }
  Double_t binSize = minDiff;

  // Apply the user-provided correction term
  binSize /= binCorr;

  // Determine the number of bins based on the bin size
  Int_t nBins = static_cast<Int_t>(TMath::CeilNint((maxPeak - minPeak) / binSize));

  // Create a histogram with calculated binning and range
  TH1D *histogram = new TH1D("histogram", (csvFilename + " - Peak Distribution; Voltage (V); Frequency").c_str(), nBins, minPeak, maxPeak);

  // Loop over the entries in the TTree and fill the histogram
  for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
    tree->GetEntry(iEntry);
    histogram->Fill(peak);
  }


  TF1 *fit = new TF1("fit", "gaus(0) + expo(3)", fitMin, fitMax);
  //-----------------------------------------
  if(fitting){

    fit->SetParameters(amp, mean, std, exp1, exp2); //map, mean, std, exp1, exp2
    histogram->Fit("fit", "", "", fitMin, fitMax);

    gStyle -> SetOptFit();
  }
  //-----------------------------------------

  // Create a canvas to draw the histogram
  TCanvas *canvas = new TCanvas("canvas", "Histogram Canvas");

  // Draw the histogram with axis labels and title
  histogram->GetXaxis()->SetTitle("Voltage (V)");
  histogram->GetYaxis()->SetTitle("Frequency");
  histogram->SetTitle(csvFilename.c_str());
  histogram->Draw();

  if(fitting){
    gStyle -> SetOptFit();
    fit->Draw("same");
  }

  // Save the histogram to an image file
  canvas->SaveAs((csvFilename + "Hist.png").c_str());

  // Clean up memory
  delete canvas;
  delete histogram;
  delete tree;
  delete file;
}

int main() {


  bool fitting = 0; //set false (no fit)
  double amp, mean, std, exp1, exp2;
  double fitMin, fitMax;

  std::string filename;
  double binCorr = 1.0;

  std::cout << "Enter the first 4/5 letters of the CSV filename: ";
  std::cin >> filename;

  if(!fitting){
    amp = 0.;
    mean = 0.;
    std = 0.;
    exp1 = 0.;
    exp2 = 0.;
    fitMin = 0.;
    fitMax = 0.; //set values (might get mad if not set but passed...)

    createHistogram(filename.c_str(), binCorr, fitting, amp, mean, std, 
		    exp1, exp2, fitMin, fitMax); //yikes

    fitting = 1; //all function calls following will attempt to fit

    std::cout << "\n";
    std::cout << "scp bhartsock@charpak.physics.wichita.edu:/home/bhartsock/Desktop/root/SiPM/" << filename << "Hist.png ~/Desktop\n";
    std::cout << "\n";

  }

  while(fitting){

    std::cout << "Bin correction term (1.): ";
    std::cin >> binCorr;

    std::cout << "Fit ampltude: ";
    std::cin >> amp;

    std::cout << "Fit mean: ";
    std::cin >> mean;

    std::cout << "Fit std: ";
    std::cin >> std;
  
    std::cout << "Fit exp1 (5.): ";
    std::cin >> exp1;

    std::cout << "Fit exp2 (-10.): ";
    std::cin >> exp2;

    std::cout << "Fit min: ";
    std::cin >> fitMin;

    std::cout << "Fit max: ";
    std::cin >> fitMax;

    createHistogram(filename.c_str(), binCorr, fitting, amp, mean, std, 
		    exp1, exp2, fitMin, fitMax); //yikes

    std::cout << "Try again (1/0)? ";
    std::cin >> fitting;

  }

  return 0;
}
