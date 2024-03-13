//rootn't
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

//root
#include <TH2.h>
#include <TCanvas.h>
#include <TPaveStats.h>


//globals
std::vector<double> energy, slopeY, slopeZ, estSlopeY, estSlopeZ;
int rows;


bool readStoreData(const std::string& filename) {


  std::ifstream inputFile(filename.c_str());

  if(!inputFile.is_open()){
    std::cerr << "Error opening file: " << filename << std::endl;
    return 0;
  }

  std::string line;

  std::getline(inputFile, line);//skip the header line

  while(std::getline(inputFile, line)){

    std::istringstream ss(line);
    std::string token;
    int columnCount = 0;

    while(std::getline(ss, token, ',')){

      ++columnCount;
      if (columnCount == 3 || columnCount == 4 || columnCount == 5 
	  || columnCount == 6 || columnCount == 7) {
	std::istringstream convert(token);
	double value;

	if (convert >> value) {
	  if (columnCount == 3) energy.push_back(value);
	  else if (columnCount == 4) slopeY.push_back(value);
	  else if (columnCount == 5) slopeZ.push_back(value);
	  else if (columnCount == 6) estSlopeY.push_back(value);
	  else if (columnCount == 7) estSlopeZ.push_back(value);
	}
      }

    }

  }
  rows = energy.size();


  return 1;
}


void displayData(){


  for(int i = 0; i < 10; ++i){
    std::cout << energy[i] << "," << slopeY[i] << "," << slopeZ[i] << ","
	      << estSlopeY[i] << "," << estSlopeZ[i] << "\n"; 
  }


}


void makeHistogram(){


  double minSlopeY = slopeY[0];
  double maxSlopeY = slopeY[0];
  double minSlopeZ = slopeZ[0];
  double maxSlopeZ = slopeZ[0];

  for(int i = 0; i < rows; ++i){
    if(slopeY[i] < minSlopeY){
      minSlopeY = slopeY[i];
    }
    if (slopeY[i] > maxSlopeY) {
      maxSlopeY = slopeY[i];
    }
        
    if (slopeZ[i] < minSlopeZ) {
      minSlopeZ = slopeZ[i];
    }
    if (slopeZ[i] > maxSlopeZ) {
      maxSlopeZ = slopeZ[i];
    }
  }

  double axisRange = 30.;
  int bins = 128;

  for(int iNo = 0; iNo < 6; ++iNo){

    double energyMin;
    double energyMax;

    if(iNo == 0){
      energyMin = 0.;
      energyMax = 10.;
    }
    else if(iNo == 1){
      energyMin = 0.;
      energyMax = 2.;
    }
    else if(iNo == 2){
      energyMin = 2.;
      energyMax = 4.;
    }
    else if(iNo == 3){
      energyMin = 4.;
      energyMax = 6.;
    }
    else if (iNo == 4){
      energyMin = 6.;
      energyMax = 8.;
    }
    else{
      energyMin = 8.;
      energyMax = 10.;
    }

    std::string mainTitle = "Linear Regression Estimation (" 
      + std::to_string((int) energyMin) 
      + " < E < " + std::to_string((int) energyMax) + " MeV)";

    TH2D *histogram = new TH2D("histogram",
			       mainTitle.c_str(),
			       bins, -axisRange / 2, axisRange / 2, bins,
			       -axisRange / 2, axisRange / 2);

    // Create a 2D histogram to keep track of the number of events in each bin
    TH2D *eventCountHistogram = new TH2D("eventCountHistogram", "Event Count",
					 bins, -axisRange / 2, axisRange / 2, bins, -axisRange / 2, axisRange / 2);

    for (int i = 0; i < rows; ++i) {
      if (energy[i] > energyMin && energy[i] < energyMax) {
	int binX = histogram->GetXaxis()->FindBin(slopeY[i]);
	int binY = histogram->GetYaxis()->FindBin(slopeZ[i]);

	// Fill the event count histogram to keep track of the number of events in each bin
	eventCountHistogram->Fill(slopeY[i], slopeZ[i]);

	// Fill the histogram with the average value (value)
	if (eventCountHistogram->GetBinContent(binX, binY) > 0) {
	  double value = histogram->GetBinContent(binX, binY);
	  value += ((estSlopeY[i] - slopeY[i]) + (estSlopeZ[i] - slopeZ[i])) / 2.;
	  histogram->SetBinContent(binX, binY, value / eventCountHistogram->GetBinContent(binX, binY));
	} //else {
	//histogram->SetBinContent(binX, binY, ((estSlopeY[i] - slopeY[i]) + (estSlopeZ[i] - slopeZ[i])) / 2.);
	//}
      }
    }

    // Create a canvas to display the histogram
    TCanvas *canvas = new TCanvas("canvas", "Canvas", 1080, 1080);
    canvas->SetLeftMargin(.15);
    canvas->SetRightMargin(.20);
    canvas->SetBottomMargin(.1);
    canvas->SetTopMargin(.15);
    canvas->cd();

    // Draw the histogram
    histogram->Draw("colz"); // "colz" option creates a heatmap-like representation
    canvas->Update(); // Update canvas to make sure statbox is generated

    // Retrieve the TPaveStats
    TPaveStats *stats = (TPaveStats *)histogram->GetListOfFunctions()->FindObject("stats");

    if (stats) {
      // Set the position and size of the statbox
      stats->SetX1NDC(0.15); // Set the x position (0.1 = 10% from the left edge)
      stats->SetX2NDC(0.45); // Set the x position for the width
      stats->SetY1NDC(0.7); // Set the y position (0.7 = 70% from the bottom edge)
      stats->SetY2NDC(0.9); // Set the y position for the height
      stats->SetOptStat(1); // Show only entries in the statbox

      // Update the canvas
      canvas->Modified();
    } else {
      std::cerr << "Statbox not found." << std::endl;
    }

    // Optional: Add labels and other decorations
    histogram->GetXaxis()->SetTitle("Slope Y");
    histogram->GetYaxis()->SetTitle("Slope Z");
    histogram->GetZaxis()->SetTitle("Mean difference (estimated - actual)");
    histogram->GetZaxis()->SetTitleOffset(1.3);

    

    std::string outfilename = "linRegHist" + std::to_string((int) energyMin) 
      + "-" + std::to_string((int) energyMax) + ".pdf";
    canvas->SaveAs(outfilename.c_str());

    //std::cout << "\n";
    //std::cout << "scp bhartsock@charpak.physics.wichita.edu:/home/bhartsock/Desktop/root/2DHists/" << outfilename << " ~/Desktop/\n";
    //std::cout << "\n";

    // Clean up memory
    delete histogram;
    delete canvas;
    delete eventCountHistogram;

  }


}


int main(){


  std::string filename = "outdata.csv";

  if(readStoreData(filename)){
    //displayData();
    makeHistogram();
  }


  return 0;
}
