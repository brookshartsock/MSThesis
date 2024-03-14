#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TF1.h>
#include <TPaveText.h>

// Declare the arrays globally
std::vector<double> energy;
std::vector<double> voltage;
std::vector<double> voltageE;


////written with the help of GPT4
void readStoreData(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }

    std::string line;
    std::getline(inputFile, line); // Skip the header line

    // Process each line of the CSV file
    while (std::getline(inputFile, line)) {
        std::stringstream ss(line);
        std::string field;

        // Skip the first field (HL)
        std::getline(ss, field, ',');

        // Get the energy value
        std::getline(ss, field, ',');
        if (!field.empty()) {
            energy.push_back(std::stod(field));
        }

        // Get the voltage value
        std::getline(ss, field, ',');
        if (!field.empty()) {
            voltage.push_back(std::stod(field));
        }

        // Get the voltageE value
        std::getline(ss, field, ',');
        if (!field.empty()) {
            voltageE.push_back(std::stod(field));
        }
    }
}


//written with the help of GPT4
void plotDataWithLinearFit(const std::string& outputFilename) {
    int nPoints = energy.size();

    // Convert std::vector to C-style arrays
    double* energyArr = new double[nPoints];
    double* voltageArr = voltage.data();
    double* voltageEArr = voltageE.data();

    // Convert energy to MeV and store in energyArr
    for (int i = 0; i < nPoints; ++i) {
        energyArr[i] = energy[i]; // Convert keV to MeV
    }

    // Create arrays for x and y errors
    double* xErrors = new double[nPoints];
    double* yErrors = voltageEArr;

    // Initialize xErrors to 0 (no error)
    for (int i = 0; i < nPoints; ++i) {
        xErrors[i] = 0.0;
    }

    // Create a TGraphErrors to hold the data points with error bars
    TGraphErrors *graph = new TGraphErrors(nPoints, energyArr, voltageArr, xErrors, yErrors);

    // Create a canvas to draw the graph
    TCanvas *canvas = new TCanvas("canvas", "Graph with Errors and Fit");

    //gPad -> DrawFrame(energyArr[0]-.1,voltageArr[0]-.01,energyArr[2]+.01, voltageArr[2]*1.15);

    // Draw the graph with error bars
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.0);
    graph->SetMarkerColor(kBlack);
    graph->SetTitle(outputFilename.c_str());
    graph->GetXaxis()->SetTitle("Energy (MeV)");
    graph->GetYaxis()->SetTitle("Voltage (V)");
    graph->SetMinimum(voltageArr[3]-.01); //min y
    graph->SetMaximum(voltageArr[2]+.1); //max y
    graph->Draw("AP"); // "A" for drawing the markers, "P" for drawing the error bars

    // Fit the data with a linear function (y = m*x + b)
    TF1 *fit = new TF1("fit","[1]*x + [0]");
    fit->SetParNames("Intercept","Slope");
    graph->Fit("fit", "Q"); // "Q" to suppress graphical output

    // Draw the fit function on top of the data points
    fit->SetLineColor(kRed);
    //fit->SetParNames("Intercept","Slope");
    fit->Draw("same");

    // Get fit parameters
    double slope = fit->GetParameter(1);
    double intercept = fit->GetParameter(0);
    double slopeError = fit->GetParError(1);
    double interceptError = fit->GetParError(0);

    gStyle->SetOptFit(1);

    // Save the plot as an image file
    std::string opFn = outputFilename + "CC.png";
    canvas->SaveAs(opFn.c_str());

    // Clean up memory
    delete[] energyArr;
    delete[] xErrors;
    delete graph;
    delete fit;
    delete canvas;
}


//written with the help of GPT4
int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename_prefix>" << std::endl;
        return 1;
    }

    std::string filenamePrefix = argv[1];
    std::string csvFilename = filenamePrefix + ".csv";

    readStoreData(csvFilename);

    // Modify the output filename
    std::string outputFilename = filenamePrefix + "CC.png";
    
    plotDataWithLinearFit(filenamePrefix);

    return 0;
}
