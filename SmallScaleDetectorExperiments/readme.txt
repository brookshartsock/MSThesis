Small-scale Detector Experiments

Lab data was trimmed using DSDR.cc which is based on landauFitPeaks.cc

Further trimming done on RStudio:
fLabData = labData[(labData$ch2Edep.MeV. + labData$ch1Edep.MeV. - sqrt(labData$ch2EdepErr.MeV.^2 + labData$ch1EdepErr.MeV.^2)) < max(fSimData$onexEdep + fSimData$twoxEdep) & (labData$tek. != "NA") & (labData$time.ns. > -10) & (labData$time.ns. < 46),]

data <- data[apply(data, 1, function(row) length(na.omit(row)) == 4), ]

Geant sim is betaTest to mimic the lab tests

Both gamma tests where analyzed with gammaTests.cc

the double SiPM crystal was analyzed with DSDR.cc again based on landauFitPeaks

Lab data for DSDR came from doubleSiPM and was plotted and fit with DSDRFit.cc
