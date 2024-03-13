Preliminary Simulations

The Neutron-Electron Angle plot was done using nuEAngle.cc

The average range plots for e- and gamams in GAGG where done with a modified TestEm1 - this is the only sim in Geant4.10.7

The 175 Gamma Range in GAGG plot was done from trueRange which I wrote based from B1. It doesn't work for e- btw

The DEP plots where made from the singleVoxel based on B1

The plate style voxel detector is from VLXe which is based on LXe

The one made of 1x1x10 pillars is hodoscope, again based on LXe

The anti-nu detector is just based on the hodoscope, but the rotation is turned off

Neutron-Electron Angle and average range plots were made on RStudio with something like this:

x = data$Energy..MeV.
y = data$Range..mm.
z = data$rms..mm.

plot(x, y, pch = 20, main = "Average Range for Electrons GAGG:Ce",
     xlab = "Energy (MeV)", ylab = "Range (mm)", ylim = c(0, 18), log = 'x')
arrows(x, y + z, x, y - z, length = 0.02, angle = 90, code = 3)
