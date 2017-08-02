# TCSPC-image-simulation
This code provides an easy way to simulate TCSPC image data which includes effect of the IRF, background after-pulsing and laser repetition rate. Images are saved as OME.tiff using Bioformat.

The method uses a Monte-Carlo simulation of photon arrival times based on probability density functions of emission (exponential) and excitation (Gaussian). 
The OME.tiff save was kindly provided by Ian Munro <i.munro@imperial.ac.uk> from the Photonics Group at Imperial College, London and Sebastien Besson <s.besson@dundee.ac.uk> from Dundee, who write the BioFormat MATLAB Utilities.

You need to download the Bio-formats Matlab plugin from 
http://downloads.openmicroscopy.org/bio-formats/5.5.0/
