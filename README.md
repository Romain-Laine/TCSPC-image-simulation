# TCSPC-image-simulation
This MATLAB code provides an easy way to simulate TCSPC image data which includes effect of the IRF, background after-pulsing and laser repetition rate. Images are saved as OME.tiff using Bioformat.

The method uses a Monte-Carlo simulation of photon arrival times based on probability density functions of emission (exponential) and excitation (Gaussian). 

The OME.tiff save was kindly provided by Ian Munro <i.munro@imperial.ac.uk> from the Photonics Group at Imperial College, London. Thanks to the OME team OME team who maintains the Bio-Formats library. Contact: Sebastien Besson <s.besson@dundee.ac.uk> from Dundee.

You need to download the Bio-formats Matlab plugin from 
http://downloads.openmicroscopy.org/bio-formats/
and set the path to it as shown in the file.
