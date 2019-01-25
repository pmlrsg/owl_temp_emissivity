# Owl Temperature and Emmisivity Script #

Calculates the temperature of each pixel of a thermal hyperspectral image based on the
measured spectra. Also calculates the error in the temperature as well as the emissivity
spectra. The method assumes that the emissivity approaches 1 somewhere in the spectrum,
which is not true for metals, so the temperature will be underestimated. A corresponding
higher error will be given for these low emissivity materials. This script is a research
tool and as such no guarentee is given for the accuracy of the output data.

The simple method is based on inverting Planck's Law for black body radiation to find the
temperature given the radiance. The theory is well covered in physics textbooks and many
articles may be found online for further reading. The error used is 2.5 times the
standard deviation of the different temperatures calculated from the different bands.
Materials that are closer to black bodies will have a more accurate temperature and there
will be less range in the temperatures calculated from the radiance, whereas materials
that do not emit as a black body and thus have the most error have a larger range of
temperatures calculated from the different bands. Thus this gives a realistic indication
of the temperature error.

Authors: Laura Harris (NERC-ARF-DAN), Hannu Holma (Specim)

## Installation ##

The script was developed under Linux but should work under Windows and macOS.

### Python and NumPy ##

The script is written in Python and uses [numpy](http://www.numpy.org/).

To install these under Windows / macOS it is recommended to use the following steps:

1. Download minconda from http://conda.pydata.org/miniconda.html#miniconda and follow the instructions to install.

2. Open a 'Command Prompt' or Terminal window and install numpy by typing:
```
conda install numpy
```

For Linux numpy can be installed via the package manager.

### NERC-ARF Tools ###

To read data in ENVI format the `arsf_envi_reader` is used, which is part of the NERC-ARF tools repository (https://github.com/pmlrsg/arsf_tools).
Download and unzip the tools and navigate to this location in a command prompt Window (e.g., `cd Downloads`) then run the following:

```
python setup.py install
```

## Usage ##

```
python owl_temp_emissivity.py -f o249011b.bil -o owl_outputs/
```

## Sample Data ##

NERC-ARF processed Owl data can be downloaded from [CEDA](http://www.ceda.ac.uk/):

* [EUFAR15/48 - Hidhaz project, Iceland](http://data.ceda.ac.uk/neodc/arsf/2015/EUFAR15_48/EUFAR15_48-2015_249_Hidhaz/hyperspectral-owl/)
