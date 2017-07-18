#! /usr/bin/env python
#Description: Script to calculate temperature and emissivity from Owl data.

"""
Function owl_temp_emissivity
Calculates the temperature of each pixel of a thermal hyperspectral image based on the
measured spectra. Also calculates the error in the temperature as well as the emissivity
spectra. The method assumes that the emissivity approaches 1 somewhere in the spectrum,
which is not true for metals, so the temperature will be underestimated. A corresponding
higher error will be given for these low emissivity materials. This script is a research
tool and as such no guarentee is given for the accuracy of the output data.

Authors: Laura Harris (NERC-ARF-DAN), Hannu Holma (Specim)

Change history
30/06/2017: Created

Known issues: Assumes that the emissivity approaches 1 somewhere in the spectrum. Not so
              true for metals, so will under estimate temperature. This will be indicated
              in the error file.

"""
###########################################################
# This file has been created by the NERC-ARF Data Analysis Node and
# is licensed under the GPL v3 Licence. A copy of this
# licence is available to download with this file.
###########################################################

from __future__ import print_function
import argparse
import numpy
from arsf_envi_reader import envi_header
from arsf_envi_reader import numpy_bin_reader
import sys
import os

# stop print out of division warnings
numpy.seterr(invalid="ignore", divide="ignore")

######################
def planck_T(wavelengths, radiance):
    """
    Function planck_T
    Calculates the temperature given the radiance using Planck's Law.

    Arguments
    wavelengths: numpy array of wavelengths, with units of nm
    radiance: measured radiance in units of W/(m**2 sr um)

    Returns the temperature for each band in degrees C
    """

    c = 299792458 # speed of light
    h = 6.626e-34 # Planck's constant
    k = 1.38e-23 # Boltzmann constant

    wavelengths = wavelengths / 1e+6 # convert to m
    radiance = radiance * 1e+3 # convert to W/(m**2 sr m)

    T = (h*c /(wavelengths[:, None] * k)) / (numpy.log 
                     (1/ (radiance * wavelengths[:, None]**5 /(2*h*c**2)) + 1))

    TC = T - 273.15 # temp in degrees C

    return TC

######################
def calc_em(wavelengths, radiance, temp):
    """
    Function calc_em
    Calculates the emissivity with the estimated temperature.

    Arguments
    wavelengths: numpy array of wavelengths, with units of nm
    radiance: measured radiance in units of W/(m**2 sr um)
    temp: temperature in degrees C

    Returns the spectral emissivity
    """

    c = 3.0e+8 # speed of light
    h = 6.626e-34 # Planck's constant
    k = 1.38e-23 # Boltzmann constant

    wavelengths = wavelengths / 1e+6 # convert to m
    T = temp + 273.15 # convert to K

    # Predicted radiance for emissivity = 1 (black body)
    L = ((2*h*c**2/wavelengths[:, None]**5) * 
                          1/(numpy.exp(h*c/(wavelengths[:, None]*k*T))-1)*1e-3)

    E = numpy.divide(radiance, L)

    return E

######################

def _run():
    parser = argparse.ArgumentParser(description="Calculates temperature and "
                                     "emissivity from owl data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--file", help="Full path to processed file", required = True)
    parser.add_argument("-o", "--outdir", help="Full path to output directory", required = True)
    parser.add_argument("-p", "--percentile", help="The percentile of the "
                        "spectra used to determine the temperature. Should be set "
                        "to ignore spikes, but not too low that good data are ignored. "
                        "Ideally spikes should be masked and percentile set high.",
                        default = 99)
    parser.add_argument("-s", "--start_band", help="Start band to use. Bands at "
                        "the edges of the detector are more noisy.", default = 5,
                        type = int)
    parser.add_argument("-e", "--end_band", help="End band to use. Bands at "
                        "the edges of the detector are more noisy.", default = 95,
                        type = int)

    args = parser.parse_args()

    # get filename without path and extension
    infile = os.path.splitext(os.path.split(args.file)[1])[0]

    # set output files
    temp_outfile = os.path.join(args.outdir, infile + "_temperature.bil")
    em_outfile = os.path.join(args.outdir, infile + "_emissivity.bil")
    error_outfile = os.path.join(args.outdir, infile + "_temp_error.bil")

    temp_hdr = temp_outfile[:-4] + ".hdr"
    em_hdr = em_outfile[:-4] + ".hdr"
    error_hdr = error_outfile[:-4] + ".hdr"

    # Read in wavelengths
    hdr = envi_header.find_hdr_file(args.file)
    hdrData = envi_header.read_hdr_file(hdr)
    wavelengths = numpy.asarray(hdrData["wavelength"].split(","), dtype = float)/1000

    # check bands are correct
    if args.start_band < 0 or args.start_band > len(wavelengths):
        print("Start band specified is out of range ({} bands)".format(hdrData["bands"]))
        sys.exit()
    if args.end_band < 0 or args.end_band > len(wavelengths):
        print("End band specified is out of range ({} bands)".format(hdrData["bands"]))
        sys.exit()

    # Open bil file reader
    in_file = numpy_bin_reader.BilReader(args.file)
    # open output files
    temp_out = open(temp_outfile,"wb")
    em_out = open(em_outfile,"wb")
    error_out = open(error_outfile,"wb")

    # Read in data
    for line, dataline in enumerate(in_file):
        sys.stdout.write("\r working on line {}".format(line))

        # calculate temperature by line
        spectral_T = planck_T(wavelengths, dataline)

        # calculate a single temperture from specified bands, ignoring noise
        line_temp = numpy.percentile(spectral_T[args.start_band:args.end_band],
                                  args.percentile, axis = 0).astype(numpy.float32)

        # estimate the error in temperature
        line_temp_error = 2.5*numpy.std(spectral_T[args.start_band:args.end_band]
                                     , axis = 0).astype(numpy.float32)
        # if error is more than 3 degrees, black and grey bodies are bad approximations
        # and the errors are also underestimates

        # calculate emissivity
        emissivity = calc_em(wavelengths, dataline, line_temp)

        # convert NaN to 0s (may be in mapped files)
        emissivity = numpy.nan_to_num(emissivity).astype(numpy.float32)

        # write out temp and emissivity bil files
        line_temp.tofile(temp_out)
        emissivity.tofile(em_out)
        line_temp_error.tofile(error_out)

    sys.stdout.write("\r all lines complete                   \n")

    #close input file
    dataline = None
    #close output file
    temp_out.close()
    em_out.close()
    error_out.close()

    # write headers
    hdr_dict = {"samples": int(hdrData["samples"]), "lines": int(hdrData["lines"]),
                     "bands": 1, "data type": 4}
    hdr_dict["_comments"] = "file written by owl_temp_emissivity.py"

    # if mapped file, map info might be useful to keep
    if "map info" in hdrData:
        hdr_dict["map info"] = hdrData["map info"]

    envi_header.write_envi_header(temp_hdr, hdr_dict)

    envi_header.write_envi_header(error_hdr, hdr_dict)

    hdr_dict["bands"] = int(hdrData["bands"])
    hdr_dict["interleave"] = "bil"
    hdr_dict["wavelength"] = hdrData["wavelength"]
    hdr_dict["band names"] = hdrData["wavelength"]

    envi_header.write_envi_header(em_hdr, hdr_dict)

if __name__ == "__main__":
   _run()

