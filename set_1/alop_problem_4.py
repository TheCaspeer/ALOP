import numpy as np
import alopphot as alop
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import os
import sys

class filter:
    def __init__(self,band,exp_time,mag_0point):
        self.band = band
        self.exp_time = exp_time
        self.mag_0point = mag_0point


B_band = filter("B", 30, 16.43)
V_band = filter('V', 10, 17.40)

''' 
DATA ASSUMES FOLLOWING FILE STRUCTURE
CWD (ALOP)/
set_1/
    Data/
        spectra_list.csv
        spectra/
         <various .csv files with spectra data>
    alop_set_1.py
'''
cwd = os.getcwd()
folder = cwd + "/set_1/Data/"
data_B = fits.getdata(f"{folder}exposure_B_alop_2026_2027.fits")
data_V = fits.getdata(f"{folder}exposure_V_alop_2026_2027.fits")

B_band.data = data_B
V_band.data = data_V

fwhm = 5.0 # Verwachte sterbreedte (in pixels)
B_band.sky_median = np.median(B_band.data)
B_band.sky_std = np.std(B_band.data)
V_band.sky_median = np.median(V_band.data)
V_band.sky_std = np.std(V_band.data)

N = 5 
B_band.threshold = B_band.sky_median + N * B_band.sky_std
V_band.threshold = V_band.sky_median + N * V_band.sky_std

B_band.starpositions = alop.findstars(B_band.data, B_band.threshold, fwhm)
V_band.starpositions = alop.findstars(V_band.data, V_band.threshold, fwhm)

radius = 5 # pixels
B_band.flux, B_band.sky = alop.aperphot(B_band.data, B_band.starpositions['x'], B_band.starpositions['y'], apradius=radius, skyradii=[10, 15])
V_band.flux, V_band.sky = alop.aperphot(V_band.data, V_band.starpositions['x'], V_band.starpositions['y'], apradius=radius, skyradii=[10, 15])

B_band.rate = B_band.flux / B_band.exp_time
V_band.rate = V_band.flux / V_band.exp_time
B_band.magnitude = B_band.mag_0point - 2.5 * np.log10(B_band.rate)
V_band.magnitude = V_band.mag_0point - 2.5 * np.log10(V_band.rate)


BV_color = B_band.magnitude - V_band.magnitude

# Opdracht 6

gain = 5 # [electrons/ADU]