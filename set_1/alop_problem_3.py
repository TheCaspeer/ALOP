from astropy . table import Table
import sys
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.constants as const



# PROBLEM 3

# reading data
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
spectra_list = Table.read(f"{folder}spectra_list.csv", format="ascii.csv") # row0=filename, row1=star_type
spectra_folder = f"{folder}spectra/" # files in spectra folder have row0=wavelength [Å], row1= flux density [erg/s/cm^2/Å]

# λ_0 ,  ∆λ
bands = {
    "U": (3659, 660),
    "B": (4582, 940),
    "V": (5448,880)
}
def get_wavelength_range(band):
    lambda_0 = bands[band][0]
    width = bands[band][1]

    return (lambda_0 - width/2, lambda_0 + width/2)

# reading vega data
vega_info = Table.read(f"{folder}vega_fnu.csv", format="ascii.csv") # row0=wavelength [nm], row1= flux density [erg cm-2 s-1 Hz-1]
vega_info['wavelength'] = vega_info['wavelength'] * 10 # convert wavelength from nm to Å 
vega_info['flux'] = vega_info['flux'] *  (const.c * 1e10) / (vega_info['wavelength']**2) # convert flux density from erg cm-2 s-1 Hz-1 to erg cm-2 s-1 Å-1


u_band_range = get_wavelength_range("U")
b_band_range = get_wavelength_range("B")
v_band_range = get_wavelength_range("V")

u_mask = (vega_info['wavelength'] >= u_band_range[0]) & (vega_info['wavelength'] <= u_band_range[1])
vega_u_band = vega_info[u_mask]
b_mask = (vega_info['wavelength'] >= b_band_range[0]) & (vega_info['wavelength'] <= b_band_range[1])
vega_b_band = vega_info[b_mask]
v_mask = (vega_info['wavelength'] >= v_band_range[0]) & (vega_info['wavelength'] <= v_band_range[1])
vega_v_band = vega_info[v_mask]


class spectra:
    def __init__(self, filename, star_type):
        self.filename = filename
        self.star_type = star_type
        self.data = Table.read(f"{spectra_folder}{filename}", format="ascii.csv")
        u_band = self.get_band_flux("U")
        b_band = self.get_band_flux("B")
        v_band = self.get_band_flux("V")

    def get_band_flux(self, band):
        range = get_wavelength_range(band)
        mask  = (self.data['wavelength'] >= range[0]) & (self.data['wavelength'] <= range[1])

        return  self.data[mask]

spectras_computed = []
for csv in spectra_list:
    s = spectra(csv[0], csv[1])
    spectras_computed.append(s)
