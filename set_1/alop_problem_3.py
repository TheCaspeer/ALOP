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

vega_bands = { 
    "U": vega_u_band,
    "B": vega_b_band,
    "V": vega_v_band
}

class spectra:
    def __init__(self, filename, star_type):
        self.filename = filename
        self.star_type = star_type
        self.data = Table.read(f"{spectra_folder}{filename}", format="ascii.csv")

    def compute_magnitude(self, band, C_x):
        wavelenghts = np.arange(get_wavelength_range(band)[0], get_wavelength_range(band)[1], 100)
        interpolated_self = interpolate_flux(self.data, wavelenghts)
        vega_int_band = interpolated_vega[band]

        # Assuming S_x(λ) = 1 for all relevant wavelengths in bound
        numerator = np.trapezoid(interpolated_self/wavelenghts, wavelenghts)
        denominator = np.trapezoid(vega_int_band/wavelenghts, wavelenghts)

        return -2.5 * np.log10(numerator / denominator) + C_x

spectras_computed = []
for csv in spectra_list:
    s = spectra(csv[0], csv[1])
    spectras_computed.append(s)
 # Determinining C constant for each band
'''
Using formula:
M_x = -2.5 log ( ∫ F_λ λ  S_x(λ) dλ / ∫ F_λ^0 S_x(λ) dλ ) + C_x
Given the problem at hand we can assume that S_x(λ) = 1 for all relevant wavelenghts in bound.
For vega it is known that M_x = 0  for all bands, thusly:
C_x = 2.5 log ( ∫ F_λ(λ) λ  dλ / λ dλ )
'''

def compute_C(vega_band):
    numerator = np.trapezoid(vega_band['flux']*vega_band['wavelength'], vega_band['wavelength'])
    denominator = np.trapezoid(vega_band['wavelength'], vega_band['wavelength'])

    C = 2.5 * np.log10(numerator / denominator)
    return C

C_x = {}
for C in ["U", "B", "V"]:
    vega_band = vega_bands[C]
    C_x[C] = compute_C(vega_band)
    print(f"The C constant for {C} band is: {C_x[C]:.4f}")


# Determining magnitudes for each star

# Interpolating arrays to have proper values
def interpolate_flux(spectra, band_range):
    return np.interp(band_range, spectra['wavelength'], spectra['flux'])

interpolated_vega = {
    "U": interpolate_flux(vega_info, np.arange(u_band_range[0], u_band_range[1], 100)),
    "B": interpolate_flux(vega_info, np.arange(b_band_range[0], b_band_range[1], 100)),
    "V": interpolate_flux(vega_info, np.arange(v_band_range[0], v_band_range[1], 100))
}

# Computing magnitudes for each star
for spectrum in spectras_computed:
    spectrum.m_u = spectrum.compute_magnitude("U", C_x["U"])
    spectrum.m_b = spectrum.compute_magnitude("B", C_x["B"])
    spectrum.m_v = spectrum.compute_magnitude("V", C_x["V"])

    spectrum.m_ub = spectrum.m_u - spectrum.m_b
    spectrum.m_bv = spectrum.m_b - spectrum.m_v



