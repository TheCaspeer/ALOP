from unittest import case

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

# Settings
interpolation_step = 5

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
        wavelenghts = np.arange(get_wavelength_range(band)[0], get_wavelength_range(band)[1], interpolation_step)
        interpolated_self = interpolate_flux(self.data, wavelenghts)
        vega_int_band = interpolated_vega[band]

        # Assuming S_x(λ) = 1 for all relevant wavelengths in bound
        numerator = np.trapezoid(interpolated_self*wavelenghts, wavelenghts)
        denominator = np.trapezoid(wavelenghts, wavelenghts)

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
    # print(f"The C constant for {C} band is: {C_x[C]:.4f}")


# Determining magnitudes for each star

# Interpolating arrays to have proper values
def interpolate_flux(spectra, band_range):
    return np.interp(band_range, spectra['wavelength'], spectra['flux'])

interpolated_vega = {
    "U": interpolate_flux(vega_info, np.arange(u_band_range[0], u_band_range[1], interpolation_step)),
    "B": interpolate_flux(vega_info, np.arange(b_band_range[0], b_band_range[1], interpolation_step)),
    "V": interpolate_flux(vega_info, np.arange(v_band_range[0], v_band_range[1], interpolation_step))
}

# Computing magnitudes for each star
for spectrum in spectras_computed:
    spectrum.m_u = spectrum.compute_magnitude("U", C_x["U"])
    spectrum.m_b = spectrum.compute_magnitude("B", C_x["B"])
    spectrum.m_v = spectrum.compute_magnitude("V", C_x["V"])

    spectrum.m_ub = spectrum.m_u - spectrum.m_b
    spectrum.m_bv = spectrum.m_b - spectrum.m_v


def planck(wavelength,T):
    # Wavelength assumed to be in Ångström, T in Kelvin

    wavelength_m = wavelength * 1e-10 # convert to meters
    exponent = (const.h * const.c) / (wavelength_m * const.k * T)
    B = 2 * const.h * const.c**2 / (wavelength_m**5) * 1/(np.expm1(exponent)) 

    # convert B to erg/s/cm^2/Å
    # 1 W = 1 J/s = 10^7 erg/s, 1 m^2 = 10^4 cm², 1 m = 10^10 Å -> 10^-7 erg/s/cm²/Å
    B_converted = B * 10**-7 

    return B_converted

class  blackbody:
    def __init__(self, T):
        self.T = T
        self.fluxs = {}

    def compute_magnitude(self,band,C_x):
            wavelenghts = np.arange(get_wavelength_range(band)[0], get_wavelength_range(band)[1], interpolation_step)
            vega_int_band = interpolated_vega[band]
            flux = planck(wavelenghts, self.T)
            self.fluxs[band] = flux

            # Assuming S_x(λ) = 1 for all relevant wavelengths in bound
            numerator = np.trapezoid(flux*wavelenghts, wavelenghts)
            denominator = np.trapezoid(wavelenghts, wavelenghts)

            return -2.5 * np.log10(numerator / denominator) + C_x

blackbodies = []
for T in np.logspace(np.log10(2000),np.log10(30000),num=500,base=10):
    bb = blackbody(T)
    blackbodies.append(bb)

for body in blackbodies:
    body.m_u = body.compute_magnitude("U", C_x["U"])
    body.m_b = body.compute_magnitude("B", C_x["B"])
    body.m_v = body.compute_magnitude("V", C_x["V"])

    body.m_ub = body.m_u - body.m_b
    body.m_bv = body.m_b - body.m_v
    
spectral_colors = {
    "O": "aqua",
    "B": "blue",
    "A": "skyblue",
    "F": "snow",
    "G": "yellow",
    "K": "orange",
    "M": "red"
} # not fully as required by assignment,  TODO: FIX


fig, ax = plt.subplots(figsize=(10,6))
for s in spectras_computed:
    if s.star_type[0:2] == "A0":
        ax.scatter(s.m_bv, s.m_ub, label=f"{s.star_type} (A0)", marker="s", color="pink", s=6)
        print(f"Found A0 star: {s.filename} with B-V color: {s.m_bv:.3f} and U-B color: {s.m_ub:.3f}")
    else:
        ax.scatter(s.m_bv, s.m_ub, label=f"{s.star_type}", marker="*", color=spectral_colors[s.star_type[0]], s=6)
# ax.scatter([s.m_bv for s in spectras_computed], [s.m_ub for s in spectras_computed], label="Stars",marker="*", color=[spectral_colors[s.star_type[0]] for s in spectras_computed], s=6)
ax.plot([bb.m_bv for bb in blackbodies], [bb.m_ub for bb in blackbodies], label="Blackbodies", color="black", linewidth=1)
ax.set_xlabel("B-V")
ax.set_ylabel("U-B")
ax.invert_yaxis()


legend_elements = [
    plt.Line2D([0], [0], marker='*', color='w', label=s_type, markerfacecolor=color, markersize=10) 
    for s_type, color in spectral_colors.items()
] + [
    plt.Line2D([0], [0], marker='o', color='w', label='Blackbodies', markerfacecolor='black', markersize=5)
]

# TODO: maybe implement for blackbodies: cmap = plt.colormaps["coolwarm"], https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html

fig.legend(handles=legend_elements, loc='upper right', title="Spectral Types")
ax.set_facecolor('lavenderblush')


# Exercise 4

def find_blackbody_temp_bv(spectrum):
    # Find the blackbody temperature that corresponds to a given B-V color index
    closest_bb = min(blackbodies, key=lambda bb: abs(bb.m_bv - spectrum.m_bv)) # computationally expensive, but easy :)
    spectrum.closest_bb = closest_bb
    


specific_types = ["B0", "A0", "F0", "G0"]
specific_stars = []

for type in specific_types: 
    for spectrum in spectras_computed:
        if spectrum.star_type[0:2] == type:
            specific_stars.append(spectrum)
            find_blackbody_temp_bv(spectrum)

            wavelenghts = np.arange(
                3000, 9000, interpolation_step
            )  
            interpolated_self = interpolate_flux(spectrum.data, wavelenghts) 
            interpolated_blackbody = planck(wavelenghts, spectrum.closest_bb.T)

            reference_wavelength = bands["V"][0] # reference wavelength for normalization, using V band central wavelength
            closest_index = np.argmin(np.abs(wavelenghts - reference_wavelength))
            normalization_factor = interpolated_self[closest_index] / interpolated_blackbody[closest_index]
            normalized_blackbody = interpolated_blackbody * normalization_factor
            
            spectrum.normalized_flux = interpolated_self / normalized_blackbody
            spectrum.normalized_wavelengths = wavelenghts
            break # only want one star of each type, so break after finding the first match



fig2, ax2 = plt.subplots(figsize=(10,6))

for star in specific_stars:
    ax2.plot(star.normalized_wavelengths, star.normalized_flux, label=f"{star.star_type} (normalized)", color=spectral_colors[star.star_type[0]], alpha=0.5)

ax2.set_xlabel("Wavelength [Å]")
ax2.set_ylabel("Normalized Flux [F_λ / F_λ,bb]")

plt.show() 








    