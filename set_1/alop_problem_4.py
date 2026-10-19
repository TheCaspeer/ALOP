import numpy as np
import alopphot as alop
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import os
import sys
import alop_problem_3 as alop3

# Settings
gain = 5 # [electrons/ADU]
radius = 5 # pixels
k = np.pi * radius**2 # Aantal pixels in de aperture
RON = 250 # Read-out noise [electrons]
fwhm = 5.0 # Verwachte sterbreedte (in pixels)
N = 5 

class filter:
    def __init__(self,band,exp_time,mag_0point):
        self.band = band
        self.exp_time = exp_time
        self.mag_0point = mag_0point

    def compute_snr(self):
        self.snr = gain * self.fluxs / np.sqrt(gain * self.fluxs + k * (gain * self.sky + RON**2))

    def sigma_magnitude(self):
        self.sigma_m = 2.5 / (np.log(10) * self.snr)



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

B_band.sky_median = np.median(B_band.data)
B_band.sky_std = np.std(B_band.data)
V_band.sky_median = np.median(V_band.data)
V_band.sky_std = np.std(V_band.data)

B_band.threshold = B_band.sky_median + N * B_band.sky_std
V_band.threshold = V_band.sky_median + N * V_band.sky_std

B_band.starpositions = alop.findstars(B_band.data, B_band.threshold, fwhm)
V_band.starpositions = alop.findstars(V_band.data, V_band.threshold, fwhm)

print(f"Aantal gedetecteerde sterren in B-band: {len(B_band.starpositions)}")
print(f"Aantal gedetecteerde sterren in V-band: {len(V_band.starpositions)}")

B_band.fluxs, B_band.sky = alop.aperphot(B_band.data, B_band.starpositions['x'], B_band.starpositions['y'], apradius=radius, skyradii=[10, 15])
V_band.fluxs, V_band.sky = alop.aperphot(V_band.data, V_band.starpositions['x'], V_band.starpositions['y'], apradius=radius, skyradii=[10, 15])

B_band.rates = B_band.fluxs / B_band.exp_time
V_band.rates = V_band.fluxs / V_band.exp_time
B_band.magnitudes = B_band.mag_0point - 2.5 * np.log10(B_band.rates)
V_band.magnitudes = V_band.mag_0point - 2.5 * np.log10(V_band.rates)


BV_colors = B_band.magnitudes - V_band.magnitudes

print("--- Overzicht van alle sterren ---")
for i in range(len(B_band.fluxs)):
    print(f"Ster {i+1}: Flux_B = {B_band.fluxs[i]:.1f}, Magnitude_B = {B_band.magnitudes[i]:.2f}, B-V = {BV_colors[i]:.3f}")

# Opdracht 6

B_band.compute_snr()
V_band.compute_snr()

B_band.sigma_magnitude()
V_band.sigma_magnitude()

BV_sigma = np.sqrt(B_band.sigma_m**2 + V_band.sigma_m**2)

# Opdracht 7 
for star in BV_colors:
    closest_star = min(alop3.spectras_computed, key=lambda s: abs(s.m_bv - star))
    print(f"Closest star: {closest_star.filename} (Type: {closest_star.star_type}) with B-V color: {closest_star.m_bv:.3f}")

# Opdracht 8 
# Star 3 is clasified as a B9V, but probably an A0 star. Given it has a B-V of -0.001, thusly we use star 3 instead of star 5 (an A0V star)
star_index = 2 #(3-1)

'''
Vega has a B-vand and V-band magnitude of 0. 
Given both are A0V stars, they should have the same intensity, and only differ by distance. Thusly we use the distance modulus:
d_star = d_vega * 10**(m_star/5)
'''
d_vega = 7.75

d_B = d_vega * 10**(B_band.magnitudes[star_index]/5 ) # pc
d_V = d_vega * 10**(V_band.magnitudes[star_index]/5 ) # pc
err_d_B = d_B * np.log(10) * B_band.sigma_m[star_index] / 5
err_d_V = d_V * np.log(10) * V_band.sigma_m[star_index] / 5

d_mean = (d_B + d_V) / 2
d_err_mean = np.sqrt(err_d_B**2 + err_d_V**2) / 2

print(f"--- Distance estimation for star {star_index+1} (Type: {alop3.spectras_computed[star_index].star_type}) ---")
print(f"Distance (B-band): {d_B:.2f} ± {err_d_B:.2f} pc")
print(f"Distance (V-band): {d_V:.2f} ± {err_d_V:.2f} pc")
print(f"Mean distance: {d_mean:.2f} ± {d_err_mean:.2f} pc")