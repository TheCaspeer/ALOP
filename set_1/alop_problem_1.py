from astropy . table import Table
import sys
import matplotlib.pyplot as plt
import numpy as np
import os




# PROBLEM 1

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


# selecting stars to plot
number_of_stars = 3
selected_stars = np.random.randint(0, len(spectra_list), size=number_of_stars)

# reading star data from files
stars = {
    "name": [],
    "stellar_type": [],
    "file": [],
    "data": []
}
for i in selected_stars:
    star_name = spectra_list[i][0][:-4] # removing the .csv extension
    star_file = spectra_list[i][0]
    star_type = spectra_list[i][1]
    star_data = Table.read(f"{spectra_folder}{star_file}", format="ascii.csv")

    stars['name'].append(star_name)
    stars['stellar_type'].append(star_type)
    stars['file'].append(star_file)
    stars['data'].append(star_data) # star_data has row0=wavelength [Å], row1= flux density [erg/s/cm^2/Å]
    print(f"Star {star_name} of type {star_type} read from file {star_file}")


# plotting the spectra of the selected stars
fig, ax = plt.subplots(figsize=(13, 9))
spectra_range = [3000,9000] # wavelength range [Å]
for i in range(number_of_stars):
    data = stars['data'][i]
    wavelength = data['wavelength']
    flux_density = data['flux']
    mask = (wavelength >= spectra_range[0]) & (wavelength <= spectra_range[1]) # masking the wavelength range
    wavelength = wavelength[mask]
    flux_density = flux_density[mask]

    ax.plot(wavelength, flux_density, label=f"{stars['name'][i]} [{stars['stellar_type'][i]}] ", linewidth=1, zorder=-1, ls='-', alpha=0.7)
    ax.scatter(wavelength, flux_density, s=1, zorder=0)

# figure configuration
ax.legend()
ax.set_title(f"Spectra of {number_of_stars} random stars")
ax.set_xlabel("Wavelength [Å]")
ax.set_ylabel("Flux density [erg/s/cm^2/Å]")

ax.set_xlim(spectra_range)
ax.set_yscale('log')
ax.grid(True, zorder=10)
fig.tight_layout()


# saving & showing plt
plt.savefig(f"{folder}random_stars_spectra.svg")
plt.show()









