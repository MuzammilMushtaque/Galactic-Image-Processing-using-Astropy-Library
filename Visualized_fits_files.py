# Astrophysics Code for Galaxy Analysis

import numpy as np
import matplotlib.pyplot as plt
import time
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gzip

# Record the start time of the code execution
start = time.time()

# Function to calculate various column densities
def column_densities():
    mu = 2.0
    con_kpc = 3.0856776e+19  # Kilo parsec into meters

    # Read the FITS file containing 3D data of the galaxy
    inp_list = fits.open('input_midplane_3d.fits.gz')
    inp_list.info()

    # Extract necessary parameters from the FITS header
    sidelength = np.abs(inp_list[0].header['CRVAL1']) / con_kpc
    distance = np.abs(inp_list[0].header['CDELT1'])

    # Extract data for gas and dust densities
    cube_gas_dens = inp_list[0].data[1, :, :, :]
    cube_dust_dens = inp_list[0].data[0, :, :, :]
    cube_both_dens = cube_dust_dens

    # Calculate the total mass of dust in the galaxy
    cube_both_dens = cube_both_dens * (distance)**3 * 5.0279e-31
    dust_mass = np.sum(cube_both_dens)
    print("Dust Mass in solar masses =", dust_mass)

    # Calculate optical depth for different wavelengths
    kappa = [14673.0, 13399.9, 2567.5, 234.6, 9.8]
    legends = [r'$(\log_{\mathrm{10}}(\tau(0.09 \mu m)))$', r'$(\log_{\mathrm{10}}(\tau(0.15 \mu m)))$',
               r'$(\log_{\mathrm{10}}(\tau(0.97 \mu m)))$', r'$(\log_{\mathrm{10}}(\tau(10 \mu m)))$',
               r'$(\log_{\mathrm{10}}(\tau(102 \mu m)))$']

    for k in range(len(kappa)):
        # Calculate optical depth and plot the distributions
        cube_both_dens = cube_dust_dens * kappa[k] * distance
        dust_dens_axis0 = np.sum(cube_both_dens, axis=0)

        A = plt.subplot(3, 2, k + 1)
        A = plt.imshow(np.log10(dust_dens_axis0), extent=[-sidelength, sidelength, -sidelength, sidelength],
                       origin='lower', cmap=plt.cm.Blues)
        ax = plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = plt.colorbar(A, cax=cax)
        cb.set_label(str(legends[k]), labelpad=20, rotation=270, fontsize=6)
        ax.set_ylabel(r'Y [kpc]', fontsize=10)
        ax.set_xlabel(r'X [kpc]', fontsize=10)

    plt.savefig('optical_depth_freq.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    # Read dust density data and plot column density distributions
    cube_dust_dens = inp_list[0].data[0, :, :, :]
    cube_dust_dens = cube_dust_dens / (mu * 1.67e-27)  # 1.67e-27 = hydrogen mass in kg
    dust_dens_axis0 = np.sum(cube_dust_dens * distance, axis=0)
    dust_dens_axis1 = np.average(cube_dust_dens, axis=1)
    dust_dens_axis2 = np.average(cube_dust_dens, axis=2)
    col_den_dust = dust_dens_axis0

    A = plt.subplot(211)
    A = plt.imshow(np.log10(col_den_dust / (1e4)), extent=[-sidelength, sidelength, -sidelength, sidelength],
                   origin='lower', cmap=plt.cm.gist_heat)
    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(A, cax=cax)
    cb.set_label(r'$\log_{\mathrm{10}}(n_{\mathrm{dust}}/\mathrm{cm}^{-2})$', labelpad=20, rotation=270, fontsize=10)
    ax.set_ylabel(r'Y [kpc]', fontsize=12)
    ax.set_xlabel(r'X [kpc]', fontsize=12)
    plt.savefig('input_column_density.pdf', dpi=300, bbox_inches='tight')

    # Read gas density data and plot column density distributions
    cube_gas_dens = inp_list[0].data[1, :, :, :]
    cube_gas_dens = cube_gas_dens / (mu * 1.67e-27)
    gas_dens_axis0 = np.sum(cube_gas_dens * distance, axis=0)
    gas_dens_axis1 = np.average(cube_gas_dens, axis=1)
    gas_dens_axis2 = np.average(cube_gas_dens, axis=2)
    col_den_gas = gas_dens_axis0

    B = plt.subplot(212)
    B = plt.imshow(np.log10(col_den_gas / (1e4)), extent=[-sidelength, sidelength, -sidelength, sidelength],
                   origin='lower', cmap=plt.cm.gist_heat)
    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(B, cax=cax)
    cb.set_label(r'$\log_{\mathrm{10}}(n_{\mathrm{gas}}/\mathrm{cm}^{-2})$', labelpad=20, rotation=270, fontsize=10)
    ax.set_ylabel(r'Y [kpc]', fontsize=12)
    ax.set_xlabel(r'X [kpc]', fontsize=12)
    plt.savefig('input_column_density.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    # Calculate metal column density and plot the distribution
    col_den_metal = col_den_dust / col_den_gas
    C = plt.imshow(np.log10(col_den_metal), extent=[-sidelength, sidelength, -sidelength, sidelength],
                   origin='lower', cmap=plt.cm.gist_heat)
    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(C, cax=cax)
    cb.set_label(r'$\log_{\mathrm{10}}(n_{\mathrm{dust}}/n_{\mathrm{gas}})$', labelpad=20, rotation=270, fontsize=10)
    ax.set_ylabel(r'Y [kpc]', fontsize=12)
    ax.set_xlabel(r'X [kpc]', fontsize=12)

    plt.savefig('column_density_metallicity.pdf', dpi=300, bbox_inches='tight')
    plt.close()

# Run the column_densities function when the script is executed
if __name__ == "__main__":
    column_densities()

# Record the end time of the code execution and calculate the elapsed time
end = time.time()
print('Time process takes in minutes =', (end - start) / 60)
