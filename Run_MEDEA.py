"""
File: medea/Run_MEDEA.py
Author: Joshua J. Hibbard
Date: January 2024

Description: Example script opening Cryo-coefficients and the Cryo-basis
			 for a particular experiment and giving them
			 to MEDEA for an analytical horizontal dipole suspended above
			 a PEC. In this case we use the flat horizon Cryo-basis,
			 and the corresponding Cryo-coefficients for this basis
			 and the analytical dipole beam.
"""
import numpy as np
from medea import BeamEmulator
import h5py
import os

input_file_path = os.getenv('MEDEA')
cryo_coefficients_path = input_file_path + \
	'/input/cryo_coeff_flat_horizon_horizontal_dipole_PEC.hdf5'
cryo_basis_path = input_file_path + \
	'/input/cryo_basis_flat_horizon_nside32.hdf5'
	
with h5py.File(cryo_basis_path, 'r') as cfile:
	cryo_basis = cfile['Basis']['Transform_to_map'][()]
	unmasked_indices = cfile['Basis']['nonmasked_indices'][()]
	
frequencies = np.linspace(50,99,50)
hyper_parameter_array = np.round(np.linspace(1,3,21),3) #in meters between 1 and 3
nside = 32

cryo_coefficients = []
with h5py.File(cryo_coefficients_path, 'r') as cofile:
	for par in hyper_parameter_array:
		par_by_frequency = []
		for frequency in frequencies:
			par_by_frequency.append(cofile[str(par)]['coefficients']\
				['Freq_'+str(int(frequency))][()])
		cryo_coefficients.append(par_by_frequency)
			
cryo_coefficients = np.array(cryo_coefficients)

print(cryo_coefficients.shape)

medea_emulator = BeamEmulator(frequencies, nside, cryo_basis, \
	cryo_coefficients, hyper_parameter_array,\
	unmasked_indices, interpolation_order=3, use_gpr=False)
	
print(medea_emulator([1.234]).shape)

try:
	import healpy as hp
	import matplotlib.pyplot as plot
	hp.orthview(medea_emulator([1.234])[0,:])
	plot.show()
except:
	raise ImportError('Install healpy to view interpolated beam maps!')
