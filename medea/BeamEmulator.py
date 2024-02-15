"""
File: medea/BeamEmulator.py
Author: Joshua J. Hibbard
Date: January 2024

Description: Object which takes an input set of beams (analytical
			 or numerical) defined at certain beam hyper-parameters,
			 a basis, and Cryo-coefficients
			 for the corresponding input set of beams,
			 and produces a model (MEDEA) capable of quickly emulating
			 beam topological patterns at arbitrary hyper-parameters
			 between values within the input set. 
"""
import numpy as np
import h5py
from scipy.interpolate import make_interp_spline, interp1d
import healpy as hp
try:
	from sklearn.gaussian_process import GaussianProcessRegressor
	from sklearn.gaussian_process.kernels import RBF, Matern, RationalQuadratic
	from sklearn.multioutput import MultiOutputRegressor
	from sklearn.multioutput import RegressorChain
except:
	print('Install sklearn in order to use GPR for interpolation.')

class BeamEmulator(object):
	def __init__(self, frequencies, nside, cryobasis, \
		cryo_coefficients, hyper_parameter_array,\
		unmasked_indices, interpolation_order=3, use_gpr=False,\
		coefficient_order=None, gpr_kwargs={},\
		return_beam_above_horizon=False,unmasked_indices_in_julia_ordering=True):
		'''
		frequencies: a 1D array of the frequencies which apply to the beam to be emulated. 
					 Note: Must have corresponding Cryo-coefficients for all frequencies
					 in this array.
		nside: the healpy nside parameter of the final emulated beams.
		cryobasis: the 2D numpy array containing the basis to transform from coefficient space
			 	   to pixel space. That is, this is the basis Y_il which takes
				   cryo-coefficients k_l and produces the beam pattern above the horizon: 
				   B_i = Y_il k_l.
		cryo_coefficients: the Cryo-coefficients as a multi-dimensional array
						   of shape: [num_hyper_par_1, num_hyper_par_2, ... , num_frequencies,\
						   num_pixels]. Note that the num_pixels corresponds to the pixels
						   above the horizon i, or Y_il in the basis, and num_frequencies
						   is len(self.frequencies). The first dimensions correspond
						   to the number of hyper_parameters of a given type to interpolate
						   over. For instance, num_hyper_par_1 may correspond to 
						   soil dielectrics, while num_hyper_par_2 could be height of 
						   dipole above a ground plane.
						   NOTE: only spline interpolation is supported for multi-dimensional
						   hyper-parameter interpolation.
		hyper_parameter_array: An ND-array containing the hyper-parameters of the input beam set,
							   where N denotes the number of hyper-parameters to interpolate over.
							   Each axis corresponds to a separate set of hyper-parameters,
							   increasing monotonically. If GPR interpolation is used, 
							   hyper_parameter_array must be a 1D array.
		unmasked_indices: the unmasked pixel indices (i.e. above the horizon) in healpy index RING format.
		interpolation_order: the order of the spline to be used for interpolation,
							 with a Default of 3 (cubic order).
		use_gpr: Boolean denoting whether to use Gaussian Process Regression (GPR)
				 instead of the splines for coefficient interpolation.
		coefficient_order: the truncation order of the Cryo-coefficients to use when generating
						   the Cryobeam. Default None (all coefficients are used).
		gpr_kwargs: A dictionary of keys to give to the GPR calculator, including:
					'kernel': with possible values of 'RBF', 'Matern_2_5', 'Matern_1_5',
							  and 'RQ'.
					'normalize_y': whether to normalize the coefficients before finding
								   best-fit kernel hyper-parameters. In most cases this should
								   be True.
					'regressor': which kind of regressor to implement, either
								 'MultiOutput' which fits one regressor per coefficient per frequency, or
								 'RegressorChain' which fits one chain of regressors to all
								 coefficients at all frequencies. Note that this regressor method is
								 still an active area of research and so should be used with caution.
		return_beam_above_horizon: Boolean determining whether BeamEmulator returns only the
								   beam pixels above the horizon (True) or the entire beam for all
								   pixels in the 4pi steradian sky (False). Default: False.
		unmasked_indices_in_julia_ordering: Boolean determining whether the self.unmasked_indices
											are in Julia array ordering (i.e. starting from index 1
											in stead of 0 in python).
		'''
		
		self.frequencies = frequencies
		self.nside = nside
		self.cryobasis = cryobasis
		self.cryo_coefficients = cryo_coefficients
		self.hyper_parameter_array = hyper_parameter_array
		self.unmasked_indices = unmasked_indices
		self.interpolation_order = interpolation_order
		self.use_gpr = use_gpr
		self.coefficient_order = coefficient_order
		self.gpr_kwargs = gpr_kwargs
		self.return_beam_above_horizon = return_beam_above_horizon
		self.unmasked_indices_in_julia_ordering = unmasked_indices_in_julia_ordering
		
		if self.unmasked_indices_in_julia_ordering:
			self.unmasked_indices = self.unmasked_indices - 1
		
	@property
	def hyper_parameter_interpolater(self):
		if not hasattr(self, '_hyper_parameter_interpolater'):
			reshaped_coeff = np.reshape(self.cryo_coefficients,\
				(self.hyper_parameter_array.shape[0],-1))
			if self.use_gpr:
				print('Using Gaussian Process Regression for interpolation...')
				desired_kernel = self.gpr_kwargs['kernel']
				if desired_kernel == 'RBF':
					kernel = 1 * RBF()
				elif desired_kernel == 'Matern_2_5':
					kernel = 1.0 * Matern(length_scale=1.0, nu=2.5)
				elif desired_kernel == 'Matern_1_5':
					kernel = 1.0 * Matern(length_scale=1.0, nu=1.5)
				elif desired_kernel == 'RQ':
					kernel = 1.0 * RationalQuadratic(length_scale=1.0, alpha=1.5)
					
				estimator = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=20,\
					normalize_y=self.gpr_kwargs['normalize_y'])
				
				if self.gpr_kwargs['regressor'] == 'MultiOutput':
					regr = MultiOutputRegressor(estimator,\
						n_jobs=-1).fit(self.hyper_parameter_array.reshape(-1, 1), \
						reshaped_coeff)
						
				elif self.gpr_kwargs['regressor'] == 'RegressorChain':
					regr = RegressorChain(estimator).fit(self.hyper_parameter_array.reshape(-1, 1), \
						reshaped_coeff)

				all_estimator_hyper_parameters = []
				all_estimator_log_marginal_likelihoods = []
				
				for estimator in regr.estimators_:
					all_estimator_hyper_parameters.append(estimator.kernel_.get_params())
					all_estimator_log_marginal_likelihoods.append(estimator.log_marginal_likelihood())
					
				self.estimator_log_marginal_likelihoods = np.array(all_estimator_log_marginal_likelihoods)
				self.estimator_hyper_parameters = all_estimator_hyper_parameters
					
				self.gpr_score_of_training_set = regr.score(self.hyper_parameter_array.reshape(-1, 1), \
					reshaped_coeff)
				
				self._hyper_parameter_interpolater = regr
		
			else:
				print('Using splines for interpolation...')
				self._hyper_parameter_interpolater = \
					make_interp_spline(self.hyper_parameter_array, self.cryo_coefficients, \
					k=self.interpolation_order)
				predicted_outputs = \
					self._hyper_parameter_interpolater(self.hyper_parameter_array)
				predicted_outputs = np.reshape(predicted_outputs, (self.hyper_parameter_array.shape[0],-1))
				try:
					from sklearn.metrics import r2_score
					self.spline_score_of_training_set = r2_score(reshaped_coeff, predicted_outputs,\
						multioutput='uniform_average')
		return self._hyper_parameter_interpolater
		
	def coefficient_predicter(self, inputs):
		if self.use_gpr:
			return self._hyper_parameter_interpolater.predict(inputs)
		else:
			return self._hyper_parameter_interpolater(inputs[0])
		
	@hyper_parameter_interpolater.setter
	def hyper_parameter_interpolater(self, value):
		self._hyper_parameter_interpolater = value
		
	@property
	def estimator_hyper_parameters(self):
		if not hasattr(self, '_estimator_hyper_parameters'):
			raise AttributeError('Estimator Hyper Parameters was referenced before it was set!')
		return self._estimator_hyper_parameters
		
	@estimator_hyper_parameters.setter
	def estimator_hyper_parameters(self, value):
		self._estimator_hyper_parameters = value
		
	@property
	def estimator_log_marginal_likelihoods(self):
		if not hasattr(self, '_estimator_log_marginal_likelihoods'):
			raise AttributeError('estimator_log_marginal_likelihoods was referenced before it was set!')
		return self._estimator_log_marginal_likelihoods
		
	@estimator_log_marginal_likelihoods.setter
	def estimator_log_marginal_likelihoods(self, value):
		self._estimator_log_marginal_likelihoods = value
		
	@property
	def gpr_score_of_training_set(self):
		if not hasattr(self, '_gpr_score_of_training_set'):
			raise AttributeError('Score was referenced before it was set!')
		return self._gpr_score_of_training_set
		
	@gpr_score_of_training_set.setter
	def gpr_score_of_training_set(self, value):
		self._gpr_score_of_training_set = value
		
	@property
	def spline_score_of_training_set(self):
		if not hasattr(self, '_spline_score_of_training_set'):
			raise AttributeError('Score was referenced before it was set!')
		return self._spline_score_of_training_set
		
	@spline_score_of_training_set.setter
	def spline_score_of_training_set(self, value):
		self._spline_score_of_training_set = value
	
	@property
	def cryo_basis_kl_to_beam(self):
		if not hasattr(self, '_cryo_basis_kl_to_beam'):
			self._cryo_basis_kl_to_beam = self.cryo_basis.T
			self._cryo_basis_kl_to_beam = self._cryo_basis_kl_to_beam[:,:self.coefficient_order]
		return self._cryo_basis_kl_to_beam
		
	@cryo_basis_kl_to_beam.setter
	def cryo_basis_kl_to_beam(self, value):
		self._cryo_basis_kl_to_beam = value
		
	@cryo_basis_beam_to_kl(self):
		if not hasattr(self, '_cryo_basis_beam_to_kl'):
			no_truncation_basis = self.cryo_basis.T
			self._cryo_basis_beam_to_kl = \
				np.linalg.inv(no_truncation_basis)
			self._cryo_basis_beam_to_kl = \
				self._cryo_basis_beam_to_kl[:,:self.coefficient_order]
		return self._cryo_basis_beam_to_kl
		
	@cryo_basis_beam_to_kl.setter:
	def cryo_basis_beam_to_kl(self,value):
		self._cryo_basis_beam_to_kl = value		
	
	@property
	def hyper_parameter_array(self):
		if not hasattr(self, '_hyper_parameter_array'):
			raise AttributeError('hyper_parameter_array referenced before set!')
		return self._hyper_parameter_array
		
	@hyper_parameter_array.setter
	def hyper_parameter_array(self,value):
		self._hyper_parameter_array = value
		
	@property
	def interpolation_order(self):
		if not hasattr(self,'_interpolation_order'):
			raise AttributeError('interpolation_order referenced before set!')
		return self._interpolation_order
	
	@interpolation_order.setter
	def interpolation_order(self,value):
		self._interpolation_order = value

	def __call__(self, parameter):
		if self.use_gpr:
			interp_coeff = \
				self.hyper_parameter_interpolater.predict(np.array(parameter[0]).reshape(1,-1))
			interp_coeff = np.reshape(interp_coeff, (len(self.frequencies), -1))
		else:
			interp_coeff = \
				self.hyper_parameter_interpolater(parameter[0])
		if self.return_beam_above_horizon:
			funks = np.matmul(self.cryo_basis_kl_to_beam,\
				interp_coeff.T)
			return np.array(funks)
		else:
			cryo_beam_list = []
			for ifreq, freq in enumerate(self.frequencies):
				beam_map = np.zeros(hp.nside2npix(self.nside))
				funks = np.matmul(self.cryo_basis_kl_to_beam,\
					interp_coeff[ifreq,:])
				beam_map[self.unmasked_indices] = funks
				cryo_beam_list.append(beam_map)
			return np.array(cryo_beam_list)	
		
