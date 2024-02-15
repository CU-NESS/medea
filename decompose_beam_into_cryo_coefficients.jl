using Healpix
using CryoFaBs
using HDF5

nside = 32
prenormalize=true
save_basis_in_hdf5=true
medea_env_var = get(ENV, "MEDEA", "default_value_if_not_set")
horizon_name_str = "flat_horizon"
beam_name = "horizontal_dipole_PEC"
coeff_save_filepath = joinpath(medea_env_var, "input", "cryo_coeff_"*horizon_name_str*beam_name*".hdf5")
healpy_beam_maps_filepath = joinpath(medea_env_var, "input", beam_name*"_beam_maps.hdf5")
basis_filepath_cfb = joinpath(medea_env_var, "input","cryo_basis_"*horizon_name_str*".cfb")
basis_filepath_hdf5 = joinpath(medea_env_var,"input","cryo_basis_"*horizon_name_str*".hdf5")

horizon_hdf5_key = "/"*horizon_name_str*"_healpy_map_beam_frame_nside_"*string(nside)

frequencies = collect(50:100)
frequency_indices = collect(1:50)
hyper_parameter_array = collect(1.0:0.025:3.0)

horizon_mask = h5read(joinpath(medea_env_var,"input","horizon_files.hdf5"),horizon_hdf5_key)
float_mask = convert(AbstractVector{Float64},horizon_mask)
horizon = HealpixMap{Float64, Healpix.RingOrder}(float_mask)
horizon = Healpix.udgrade(horizon,nside)

nonempty_indices = findall(>(0), horizon)

horizon[nonempty_indices] .= 1.

if isfile(basis_filepath_cfb)
	cfb = AngularCryoFaB(basis_filepath_cfb)
else
	cfb = AngularCryoFaB(horizon)
	write(basis_filepath_cfb, cfb)
end

if save_basis_in_hdf5
	if isfile(basis_filepath_hdf5)
		cfb_file = h5read(basis_filepath_hdf5, "/Basis")
		cfb_beam_to_kl = cfb_file["Transform_to_kl"]
		cfb_kl_to_beam = cfb_file["Transform_to_map"]
	else
		cfb_filename_hdf5 = h5open(basis_filepath_hdf5, "cw")
		cfb_group = create_group(cfb_filename_hdf5, "Basis")
		write(cfb_group, "Transform_to_kl", cfb.TransInv)
		write(cfb_group, "Transform_to_map", cfb.Trans)
		write(cfb_group, "nonmasked_indices", nonempty_indices)
	end
end

for hyper_par in hyper_parameter_array
	println("Decomposing hyper-parameter "*string(hyper_par))
	coeff_file =
		h5open(coeff_save_filepath,"cw")
				
	parameter_group = create_group(coeff_file, string(hyper_par))
	coeff_group = create_group(parameter_group, "coefficients")
	
	beam_maps = h5read(healpy_beam_maps_filepath, string(hyper_par))
	beams = beam_maps["beam_maps"]

	if prenormalize
		(beams = [HealpixMap{Float64, Healpix.RingOrder}(beams[:,ibeam]/(sum(beams[:,ibeam]))) 
			for ibeam in frequency_indices])
	else
		beams = [HealpixMap{Float64, Healpix.RingOrder}(beams[:,ibeam]) for ibeam in frequency_indices]
	end

	for (ibeam, beam) in enumerate(beams)
		beam_map = beam[horizon.>0]
		dlnm = cfb_beam_to_kl * beam_map
		write(coeff_group, "Freq_"*string(frequencies[ibeam]), dlnm)
	end

end
	
	
	
	
	
	
	
	
	
	
