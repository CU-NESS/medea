using Healpix
using CryoFaBs
using HDF5

nside = 32
prenormalize=true
save_basis_in_hdf5=true
coeff_save_filepath = "cryo_coeff_test_flat_horizon"
healpy_beam_maps_filepath = "test_beam_maps"
basis_filepath = "cryo_basis_test_flat_horizon"
basis_filepath_hdf5 = "cryo_basis_test_flat_horizon.hdf5"

if prenormalize
	prenorm_str = "_prenormalized"
else
	prenorm_str = ""
end

frequencies = collect(50:100)
frequency_indices = collect(1:50)
hyper_parameter_array = collect(1.0:0.05:3.0)

npix=nside2npix(nside)

horizon_mask = h5read("/medea/input/flat_horizon_nside32.hdf5","/map")
float_mask = convert(AbstractVector{Float64},horizon_mask)
horizon = HealpixMap{Float64, Healpix.RingOrder}(float_mask)
horizon = Healpix.udgrade(horizon,nside)

nonempty_indices = findall(>(0), horizon)

horizon[nonempty_indices] .= 1.

if isfile(basis_filepath".cfb")
	cfb = AngularCryoFaB(cfb_filename*".cfb")
else
	cfb = AngularCryoFaB(horizon)
	write(basis_filepath*cfb_filename*".cfb", cfb)
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
		h5open(coeff_save_filepath*"_nside"*string(nside)*"_50_100MHz_"*prenorm_str,"cw")
				
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
	
	
	
	
	
	
	
	
	
	
