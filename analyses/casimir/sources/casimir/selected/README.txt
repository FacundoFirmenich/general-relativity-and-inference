This Readme is intended to explain the contents of the dataset accompanying the publication

"Measurement of the Casimir force between superconductors"

by

Matthijs H.J. de Jong, Evren Korkmazgil, Louise Banniard, Mika A. Sillanpää, and Laure Mercier de Lépinay.

This dataset was updated on 2026-02-17 by Matthijs H.J. de Jong.
For more information and/or questions, please don't hesitate to reach out to
matthijsdejong1994@gmail.com
laure.mercierdelepinay@aalto.fi

----

This dataset contains the following files (organized by subfolder).

in 'Additional information':
	2023_Vistec_W1P3_2cav2drum.gds		CAD file containing the design of the drums, resonator and chip. 

	Casimir_01.png - Casimir_05.png		Microscope images, taken postmortem after the measurements. 
						This design contained two drums, but the bigger drum (right) collapsed during 
						fabrication. The smaller drum (left) survived and is the one measured in the 
						experiments reported in the paper. It had a separate connection to the coupling 
						capacitor (visible in Casimir_03.png), which was destroyed by a piece of dirt 
						during fabrication. The collapse of the large drum shorted one side of the
						intended double-cavity, and the breaking of the separate coupling connection
						'restored' the resulting frequency of the cavity to the 5.6 GHz mentioned in the
						paper.
	E2-100X.png				Microscope image, taken before the measurements.
	2c2m-e3-bigdrum.tiff			AFM measurement of the larger of the two drums, which is collapsed (small drum in
						'Figure 7' folder.
	2c2m-e3-bigdrum_cleaned.png		Image of the above AFM measurement after the same data processing as in Figure 7 small drum.

	in 'Fit_of_Casimir_data_for_other_n':

	Fit_Casimir_weirdn_v3.py		Python script wherein I calculate and fit some of our data sets assuming different
						values of the Casimir force coefficient P_c and force-distance scaling n. The main
						conclusion of this script is that - yes - we can fit our data for different 
						combinations of P_c and n, which yield different values for d than used in the paper.
						However, the combinations require the existence of some alternative force which is 
						comparable in strength and nonlinearity to the Casimir force, while simultaneously 
						requiring the complete absence of the Casimir force. This was not included in the main
						paper to avoid confusion for the readers, but it was included in internal discussions
						between the first author and the last author.
	Optimized_mech_freqs_for*.txt		Text file containing a list of separations d and mechanical frequencies where the final
						frequency after taking into account the Casimir force correspond to the measured 
						mechanical frequency. 
	Plot_of_Casimir_fits_*.png		Image showing fits to (some of) the data for different values of n and P_c, and the 
						separation d that results in the fit with the correct nonlinearity.
	Plot_of_Casimir_pressure_*.png		Images showing the pressure-distance curves assumed for the various fits in this folder.

	in 'Casimir-VdW_crossover':

	Calculate_Casimir_Pressure_distance_wplasma.py	Python script used to calculate the Casimir pressure at various distances.
						These calculations are based on Bimonte, PRA 99, 052507 (2019).
	Casimir_pressure_dist_v4.h5		hdf5 file containing the pressures and distances as evaluated by the above script
						 (factor pi may be missing from pressure values, known oversight in calculation script).
	Casimir_pressure_dist_v5_closedist.h5	hdf5 file idem to the above, but with more values at closer distances to show the VdW regime.
	Casimir_vdW_crossover.py		Plot of the dispersion force calculated by the above script, in the VdW and Casimir regimes.
						Some of the assumptions on the distribution of matter as a continuum start to break down in this regime,
						so it is likely the actual curves differ somewhat (others found experimentally the regime
						goes from Casimir to VdW scaling around 10 nm for Al).
	Plot_Casimir_closedist.py		Python script to turn the above data files into the plot.

	in 'Drive_mechanism':
	
	Blue_sideband_drive_vs_our_drive.png	Plot showing the simulated mechanical amplitude that can be reached for different drive
						powers (on the red sideband or our two-tone drive, with fixed cavity power, and for the 
						blue sideband without cavity power). The main conclusion is that our drive mechanism results
						in more control over the mechanical amplitude: the mechanical amplitude is linear with
						drive power over a very large regime. Conversely, with blue-sideband driving the transition
						from low-amplitude to high-amplitude is very sharp at the onset of the self-oscillation regime.
						This figure/simulation was done in response to the referee questions regarding the feasibility
						of using a single blue-sideband tone as drive instead of our (slightly more complicated) drive setup.
	Estimate_bluesidebanddrivingamplitude.py	Python script containing the calculations of the above plot.

	in 'Patch_measurements':
	Analyze_v2.py				Python script used to analyze the AFM files in hdf5 format
	Compare_patch_statistics.py		Script where patch statistics were compared
	Determine_featuresize*.py		Python script to determine the feature size in the various patch measurements
	W1P3_E4_drums.png			Screenshot of the sample through the AFM microscope
	W1P3_E4_Measurement_area_*.png		Screenshot of the sample areas that were measured through the AFM microscope
	W1P3_E4_measarea*.h5			AFM measurement data (various measurement areas, various scan sizes)
	W1P3_E4_measarea*.png			AFM and KPFM data, analyzed and post-processed.

in 'Figure 1':
	--- I've added a lot of the data/simulation/calculation files to this folder that may be duplicated in other folders, as these files are used in the
	calibration of the amplitude and the fit procedure. For completion, all of them have been added to this folder but that has increased the complexity.
	In the details below, I've first described the files that are necessary for the parts of Figure 1 that are visible. The individual calibration and
	characterization files/scripts each get their own little subsection. ---
	Figure1.pdf				Figure 1 of the paper.
	Casimir_-20dBm_inset_v4.png		Plot of the measured data showing the Casimir nonlinearity, right panel of Figure 1c (inset)
	Casimir_-20dBm_main_v4.png		Plot of the measured data showing the Casimir nonlinearity, left panel of Figure 1c

	Plot_potentials.py			Script where I calculate the shape of the potential due to the Casimir force and plot it.
	FIt_Casimir_model_v27_fit_d.py		This is the big simulation script where loads of things come together: It is what I used to manually
						fit the distance d to our data set. Many of the other data/script files go into this script, and listing
						everything it does is rather tedious: Please open the script, read it, and email
						matthijsdejong1994@gmail.com if you have any questions.
	Casimir_single.m			This file contains the matlab code that can be called from python to generate a single frequency-sweep
						of a resonator with a Casimir nonlinearity. This is in essence the lines we plot throughout most of this
						paper (here in panel c, black lines).
	casimir_.m				This is the MATCONT function where we defined the Casimir nonlinearity.

	--- This subsection has the files dedicated to the frequency shift due to the Casimir force. We measure the final experimental frequency and
	don't know the exact position in absence of the Casimir force, nor the initial frequency that follows solely from the resonator stiffness.
	We also treat the final equilibrium position as a fit parameter, so we need to be flexible with that as well. The files mentioned in this
	subsection deal with that. ---
	Comparison_precalculated_recalculated_stable_positions.png 	Plot showing how the stable positions changed after optimizing the stable 
									points with MatCont
	Optimized_mech_freqs.txt		Simulation results where for various initial positions d, I calculated how much the Casimir force
						shifts the frequency and where the equilibrium position is. Then I optimize the initial frequency
						(i.e. resonator stiffness) and simulation start position such that the simulation converges to 
						the equilibrium position in a stable fashion, and that the final frequency is equal to the 
						experimentally measured frequency
	Optimized_mech_freqs_v2.txt		As above, with the results refined based on MatCont calculation rather than just Python script
						where I take a simple fit to the potential well.
	Precalculate_mechfreqs_v2.py		Python script where I simulate the change in frequency and position due to the Casimir force
						for given resonator initial position, elastic potential offset (d), and resonator stiffness. 
						This script outputs Optimized_mech_freqs.txt
	Recalculate_mechfreqs.py		Python script where I re-evaluate the output of the above script using MatCont.

	--- This subsection has the files dedicated to figuring out the correct mechanical amplitude for the different MW and VNA driving powers. ---
	Fit_Casimir_model_v27_Peakpowers.txt	This data file contains the experimentally measured maximum powers for all 6 sidebands + the
						off-resonant power.
	Fit_Casimir_model_v27_Simdata_*.txt	These files contain the simulated maximum powers for all 6 sidebands + the off-resonant power, and 
						the mechanical displacement amplitude from those simulations (which are matched to the experimental data).
						The numbers in the file name give the microwave drive power (resonant with the Cavity) for which the simulation applies.
	Fit_Casimir_model_v27_sim_data_match.py	This python script contains the code where I matched the experimentally observed peak powers to the 
						simulated peak powers.
	Fit_Casimir_model_v27_find_max_ampl.py	This python script contains some initial calculations/estimations for good phonon numbers to quickly
						reach the steady state in the simulations of the script above.


	--- This subsection has the files dedicated to simulating the power factor, where the drive becomes less efficient at large displacement
	amplitudes because the cavity frequency shifts so much that the drives become off-resonant (i.e. basic dispersive optomechanical coupling). ---
	Power_factors_9_12e-11.txt		This data file contains a list of mechanical amplitudes and the corresponding steady state cavity
						photon number under a given drive power (drive resonant with cavity frequency).
	Precalculate_powerfactor.py		This python script contains the simulation code to calculate the power factor, and it outputs the 
						above file.

in 'Figure 2':
	Figure2.pdf				This is figure 2 from the paper.
	Fit_Casimir_model_v27_fit_d_*.py	These python scripts have the calculation and plots to match the theory curves to the experimental
						data, like Fit_Casimir_model_v27_fit_d.py in the Figure 1 folder. They are based off the same calibrations
						but I have not copied those here.
	Casimir_*dBm_main_2.png			These plots are the individual panels of the figure, and the number in the name denotes the on-resonant
						cavity drive power.
	Power_factors_*.txt			These data files contain the various computed power factors (see 'Precalculate_powerfactor.py' in Figure 1)

in 'Figure 3':
	Figure3.pdf				This is figure 3 from the paper.
	Fit_Casimir_model_v27_fit_d.py		This script is a copy of the one in 'Figure 1'. It outputs the file below
	Casimir_fit_d_plot.png			This shows the fitted curve for d = 18 nm, and two ones with slightly different values.
	Precalculate_powerfactor.py		This script calculates the power factor (as on 'Figure 1')
	Casimir_Powerfactor_plot.png		This plot shows the effect of the taking the power factor into account on one of the experimentally measured
						data traces.
	Powerfactor_optomechanics_v10.png	This plot shows the power factor, the steady-state photon occupation of the cavity (normalized) for
						different mechanical amplitudes.	

in 'Figure 4':
	Figure4.pdf				This is figure 4 from the paper.
	Plot_allsidebands.py			Python script that takes in the data (dataset 2024-03-21-1858 in 'Raw and filtered data') and plots the
						sidebands, all separately and all overlapped.
	All_sidebands.png			Plot showing all sidebands as measured
	Overlapping_sidebands_v2.png		Plot showing the sidebands all overlap when suitably scaled and offset.

in 'Figure 5':
	Figure5.pdf				This is figure 5 from the paper.
	Calculate_Casimir_Pressure_distance_wplasma.py	This python script calculates the Casimir pressure according to the model developed by 
						Bimonte in PRA 99, 052507 (2019).
	Casimir_pressure_dist_v4.h5		This data file (hdf5) contains the output from the above script, distance, pressures (various models). The
						pressures may be too low by a factor of pi, which was an error in the script (can be validated by comparing the
						output values to those of Bimonte's paper).
	Fit_Casimir_scaling.py			This script fits the Casimir pressure-distance relation with a simple 1/d^n fit. 
	Plot_Casimir_pressure_distance.py	This script generates the plot for this figure.
	Plot_casimir_pressure_full_v2.png	Image showing the png version of this figure.	

in 'Figure 6':
	Figure6.pdf				This is figure 6 from the paper
	Setup_schematic_full.png		This is the png version of the same figure.

in 'Figure 7':
	Figure7.pdf				This is figure 7 from the paper
	Drum_simulations_full.png		This is the png version of the same figure
	2c2m_e3_small_drum.tiff			This is the AFM measurement (raw data) of the small drum (there are 2 drums and 2 cavities on this chip,
						but the larger drum collapsed and shorted one of the cavities).
	2c2m-e3-smalldrum_cleaned.png		This is the processed AFM measurement data
	Analyze_AFMdata_v2.py			Python script used to analyze/process the AFM measurement data
	Linecuts_through_smalldrum.png		Panel showing linecuts through the processed AFM measurement.

in 'Figure 8':
	Figure8.pdf				This is figure 8 from the paper
	Drum_simulations_2_full.png		Idem, but the png version
	Smalldrum_Louise_v5_out.mph		COMSOL model where the eigenmodes of the drum device were calculated, both at room 
						temperature and at 10 mK (including the thermal expansion)
	Smalldrum_Louise_mK_Mises_2.png		Linecut through the drum showing the static deflection at 10 mK due to thermal
						expansion
	Smalldrum_Louise_milliKelvin_inset_2	Inset of the simulation around the center of the drum
	Smalldrum_Louise_RT_Mises_2		Linecut through the drum showing the static deflection at room temperature
	Smalldrum_Louise_Roomtemperature_inset_2	Inset of the simulation around the center of the drum
	Smalldrum_Louise_Mode.png		Simulated mode shape of the first eigenmode of the drum. Panel B of the figure.

in 'Figure S1':
	FigureS1.pdf				Figure S1, showing the sideband amplitudes calculated via the repeated-scattering model
	Sideband_amps_MW-XX.npy			Numpy data files (read via numpy.load) containing the calculated sideband amplitudes
						for different powers of the cavity drive (file names).
	Calculate_scattering_rates_v4_messed.py	Python script used to calculate the scattering rates, sideband amplitudes, fit those
						to the data and plot the subpanels of this figure.
	Sidebands_fitted_scattering.png		Output from the python script.

in 'Figure S2 S3 S4':
	Figure S2.pdf				This is figure S2 from the paper
	Absolute_amplitude_calibration_full.png	Idem, but the PNG version.
	Figure S3.pdf				This is figure S3 from the paper
	Absolute_amplitude_calibration_offresonant.png	Idem, but the PNG version
	Figure S4.pdf				This is figure S4 from the paper
	Absolute_amplitude_calibration_amplitude.png	Idem, but the PNG version
	Fit_Casimir_model_v27_sim_data_match.py	Python script used to find the match between simulated maximum powers in the various sidebands,
						and the measured maximum powers in the sidebands (see also 'Figure 1').
	Fit_Casimir_model_v27_Simdata_*.txt	Data file with the simulated maximum powers of the various sidebands, for different
						resonant MW drive powers.
	Simulated_vs_measured_blue_1_v2.png	Output of the python file, one of the panels of Figure S2
	Simulated_vs_measured_red_1_v2.png	Output of the python file, the other panel of Figure S2
	Simulated_vs_measured_pedestal_v2.png	Output of the python file, the panel of Figure S3
	Simulated_vs_measured_mech_amplitude.png	Output of the python file, the panel of Figure S4.

in 'Figure S5':
	Figure S5.pdf				This is figure S5 from the paper
	Casimir_n_scaling_full.png		The PNG version of the above.
	Calculate_Casimir_Pressure_distance_wplasma.py	Python script used to calculate the Casimir pressure at various distances between
						two plates, with various material properties. Based on the model of Bimonte in PRA
						99, 052507 (2019). The script output may be too low by a factor of pi, which is a known
						issue in various versions of this script.
	Casimir_pressure_dist_v*.h5		Data files containing the Casimir pressures calculated at different distances for some
						various sets of material parameters. All are outputs from the above python script.
	Fit_Casimir_scaling_v2.py		Python script used to plot figure S5
	Plot_casimir_pressure_epsilon.png	One of the panels of figure S5, output from the above script
	Plot_casimir_pressure_omegap.png	Idem
	Plot_casimir_pressure_RRR.png		Idem

in 'Figure S6':
	Figure S6.pdf				This is figure S6 from the paper
	Alternative_fits_full.png		The PNG version of the above
	Casimir_fit_for_omegap_*.png		These are the individual panels of the above
	Casimir_pressure_dist_*.h5		These data files contain the Casimir pressures for various distances for the different
						material parameters used in this figure.
	Fit_Casimir_alternative_dn.py		Python script that fits the Casimir force to the measured data assuming the plasma
						frequency is 15 eV. It works similarly to Fit_Casimir_model_v27 which is explained in
						the Figure 1 folder, and some of the data/other files can be found there as well
	Fit_Casimir_alternative_dn2.py		Idem, but assuming the plasma frequency is 10 eV.
	Fit_Casimir_model_v27_Simdata_-20dBm_MW_altfreq.txt	Data file containing the simulated maximum amplitudes of the sidebands	
						which are calibrated as described in the Figure 1 folder/paper text. 
	Optimized_mech_freqs_altCas*.txt	Data file containing the calculated start positions/resonator stiffnesses that result in 
						a stable resonance and correct final frequency for various start positions of d.

in 'Figure S7 S8':
	Figure S7.pdf				This is figure S7 from the paper
	Figure S8.pdf				This is figure S8 from the paper
	Optimized_mech_freqs.txt		Simulation results where for various initial positions d, I calculated how much the Casimir force
						shifts the frequency and where the equilibrium position is. Then I optimize the initial frequency
						(i.e. resonator stiffness) and simulation start position such that the simulation converges to 
						the equilibrium position in a stable fashion, and that the final frequency is equal to the 
						experimentally measured frequency
	Optimized_mech_freqs_v2.txt		As above, with the results refined based on MatCont calculation rather than just Python script
						where I take a simple fit to the potential well.
	Precalculate_mechfreqs_v2.py		Python script where I simulate the change in frequency and position due to the Casimir force
						for given resonator initial position, elastic potential offset (d), and resonator stiffness. 
						This script outputs Optimized_mech_freqs.txt
	Recalculate_mechfreqs.py		Python script where I re-evaluate the output of the above script using MatCont.
	Perturbed_versus_unperturbed_frequency_distance_2.png	PNG version of figure S7
	Stable_position_calculated_for_final_frequency.png	PNG version of figure S8

in 'Figure S9':
	Figure S9.pdf				This is figure S9 from the paper
	cas_parameters_v3.png			Idem, but the PNG version
	casimir_.m				MATCONT function where we defined the Casimir nonlinearity
	Casimir_MJ_*.m				MATLAB script where we ran the simulations shown in this figure
	Casimir_vary_n.py			Python script ran for the revised figure, where the n-varying simulation was re-done
	Casimir_for_different_n.png		Results from the above python script
	dmapfig.png				Results from the simulation where the damping rate was varied
	distfig.png				Results from the simulation where the distance was varied
	drivefig.png				Results from the simulation where the drive strength was varied

in 'Figure S10':
	Figure S10.pdf				This is Figure S10 from the paper
	Casimir_01.png				Microscope image, post-experiment.
	E2-100X.png				Microscope image, pre-experiment.
	Sonnet_Casimir_drum.png			Screenshot of Sonnet simulation showing the cavity design.
	MMCC_*.csv/son				Sonnet files used in simulating the microwave cavity, for 15.0 and 15.1 nm gap distances
	in sondata folder			Specific Sonnet files used in simulating the microwave cavity.

in 'Figure S11':
	Figure S11.pdf				This is figure S11 from the paper
	Cavity_calibration_full.png		Idem, but the png version
	2024-03-13-1405_MeasCNEMScool_SMA100B_test.m	Matlab script used to measure the cavity response
	2024-03-13-1405-W1P3-E3-SidebandCooling-[-40_-10]_XX.fig	Matlab figures showing the cavity response and spectra
						(Note: The data for these figures is stored in the raw data folder)
	Cavity fit try X.png			Intermediate results of the cavity calibration.
	Sideband_cooling_v2.py			Python script used to calibrate the microwave cavity (also does more, used for other
						figures involving the calibration of the microwave optomechanical system).

in 'Figure S12':
	Figure S12.pdf				This is figure S12 from the paper
	2024-03-12-1830.mat			Data file containing the MATLAB working environment after the thermalization measurement
	2024-03-12-1830_MeasCNEMScool_SMA100B_test.m	MATLAB script file of the measurement of the thermalization.
	Analyze_mode_thermalization.py		Python script containing the analysis of the above measurement data
	Thermal *.png				Individual thermalization analysis results of amplitude, frequency, linewidth and peak amplitude.

in 'Figure S13':
	FigureS13.pdf				Figure S13, showing the scattering amplitude used to calibrate single-photon opto-
						mechanical coupling rate g_0
	Calibrated_g0_full.png			Idem, but the png version
	Fit_g0_QNOF_v3.py			Python script used to fit and analyze the g_0 data by comparing the scattering rate
						from the sideband peak of known amplitude.
	W1P3-E3-mesureg0_XXHz_Pump[XX]_extractedworkspace.mat	Data measured showing the drive and scattered peaks for various drive
						powers.

in 'Figure S14':
	FigureS14.pdf				This is figure S14 from the paper, showing the calibrated frequency, linewidth, cavity- and mechanical 
						occupation numbers. 
	Sideband_calibration_full.png		Idem, but the png version
	2024-03-13-1405.mat			Raw data from the red-sideband drive power sweep used to calibrate the system
	2024-03-13-1405_MeasCNEMScool_SMA100B_test.m	Measurement script used to run the red-sideband calibration
	Final mechanical parameters fitted_forpaper.png	Panel A of Figure S11, showing the frequency and linewidth for various
						powers of the red sideband drive.
	
	Final photon and phonon number.png	Panel B of Figure S11, showing the mechanical and cavity occupation numbers for
						various powers of the red sideband drive.
	Sideband_cooling_v2.py			Python script used to analyze and plot the red-sideband calibration measurement.
	
in 'Figure S15':
	FigureS15.pdf				This is figure S15 from the paper, showing the before and after filtering.
	Data_filtering_full_v2.png		Idem, but the png version.
	2024-03-19-1759.mat			Matlab workspace file containing raw data output from the measurement script		
	2024-03-19-1759_filtdata_down.npy	Numpy data file containing the filtered powers of the drive frequency sweep in 
						decreasing frequency direction. Read with numpy.load()
	2024-03-19-1759_filtdata_up.npy		Numpy data file containing the filtered powers of the drive frequency sweep in 
						increasing frequency direction.
	2024-03-19-1759_filtfreqs.npy		Numpy data file containing the frequencies at which the red-sideband drive was parked
						and at which the data is filtered.
	2024-03-19-1759_QNOF_VNApwr_sweep.m	Matlab script used to run the measurement that generated the raw data.
	Blue sideband X.png			Plots of the raw and filtered data for the different sidebands on the blue side
						(higher frequency than the cavity)
	Plot_Matlabdata.py			Python script used to analyze and filter the data.
	Red sideband X.png			Plots of the raw and filtered data for the differend sidebands on the red side
						(lower frequency than the cavity)
	Unfiltered data.png / zoom.png		Plots of the filtered and unfiltered data, panel A of Figure S12.

in 'Figure S16':
	FigureS16.pdf				This is figure S16 from the paper, showing the electrostatic force as a function of distance for a 1V	
						static potential difference between the drum top and bottom plates.
	Electrostatic_force_full.png		Idem, but the png version
	Electrostatic_energy_1V.txt		Text file containing the numerical data of electrostatic energy as a function of
						distance between the plates
	Plot_electrostatic_nonlinearity.py	Python script for plotting Figure S13
	Smalldrum_Louise_v5_electrostatics.mph	COMSOL file containing the simulation setup for calculation of the electrostatic
						force as a function of energy. (run with COMSOL 6.3)

in 'Figure S17':
	FigureS17.png				This is figure S17 from the paper, showing the pressure from potential patches in the plate-plate
						geometry (for various patch sizes). The patch pressure calculations were done using the FEM
						model developed in M.H.J. de Jong et al., https://arxiv.org/abs/2408.16323, and the calculation
						scripts can be found in the data set associated with that paper.
	Casimir_pressure_dist_v6_largedist*.h5	Data file containing the calculated Casimir pressure for various distances (calculation based on 
						model by Bimonte, PRA 99, 052507 (2019), again may be missing the factor of pi mentioned elsewhere
						in this readme).
	Casimir_vs_patch_pressure.png		Raw version of FigureS17.png
	Simulate_Speak_Trenkel_plate_plate_v3.py	Python script used to simulate (numerically) the plate-plate patch force.
	Sim_results_plateplate_v6.txt		Data file containing the output of the above script.

in 'Figure S18 S19 S20':
	FigureS18.png				This is figure S18 from the paper
	FigureS19.png				This is figure S19 from the paper
	FigureS20.png				This is figure S20 from the paper
	W1P3_E4_measarea3 0003_corrected_data.png	Figure S14 showing the AFM micrograph of the potential and height, showing
						a part of the sample from the same batch as the one used for the measurements in 
						this paper (i.e. same material condition/treatment).
	W1P3_E4_measarea3 0003_Potential_height_distribution.png	Figure S15, showing the distribution of potentials and 
						heights for the data shown in Figure S14.
	W1P3_E4_measarea3 0003_featuresize_uncorrected.png	Figure S16, showing the autocorrelation of the potential shown in
						Figure S14.
	Analyze_v2.py				Python script used to analyze and plot the data
	Determine_featuresize_2.py		Python script used to analyze and plot the autocorrelation/feature size of potential
	Determine_featuresize_2h.py		Python script used to analyze and plot the autocorrelation/feature size of height
	W1P3_E4_drums.py			Microscope image from the AFM showing the area where the measurements were taken
	W1P3_E4_measarea3 0003.h5		Raw measured data output from the Olympus Jupiter AFM used for measuring.
	W1P3_E4_measarea3 0003_corrected_height.png	Image of the sample height
	W1P3_E4_measarea3 0003_corrected_potential.png	Image of the sample potential
	W1P3_E4_measarea3 0003_FFT_potential.png	FFT of the image of sample potential
	W1P3_E4_measarea3 0003_FFT_potential_windowed.png	Windowed version of the above showing the center only.
	W1P3_E4_measarea3 0003_height_autocorr.csv	Comma-separated value file for height autocorrelations
	W1P3_E4_measarea3 0003_potential_autocorr.csv	Comma-separated value file for potential autocorrelations
	W1P3_E4_measarea3 0003_Potential_autocorrelation.png	Plot of the autocorrelation of potential
	W1P3_E4_measarea3 0003_potential_color_centered.png	Color version of the image showing potential
	W1P3_E4_measarea3 0003_Potential_correlation_linecuts.png	Line-cut along x- and y-axis of the dataset showing the
						Autocorrelation of potential.
	W1P3_E4_measarea3 0003_Potential_height_correlation.png	Plot showing the cross-correlation between height and potential
	W1P3_E4_measarea3 0003_potential_height_scatterplot.png	Scatterplot showing the cross-correlation between height and
						potential
	W1P3_E4_measarea3 0003_potheight_crosscorr.csv	Comma-separated value file for the height-potential cross correlation.
	W1P3_E4_measarea3 0003_raw_data.png	Plot of the raw, untreated data of height and potential.

in 'Figure S21':
	FigureS21.png				This is figure S21 of the paper.
	rsznb8vna_and_mg3692c_and_ms2682a_segm_sweep.py	Measurement script used to probe the second device (in new setup, different fridge)
	Readme.txt				Note for associating the dataset within the larger data management system used in the group
	Analyze_Casimir_segm_sweep.py		Python script for analyzing the measurements
	Analyze_Casimir_v1.py			Second Python script for analyzing the measurements
	Compute_maximum_amplitudes.py		Python script to calculate the maximum amplitudes reached in this second device
	Data.npz				Zipped numpy data file containing the output of the measurement script.
	Extracted_maxima.txt			Maximum sideband powers extracted from the measurement
	Filter_SAdata.py			Python script to filter the output of the rsznb8vna measurement script in a similar
						fashion to the script included in Figure S15 folder.
	Filtered_data.npz			Zipped numpy data file containing the filtered measurement data	
	Fitted_Casimir_curve_19_1nm.png		PNG version of Figure S21
	Optimized_mech_freqs_dev13MHz.txt	Data file containing the calculated initial position, resonator stiffness and start point 
						such that the final resonator including the Casimir force comes to the correct frequency and
						with a stable initial position so we arrive in the locally stable position.
	Power_factors_*.txt			Data file containing the calculated power factors for the new cavity linewidth and frequency, 
						including estimation of the cable losses and new setup attenuation.
	Precalculate_mechfreqs_v2.py		Python script to obtain the optimized mechanical frequencies for 13 MHz
	Precalculate_powerfactor.py		Python script to obtain the power factors for this device

in 'Raw and filtered data':
	XXXX-XX-XX-XXXX.mat 			MATLAB workspace file containing the data as measured.
	XXXX-XX-XX-XXXX_filtdata_down.npy	Filtered data, power as measured in the downward frequency direction
	XXXX-XX-XX-XXXX_filtdata_up.npy		Filtered data, power as measured in the upward frequency direction
	XXXX-XX-XX-XXXX_filtfreqs.npy		Filtered data, frequencies
	XXXX-XX-XX-XXXX_QNOF_VNApwr_sweep.m	Matlab script for running the measurement.

	