from Modules import Execution_Functions
from Parameters import Library_Parameters, Analysis_Parameters
import sys
import datetime


custom_glycans_list = Library_Parameters.custom_glycans_list
min_max_monos = Library_Parameters.min_max_monos
min_max_hex = Library_Parameters.min_max_hex
min_max_hexnac = Library_Parameters.min_max_hexnac
min_max_sia = Library_Parameters.min_max_sia
min_max_fuc = Library_Parameters.min_max_fuc
min_max_ac = Library_Parameters.min_max_ac
min_max_gc = Library_Parameters.min_max_gc
force_nglycan = Library_Parameters.force_nglycan
max_adducts = Library_Parameters.max_adducts
max_charges = Library_Parameters.max_charges
tag_mass = Library_Parameters.tag_mass
fast_iso = Library_Parameters.fast_iso
high_res = Library_Parameters.high_resolution_isotopic_dist
imp_exp_library = (Library_Parameters.imp_library, Library_Parameters.exp_library)
only_gen_lib = Library_Parameters.only_gen_lib

multithreaded_analysis = Analysis_Parameters.multithreaded_analysis
analyze_ms2 = Analysis_Parameters.analyze_ms2
accuracy_unit = Analysis_Parameters.accuracy_unit
accuracy_value = Analysis_Parameters.accuracy_value
ret_time_interval = Analysis_Parameters.ret_time_interval
min_isotopologue_peaks = Analysis_Parameters.min_isotopologue_peaks
min_ppp = Analysis_Parameters.min_points_per_peak
max_ppm = Analysis_Parameters.max_ppm
iso_fit_score = Analysis_Parameters.isotopic_fitting_score
curve_fit_score = Analysis_Parameters.curve_fitting_score
s_to_n = Analysis_Parameters.signal_to_noise
custom_noise = Analysis_Parameters.custom_noise_level
samples_list = Analysis_Parameters.samples_list
save_path = Analysis_Parameters.results_save_path
reanalysis = Analysis_Parameters.reanalysis

verbose = False
multithreaded_execution = (False, 0)

##-----------------------------------------------------------------------------

begin_time = datetime.datetime.now()

if reanalysis[0]:
    Execution_Functions.output_filtered_data(curve_fit_score,
                                             iso_fit_score,
                                             s_to_n,
                                             max_ppm,
                                             reanalysis,
                                             save_path,
                                             multithreaded_analysis,
                                             analyze_ms2[0])

else:
    if multithreaded_execution[0]:
        print('Multithreaded Execution: '+str(multithreaded_execution[1]))
    samples_names = Execution_Functions.sample_names(samples_list)
    library = Execution_Functions.imp_exp_gen_library(multithreaded_analysis,
                                                      multithreaded_execution,
                                                      samples_names,
                                                      custom_glycans_list,
                                                      min_max_monos,
                                                      min_max_hex,
                                                      min_max_hexnac,
                                                      min_max_sia,
                                                      min_max_fuc,
                                                      min_max_ac,
                                                      min_max_gc,
                                                      force_nglycan,
                                                      max_adducts,
                                                      max_charges,
                                                      tag_mass,
                                                      fast_iso,
                                                      high_res,
                                                      imp_exp_library,
                                                      only_gen_lib)
    print('Library length: '+str(len(library)))
    Execution_Functions.print_sep()
    tolerance = Execution_Functions.tolerance(accuracy_unit,
                                              accuracy_value)
    data = Execution_Functions.list_of_data(samples_list)
    print("Indexing spectra...", end = "", flush = True)
    ms1_index = Execution_Functions.index_ms1_from_file(data)
    if analyze_ms2[0]:
        data = Execution_Functions.list_of_data(samples_list)
        ms2_index = Execution_Functions.index_ms2_from_file(data)
    print("Done!")
    Execution_Functions.print_sep()
    lib_size = len(library)
    analyzed_data = Execution_Functions.analyze_files(library,
                                                      lib_size,
                                                      data,
                                                      ms1_index,
                                                      tolerance,
                                                      ret_time_interval,
                                                      min_isotopologue_peaks,
                                                      min_ppp,
                                                      max_charges,
                                                      custom_noise,
                                                      verbose)
    if analyze_ms2[0]:
        Execution_Functions.print_sep()
        analyzed_data = Execution_Functions.analyze_ms2(ms2_index, 
                                                        data, 
                                                        analyzed_data, 
                                                        ret_time_interval,
                                                        tolerance,
                                                        min_max_monos,
                                                        min_max_hex,
                                                        min_max_hexnac,
                                                        min_max_sia,
                                                        min_max_fuc,
                                                        min_max_ac,
                                                        min_max_gc,
                                                        max_charges,
                                                        tag_mass,
                                                        analyze_ms2[1])
    Execution_Functions.print_sep()
    Execution_Functions.arrange_raw_data(analyzed_data,
                                         samples_names,
                                         multithreaded_analysis,
                                         multithreaded_execution,
                                         analyze_ms2[0])
    Execution_Functions.print_sep()
    Execution_Functions.output_filtered_data(curve_fit_score,
                                             iso_fit_score,
                                             s_to_n,
                                             max_ppm,
                                             reanalysis,
                                             save_path,
                                             multithreaded_analysis,
                                             analyze_ms2[0])
                                             
input('Execution complete. Time elapsed: '+str(datetime.datetime.now() - begin_time)+'\nPress Enter to exit.')
