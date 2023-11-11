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
intensoids_clumping = Library_Parameters.intensoids_clumping
tag_mass = Library_Parameters.tag_mass
fast_iso = Library_Parameters.fast_iso
imp_exp_library = Library_Parameters.imp_exp_library
only_gen_lib = Library_Parameters.only_gen_lib

multithreaded_analysis = Analysis_Parameters.multithreaded_analysis
accuracy_unit = Analysis_Parameters.accuracy_unit
accuracy_value = Analysis_Parameters.accuracy_value
threshold = Analysis_Parameters.threshold
ret_time_interval = Analysis_Parameters.ret_time_interval
peak_width = Analysis_Parameters.peak_width
peak_area_threshold = Analysis_Parameters.peak_area_threshold
max_peaks = Analysis_Parameters.max_peaks
peaks_relative_intensity = Analysis_Parameters.peaks_relative_intensity
samples_list = Analysis_Parameters.samples_list
analyze_ms2 = Analysis_Parameters.analyze_ms2
output_file = Analysis_Parameters.output_file
multithreaded_execution = False

##-----------------------------------------------------------------------------

begin_time = datetime.datetime.now()
samples_names = Execution_Functions.sample_names(samples_list)

if imp_exp_library[0]:
    from glycans_library import full_library
else:
    full_library = 'Empty' 

library = Execution_Functions.imp_exp_gen_library(full_library, multithreaded_analysis, samples_names, custom_glycans_list, min_max_monos, min_max_hex, min_max_hexnac, min_max_sia, min_max_fuc, min_max_ac, min_max_gc, force_nglycan, max_adducts, max_charges, intensoids_clumping, tag_mass, fast_iso, imp_exp_library)
print('Library length: '+str(len(full_library)))
if only_gen_lib:
    input("Library generated/imported in "+str(datetime.datetime.now()-begin_time)+". "+
          "Check it in glycans_library.py. If you wish to analyze files, set 'only_gen"+
          "_lib' to False and input remaining parameters.\n Press Enter to exit.")
    sys.exit()
Execution_Functions.print_sep()
tolerance = Execution_Functions.tolerance(accuracy_unit, accuracy_value)
data = Execution_Functions.list_of_data(samples_list)
ms1_index = Execution_Functions.index_ms1_from_file(data)
if analyze_ms2:
    ms2_index = Execution_Functions.index_ms2_from_file(data)
lib_size = len(library)
analyzed_data = Execution_Functions.analyze_files(library, lib_size, data, ms1_index, tolerance, threshold, peak_width, ret_time_interval, max_peaks, peaks_relative_intensity, peak_area_threshold)
Execution_Functions.print_sep()
lines_to_print = Execution_Functions.generate_lines_to_print(analyzed_data, samples_names, peak_width, multithreaded_execution)
Execution_Functions.print_sep()
Execution_Functions.print_to_file(lines_to_print, output_file)
print(datetime.datetime.now() - begin_time)
input('Execution complete. Press Enter to exit.')
