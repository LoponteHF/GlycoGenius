# GlycoGenius: Glycomics Data Analysis Tool
# Copyright (C) 2023 by Hector Franco Loponte
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or 
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. It is accessible within the program files
# or by typing 'license' after running it stand-alone in the terminal
# by typing 'glycogenius'. If not, see <https://www.gnu.org/licenses/>.

from . import Execution_Functions
from . import General_Functions
from . import CLI
from pathlib import Path
from pyteomics import mass
import time
import sys
import os
import datetime
import configparser
import shutil
import tempfile
import pathlib
import importlib

# date = datetime.datetime.now()
# begin_time = str(date)[2:4]+str(date)[5:7]+str(date)[8:10]+"_"+str(date)[11:13]+str(date)[14:16]+str(date)[17:19]
# temp_folder = os.path.join(tempfile.gettempdir(), "gg_"+begin_time)

def config_handler(from_GUI = False, param_file_path = ''):
    '''A function that handles the input of configs through a pipelined
    parameters .ini file and turns them into arguments for the functions
    used in core module.
    
    Parameters
    ----------
    none
    
    Uses
    ----
    pyteomics.mass
        Used for calculating the molecular weight of molecules.
    pathlib.Path
        Class of paths that can be used to create, change and access 
        directories in the computer.
    
    Returns
    -------
    functions arguments
        Organized arguments for the various functions used in core module.
    '''

    date = datetime.datetime.now()
    begin_time = str(date)[2:4]+str(date)[5:7]+str(date)[8:10]+"_"+str(date)[11:13]+str(date)[14:16]+str(date)[17:19]
    temp_folder = os.path.join(tempfile.gettempdir(), "gg_"+begin_time)
    
    custom_glycans_list = [False, []]
    min_max_monos = [0, 0]
    min_max_hex = [0, 0]
    min_max_hexnac = [0, 0]
    min_max_sia = [0, 0]
    min_max_fuc = [0, 0]
    min_max_ac = [0, 0]
    min_max_gc = [0, 0]
    force_nglycan = False
    max_adducts = {}
    adducts_exclusion = []
    max_charges = 0
    reducing_end_tag = 0.0
    permethylated = False
    reduced = False
    lactonized_ethyl_esterified = False
    fast_iso = True
    high_res = False
    internal_standard = 0.0
    imp_exp_library = [False, False]
    exp_lib_name = ''
    library_path = ''
    only_gen_lib = False

    multithreaded_analysis = True
    number_cores = 'all'
    analyze_ms2 = [False, False, False]
    reporter_ions = []
    tolerance = ['mz', 0.01]
    ret_time_interval = [0, 99999, 0.2]
    rt_tolerance_frag = 0.2
    min_isotopologue_peaks = 2
    min_ppp = [False, 0]
    close_peaks = [False, 3]
    align_chromatograms = True
    percentage_auc = 0.1
    max_ppm = 10
    iso_fit_score = 0.9
    curve_fit_score = 0.9
    s_to_n = 3
    custom_noise = [False, []]
    samples_path = ''
    save_path = ''
    plot_metaboanalyst = [False, []]
    compositions = True
    iso_fittings = False
    reanalysis = False
    reanalysis_path = ''
    output_plot_data = False

    samples_list = []
    samples_names = []
    
    CLI.print_header(False)
    config = configparser.ConfigParser()
    configs = ""
    if from_GUI:
        with open(param_file_path, "r") as f:
            for line in f:
                configs+=line
            f.close()
    else:
        for line in sys.stdin:
            configs+=line
    config.read_string(configs)
    
    #working directory
    save_path = config['running_modes']['working_directory']
    save_path = save_path.strip()
    save_path = save_path.strip("'")
    save_path = save_path.strip("\"")
    for i_i, i in enumerate(save_path):
        if i == "\\":
            save_path = save_path[:i_i]+"/"+save_path[i_i+1:]
    if save_path[-1] != "/":
        save_path = save_path+"/"
    Path(save_path).mkdir(exist_ok = True, parents = True)
    
    #multithreading settings
    multithreaded_analysis = config['running_modes'].getboolean('use_multiple_CPU_cores')
    number_cores = (config['running_modes']['number_cores']).strip()
    
    #running mode
    general_mode = config['running_modes']['mode']
    
    if general_mode == 'reanalysis':
        reanalysis = True
        reanalysis_path = config['running_modes']['file_for_reanalysis']
        reanalysis_path = reanalysis_path.strip()
        reanalysis_path = reanalysis_path.strip("'")
        reanalysis_path = reanalysis_path.strip("\"")
        for i_i, i in enumerate(reanalysis_path):
            if i == "\\":
                reanalysis_path = reanalysis_path[:i_i]+"/"+reanalysis_path[i_i+1:]
        
    elif general_mode == 'library':
        only_gen_lib = True
        exp_lib_name = config['library_building_modes']['exported_library_name'].strip()
        
    elif general_mode != 'reanalysis' and general_mode != 'library' and general_mode != 'analysis':
        print("Wrong mode selected. Use 'analysis',\n'reanalysis' or 'library' in mode under\n'running_modes', in the parameters file.")
        print("Close the window or press CTRL+C to exit.")
        try:
            while True:
                time.sleep(3600)
        except KeyboardInterrupt:
            os._exit(1)
        
    
    #library building mode
    library_mode = config['library_building_modes']['mode']
    
    if library_mode == 'import_library':
        imp_exp_library[0] = True
        
        library_path = config['library_building_modes']['import_library_path']
        library_path = library_path.strip()
        library_path = library_path.strip("'")
        library_path = library_path.strip("\"")
        for i_i, i in enumerate(library_path):
            if i == "\\":
                library_path = library_path[:i_i]+"/"+library_path[i_i+1:]
                
    elif library_mode == 'custom_library':
        custom_glycans_list[0] = True
        
        if not reanalysis:
            if (len(config['library_building_modes']['custom_glycans_list'].split("\\")) > 1 or len(config['library_building_modes']['custom_glycans_list'].split("/")) > 1):
                try:
                    temp_path_custom_glycans = config['library_building_modes']['custom_glycans_list'].strip("\"").strip("'")
                    for i_i, i in enumerate(temp_path_custom_glycans):
                        if i == "\\":
                            temp_path_custom_glycans = temp_path_custom_glycans[:i_i]+"/"+temp_path_custom_glycans[i_i+1:]
                    temp_custom_glycans_list = ""
                    with open(temp_path_custom_glycans.strip("\""), 'r') as f:
                        for i in f:
                            temp_custom_glycans_list += i
                    f.close()
                    custom_glycans = temp_custom_glycans_list.split(",")
                    if len(custom_glycans) == 1:
                        custom_glycans = temp_custom_glycans_list.split("\n")
                    to_remove = []
                    if len(custom_glycans) > 1:
                        for i_i, i in enumerate(custom_glycans):
                            custom_glycans[i_i] = i.strip()
                            if len(i) == 0:
                                to_remove.append(i_i)
                    for i in sorted(to_remove, reverse = True):
                        del custom_glycans[i]
                except:
                    print("Custom glycans file not found. Check the path to\nthe file you've inputted in the parameters file.")
                    print("Close the window or press CTRL+C to exit.")
                    try:
                        while True:
                            time.sleep(3600)
                    except KeyboardInterrupt:
                        os._exit(1)
            else:
                custom_glycans = config['library_building_modes']['custom_glycans_list'].split(",")
                for i_i, i in enumerate(custom_glycans):
                    custom_glycans[i_i] = i.strip()
                    if len(i) == 0:
                        custom_glycans = custom_glycans[:i_i]+custom_glycans[i_i+1:]
            custom_glycans_list[1] = custom_glycans
    elif library_mode != 'custom_library' and library_mode != 'import_library' and library_mode != 'generate_library':
        print("Wrong mode selected. Use 'import_library',\n'custom_library' or 'generate_library' in mode\nunder 'library_building_modes', in the\nparameters file.")
        print("Close the window or press CTRL+C to exit.")
        try:
            while True:
                time.sleep(3600)
        except KeyboardInterrupt:
            os._exit(1)
    
    elif library_mode == 'generate_library': 
        total_monosaccharides = config['library_building_modes']['total_monosaccharides'].split(',')
        min_max_monos = [int(total_monosaccharides[0]), int(total_monosaccharides[1])]
        
        total_hexoses = config['library_building_modes']['hexoses'].split(',')
        min_max_hex = [int(total_hexoses[0]), int(total_hexoses[1])]
        
        total_hexnacs = config['library_building_modes']['hexnacs'].split(',')
        min_max_hexnac = [int(total_hexnacs[0]), int(total_hexnacs[1])]
        
        total_sialics = config['library_building_modes']['sialic_acids'].split(',')
        min_max_sia = [int(total_sialics[0]), int(total_sialics[1])]
        
        total_fucoses = config['library_building_modes']['fucoses'].split(',')
        min_max_fuc = [int(total_fucoses[0]), int(total_fucoses[1])]
        
        total_neu5ac = config['library_building_modes']['neu5ac'].split(',')
        min_max_ac = [int(total_neu5ac[0]), int(total_neu5ac[1])]
        
        total_neu5gc = config['library_building_modes']['neu5gc'].split(',')
        min_max_gc = [int(total_neu5gc[0]), int(total_neu5gc[1])]
        
    imp_exp_library[1] = config['library_building_modes'].getboolean('export_library')
    if imp_exp_library[1]:
        exp_lib_name = config['library_building_modes']['exported_library_name'].strip()
    force_nglycan = config['common_library_building_settings'].getboolean('force_nglycan')
    max_adducts = General_Functions.form_to_comp(config['common_library_building_settings']['max_adducts'])
    adducts_exc = config['common_library_building_settings']['adducts_exclusion'].split(",")
    if len(adducts_exc) > 0:
        for i in adducts_exc:
            adducts_exclusion.append(General_Functions.form_to_comp(i.strip()))
    max_charges = int(config['common_library_building_settings']['max_charges'])
    reducing_end_tag = config['common_library_building_settings']['reducing_end_tag']
    permethylated = config['common_library_building_settings'].getboolean('permethylated')
    reduced = config['common_library_building_settings'].getboolean('reduced')
    lactonized_ethyl_esterified = config['common_library_building_settings'].getboolean('aminated_ethyl_esterified')
    fast_iso = config['common_library_building_settings'].getboolean('fast_iso')
    high_res = config['common_library_building_settings'].getboolean('high_resolution_isotopic_dist')
    internal_standard = config['common_library_building_settings']['internal_standard_mass']
    if len(internal_standard.strip()) == 0:
        internal_standard = 0.0
    else:
        internal_standard = float(internal_standard)
        
        
    #analysis running mode
    if not reanalysis and not only_gen_lib:
        samples_path = config['running_modes']['samples_directory']
        samples_path = samples_path.strip()
        samples_path = samples_path.strip("'")
        samples_path = samples_path.strip("\"")
        for i_i, i in enumerate(samples_path):
            if i == "\\":
                samples_path = samples_path[:i_i]+"/"+samples_path[i_i+1:]
        if samples_path[-1] != "/":
            samples_path = samples_path+"/"
        samples_list = Execution_Functions.samples_path_to_list(samples_path)
        if len(samples_list) == 0 and not from_GUI:                                             
            if not os.isatty(0):
                Execution_Functions.print_sep()
                print("No sample files to analyze.")
                os._exit(1)
            else:
                Execution_Functions.print_sep()
                input("No sample files to analyze. Press Enter\nto exit.")
                os._exit(1)
        
        samples_names = Execution_Functions.sample_names(samples_list)
        print("Sample files detected: "+str(len(samples_names)))
        for i in samples_names:
            print("--> "+i)
        Execution_Functions.print_sep()
            
        ms2_analysis = config['analysis_parameters'].getboolean('analyze_ms2')
        if ms2_analysis:
            analyze_ms2 = (ms2_analysis, config['analysis_parameters'].getboolean('force_fragments_to_glycans'), config['analysis_parameters'].getboolean('unrestricted_fragments'))
            if 'ret_time_tolerance_ms2' in config['analysis_parameters']:
                rt_tolerance_frag = float(config['analysis_parameters']['ret_time_tolerance_ms2'])
        accuracy_unit = config['analysis_parameters']['accuracy_unit']
        accuracy_value = float(config['analysis_parameters']['accuracy_value'])
        tolerance = [accuracy_unit, accuracy_value]
        ret_time_config = config['analysis_parameters']['ret_time_interval'].split(",")
        ret_time_interval[0] = float(ret_time_config[0])
        ret_time_interval[1] = float(ret_time_config[1])
        if 'ret_time_tolerance' in config['analysis_parameters']:
            ret_time_interval[2] = float(config['analysis_parameters']['ret_time_tolerance'])
        custom_ppp = config['analysis_parameters'].getboolean('custom_min_points_per_peak')
        if custom_ppp:
            min_ppp = (config['analysis_parameters'].getboolean('custom_min_points_per_peak'), int(config['analysis_parameters']['number_points_per_peak']))
        limit_peaks = config['analysis_parameters'].getboolean('limit_peaks_picked')
        if limit_peaks:
            close_peaks = (config['analysis_parameters'].getboolean('limit_peaks_picked'), int(config['analysis_parameters']['max_number_peaks']))
        if 'noise_levels' in config['analysis_parameters']:
            noise_levels = config['analysis_parameters']['noise_levels'].split(",")
            for i_i, i in enumerate(noise_levels):
                noise_levels[i_i] = int(i.strip())
            custom_noise = (config['analysis_parameters'].getboolean('custom_noise_level'), noise_levels)
        
    #post-analysis for analysis running mode and reanalysis running mode
    if not only_gen_lib:
        if analyze_ms2[0]:
            reporter_ions = config['post-analysis/reanalysis']['filter_ms2_by_reporter_ions'].split(",")
            for i_i, i in enumerate(reporter_ions):
                reporter_ions[i_i] = i.strip()
                if len(i) == 0:
                    reporter_ions = reporter_ions[:i_i]+reporter_ions[i_i+1:]
        align_chromatograms = config['post-analysis/reanalysis'].getboolean('align_chromatograms')
        percentage_auc = float(config['post-analysis/reanalysis']['auc_percentage_threshold'])/100
        max_ppm = int(config['post-analysis/reanalysis']['max_ppm_threshold'])
        iso_fit_score = float(config['post-analysis/reanalysis']['isotopic_fitting_score_threshold'])
        curve_fit_score = float(config['post-analysis/reanalysis']['curve_fitting_score_threshold'])
        s_to_n = float(config['post-analysis/reanalysis']['signal_to_noise_threshold'])
        compositions = config['post-analysis/reanalysis'].getboolean('analyze_compositions')
        plot_metaboanalyst_file = config['post-analysis/reanalysis'].getboolean('output_metaboanalyst_file')
        if plot_metaboanalyst_file:
            metaboanalyst_groups = config['post-analysis/reanalysis']['metaboanalyst_groups'].split(",")
            for i_i in range(len(metaboanalyst_groups)-1, -1, -1):
                metaboanalyst_groups[i_i] = metaboanalyst_groups[i_i].strip()
                metaboanalyst_groups[i_i] = metaboanalyst_groups[i_i].strip("'")
                metaboanalyst_groups[i_i] = metaboanalyst_groups[i_i].strip("\"")
                if len(metaboanalyst_groups[i_i]) == 0:
                    del metaboanalyst_groups[i_i]
            plot_metaboanalyst = (plot_metaboanalyst_file, metaboanalyst_groups)
        iso_fittings = config['post-analysis/reanalysis'].getboolean('output_fittings_data')
        output_plot_data = config['post-analysis/reanalysis'].getboolean('output_plot_data')
        
    if from_GUI:
        return (custom_glycans_list, min_max_monos, min_max_hex, min_max_hexnac, min_max_sia, min_max_fuc, min_max_ac, min_max_gc, force_nglycan, max_adducts, adducts_exclusion, max_charges, reducing_end_tag, permethylated, reduced, lactonized_ethyl_esterified, fast_iso, high_res, internal_standard, imp_exp_library, exp_lib_name, library_path, only_gen_lib), (multithreaded_analysis, number_cores, analyze_ms2, reporter_ions, tolerance, ret_time_interval, rt_tolerance_frag, min_isotopologue_peaks, min_ppp, close_peaks, align_chromatograms, percentage_auc, max_ppm, iso_fit_score, curve_fit_score, s_to_n, custom_noise, samples_path, save_path, plot_metaboanalyst, compositions, iso_fittings, reanalysis, reanalysis_path, output_plot_data)
    
    
    shutil.copy(library_path, os.path.join(temp_folder, 'glycans_library.py'))
    spec = importlib.util.spec_from_file_location("glycans_library", temp_folder+"/glycans_library.py")
    lib_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lib_module)
    try:
        library_metadata = lib_module.metadata
    except:
        library_metadata = []
    if len(library_metadata) > 0:
        min_max_monos = library_metadata[0]
        min_max_hex = library_metadata[1]
        min_max_hexnac = library_metadata[2]
        min_max_fuc = library_metadata[3]
        min_max_sia = library_metadata[4]
        min_max_ac = library_metadata[5]
        min_max_gc = library_metadata[6]
        force_nglycan = library_metadata[7]
        max_adducts = library_metadata[8]
        max_charges = library_metadata[9]
        tag_mass = library_metadata[10]
        internal_standard = library_metadata[11]
        permethylated = library_metadata[12]
        lactonized_ethyl_esterified = library_metadata[13]
        reduced = library_metadata[14]
        fast_iso = library_metadata[15]
        high_res = library_metadata[16]
                
    #args to execution functions:
    output_filtered_data_args = [curve_fit_score, iso_fit_score, s_to_n, max_ppm, percentage_auc, reanalysis, reanalysis_path, save_path, analyze_ms2[0], analyze_ms2[2], reporter_ions, plot_metaboanalyst, compositions, align_chromatograms, force_nglycan, ret_time_interval[2], rt_tolerance_frag, iso_fittings, output_plot_data, multithreaded_analysis, number_cores, 0.0]

    imp_exp_gen_library_args = [custom_glycans_list, min_max_monos, min_max_hex, min_max_hexnac, min_max_sia, min_max_fuc, min_max_ac, min_max_gc, force_nglycan, max_adducts, adducts_exclusion, max_charges, reducing_end_tag, fast_iso, high_res, imp_exp_library, library_path, exp_lib_name, only_gen_lib, save_path, internal_standard, permethylated, lactonized_ethyl_esterified, reduced]

    list_of_data_args = [samples_list]

    index_spectra_from_file_ms1_args = [None, 1, multithreaded_analysis, number_cores]

    index_spectra_from_file_ms2_args = [None, 2, multithreaded_analysis, number_cores]

    analyze_files_args = [None, None, None, None, tolerance, ret_time_interval, min_isotopologue_peaks, min_ppp, max_charges, custom_noise, close_peaks, multithreaded_analysis, number_cores]

    analyze_ms2_args = [None, None, None, ret_time_interval, tolerance, min_max_monos, min_max_hex, min_max_hexnac,  min_max_sia, min_max_fuc, min_max_ac, min_max_gc, max_charges, reducing_end_tag, force_nglycan, permethylated, reduced, lactonized_ethyl_esterified, analyze_ms2[1], analyze_ms2[2], ret_time_interval[2], multithreaded_analysis, number_cores]

    arrange_raw_data_args = [None, samples_names, analyze_ms2[0], save_path, [(custom_glycans_list, min_max_monos, min_max_hex, min_max_hexnac, min_max_sia, min_max_fuc, min_max_ac, min_max_gc, force_nglycan, max_adducts, adducts_exclusion, max_charges, reducing_end_tag, permethylated, reduced, lactonized_ethyl_esterified, fast_iso, high_res, internal_standard, imp_exp_library, exp_lib_name, library_path, only_gen_lib), (multithreaded_analysis, number_cores, analyze_ms2, reporter_ions, tolerance, ret_time_interval, rt_tolerance_frag, min_isotopologue_peaks, min_ppp, close_peaks, align_chromatograms, percentage_auc, max_ppm, iso_fit_score, curve_fit_score, s_to_n, custom_noise, samples_path, save_path, plot_metaboanalyst, compositions, iso_fittings, reanalysis, reanalysis_path, output_plot_data)]]

    return output_filtered_data_args, imp_exp_gen_library_args, list_of_data_args, index_spectra_from_file_ms1_args, index_spectra_from_file_ms2_args, analyze_files_args, analyze_ms2_args, arrange_raw_data_args, samples_names, reanalysis, analyze_ms2[0]
