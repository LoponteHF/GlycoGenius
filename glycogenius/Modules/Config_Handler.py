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
import dill

forced_structures = ['none', 'n_glycans', 'o_glycans', 'gags']

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
    
    custom_glycans_list = [False, []]
    min_max_monos = [0, 0]
    min_max_hex = [0, 0]
    min_max_hn = [0, 0]
    min_max_hexnac = [0, 0]
    min_max_xyl = [0, 0]
    min_max_sia = [0, 0]
    min_max_fuc = [0, 0]
    min_max_ac = [0, 0]
    min_max_gc = [0, 0]
    min_max_ua = [0, 0]
    custom_monosaccharides = []
    forced = 'none'
    min_max_proton = [0, 0]
    custom_adducts = []
    max_charges = 0
    reducing_end_tag = 0.0
    permethylated = False
    reduced = False
    lactonized_ethyl_esterified = False
    min_max_sulfation = [0, 0]
    min_max_phosphorylation = [0, 0]
    lyase_digested = False
    fast_iso = True
    high_res = False
    internal_standard = '0.0'
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
    min_samples = 0
    max_ppm = 10
    iso_fit_score = 0.9
    curve_fit_score = 0.9
    s_to_n = 3
    fill_gaps = (False, 50, 0.2, False)
    custom_noise = [False, []]
    samples_path = ''
    save_path = ''
    plot_metaboanalyst = [False, ""]
    compositions = True
    iso_fittings = False
    reanalysis = False
    reanalysis_path = ''
    output_plot_data = False

    samples_list = []
    samples_names = []
    
    if not from_GUI:
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
    save_path = save_path.strip('"')
    save_path = Path(save_path)
    save_path.mkdir(exist_ok = True, parents = True)
    
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
        reanalysis_path = reanalysis_path.strip('"')
        reanalysis_path = Path(reanalysis_path)
        
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
        library_path = library_path.strip('"')
        library_path = Path(library_path)
        
        with open(library_path, 'rb') as f:
            library_data = dill.load(f)
            f.close()
        full_library = library_data[0]
        library_metadata = library_data[1]
        
        custom_glycans_list[0] = library_metadata.get('custom glycans list', [False, []])[0]
        if library_metadata.get('custom glycans list', [False, []])[0]:
            custom_glycans_list[1] = library_metadata.get('custom glycans list', [False, []])[1]
        min_max_monos = library_metadata.get('min/max monosaccharides', [0, 0])
        min_max_hex = library_metadata.get('min/max hexoses', [0, 0])
        min_max_hexnac = library_metadata.get('min/max hexnac', [0, 0])
        min_max_fuc = library_metadata.get('min/max fucoses', [0, 0])
        min_max_sia = library_metadata.get('min/max sialic acids', [0, 0])
        min_max_ac = library_metadata.get('min/max acetyl sialic acids', [0, 0])
        min_max_gc = library_metadata.get('min/max glycolyl sialic acids', [0, 0])
        min_max_xyl = library_metadata.get('min/max xyloses', [0, 0])
        min_max_hn = library_metadata.get('min/max hexosamines', [0, 0])
        min_max_ua = library_metadata.get('min/max uronic acids', [0, 0])
        min_max_sulfation = library_metadata[21]
        min_max_phosphorylation = library_metadata[22]
        forced = library_metadata.get('glycan class', None)
        min_max_proton = library_metadata.get('min/max protons', [0, 0])
        max_charges = library_metadata.get('maximum charges', 3)
        adduct_combos = library_metadata.get('adduct combos', [])
        reducing_end_tag = library_metadata.get('reducing end tag', 0.0)
        internal_standard = library_metadata.get('internal standard', '0.0')
        permethylated = library_metadata.get('permethylated', False)
        lactonized_ethyl_esterified = library_metadata.get('lactonized/ethyl esterified', False)
        reduced = library_metadata.get('reducing end reduced', False)
        fast_iso = library_metadata.get('fast isotopic pattern calculation', True)
        high_res = library_metadata.get('high resolution isotopic pattern', False)
        lyase_digested = library_metadata.get('lyase digested', False)
        custom_monosaccharides = library_metadata.get('custom monosaccharides', [])
                
    elif library_mode == 'custom_library':
        custom_glycans_list[0] = True
        
        if not reanalysis:
            if (len(config['library_building_modes']['custom_glycans_list'].split("\\")) > 1 or len(config['library_building_modes']['custom_glycans_list'].split("/")) > 1):
                try:
                    # Adjust custom_glycans path
                    temp_path_custom_glycans = config['library_building_modes']['custom_glycans_list']
                    temp_path_custom_glycans = temp_path_custom_glycans.strip()
                    temp_path_custom_glycans = temp_path_custom_glycans.strip('"')
                    temp_path_custom_glycans = temp_path_custom_glycans.strip("'")
                    temp_path_custom_glycans = Path(temp_path_custom_glycans)
                    
                    # Extract the list of glycans
                    temp_custom_glycans_list = ""
                    with open(temp_path_custom_glycans, 'r') as f:
                        for i in f:
                            temp_custom_glycans_list += i
                        f.close()
                    
                    # Separate the glycans in a list, doesn't matter if they are comma separated or line separated
                    custom_glycans = temp_custom_glycans_list.split(",")
                    if len(custom_glycans) == 1:
                        custom_glycans = temp_custom_glycans_list.split("\n")
                        
                    # Sort out the glycans in the list
                    to_remove = []
                    to_add = []
                    if len(custom_glycans) > 1:
                        for i_i, i in enumerate(custom_glycans):
                            custom_glycans[i_i] = i.strip()
                            if len(custom_glycans[i_i]) == 0:
                                to_remove.append(i_i)
                                continue
                            if len(custom_glycans[i_i].split("/")) > 1:
                                splitted_glycan = custom_glycans[i_i].split("/")
                                custom_glycans[i_i] = splitted_glycan[0]
                                to_add.append(splitted_glycan[1])
                    for i in sorted(to_remove, reverse = True):
                        del custom_glycans[i]
                    for i in to_add:
                        custom_glycans.append(i)
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
                to_add = []
                for i_i, i in enumerate(custom_glycans):
                    custom_glycans[i_i] = i.strip()
                    if len(custom_glycans[i_i]) == 0:
                        custom_glycans = custom_glycans[:i_i]+custom_glycans[i_i+1:]
                        continue
                    if len(custom_glycans[i_i].split("/")) > 1:
                        splitted_glycan = custom_glycans[i_i].split("/")
                        custom_glycans[i_i] = splitted_glycan[0]
                        to_add.append(splitted_glycan[1])
                for i in to_add:
                    custom_glycans.append(i)
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
        
        total_xyl = config['library_building_modes']['xyloses'].split(',')
        min_max_xyl = [int(total_xyl[0]), int(total_xyl[1])]
        
        total_sialics = config['library_building_modes']['sialic_acids'].split(',')
        min_max_sia = [int(total_sialics[0]), int(total_sialics[1])]
        
        total_fucoses = config['library_building_modes']['fucoses'].split(',')
        min_max_fuc = [int(total_fucoses[0]), int(total_fucoses[1])]
        
        total_neu5ac = config['library_building_modes']['neu5ac'].split(',')
        min_max_ac = [int(total_neu5ac[0]), int(total_neu5ac[1])]
        
        total_neu5gc = config['library_building_modes']['neu5gc'].split(',')
        min_max_gc = [int(total_neu5gc[0]), int(total_neu5gc[1])]
        
        total_hn = config['library_building_modes']['hexosamines'].split(',')
        min_max_hn = [int(total_hn[0]), int(total_hn[1])]
        
        total_ua = config['library_building_modes']['uronic_acids'].split(',')
        min_max_ua = [int(total_ua[0]), int(total_ua[1])]
        
    imp_exp_library[1] = config['library_building_modes'].getboolean('export_library')
    if imp_exp_library[1]:
        exp_lib_name = config['library_building_modes']['exported_library_name'].strip()
    if 'force_nglycan' in config['common_library_building_settings']:
        if config['common_library_building_settings'].getboolean('force_nglycan'):
            forced = 'n_glycans'
        else:
            forced = 'none'
    else:
        if config['common_library_building_settings']['force_class_structure'].lower() in forced_structures:
            forced = config['common_library_building_settings']['force_class_structure'].lower()
        else:
            print("Wrong 'force class structure' chosen.\nYou must choose between 'n_glycan', 'o_glycan'\nor 'gags'.")
            print("Close the window or press CTRL+C to exit.")
            try:
                while True:
                    time.sleep(3600)
            except KeyboardInterrupt:
                os._exit(1)
                
    custom_monosaccharides_list = config['library_building_modes'].get('custom_monosaccharides', [])
                    
    single_letters_in_use = set([x[-1] for x in General_Functions.monosaccharides.values()]+['T'])
    
    short_codes_in_use = set([x.upper() for x in General_Functions.monosaccharides.keys()]+['T'])
    
    if len(custom_monosaccharides_list) > 0:
        for cm in custom_monosaccharides_list.split("),"):
            cm = cm.strip().strip("(").strip(")")
            cm = [element.strip() for element in cm.split(",")]
            
            # If short-code in use
            if cm[1] in short_codes_in_use:
                print(f"Custom monosaccharide {cm[0]} short-code already in use. Use one not listed below:")
                for short_code in [x.upper() for x in General_Functions.monosaccharides.keys()]+['T']:
                    print(short_code)
                print("Close the window or press CTRL+C to exit.")
                try:
                    while True:
                        time.sleep(3600)
                except KeyboardInterrupt:
                    os._exit(1)
            
            # Generate new one-letter code
            alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            for letter in alphabet:
                if letter not in single_letters_in_use:
                    single_letter_code = letter
                    single_letters_in_use.add(single_letter_code)
                    break
                    
            custom_monosaccharides.append({
                'cm_name' : cm[0], 
                'cm_short_code' : cm[1],
                'cm_single_letter_code' : single_letter_code,
                'cm_chem_comp' : cm[2],
                'cm_min' : int(cm[3]),
                'cm_max' : int(cm[4]),
                'sialic' : True if cm[5] == 'yes' else False
                })
    
        
    total_protons = config['common_library_building_settings']['min_max_proton_adducts'].split(',')
    min_max_proton = [int(total_protons[0]), int(total_protons[1])]
    
    # Fill in details to use custom adducts
    custom_adducts_list = config['common_library_building_settings'].get('custom_adducts', [])
    if len(custom_adducts_list) > 0:
        for ca in custom_adducts_list.split("),"):
            ca = ca.strip().strip("(").strip(")")
            ca = [element.strip() for element in ca.split(",")]
            custom_adducts.append(ca)
    
    max_charges = int(config['common_library_building_settings']['max_charges'])
    reducing_end_tag = config['common_library_building_settings']['reducing_end_tag']
    permethylated = config['common_library_building_settings'].getboolean('permethylated')
    reduced = config['common_library_building_settings'].getboolean('reduced')
    lactonized_ethyl_esterified = config['common_library_building_settings'].getboolean('aminated_ethyl_esterified')
        
    total_sulfation = config['common_library_building_settings']['min_max_sulfation_per_glycan'].split(',')
    min_max_sulfation = [int(total_sulfation[0]), int(total_sulfation[1])]
        
    total_phosphorylation = config['common_library_building_settings']['min_max_phosphorylation_per_glycan'].split(',')
    min_max_phosphorylation = [int(total_phosphorylation[0]), int(total_phosphorylation[1])]
    
    if forced == 'gags':
        lyase_digested = config['common_library_building_settings'].getboolean('lyase_digested')
        
    fast_iso = config['common_library_building_settings'].getboolean('fast_iso')
    high_res = config['common_library_building_settings'].getboolean('high_resolution_isotopic_dist')
    internal_standard = config['common_library_building_settings']['internal_standard_mass']
    if len(internal_standard.strip()) == 0:
        internal_standard = '0.0'
    else:
        internal_standard = internal_standard
        
        
    #analysis running mode
    if not reanalysis and not only_gen_lib:
        samples_path = config['running_modes']['samples_directory']
        samples_path = samples_path.strip()
        samples_path = samples_path.strip("'")
        samples_path = samples_path.strip('"')
        samples_path = Path(samples_path)
        
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
        if not from_GUI:
            print("Sample files detected: "+str(len(samples_names)))
            for i in samples_names:
                print("--> "+i)
            Execution_Functions.print_sep()
            
        ms2_analysis = config['analysis_parameters'].getboolean('analyze_ms2')
        if ms2_analysis:
            analyze_ms2 = [ms2_analysis, config['analysis_parameters'].getboolean('force_fragments_to_glycans'), config['analysis_parameters'].getboolean('unrestricted_fragments')]
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
            min_ppp = [config['analysis_parameters'].getboolean('custom_min_points_per_peak'), int(config['analysis_parameters']['number_points_per_peak'])]
        limit_peaks = config['analysis_parameters'].getboolean('limit_peaks_picked')
        if limit_peaks:
            close_peaks = [config['analysis_parameters'].getboolean('limit_peaks_picked'), int(config['analysis_parameters']['max_number_peaks'])]
        if 'noise_levels' in config['analysis_parameters']:
            noise_levels = config['analysis_parameters']['noise_levels'].split(",")
            for i_i, i in enumerate(noise_levels):
                noise_levels[i_i] = int(i.strip())
            custom_noise = [config['analysis_parameters'].getboolean('custom_noise_level'), noise_levels]
        
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
        min_samples = int(config['post-analysis/reanalysis']['minimum_samples'])
        ppm_setting = config['post-analysis/reanalysis']['max_ppm_threshold'].split(",")
        if len(ppm_setting) > 1:
            max_ppm = [float(ppm_setting[0]), float(ppm_setting[1])]
        else:
            max_ppm = float(config['post-analysis/reanalysis']['max_ppm_threshold'])
        iso_fit_score = float(config['post-analysis/reanalysis']['isotopic_fitting_score_threshold'])
        curve_fit_score = float(config['post-analysis/reanalysis']['curve_fitting_score_threshold'])
        s_to_n = float(config['post-analysis/reanalysis']['signal_to_noise_threshold'])
        
        fill_gaps_on_off = config['post-analysis/reanalysis'].get('fill_data_gaps', 'no')
        if fill_gaps_on_off.lower() == 'no' or fill_gaps_on_off.lower() == 'false':
            fill_gaps_on_off = False
        else:
            fill_gaps_on_off = True
        fill_gaps_percentage = float(config['post-analysis/reanalysis'].get('fill_data_gaps_min_samples', 50))
        fill_gaps_rt_tol = float(config['post-analysis/reanalysis'].get('fill_data_gaps_rt_tolerance', 0.2))
        fill_gaps_noise = config['post-analysis/reanalysis'].get('fill_remaining_gaps_with_noise', 'no')
        if fill_gaps_noise.lower() == 'no' or fill_gaps_noise.lower() == 'false':
            fill_gaps_noise = False
        else:
            fill_gaps_noise = True
            
        fill_gaps = (fill_gaps_on_off, fill_gaps_percentage, fill_gaps_rt_tol, fill_gaps_noise)
        
        compositions = config['post-analysis/reanalysis'].getboolean('analyze_compositions')
        plot_metaboanalyst_file = config['post-analysis/reanalysis'].getboolean('output_metaboanalyst_file')
        if plot_metaboanalyst_file:
            metaboanalyst_groups_path = config['post-analysis/reanalysis']['sample_groups'].strip().strip("'").strip('"')
            metaboanalyst_groups_path = Path(metaboanalyst_groups_path)
            plot_metaboanalyst = (plot_metaboanalyst_file, metaboanalyst_groups_path)
        iso_fittings = config['post-analysis/reanalysis'].getboolean('output_fittings_data')
        output_plot_data = config['post-analysis/reanalysis'].getboolean('output_plot_data')
        
        
    
    analysis_parameters = {
                            'library settings': {
                                'min/max monosaccharides': min_max_monos,
                                'min/max hexoses': min_max_hex,
                                'min/max hexnac': min_max_hexnac,
                                'min/max fucoses': min_max_fuc,
                                'min/max sialic acids': min_max_sia,
                                'min/max acetyl sialic acids': min_max_ac,
                                'min/max glycolyl sialic acids': min_max_gc,
                                'min/max xyloses': min_max_xyl,
                                'min/max hexosamines': min_max_hn,
                                'min/max uronic acids': min_max_ua,
                                'glycan class': forced,
                                'min/max protons': min_max_proton,
                                'custom adducts': custom_adducts,
                                'maximum charges': max_charges,
                                'import/export library': imp_exp_library,
                                'library path': library_path,
                                'exported library name': exp_lib_name,
                                'only generate library': only_gen_lib,
                                'reducing end tag': reducing_end_tag,
                                'internal standard': internal_standard,
                                'permethylated': permethylated,
                                'lactonized/ethyl esterified': lactonized_ethyl_esterified,
                                'reducing end reduced': reduced,
                                'fast isotopic pattern calculation': fast_iso,
                                'high resolution isotopic pattern': high_res,
                                'custom glycans list': custom_glycans_list,
                                'min/max sulfation': min_max_sulfation,
                                'min/max phosphorylation': min_max_phosphorylation,
                                'lyase digested': lyase_digested,
                                'custom monosaccharides': custom_monosaccharides
                                },
                            'analysis settings': {
                                'multithreaded': multithreaded_analysis,
                                'number of cpu cores': number_cores,
                                'analyze ms2': analyze_ms2,
                                'reporter ions': reporter_ions,
                                'mass tolerance': tolerance,
                                'retention time interval': ret_time_interval,
                                'fragmentation retention time tolerance': rt_tolerance_frag,
                                'min isotopologue peaks': min_isotopologue_peaks,
                                'min points per peak': min_ppp,
                                'only x most intense peaks': close_peaks,
                                'align chromatograms': align_chromatograms,
                                'minimum percentage auc threshold': percentage_auc,
                                'maximum ppm threshold': max_ppm,
                                'isotopic fitting score threshold': iso_fit_score,
                                'curve fitting score threshold': curve_fit_score,
                                'signal to noise ratio threshold': s_to_n,
                                'custom noise level': custom_noise,
                                'samples path': samples_path,
                                'save path': save_path,
                                'output abundance table': plot_metaboanalyst,
                                'output composition-separated data': compositions,
                                'output isotopic fitting data': iso_fittings,
                                'output chromatogram plotting data': output_plot_data,
                                'reanalysis': reanalysis,
                                'reanalysis_path': reanalysis_path,
                                'fill data gaps': fill_gaps
                                }
                            }
        
    if from_GUI:
        return analysis_parameters
                
    #args to execution functions:
    output_filtered_data_args = {
                                'curve fitting score threshold': curve_fit_score,
                                'isotopic fitting score threshold': iso_fit_score,
                                'signal to noise ratio threshold': s_to_n,
                                'maximum ppm threshold': max_ppm,
                                'minimum percentage auc threshold': percentage_auc,
                                'reanalysis': reanalysis,
                                'reanalysis_path': reanalysis_path,
                                'save path': save_path,
                                'analyze ms2': analyze_ms2[0],
                                'unrestricted fragments': analyze_ms2[2],
                                'reporter ions': reporter_ions,
                                'output abundance table': plot_metaboanalyst,
                                'output composition-separated data': compositions,
                                'align chromatograms': align_chromatograms,
                                'glycan class': forced,
                                'retention time tolerance': ret_time_interval[2],
                                'fragmentation retention time tolerance': rt_tolerance_frag,
                                'output isotopic fitting data': iso_fittings,
                                'output chromatogram plotting data': output_plot_data,
                                'multithreaded': multithreaded_analysis,
                                'number of cpu cores': number_cores,
                                'analysis done time': None,
                                'minimum percentage of samples': min_samples,
                                'sample groups': None,
                                'temporary folder': None,
                                'fill data gaps': fill_gaps,
                                'from GUI': False
                                }

    imp_exp_gen_library_args = {
                                'custom glycans list': custom_glycans_list,
                                'min/max monosaccharides': min_max_monos,
                                'min/max hexoses': min_max_hex,
                                'min/max hexnac': min_max_hexnac,
                                'min/max xyloses': min_max_xyl,
                                'min/max sialic acids': min_max_sia,
                                'min/max fucoses': min_max_fuc,
                                'min/max acetyl sialic acids': min_max_ac,
                                'min/max glycolyl sialic acids': min_max_gc,
                                'min/max hexosamines': min_max_hn,
                                'min/max uronic acids': min_max_ua,
                                'glycan class': forced,
                                'min/max protons': min_max_proton,
                                'custom adducts': custom_adducts,
                                'maximum charges': max_charges,
                                'reducing end tag': reducing_end_tag,
                                'fast isotopic pattern calculation': fast_iso,
                                'high resolution isotopic pattern': high_res,
                                'import/export library': imp_exp_library,
                                'library path': library_path,
                                'exported library name': exp_lib_name,
                                'only generate library': only_gen_lib,
                                'save path': save_path,
                                'internal standard': internal_standard,
                                'permethylated': permethylated,
                                'lactonized/ethyl-esterified': lactonized_ethyl_esterified,
                                'reducing end reduced': reduced,
                                'min/max sulfation': min_max_sulfation,
                                'min/max phosphorylation': min_max_phosphorylation,
                                'lyase digested': lyase_digested,
                                'temporary folder': None,
                                'custom monosaccharides': custom_monosaccharides,
                                'from GUI': False
                                }

    list_of_data_args = {
                        'samples list': samples_list
                        }

    index_spectra_from_file_ms1_args = {
                                        'raw data': None,
                                        'ms level': 1,
                                        'multithreaded': multithreaded_analysis,
                                        'number of cpu cores': number_cores
                                        }
                                        
    index_spectra_from_file_ms2_args = {
                                        'raw data': None,
                                        'ms level': 2,
                                        'multithreaded': multithreaded_analysis,
                                        'number of cpu cores': number_cores
                                        }

    analyze_files_args = {
                        'library': None,
                        'library size': None,
                        'raw data': None,
                        'ms1 index': None,
                        'mass tolerance': tolerance,
                        'retention time interval': ret_time_interval,
                        'min isotopologue peaks': min_isotopologue_peaks,
                        'min points per peak': min_ppp,
                        'adduct combos': None,
                        'maximum charges': max_charges,
                        'custom noise level': custom_noise,
                        'only x most intense peaks': close_peaks,
                        'multithreaded': multithreaded_analysis,
                        'number of cpu cores': number_cores,
                        'analysis start time': None,
                        'temporary folder': None,
                        'from GUI': False
                        }

    analyze_ms2_args = {
                        'ms2 index': None,
                        'raw data': None,
                        'analyzed data': None,
                        'retention time interval': ret_time_interval,
                        'mass tolerance': tolerance,
                        'min/max monosaccharides': min_max_monos,
                        'min/max hexoses': min_max_hex,
                        'min/max hexnac': min_max_hexnac,
                        'min/max xyloses': min_max_xyl,
                        'min/max sialic acids': min_max_sia,
                        'min/max fucoses': min_max_fuc,
                        'min/max acetyl sialic acids': min_max_ac,
                        'min/max glycolyl sialic acids': min_max_gc,
                        'min/max hexosamines': min_max_hn,
                        'min/max uronic acids': min_max_ua,
                        'adduct combos': None,
                        'maximum charges': max_charges,
                        'reducing end tag': reducing_end_tag,
                        'glycan class': forced,
                        'permethylated': permethylated,
                        'reducing end reduced': reduced,
                        'lactonized/ethyl-esterified': lactonized_ethyl_esterified,
                        'assign fragments compatible with precursor': analyze_ms2[1],
                        'annotate all glycans in library': analyze_ms2[2],
                        'fragmentation retention time tolerance': ret_time_interval[2],
                        'multithreaded': multithreaded_analysis,
                        'number of cpu cores': number_cores,
                        'library': None,
                        'temporary folder': None,
                        'from GUI': False,
                        'custom monosaccharides': custom_monosaccharides
                        }

    arrange_raw_data_args = {
                            'analyzed data': None,
                            'sample names': samples_names,
                            'ms2 analyzed': analyze_ms2[0],
                            'save path': save_path,
                            'analysis parameters': analysis_parameters,
                            'adduct combos': None,
                            'library': None,
                            'temporary folder': None,
                            'results file name': None,
                            'from GUI': False,
                            'erase files': True,
                            'calculated noise levels': []
                            }

    return {
            'output_filtered_data_args': output_filtered_data_args, 
            'imp_exp_gen_library_args': imp_exp_gen_library_args, 
            'list_of_data_args': list_of_data_args, 
            'index_spectra_from_file_ms1_args': index_spectra_from_file_ms1_args,
            'index_spectra_from_file_ms2_args': index_spectra_from_file_ms2_args, 
            'analyze_files_args': analyze_files_args, 
            'analyze_ms2_args': analyze_ms2_args, 
            'arrange_raw_data_args':arrange_raw_data_args
            }