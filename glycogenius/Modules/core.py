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
from pathlib import Path
import sys
import os
import datetime
import configparser

#-----------------------------------------------------------------------------

def main():
    custom_glycans_list = (False, [])
    min_max_monos = (0, 0)
    min_max_hex = (0, 0)
    min_max_hexnac = (0, 0)
    min_max_sia = (0, 0)
    min_max_fuc = (0, 0)
    min_max_ac = (0, 0)
    min_max_gc = (0, 0)
    force_nglycan = False
    max_adducts = {}
    max_charges = 0
    tag_mass = 0.0
    permethylated = False
    reduced = False
    fast_iso = True
    high_res = False
    internal_standard = 0.0
    imp_exp_library = (False, False)
    only_gen_lib = False

    multithreaded_analysis = (False, 1)
    analyze_ms2 = (False, False, False)
    accuracy_unit = 'pw'
    accuracy_value = 0.0
    ret_time_interval = (0, 99999, 0.2)
    min_isotopologue_peaks = 0
    min_ppp = (False, 0)
    close_peaks = (False, 3)
    max_ppm = 0
    iso_fit_score = 1.0
    curve_fit_score = 1.0
    s_to_n = 0
    custom_noise = (False, [])
    samples_path = ''
    save_path = ''
    plot_metaboanalyst = (False, [])
    reanalysis = (False, False)

    multithreaded_execution = (False, 0, 0) #editted by multithreaded 1
    sneakpeek = (False, 0)
    verbose = False

    if not os.isatty(0) or multithreaded_execution[0]: #If multithreaded execution or pipelined parameters file
        Execution_Functions.print_header(False)
        config = configparser.ConfigParser()
        configs = ""
        for line in sys.stdin:              #editted by multithreaded 2
            configs+=line                   #editted by multithreaded 3
        config.read_string(configs)
        custom_glycans = config['library_building']['custom_glycans_list'].split(",")
        for i_i, i in enumerate(custom_glycans):
            custom_glycans[i_i] = i.strip()
            if len(i) == 0:
                custom_glycans = custom_glycans[:i_i]+custom_glycans[i_i+1:]
        custom_glycans_list = (config['library_building'].getboolean('use_custom_glycans_list'), custom_glycans)
        min_max_monos = (int(config['library_building']['min_monos']), int(config['library_building']['max_monos']))
        min_max_hex = (int(config['library_building']['min_hex']), int(config['library_building']['max_hex']))
        min_max_hexnac = (int(config['library_building']['min_hexnac']), int(config['library_building']['max_hexnac']))
        min_max_sia = (int(config['library_building']['min_sia']), int(config['library_building']['max_sia']))
        min_max_fuc = (int(config['library_building']['min_fuc']), int(config['library_building']['max_fuc']))
        min_max_ac = (int(config['library_building']['min_ac']), int(config['library_building']['max_ac']))
        min_max_gc = (int(config['library_building']['min_gc']), int(config['library_building']['max_gc']))
        force_nglycan = config['library_building'].getboolean('force_nglycan')
        max_adducts = General_Functions.form_to_comp(config['library_building']['max_adducts'])
        max_charges = int(config['library_building']['max_charges'])
        tag_mass = float(config['library_building']['tag_mass'])
        permethylated = config['library_building'].getboolean('permethylated')
        reduced = config['library_building'].getboolean('reduced')
        fast_iso = config['library_building'].getboolean('fast_iso')
        high_res = config['library_building'].getboolean('high_resolution_isotopic_dist')
        internal_standard = float(config['library_building']['internal_standard_mass'])
        imp_exp_library = (config['library_building'].getboolean('imp_library'), config['library_building'].getboolean('exp_library'))
        only_gen_lib = config['library_building'].getboolean('only_gen_lib')
        multithreaded_analysis = (config['analysis_parameters'].getboolean('multithreaded_analysis'), int(config['analysis_parameters']['threads_number']))
        analyze_ms2 = (config['analysis_parameters'].getboolean('analyze_ms2'), config['analysis_parameters'].getboolean('force_fragments_to_glycans'), config['analysis_parameters'].getboolean('unrestricted_fragments'))
        accuracy_unit = config['analysis_parameters']['accuracy_unit']
        accuracy_value = float(config['analysis_parameters']['accuracy_value'])
        ret_time_interval = (float(config['analysis_parameters']['ret_time_begin']), float(config['analysis_parameters']['ret_time_end']), float(config['analysis_parameters']['ret_time_tolerance']))
        min_isotopologue_peaks = int(config['analysis_parameters']['min_isotopologue_peaks'])
        min_ppp = (config['analysis_parameters'].getboolean('custom_min_points_per_peak'), int(config['analysis_parameters']['number_points_per_peak']))
        close_peaks = (config['analysis_parameters'].getboolean('limit_peaks_picked'), int(config['analysis_parameters']['max_number_peaks']))
        max_ppm = int(config['analysis_parameters']['max_ppm'])
        iso_fit_score = float(config['analysis_parameters']['isotopic_fitting_score'])
        curve_fit_score = float(config['analysis_parameters']['curve_fitting_score'])
        s_to_n = int(config['analysis_parameters']['signal_to_noise'])
        noise_levels = config['analysis_parameters']['noise_levels'].split(",")
        for i_i, i in enumerate(noise_levels):
            noise_levels[i_i] = int(i.strip())
        custom_noise = (config['analysis_parameters'].getboolean('custom_noise_level'), noise_levels)
        samples_path = config['analysis_parameters']['samples_path']
        samples_path = samples_path.strip()
        samples_path = samples_path.strip("'")
        samples_path = samples_path.strip("\"")
        for i_i, i in enumerate(samples_path):
            if i == "\\":
                samples_path = samples_path[:i_i]+"/"+samples_path[i_i+1:]
        save_path = config['analysis_parameters']['working_path']
        save_path = save_path.strip()
        save_path = save_path.strip("'")
        save_path = save_path.strip("\"")
        for i_i, i in enumerate(save_path):
            if i == "\\":
                save_path = save_path[:i_i]+"/"+save_path[i_i+1:]
        metaboanalyst_groups = config['analysis_parameters']['metaboanalyst_groups'].split(",")
        for i_i in range(len(metaboanalyst_groups)-1, -1, -1):
            metaboanalyst_groups[i_i] = metaboanalyst_groups[i_i].strip()
            metaboanalyst_groups[i_i] = metaboanalyst_groups[i_i].strip("'")
            metaboanalyst_groups[i_i] = metaboanalyst_groups[i_i].strip("\"")
            if len(metaboanalyst_groups[i_i]) == 0:
                del metaboanalyst_groups[i_i]
        plot_metaboanalyst = (config['analysis_parameters'].getboolean('plot_metaboanalyst'), metaboanalyst_groups)
        reanalysis = (config['analysis_parameters'].getboolean('reanalysis'), config['analysis_parameters'].getboolean('output_plot_data'))
    else: #If no parameters file pipelines, run CLI
        parameters = Execution_Functions.interactive_terminal()
        Execution_Functions.print_sep()
        if parameters[0][0] == 69: #sneakpeek feature
            sneakpeek = (True, parameters[2])
            config = configparser.ConfigParser()
            configs = ""
            with open(parameters[1]+'glycogenius_parameters.ini', 'r') as f:
                for line in f:
                    configs+=line
            config.read_string(configs)
            max_ppm = int(config['analysis_parameters']['max_ppm'])
            iso_fit_score = float(config['analysis_parameters']['isotopic_fitting_score'])
            curve_fit_score = float(config['analysis_parameters']['curve_fitting_score'])
            s_to_n = int(config['analysis_parameters']['signal_to_noise'])
            analyze_ms2 = (config['analysis_parameters'].getboolean('analyze_ms2'), config['analysis_parameters'].getboolean('force_fragments_to_glycans'), config['analysis_parameters'].getboolean('unrestricted_fragments'))
            save_path = parameters[1]
            reanalysis = (True, True)
        if parameters[0][0] == 1 or parameters[0][0] == 2:
            if parameters[0][1] == 1:
                custom_glycans_list = (True, parameters[1])
                max_adducts = parameters[2]
                max_charges = parameters[3]
                tag_mass = parameters[4]
                fast_iso = parameters[5]
                high_res = parameters[6]
                permethylated = parameters[8]
                reduced = parameters[9]
            if parameters[0][1] == 2:
                min_max_monos = (parameters[1][0], parameters[1][1])
                min_max_hex = (parameters[1][2], parameters[1][3])
                min_max_hexnac = (parameters[1][4], parameters[1][5])
                min_max_fuc = (parameters[1][6], parameters[1][7])
                min_max_ac = (parameters[1][8], parameters[1][9])
                min_max_gc = (parameters[1][10], parameters[1][11])
                min_max_sia = (parameters[1][12], parameters[1][13])
                force_nglycan = parameters[1][14]
                max_adducts = parameters[2]
                max_charges = parameters[3]
                tag_mass = parameters[4]
                fast_iso = parameters[5]
                high_res = parameters[6]
                permethylated = parameters[8]
                reduced = parameters[9]
            if parameters[0][0] == 1:
                save_path = parameters[7]
                if save_path[-1] != "/":
                    save_path+= "/"
                only_gen_lib = True
            if parameters[0][0] == 2:
                analyze_ms2 = parameters[7]
                accuracy_unit = parameters[8]
                accuracy_value = parameters[9]
                ret_time_interval = (parameters[10][0], parameters[10][1], 0.2)
                min_isotopologue_peaks = parameters[11]
                max_ppm = parameters[12]
                iso_fit_score = parameters[13]
                curve_fit_score = parameters[14]
                s_to_n = parameters[15]
                samples_path = parameters[16]
                save_path = parameters[17]
                permethylated = parameters[18]
                reduced = parameters[19]
        if parameters[0][0] == 3:
            save_path = parameters[1]
            max_ppm = parameters[2]
            iso_fit_score = parameters[3]
            curve_fit_score = parameters[4]
            s_to_n = parameters[5]
            reanalysis = (True, True)
        if parameters[0][0] == 4:
            comments = parameters[1]
            save_path = parameters[2]
            Path(save_path).mkdir(exist_ok = True, parents = True)
            Execution_Functions.generate_cfg_file(save_path, comments)

    Path(save_path).mkdir(exist_ok = True, parents = True)
    samples_list = Execution_Functions.samples_path_to_list(samples_path)

#-----------------------------------------------------------------------------

    begin_time = datetime.datetime.now()

    if reanalysis[0]:
        Execution_Functions.output_filtered_data(curve_fit_score,
                                                 iso_fit_score,
                                                 s_to_n,
                                                 max_ppm,
                                                 reanalysis,
                                                 save_path,
                                                 multithreaded_analysis,
                                                 multithreaded_execution,
                                                 analyze_ms2[0],
                                                 analyze_ms2[2],
                                                 plot_metaboanalyst,
                                                 ret_time_interval[2],
                                                 sneakpeek)

    else:
        if multithreaded_execution[0]:
            print('Multithreaded Execution: '+str(multithreaded_execution[1]))
        samples_names = Execution_Functions.sample_names(samples_list)
        if not only_gen_lib:
            print("Sample files detected: "+str(len(samples_names)))
            for i in samples_names:
                print("--> "+i)
            Execution_Functions.print_sep()
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
                                                          only_gen_lib,
                                                          save_path,
                                                          internal_standard,
                                                          permethylated,
                                                          reduced)
#line to add multithreaded library
        if multithreaded_execution[0]:
            library = full_library
        print('Library length: '+str(len(library)))
        Execution_Functions.print_sep()
        tolerance = (accuracy_unit, accuracy_value)
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
                                                          close_peaks,
                                                          fast_iso,
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
                                                            force_nglycan,
                                                            permethylated,
                                                            reduced,
                                                            analyze_ms2[1],
                                                            analyze_ms2[2],
                                                            ret_time_interval[2])
        Execution_Functions.print_sep()
        Execution_Functions.arrange_raw_data(analyzed_data,
                                             samples_names,
                                             multithreaded_analysis,
                                             multithreaded_execution,
                                             analyze_ms2[0],
                                             save_path)
        Execution_Functions.print_sep()
        Execution_Functions.output_filtered_data(curve_fit_score,
                                                 iso_fit_score,
                                                 s_to_n,
                                                 max_ppm,
                                                 reanalysis,
                                                 save_path,
                                                 multithreaded_analysis,
                                                 multithreaded_execution,
                                                 analyze_ms2[0],
                                                 analyze_ms2[2],
                                                 plot_metaboanalyst,
                                                 ret_time_interval[2],
                                                 sneakpeek)
                                                 
    if not os.isatty(0):
        print('Execution complete. Time elapsed: '+str(datetime.datetime.now() - begin_time))
    else:
        input('Execution complete. Time elapsed: '+str(datetime.datetime.now() - begin_time)+'\nPress Enter to exit.')
#here multithreaded prints execution of main()