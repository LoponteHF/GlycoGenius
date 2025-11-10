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

from . import Library_Tools
from . import General_Functions
from . import File_Accessing
from pyteomics import mzxml, mzml, mass
from itertools import islice, product
from pandas import DataFrame, ExcelWriter
from numpy import percentile
from re import split
from math import inf, isnan
from statistics import mean, median
import multiprocessing
import numpy as np
import concurrent.futures
import zipfile
import time
import os
import dill
import sys
import datetime
import platform
import copy
import pathlib
import shutil

version = '1.2.15'

##---------------------------------------------------------------------------------------
##Functions to be used for execution and organizing results data
    
def generate_cfg_file(path, comments):
    '''Generates a .ini file and a .bat file at the indicated directory, with or without comments.
    
    Parameters
    ----------
    path : string
        The path to where the file will be saved.
        
    comments : boolean
        Whether the output file will contain comments or not.
        
    Uses
    ----
    nothing
    
    Returns
    -------
    nothing
        Creates 2 files: a .ini and a .bat file.
    '''
    print("Creating parameters file...")
    glycogenius_path = pathlib.Path(__file__).parent.parent.resolve()
    curr_os = platform.system()
    pathlib.Path(path).mkdir(exist_ok = True, parents = True)
    pathlib.Path(os.path.join(path, "Sample Files")).mkdir(exist_ok = True, parents = True)
    with open(os.path.join(path, 'glycogenius_parameters.ini'), 'w') as g:
        with open(os.path.join(glycogenius_path, 'Parameters_Template.py'), 'r') as f:
            for line in f:
                if "samples_directory =" in line:
                    samples_directory = os.path.join(path, "Sample Files")
                    g.write(f"samples_directory = {samples_directory}\n")
                    continue
                if "working_directory =" in line:
                    g.write(f"working_directory = {path}\n")
                    continue
                if not comments and line[0] == ';':
                    continue
                g.write(line)
        f.close()
    g.close()
    if curr_os == "Windows":
        with open(os.path.join(path, 'Run Glycogenius.bat'), 'w') as f:
            f.write("@echo off\n")
            f.write("cd %~dp0\n")
            f.write("type .\\glycogenius_parameters.ini | glycogenius")
        f.close()
        print("Done!")
        print("Set your parameters in the file\n'glycogenius_parameters.ini' and\nrun 'Run Glycogenius.bat' to run Glycogenius\nwith the set parameters.")
    else:
        print("Done!")
        print("Set your parameters in the file\n'glycogenius_parameters.ini' and\ncat-pipeline it to glycogenius.")
    input("Press Enter to exit.")
    os._exit(1)
                
def samples_path_to_list(path):
    '''Detects all files in samples path and makes a list of the path to each file for
    internal usage.
    
    Parameters
    ----------
    path : string
        Path to the folder containing mzML and mzXML files to analyze.
        
    Uses
    ----
    itertools.product : 
        Cartesian product of input iterables.
        
    Returns
    -------
    samples_list : list
        A list containing the path to each file to be analyzed.
    '''
    if len(str(path)) == 0:
        return []
        
    mzml_possibilities = list(map(''.join, product(*zip("mzml".upper(), "mzml".lower()))))
    mzxml_possibilities = list(map(''.join, product(*zip("mzxml".upper(), "mzxml".lower()))))
    file_extensions = mzml_possibilities+mzxml_possibilities
    try:
        dir_list = os.listdir(path)
    except:
        return []
    samples_list = []
    for i_i, i in enumerate(dir_list):
        if i.split('.')[-1] in file_extensions:
            samples_list.append(os.path.join(path, i))
    return samples_list
    
def list_of_data(args_dict):
    '''Detects if a file is mzXML or mzML and processes it into a generator using
    pyteomics.

    Parameters
    ----------
    samples_list : list
        A list containing the path to each file to be analyzed.

    Uses
    ----
    pyteomics.mzxml.MzXML() : generator
        Indexes the mzXML file into a generator, allowing you to parse the file for
        analysis.

    File_Accessing.make_mzxml : class
        A wrapper that takes the output of pyteomics mzML parser and converts it to
        the mzXML pyteomics parser standard to be used within the script. Allows for 
        full support of mzML.

    Returns
    -------
    data : list
        A list containing the generator of each file name at each index.
    '''
    # Grab the args
    samples_list = args_dict.get('samples list', [])
    
    data = []
    mzml_possibilities = list(map(''.join, product(*zip("mzml".upper(), "mzml".lower()))))
    mzxml_possibilities = list(map(''.join, product(*zip("mzxml".upper(), "mzxml".lower()))))
    for i in samples_list:
        if i.split('.')[-1] in mzxml_possibilities:
            mzxml_data = mzxml.MzXML(i)
            data.append(mzxml_data)
        else:
            mzml_data = File_Accessing.make_mzxml(i)
            data.append(mzml_data)
    return data

def index_spectra_from_file(args_dict):
    '''Scans the mz(X)ML file and indexes all the MS1 scans, so that you don't have to
    go through the MSn scans as well when doing something only with the MS1.

    Parameters
    ----------
    files : list
        List with each index containing a generator created by the pyteomics function
        pyteomics.mzxml.MzXML() or File_Accessing.make_mzxml.
        
    ms_level : int
        The MS level to index (ie. MS1 or MS2, inputted as 1 or 2).
        
    multithreaded : boolean
        Whether or not to multithread the indexing.
        
    number_cores : string or int
        Number of cores to be used.

    Returns
    -------
    indexes : dict
        Returns a dictionary with each key pointing to a file and each key containing a
        list of indexes of the MS1 spectra.
    '''
    # Grab the args
    files = args_dict.get('raw data', None)
    ms_level = args_dict.get('ms level', 1)
    multithreaded = args_dict.get('multithreaded', True)
    number_cores = args_dict.get('number of cpu cores', 99)
    
    indexes = {}
    
    results = []
    if multithreaded:
        if number_cores == 'all':
            cpu_count = (os.cpu_count())-2
            if cpu_count <= 0:
                cpu_count = 1
        else:
            number_cores = int(number_cores)
            if number_cores > (os.cpu_count())-2:
                cpu_count = (os.cpu_count())-2
                if cpu_count <= 0:
                    cpu_count = 1
            else:
                cpu_count = number_cores
    else:
        cpu_count = 1
    
    with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_count if cpu_count < 60 else 60) as executor:
        for i_i, i in enumerate(files):
            result = executor.submit(get_indexes,
                                     i,
                                     ms_level,
                                     i_i)
            results.append(result)
    
        for index, i in enumerate(results):
            result_data = i.result()
            indexes[result_data[1]] = result_data[0]
            results[index] = None
        
    return indexes
    
def get_indexes(file,
                ms_level,
                file_id):
    '''Core function of index_spectra_from_file. Indexes the desired MS level
    for one file.
    
    Parameters
    ----------
    files : list
        List with each index containing a generator created by the pyteomics function
        pyteomics.mzxml.MzXML() or File_Accessing.make_mzxml.
        
    ms_level : int
        The MS level to index (ie. MS1 or MS2, inputted as 1 or 2).
        
    file_id : int
        The file of which the indexes belong to, for further use by the multithreading
        algorithm.

    Returns
    -------
    temp_indexes : list
        Returns a list containing the indexes of the MS1 spectra for a given file.
        
    file_id : int
        The file of which the indexes belong to, for further use by the multithreading
        algorithm.
    '''
    temp_indexes = [
        j_j for j_j, j in enumerate(file)
        if ('msLevel' in j and j['msLevel'] == ms_level) or ('ms level' in j and j['ms level'] == ms_level)
    ]
    return temp_indexes, file_id

def sample_names(samples_list):
    '''Extracts the sample names from the file path.
    
    Parameters
    ----------
    samples_list : list
        A list of strings containing the path to each sample file.

    Returns
    -------
    curated_samples : list
        A list of strings, each string with a sample name.
    '''
    curated_samples = []
    for i in samples_list:
        if i[-1] == "\\" or i[-1] == "/":
            i = i[:-1]
        i = i.split("/")[-1]
        i = i.split("\\")[-1]
        curated_samples.append(".".join(i.split(".")[:-1]))
    return curated_samples
    
def output_extra_library_files(full_library,
                               lactonized_ethyl_esterified,
                               custom_monos,
                               imp_exp_library,
                               library_path,
                               exp_lib_name,
                               metadata,
                               save_path):
    '''
    '''
    adduct_combos_dict = {General_Functions.comp_to_formula(adduct[0]): {'charges': adduct[1], 'comp': adduct[0]} for adduct in metadata['adduct combos']}
    
    # Here the header of the human-readable library is made
    df = {'Glycan' : [], 'Hex' : [], 'HexN' : [], 'HexNAc' : [], 'Xylose' : [], 'dHex' : []}
    if lactonized_ethyl_esterified:
        df['a2,3-Neu5Ac'] = []
        df['a2,6-Neu5Ac'] = []
        df['a2,3-Neu5Gc'] = []
        df['a2,6-Neu5Gc'] = []
    else:
        df['Neu5Ac'] = []
        df['Neu5Gc'] = []
    df['UroA'] = []
    
    if len(custom_monos) > 0:
        for cm in custom_monos:
            cm_name = cm['cm_name']
            if cm_name in df.keys():
                cm_name += '-custom'
            df[cm_name] = []
            
    df['Isotopic Distribution'] = []
    df['Neutral Mass + Tag'] = []
    
    # Here each glycan is added to the dataframe for human-readable library
    for i_i, i in enumerate(full_library):
        df['Glycan'].append(i)
        df['Hex'].append(full_library[i]['Monos_Composition']['H'])
        df['HexN'].append(full_library[i]['Monos_Composition']['HN'])
        df['HexNAc'].append(full_library[i]['Monos_Composition']['N'])
        df['Xylose'].append(full_library[i]['Monos_Composition']['X'])
        df['dHex'].append(full_library[i]['Monos_Composition']['F'])
        
        if lactonized_ethyl_esterified:
            df['a2,3-Neu5Ac'].append(full_library[i]['Monos_Composition']['Am'])
            df['a2,6-Neu5Ac'].append(full_library[i]['Monos_Composition']['E'])
            df['a2,3-Neu5Gc'].append(full_library[i]['Monos_Composition']['AmG'])
            df['a2,6-Neu5Gc'].append(full_library[i]['Monos_Composition']['EG'])
        else:
            df['Neu5Ac'].append(full_library[i]['Monos_Composition']['S'])
            df['Neu5Gc'].append(full_library[i]['Monos_Composition']['G'])
            
        df['UroA'].append(full_library[i]['Monos_Composition']['UA'])
            
        if len(custom_monos) > 0:
            for cm in custom_monos:
                cm_name = cm['cm_name']
                try:
                    df[cm_name].append(full_library[i]['Monos_Composition'].get(cm['cm_short_code'], 0))
                except:
                    df[cm_name+'-custom'].append(full_library[i]['Monos_Composition'].get(cm['cm_short_code'], 0))
        
        temp_isotopic = []
        for j in full_library[i]['Isotopic_Distribution']:
            temp_isotopic.append(float("%.3f" % round(j, 3)))
        df['Isotopic Distribution'].append(str(temp_isotopic)[1:-1])
        df['Neutral Mass + Tag'].append(float("%.4f" % round(full_library[i]['Neutral_Mass+Tag'], 4)))
        for j in full_library[i]['Adducts_mz']:
            if i_i ==0:
                df[j] = [float("%.4f" % round(full_library[i]['Adducts_mz'][j], 4))]
            else:
                df[j].append(float("%.4f" % round(full_library[i]['Adducts_mz'][j], 4)))
    
    # Here it starts to write the library to xlsx
    df = DataFrame(df)
    if imp_exp_library[0]:
        file_name = library_path.split("\\")[-1].split("/")[-1].split(".")[-2]
    else:
        file_name = exp_lib_name
    if file_name+'.xlsx' not in os.listdir(save_path):
        with ExcelWriter(os.path.join(save_path, file_name+'.xlsx')) as writer:
            df.to_excel(writer, index = False)
            General_Functions.autofit_columns_excel(df, writer.sheets['Sheet1'])
            
    # Here is the skyline transitions list creation
    if file_name+'_skyline_transitions.csv' not in os.listdir(save_path):
        with open(os.path.join(save_path, file_name+'_skyline_transitions.csv'), 'w') as f:
            f.write('Precursor Name, Precursor Formula, Precursor Adduct, Precursor Charge\n')
            for i_i, i in enumerate(full_library):
                for j_j, j in enumerate(full_library[i]['Adducts_mz']):
                    
                    adduct_comp = copy.deepcopy(adduct_combos_dict[j]['comp'])
                    adduct_charge = copy.deepcopy(adduct_combos_dict[j]['charges'])
                    
                    if len(adduct_comp) > 1 or i == "Internal Standard": #can't seem to make skyline work with mixed adducts, so have this in place for now
                        continue
                    
                    adduct = str(adduct_comp[list(adduct_comp.keys())[0]])+str(list(adduct_comp.keys())[0]) #only first adduct
                    del adduct_comp[list(adduct_comp.keys())[0]]
                    formula = General_Functions.comp_to_formula(General_Functions.sum_atoms(full_library[i]['Atoms_Glycan+Tag'], adduct_comp))
                    list_form = [i, str(formula), '[M+'+j+']', str(adduct_charge)]
                    f.write(",".join(list_form)+'\n')
            f.close()

def imp_exp_gen_library(args_dict):
    '''Imports, generates and/or exports a glycans library.

    Parameters
    ----------                
    custom_glycans_list : tuple
        A tuple with two indexes: The first one is whether or not the library should be built
        off a custom glycans list and the second one is a list containing the glycans you wish
        to analyze.
        
    min_max_monos : tuple
        Minimum and maximum amount of monosaccharides for the hypotethical glycans in the
        library. ie. (5, 20).
        
    min_max_hex : tuple
        Minimum and maximum amount of hexoses for the hypotethical glycans in the library.
        ie. (5, 20).
        
    min_max_hexnac : tuple
        Minimum and maximum amount of N-Acetyl hexosamines for the hypotethical glycans
        in the library. ie. (5, 20).

    min_max_xyl : tuple
        Minimum and maximum amount of Xyloses for the hypotethical glycans
        in the library. ie. (5, 20).
        
    min_max_sia : tuple
        Minimum and maximum amount of sialic acids for the hypotethical glycans in the
        library. ie. (5, 20).
        
    min_max_fuc : tuple
        Minimum and maximum amount of deoxyhexoses for the hypotethical glycans in the
        library. ie. (5, 20).
        
    min_max_ac : tuple
        Minimum and maximum amount of N-Acetyl Neuraminic acids for the hypotethical
        glycans in the library. ie. (5, 20).
        
    min_max_gc : tuple
        Minimum and maximum amount of N-Glicolyl Neuraminic acids for the hypotethical
        glycans in the library. ie. (5, 20).
        
    min_max_hn : tuple
        Minimum and maximum amount of Hexosamines for the hypotethical glycans in the library. ie. (5, 20).
        
    min_max_ua : tuple
        Minimum and maximum amount of Uronic Acids for the hypotethical glycans in the library. ie. (5, 20).

    forced : string
        Indicates whether the function should force strict conditions based on the
        biological knowledge of glycans in order to avoid possible false positives when
        analysing N-glycans, O-glycans or GAGs.

    max_adducts : dict
        A dictionary with keys containing each possible atomic adducts (ie. 'H', 'Na',
        'K', etc.) and the maximum amount of such adducts as the values.
        
    max_charges : int
        The maximum amount of charges to calculate per glycan.

    tag_mass : float
        The tag's added mass to the glycans, if the glycans are tagged.
        Default = 0 (No Tag).

    fast_iso : boolean
        Makes the isotopic distribution calculation fast (less accurate) or not (more
        accurate).
        Default = True
        
    high_res : boolean
        Decides whether to clump (if set to False) or not (if set to True) the neighbouring
        isotopic envelope peaks. Only works if fast is set to False.
        Default = False

    imp_exp_library : tuple
        A list with two indexes : The first one indicates whether or not the script should try
        to import the library from a glycans_library.py file and the second one indicates
        whether or not the script should export the library to a glycans_library.py file and
        an excel file for visualization.
        
    library_path : str
        Path to save library export files.
        
    exp_lib_name : str
        Name of the library file saved when exporting.
        
    only_gen_lib : boolean
        A boolean indicating whether or not this execution should stop after generating library.
        
    save_path : string
        A string containing the path to the working directory of the script.
        
    internal_standard : float
        If a internal standard is added to the sample, this allows the function
        to calculate its mass based on adducts combination, as well as adding the tag to
        it.
        
    permethylated : boolean
        Whether or not the sample was permethylated.
        
    lactonized_ethyl_esterified : boolean
        Whether the glycans were submitted to lactonization/ethyl-esterification
        derivatization, which differentiates the mass of alpha2-3 and alpha2-6 -bound 
        sialic acids.
        
    reduced : boolean
        Whether or not the sample was reduced.
        
    Uses
    ----
    Library_Tools.generate_glycans_library() : list
        A list containing dictionaries of the monosaccharides compositions of all the
        glycans generated.

    Library_Tools.full_glycans_library() : dict
        A dictionary with each key containing the glycan formula and each key containing
        a dictionary with monosaccharides composition, atoms composition with tag,
        neutral mass, neutral mass with tag, isotopic distribution and the mzs of the
        glycans with the desired adducts combination.
        
    pathlib.Path.resolve() : Path object
        Make the path absolute, resolving any symlinks. A new path object is returned

    Returns
    -------
    full_library : dict
        A dictionary with each key containing the glycan formula and each key containing
        a dictionary with monosaccharides composition, atoms composition with tag,
        neutral mass, neutral mass with tag, isotopic distribution and the mzs of the
        glycans with the desired adducts combination.
    '''
    # Grab all the args
    custom_glycans_list = args_dict.get('custom glycans list', [False, []])
    min_max_monos = args_dict.get('min/max monosaccharides', [0,0])
    min_max_hex = args_dict.get('min/max hexoses', [0,0])
    min_max_hexnac = args_dict.get('min/max hexnac', [0,0])
    min_max_xyl = args_dict.get('min/max xyloses', [0,0])
    min_max_sia = args_dict.get('min/max sialic acids', [0,0])
    min_max_fuc = args_dict.get('min/max fucoses', [0,0])
    min_max_ac = args_dict.get('min/max acetyl sialic acids', [0,0])
    min_max_gc = args_dict.get('min/max glycolyl sialic acids', [0,0])
    min_max_hn = args_dict.get('min/max hexosamines', [0,0])
    min_max_ua = args_dict.get('min/max uronic acids', [0,0])
    forced = args_dict.get('glycan class', None)
    min_max_proton = args_dict.get('min/max protons', [1,3])
    custom_adducts = args_dict.get('custom adducts', [])
    max_charges = args_dict.get('maximum charges', 3)
    tag_mass = args_dict.get('reducing end tag', 0.0)
    fast_iso = args_dict.get('fast isotopic pattern calculation', True)
    high_res = args_dict.get('high resolution isotopic pattern', False)
    imp_exp_library = args_dict.get('import/export library', [False, False])
    library_path = args_dict.get('library path', None)
    exp_lib_name = args_dict.get('exported library name', None)
    only_gen_lib = args_dict.get('only generate library', False)
    save_path = args_dict.get('save path', None)
    internal_standard = args_dict.get('internal standard', 0.0)
    permethylated = args_dict.get('permethylated', False)
    lactonized_ethyl_esterified = args_dict.get('lactonized/ethyl-esterified', False)
    reduced = args_dict.get('reducing end reduced', False)
    min_max_sulfation = args_dict.get('min/max sulfation', [0,0])
    min_max_phosphorylation = args_dict.get('min/max phosphorylation', [0,0])
    lyase_digested = args_dict.get('lyase digested', False)
    temp_folder = args_dict.get('temporary folder', None)
    custom_monos = args_dict.get('custom monosaccharides', [])
    from_GUI = args_dict.get('from GUI', False)
    
    monosaccharides = copy.deepcopy(General_Functions.monosaccharides)
    
    # Add custom monosaccharides
    if len(custom_monos) > 0:
        for cm in custom_monos:
            if cm['cm_short_code'] not in monosaccharides.keys():
                monosaccharides[cm['cm_short_code']] = (cm['cm_name'], cm['cm_chem_comp'], General_Functions.sum_atoms({"C": 0, "O": 0, "N": 0, "H": 0}, General_Functions.form_to_comp_atoms(cm['cm_chem_comp'])), cm['cm_single_letter_code'])
    
    # Define the start time
    date = datetime.datetime.now()
    begin_time = str(date)[2:4]+str(date)[5:7]+str(date)[8:10]+"_"+str(date)[11:13]+str(date)[14:16]+str(date)[17:19]
    date, time_lib = begin_time.split("_")
    is_custom = False
    
    # Generate adduct combos
    adduct_combos = General_Functions.gen_adducts_combo(min_max_proton, custom_adducts, max_charges)
    
    # If importing a library
    if imp_exp_library[0]:
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        print(time_formatted+'Importing existing library...', end = '', flush = True)
        try:
            with open(library_path, 'rb') as f:
                library_data = dill.load(f)
                f.close()
            full_library = library_data[0]
            library_metadata = library_data[1]
            
            # Load data from the metadata
            if library_metadata.get('custom glycans list', [False, []])[0]:
                is_custom = True
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
            min_max_sulfation = library_metadata.get('min/max sulfation', [0, 0])
            min_max_phosphorylation = library_metadata.get('min/max phosphorylation', [0, 0])
            forced = library_metadata.get('glycan class', None)
            min_max_proton = library_metadata.get('min/max protons', [0, 0])
            max_charges = library_metadata.get('maximum charges', 3)
            adduct_combos = library_metadata.get('adduct combos', [])
            tag_mass = library_metadata.get('reducing end tag', 0.0)
            internal_standard = library_metadata.get('internal standard', '0.0')
            permethylated = library_metadata.get('permethylated', False)
            lactonized_ethyl_esterified = library_metadata.get('lactonized/ethyl esterified', False)
            reduced = library_metadata.get('reducing end reduced', False)
            fast_iso = library_metadata.get('fast isotopic pattern calculation', True)
            high_res = library_metadata.get('high resolution isotopic pattern', False)
            lyase_digested = library_metadata.get('lyase digested', False)
            custom_monos = library_metadata.get('custom monosaccharides', [])
                
            print("Done!")
        except:
            print("\n\nGlycan library file not found. Check if the\npath used in 'library_path' is correct or set\n'imp_library' to 'no'.\n")
            if os.isatty(0):
                input("\nPress Enter to exit.")
                os._exit(1)
            else:
                print("Close the window or press CTRL+C to exit.")
                try:
                    while True:
                        time.sleep(3600)
                except KeyboardInterrupt:
                    os._exit(1)
                    
    # If not importing, but generating a library from custom glycans list
    elif custom_glycans_list[0] and not imp_exp_library[0]:
        # Here it checks if the monosaccharides inputted in the custom glycans list are valid
        custom_glycans_comp = []
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        print(time_formatted+'Building glycans library...', end = "", flush = True)
        for i in custom_glycans_list[1]:
            glycan_comp = General_Functions.form_to_comp_glycans(i)
            for i in glycan_comp:
                if i not in monosaccharides:
                    print(f"\n\nUnrecognized monosaccharide in glycan list: {i}\nCheck your custom glycans list.\n")
                    if os.isatty(0):
                        input("\nPress Enter to exit.")
                        os._exit(1)
                    else:
                        print("Close the window or press CTRL+C to exit.")
                        try:
                            while True:
                                time.sleep(3600)
                        except KeyboardInterrupt:
                            os._exit(1)
            custom_glycans_comp.append(General_Functions.sum_monos(glycan_comp))
            
        # Here it dynamically builds the constraints from the custom glycans list. This is important for fragment library building.
        monos = {}
        max_monos = 0
        min_monos = 99
        for i in monosaccharides:
            monos[i] = 0
        for i in custom_glycans_comp:
            count_monos = 0
            for j in i:
                count_monos+=i[j]
            if count_monos > max_monos:
                max_monos = count_monos
            if count_monos < min_monos:
                min_monos = count_monos
            for j in monos:
                if j in list(i.keys()) and i[j] > monos[j]:
                    monos[j] = i[j]
        if min_monos > max_monos:
            min_monos = 1
        min_max_monos = (min_monos-1, max_monos+1)
        min_max_hex = (0, (monos['H']+1) if monos['H'] != 0 else 0)
        min_max_hexnac = (0, (monos['N']+1) if monos['N'] != 0 else 0)
        min_max_fuc = (0, (monos['F']+1) if monos['F'] != 0 else 0)
        min_max_sia = (0, (max([monos['S'], monos['G'], monos['Am'], monos['E'], monos['AmG'], monos['EG']])+1) if (monos['S'] != 0 or monos['G'] != 0 or monos['Am'] != 0 or monos['E'] != 0 or monos['AmG'] != 0 or monos['EG'] != 0) else 0)
        min_max_ac = (0, (max([monos['S'], monos['Am'], monos['E']])+1) if (monos['S'] != 0 or monos['Am'] != 0 or monos['E'] != 0) else 0)
        min_max_gc = (0, (max([monos['G'], monos['AmG'], monos['EG']])+1) if (monos['G'] != 0 or monos['AmG'] != 0 or monos['EG'] != 0) else 0)
        min_max_xyl = (0, (monos['X']+1) if monos['X'] != 0 else 0)
        min_max_hn = (0, (monos['HN']+1) if monos['HN'] != 0 else 0)
        min_max_ua = (0, (monos['UA']+1) if monos['UA'] != 0 else 0)
        
        # Build custom library
        full_library = Library_Tools.full_glycans_library(custom_glycans_comp,
                                                          forced,
                                                          adduct_combos,
                                                          tag_mass,
                                                          fast_iso,
                                                          high_res,
                                                          internal_standard,
                                                          permethylated,
                                                          reduced,
                                                          min_max_sulfation,
                                                          min_max_phosphorylation,
                                                          lyase_digested,
                                                          custom_monos)
        print('Done!')
        
    # Here if it just wants to generate the whole library
    else:
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        print(time_formatted+'Building glycans library...', end = "", flush = True)
        monos_library = Library_Tools.generate_glycans_library(min_max_monos,
                                                               min_max_hex,
                                                               min_max_hexnac,
                                                               min_max_xyl,
                                                               min_max_sia,
                                                               min_max_fuc,
                                                               min_max_ac,
                                                               min_max_gc,
                                                               min_max_hn,
                                                               min_max_ua,
                                                               lactonized_ethyl_esterified,
                                                               forced,
                                                               custom_monos)
        full_library = Library_Tools.full_glycans_library(monos_library,
                                                          forced,
                                                          adduct_combos,
                                                          tag_mass,
                                                          fast_iso,
                                                          high_res,
                                                          internal_standard,
                                                          permethylated,
                                                          reduced,
                                                          min_max_sulfation,
                                                          min_max_phosphorylation,
                                                          lyase_digested,
                                                          custom_monos)
        custom_glycans_list[0] = False
        print('Done!')
        
    if is_custom:
        custom_glycans_list[0] = True
        
    if imp_exp_library[1] or only_gen_lib:
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        print(time_formatted+'Exporting glycans library...', end = '', flush = True)
        if exp_lib_name != '':
            if len(exp_lib_name.split('.')) > 1:
                exp_lib_name = exp_lib_name.split('.')[0]
            if "<date>" in exp_lib_name or "<time>" in exp_lib_name:
                temp_lib_name = []
                temp_lib_name_first = exp_lib_name.split('>')
                for i in temp_lib_name_first:
                    temp_lib_name += i.split('<')
                exp_lib_name = ''
                for word in temp_lib_name:
                    if word == 'date':
                        exp_lib_name += str(date)
                    elif word == 'time':
                        exp_lib_name += str(time_lib)
                    else:
                        exp_lib_name += word
                        
            exp_lib_name = exp_lib_name.replace("<", "_")
            exp_lib_name = exp_lib_name.replace(">", "_")
            
            counter = 0
            while True:
                if counter == 0 and os.path.isfile(os.path.join(save_path, exp_lib_name+'.ggl')):
                    counter+=1
                    continue
                elif counter != 0 and os.path.isfile(os.path.join(save_path, exp_lib_name+'('+str(counter)+').ggl')):
                    counter+=1
                    continue
                else:
                    if counter != 0:
                        exp_lib_name = exp_lib_name+'('+str(counter)+')'
                    else:
                        exp_lib_name = exp_lib_name
                    break
        else:
            exp_lib_name = begin_time+'_glycans_lib'
        if not imp_exp_library[0]:
            metadata = {
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
                        'adduct combos': adduct_combos,
                        'reducing end tag': tag_mass,
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
                        'custom monosaccharides': custom_monos
                        }
            if exp_lib_name+'.ggl' not in os.listdir(save_path):
                with open(os.path.join(save_path, exp_lib_name+'.ggl'), 'wb') as f:
                    dill.dump([full_library, metadata], f)
                    f.close()
        if not from_GUI:
            output_extra_library_files(full_library, lactonized_ethyl_esterified, custom_monos, imp_exp_library, library_path, exp_lib_name, metadata, save_path)
        print("Done!")
    if only_gen_lib:
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        print(time_formatted+'Library length: '+str(len(full_library)))
        print("File name is '"+exp_lib_name+".ggl'.")
        print("If you wish to analyze files,")
        print("set 'only_gen_lib' to False and input")
        print("remaining parameters.")
        
        shutil.rmtree(temp_folder)
        
        if os.isatty(0):
            input("\nPress Enter to exit.")
            os._exit(1)
        else:
            print("Close the window or press CTRL+C to exit.")
            try:
                while True:
                    time.sleep(3600)
            except KeyboardInterrupt:
                os._exit(1)
    time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
    print(time_formatted+'Library length: '+str(len(full_library)))
    return full_library, adduct_combos
    
def align_assignments(df, df_type, multithreaded, number_cores, temp_folder = None, gg_file = None, deltas = None, rt_tol = None, iso_fit_score = 0, curve_fit_score = 0, max_ppm = 0, s_to_n = 0):
    '''Aligns the results obtained from running the whole program and uses
    the identifications to align the chromatograms between samples for 
    visualization purposes.
    
    Parammeters
    -----------
    df : dict
        A results dataframe or chromatogram dataframe from output_filtered_data.
        
    df_type : string
        Type of the df ("total_glycans" result or "chromatograms").
        
    multithreaded : boolean
        Whether the execution is multithreaded or not.
        
    number_cores : int
        Number of cores used if execution is multithreaded.
        
    deltas : dict
        The identified delta-t per sample, for use in aligning the chromatograms.
        'deltas' obtained by running this function in "total_glycans" mode for df_type.
        
    rt_tol : float
        Retention time tolerance used for identifying glycan identification from different
        files as the same.
        
    Uses
    ----
    copy.deepcopy : same as input variable
        Creates a copy of a variable unlinked to the original variable.
        
    General_Functions.linear_regression : tuple
        Returns the slope, y-intercept and outliers of linear fit of given data.
        
    Returns
    -------
    dataframe : dict
        The aligned results or chromatograms dataframe.
        
    deltas_per_sample : dict
        A dictionary containing the deltas identified per sample.
        
    biggest_df : int
        The sample at which most identifications were found.
    '''
    dataframe = copy.deepcopy(df)
    if df_type == 'total_glycans': #this is for alignment of total_glycans dataframes... it generates delta values to be used for other alignments
        biggest_df = inf
        biggest_df_size = 0
        for i_i, i in enumerate(dataframe): #this determines the dataframe with most peaks assigned, it will be used as reference
            if len(i['Glycan']) > biggest_df_size:
                biggest_df_size = len(i['Glycan'])
                biggest_df = i_i
        if biggest_df == inf:
            biggest_df = 0
        ids_per_sample = []
        deltas_per_sample = [] #this is a list that will contain every delta of aligned glycans per sample, in a dictionary of 'rt : delta'
        skips = 0
        for i_i, i in enumerate(dataframe[biggest_df]['Glycan']): #goes through each glycan in the reference df to align the other samples with
            if skips != 0:
                skips -= 1
                continue
            peaks_number = dataframe[biggest_df]['Glycan'].count(i)
            skips = peaks_number-1
            rts_reference = dataframe[biggest_df]['RT'][i_i:i_i+peaks_number] #this collects the reference rts for the glycan in reference sample
            aucs_reference = []
            for j in range(i_i, i_i+peaks_number):
                aucs_reference.append(dataframe[biggest_df]['AUC'][j]/max(dataframe[biggest_df]['AUC'][i_i:i_i+peaks_number]))
            for j_j, j in enumerate(dataframe): #going through each target sample
                if i_i == 0:
                    ids_per_sample.append([])
                    deltas_per_sample.append({})
                if j_j == biggest_df:
                    continue
                peaks_number_target = j['Glycan'].count(i)
                if peaks_number_target == 0: #glycan not available in target sample
                    continue
                first_glycan_id = 0
                rts = []
                aucs = []
                for k_k, k in enumerate(j['Glycan']): #going through each glycan of target sample
                    if k == i:
                        for l in range(k_k, k_k+peaks_number_target):
                            ids_per_sample[j_j].append(l)
                            rts.append(j['RT'][l])
                            aucs.append(j['AUC'][l]/max(j['AUC'][k_k:k_k+peaks_number_target]))
                        used_peaks = []
                        last_assigned_rt_id = 0
                        for l_l, l in enumerate(aucs):
                            found = False
                            if len(aucs_reference[last_assigned_rt_id:]) != 0: 
                                for m_m in range(last_assigned_rt_id, len(aucs_reference)):
                                    m = aucs_reference[m_m]
                                    if m_m not in used_peaks and abs(l-m) <= 0.5*m:
                                        if m_m > aucs_reference.index(1.0) and l_l < aucs.index(1.0):
                                            break
                                        found = True
                                        last_assigned_rt_id = m_m+1
                                        deltas_per_sample[j_j][rts[l_l]] = [rts_reference[m_m]-rts[l_l], k]
                                        break
                            if not found:
                                deltas_per_sample[j_j][rts[l_l]] = [inf, k]
                        break
                          
        # print(deltas_per_sample)                         
        # print()
        
        # #pre-cleanup step to remove odd variations
        # for i_i, i in enumerate(deltas_per_sample):
            # if len(i) > 0:
                # x_list = []
                # y_list = []
                # for j_j, j in enumerate(i):
                    # if i[j][0] != inf:
                        # x_list.append(j)
                        # y_list.append(i[j][0])
                # if len(x_list) > 0:
                    # outliers = General_Functions.linear_regression(x_list, y_list)[2]
                    # if len(outliers) > 0:
                        # for j_j, j in enumerate(outliers):
                            # i[x_list[j]][0] = inf
                            
        #print(deltas_per_sample)                         
        #print()
                        
        for i_i, i in enumerate(deltas_per_sample):
            if len(i) == 0:
                continue
            negative_ids = []
            positive_ids = []
            all_deltas = []
            for j_j, j in enumerate(i):
                all_deltas.append(i[j][0])
            for j_j, j in enumerate(i):
                if i[j][0] != inf:
                    if i[j][0] > 0:
                        positive_ids.append(j)
                    if i[j][0] < 0:
                        negative_ids.append(j)
            if len(negative_ids) > len(positive_ids):
                for j_j, j in enumerate(positive_ids):
                    i[j][0] = inf
            elif len(negative_ids) < len(positive_ids):
                for j_j, j in enumerate(negative_ids):
                    i[j][0] = inf
            elif mean(all_deltas) < rt_tol:
                for j_j, j in enumerate(i):
                    i[j][0] = 0
            elif len(negative_ids) == len(positive_ids):
                for j_j, j in enumerate(i):
                    i[j][0] = mean(all_deltas)
                    
        #print(deltas_per_sample[2])                         
        #print()
                    
        for i_i, i in enumerate(deltas_per_sample):
            if len(i) > 0:
                x_list = []
                y_list = []
                for j_j, j in enumerate(i):
                    if i[j][0] != inf:
                        x_list.append(j)
                        y_list.append(i[j][0])
                if len(x_list) > 0:
                    outliers = General_Functions.linear_regression(x_list, y_list, 2)[2]
                    if len(outliers) > 0:
                        for j_j, j in enumerate(outliers):
                            i[x_list[j]][0] = inf
                           
        #print(deltas_per_sample[2])                         
        #print()
        
        for i_i, i in enumerate(dataframe): #this will apply the alignment based on ids and deltas
            rts_list_original = sorted(i['RT'])
            rts_list_adjusted = copy.deepcopy(rts_list_original)
            if len(ids_per_sample[i_i]) == 0: #this is the reference sample or blank sample
                continue
            to_fix_rt = []
            to_fix_id = []
            for j_j, j in enumerate(i['RT']):
                found = False
                for k_k, k in enumerate(deltas_per_sample[i_i]):
                    if i['RT'][j_j] == k and i['Glycan'][j_j] == deltas_per_sample[i_i][k][1] and deltas_per_sample[i_i][k][0] != inf:
                        rts_list_adjusted[rts_list_original.index(i['RT'][j_j])] = float("%.2f" % round(i['RT'][j_j] + deltas_per_sample[i_i][i['RT'][j_j]][0], 2))
                        i['RT'][j_j] = float("%.2f" % round(i['RT'][j_j] + deltas_per_sample[i_i][i['RT'][j_j]][0], 2))
                        found = True
                        break
                if not found:
                    to_fix_rt.append(j)
                    to_fix_id.append(j_j)
            if len(to_fix_rt) > 0:
                for j_j, j in enumerate(rts_list_adjusted):
                    if j in to_fix_rt:
                        before_j_j = None
                        after_j_j = None
                        if j_j > 0:
                            for k in range(j_j-1, -1, -1):
                                if rts_list_adjusted[k] not in to_fix_rt:
                                    before_j_j = rts_list_original[k]
                                    before_j_j_id = k
                        if j_j < len(rts_list_adjusted)-1:
                            for k in range(j_j+1, len(rts_list_adjusted)):
                                if rts_list_adjusted[k] not in to_fix_rt:
                                    after_j_j = rts_list_original[k]
                                    after_j_j_id = k
                        x = []
                        y = []
                        for k_k, k in enumerate(deltas_per_sample[i_i]):
                            if deltas_per_sample[i_i][k][0] != inf:
                                x.append(float(k))
                                y.append(deltas_per_sample[i_i][k][0])
                        linear_equation = General_Functions.linear_regression(x, y)
                        fixed_RT = float("%.2f" % round(j + (j*linear_equation[0]) + linear_equation[1], 2))
                        if before_j_j != None and after_j_j != None:
                            scaling_factor = (j-before_j_j)/(after_j_j-before_j_j)
                            fixed_RT = float("%.2f" % round(rts_list_adjusted[before_j_j_id]+(scaling_factor*(rts_list_adjusted[after_j_j_id]-rts_list_adjusted[before_j_j_id])), 2))
                        elif before_j_j == None and after_j_j != None:
                            scaling_factor = j/after_j_j
                            fixed_RT = float("%.2f" % round(rts_list_adjusted[after_j_j_id]*scaling_factor, 2))
                        elif before_j_j != None and after_j_j == None:
                            scaling_factor = j/before_j_j
                            fixed_RT = float("%.2f" % round(rts_list_adjusted[before_j_j_id]*scaling_factor, 2))
                        if j in list(deltas_per_sample[i_i].keys()):
                            deltas_per_sample[i_i][j][0] = fixed_RT - j
                        elif j not in list(deltas_per_sample[i_i].keys()):
                            deltas_per_sample[i_i][j] = [fixed_RT - j]
                        rts_list_adjusted[j_j] = fixed_RT
                        i['RT'][to_fix_id[to_fix_rt.index(j)]] = fixed_RT
        for i_i, i in enumerate(deltas_per_sample):
            deltas_per_sample[i_i] = dict(sorted(i.items()))
        
        # print(deltas_per_sample[2])
            
        return dataframe, deltas_per_sample, biggest_df
        
    if df_type == "chromatograms": #for when you want to align chromatograms based on existing deltas calculated previously
    
        if multithreaded:
            if number_cores == 'all':
                cpu_count = (os.cpu_count())-2
                if cpu_count <= 0:
                    cpu_count = 1
            else:
                number_cores = int(number_cores)
                if number_cores > (os.cpu_count())-2:
                    cpu_count = (os.cpu_count())-2
                    if cpu_count <= 0:
                        cpu_count = 1
                else:
                    cpu_count = number_cores
        else:
            cpu_count = 1
            
        results = []
        with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_count if cpu_count < 60 else 60) as executor:
            for i_i, i in enumerate(dataframe['File_Name']): #sample by sample
                result = executor.submit(adjust_chromatogram,
                                         i,
                                         i_i,
                                         deltas, 
                                         iso_fit_score,
                                         curve_fit_score,
                                         max_ppm,
                                         s_to_n,
                                         temp_folder,
                                         gg_file)
                results.append(result)
            
            for index, i in enumerate(results):
                result_data = i.result()
                results[index] = None

def adjust_chromatogram(i,
                        i_i,
                        deltas,
                        iso_fit_score,
                        curve_fit_score,
                        max_ppm,
                        s_to_n,
                        temp_folder,
                        gg_file):
    '''Actual alignment of chromatogram. This function does just the adjustment part of the alignment for a given chromatogram. Used in the concurrent.futures for alignment.
    
    Parameters
    ----------
    i : str
        Sample chromatogram.
        
    i_i : int
        Sample index.
        
    deltas : dict
        The identified delta-t per sample, for use in aligning the chromatograms.
        'deltas' obtained by running this function in "total_glycans" mode for df_type.
        
    Uses
    ----
    General_Functions.linear_regression : tuple
        Returns the slope, y-intercept and outliers of linear fit of given data.
        
    Returns
    -------
    i : str
        Adjusted sample chromatogram.
        
    i_i : int
        Sample index.
    '''
    eic_name = 'RTs'
    sample_RTs = General_Functions.access_chromatogram(i_i, f"{i_i}_{eic_name}", temp_folder, gg_file)
        
    if len(deltas[i_i]) > 0:
        
        if f"eics_list" not in os.listdir(temp_folder):
            General_Functions.open_gg(gg_file, temp_folder, file = f"eics_list")
        
        with open(os.path.join(temp_folder, f"eics_list"), 'rb') as f:
            loaded_eics_list = dill.load(f)
            chromatograms_list = loaded_eics_list[i_i]
            f.close()
            
        chromatogram_length_rt = sample_RTs[-1]
        chromatogram_beg_rt = sample_RTs[0]
        chromatogram_length = len(sample_RTs)
        chromatogram_interval = sample_RTs[-1]/len(sample_RTs)
        points_per_minute = int(1/chromatogram_interval)
        interval_list_rts = []
        interval_list = []
        for j in range(len(sample_RTs)-1, -1, -1):
            if sample_RTs[j] < list(deltas[i_i].keys())[0]: #this finds the zero before the peaks
                zero = True
                for k_k, k in enumerate(chromatograms_list):
                    if k_k != 0:
                        
                        # Load the target chromatogram
                        target_chromatogram = General_Functions.access_chromatogram(i_i, f"{i_i}_smoothed_{k}", temp_folder, gg_file)
                            
                        if target_chromatogram[j] > max(target_chromatogram)*0.01:
                            zero = False
                            break
                if zero:
                    interval_list_rts.append(sample_RTs[j])
                    interval_list.append(j)
                    break
        if not zero:
            interval_list_rts.append(sample_RTs[0])
            interval_list.append(0)
        for j_j, j in enumerate(sample_RTs):
            if j > list(deltas[i_i].keys())[-1]: #this finds the zero after the peaks
                zero = True
                for k_k, k in enumerate(chromatograms_list):
                    if k_k != 0:
                        
                        # Load the target chromatogram
                        target_chromatogram = General_Functions.access_chromatogram(i_i, f"{i_i}_smoothed_{k}", temp_folder, gg_file)
                            
                        if target_chromatogram[j_j] > max(target_chromatogram)*0.1:
                            zero = False
                            break
                if zero:
                    interval_list_rts.append(j)
                    interval_list.append(j_j)
                    break
        if not zero:
            interval_list_rts.append(sample_RTs[-1])
            interval_list.append(len(sample_RTs)-1)
        for j_j, j in enumerate(interval_list):
            x = []
            y = []
            if j_j == len(interval_list)-1:
                last_range = len(sample_RTs)-1
            else:
                last_range = interval_list[j_j+1]
            for k_k, k in enumerate(deltas[i_i]):
                if k > sample_RTs[j] and k < sample_RTs[last_range]:
                    x.append(float(k))
                    y.append(deltas[i_i][k][0])
            if len(x) == 0:
                continue
            linear_equation = General_Functions.linear_regression(x, y)
            lower_boundary = sample_RTs[j-1]
            upper_boundary = sample_RTs[last_range]
            lowest = inf
            lowest_id = 0
            highest = 0
            highest_id = inf
            for k in range(j, last_range):
                sample_RTs[k] = sample_RTs[k] + (sample_RTs[k]*linear_equation[0])+linear_equation[1]
                if sample_RTs[k] < lowest:
                    lowest = sample_RTs[k]
                    lowest_id = k
                if sample_RTs[k] > highest:
                    highest = sample_RTs[k]
                    highest_id = k
            if lowest != inf and lowest < lower_boundary:
                if j_j > 0:
                    interval = (lowest-sample_RTs[interval_list[j_j-1]])/(lowest_id+1-interval_list[j_j-1])
                    for l_l, l in enumerate(range(lowest_id-1, interval_list[j_j-1], -1)):
                        sample_RTs[l] = float("%.4f" % round(lowest-(interval*(l_l+1)), 4))
                else:
                    interval = (lowest)/(lowest_id+1)   
                    for l_l, l in enumerate(range(lowest_id-1, -1, -1)):
                        sample_RTs[l] = float("%.4f" % round(lowest-(interval*(l_l+1)), 4))
            if highest != 0 and highest > upper_boundary:
                if j_j != len(interval_list)-2:
                    interval = (sample_RTs[interval_list[j_j+2]]-highest)/(interval_list[j_j+2]-highest_id)
                    for l_l, l in enumerate(range(highest_id+1, interval_list[j_j+1])):
                        sample_RTs[l] = float("%.4f" % round(highest+(interval*(l_l+1)), 4))
                else:
                    interval = (chromatogram_length_rt-highest)/(chromatogram_length-1-highest_id)
                    for l_l, l in enumerate(range(highest_id+1, chromatogram_length)):
                        sample_RTs[l] = float("%.4f" % round(highest+(interval*(l_l+1)), 4))
    
    with open(os.path.join(temp_folder, f"{i_i}_aligned_{eic_name}_{iso_fit_score}_{curve_fit_score}_{max_ppm}_{s_to_n}"), 'wb') as f:
        dill.dump(sample_RTs, f)
        f.close()
    
    return None
        
def make_df1_refactor(df1,
                      df2,
                      curve_fit_score,
                      iso_fit_score,
                      sn,
                      max_ppm,
                      percentage_auc,
                      analyze_ms2,
                      unrestricted_fragments,
                      min_samples,
                      fragments_dataframes = [],
                      fill_gaps = (False, 50, 0.2, False),
                      sample_groups = {},
                      noise_levels = []):
    '''Reorganizes the raw data into a more comprehensible format and filter by the quality thresholds.

    Parameters
    ----------
    df1 : list
        List of dictionaries containing the raw results data for every sample.
        
    df2 : list
        A dictionary containing information about analyzed files.
        
    curve_fit_score : float
        The minimum curve fitting score to consider when outputting results.
        
    iso_fit_score : float
        The minimum isotopic fitting score to consider when outputting results.
        
    sn : int
        The minimum signal-to-noise ration to consider when outputting results.
        
    max_ppm : int
        The maximum amount of PPM difference to consider when outputting results.
        
    percentage_auc : float
        Percentage of the highest peak AUC that another peak AUC must have in a same chromatogram in order to be saved.
        
    analyze_ms2 : boolean
        Whether MS2 was analyzed or not.
    
    unrestricted_fragments : boolean
        Whether or not fragments are restricted to the precursor putative composition.
        
    min_samples : int
        Parameter to remove glycans not present in at least this number of samples.
        
    fragments_dataframes : list
        Dataframe containing fragments information.
        
    Returns
    -------
    df1_refactor : list
        A list of dictionaries containing organized and filtered data.
        
    fragments_dataframes : list
        Dataframe containing fragments information.
    '''
    df1_refactor = []
    remove_sub_id = []
    remove_glycan = []
    remove_adduct = []
    
    # Remove ungrouped samples before analysis
    if len(sample_groups) != 0:
        
        # Extract sample:group dictionary from path to groups or from grouping dictionary
        sample_group_dict = {}
        if type(sample_groups) != dict:
            with open(sample_groups, "r") as file:
                for index, line in enumerate(file):
                    if index == 0:
                        continue
                    sample, group = line.strip().split(",")
                    sample_group_dict[sample] = group
        else:
            for group, samples in sample_groups.items():
                for sample in samples:
                    sample_group_dict[sample] = group
        
        # Remove the samples
        if len(sample_group_dict) != 0:
            sample_to_remove = []
            for index, file_name in enumerate(df2['File_Name']):
                if file_name not in sample_group_dict.keys():
                    sample_to_remove.append(index)
                    
            sample_to_remove.reverse()
            for index in sample_to_remove:
                for key in df2:
                    del df2[key][index]
                del df1[index]
                if analyze_ms2:
                    del fragments_dataframes[index]
                if len(noise_levels) != 0:
                    del noise_levels[0][index]
                    del noise_levels[1][index]
                
            df2['Sample_Number'] = [index for index, data in enumerate(df2['File_Name'])]
            if len(noise_levels) > 0:
                noise_levels[0] = {index:values for index, values in enumerate(noise_levels[0].values())}
                noise_levels[1] = {index:values for index, values in enumerate(noise_levels[1].values())}
            sample_groups = sample_group_dict
        else:
            sample_groups = {}
    
    for i_i, i in enumerate(df2["Sample_Number"]): #QCs cutoff
        remove_sub_id.append([])
        remove_glycan.append([])
        remove_adduct.append([])
        
        df1_refactor.append({"Glycan" : [], "Adduct" : [], "mz" : [], "RT" : [], "AUC" : [], "PPM" : [], "S/N" : [], "Iso_Fitting_Score" : [], "Curve_Fitting_Score" : []})
        if analyze_ms2:
            df1_refactor[i_i]["Detected_Fragments"] = []
        for j_j, j in enumerate(df1[i_i]["Adduct"]): 
            temp_rt = df1[i_i]["RT"][j_j]
            temp_auc = df1[i_i]["AUC"][j_j]
            temp_ppm = df1[i_i]["PPM"][j_j]
            temp_sn = df1[i_i]["S/N"][j_j]
            temp_fit = df1[i_i]["Iso_Fitting_Score"][j_j]
            temp_curve = df1[i_i]["Curve_Fitting_Score"][j_j]
            
            to_remove = []
            to_remove_glycan = []
            to_remove_adduct = []
            
            # This pass selects specific glycan+adduct combos to remove if they don't pass the thresholds
            for k_k, k in enumerate(temp_sn):
                
                marked_for_removal = False
                
                if df1[i_i]["Glycan"][j_j] != "Internal Standard":
                    if float(k) < sn and not marked_for_removal:
                        marked_for_removal = True
                    if float(temp_fit[k_k]) < iso_fit_score and not marked_for_removal:
                        marked_for_removal = True
                    if float(temp_curve[k_k]) < curve_fit_score and not marked_for_removal:
                        marked_for_removal = True
                    if type(max_ppm) == float:
                        if abs(float(temp_ppm[k_k])) > max_ppm and not marked_for_removal:
                            marked_for_removal = True
                    else:
                        if (float(temp_ppm[k_k]) < max_ppm[0] or float(temp_ppm[k_k]) > max_ppm[1]) and not marked_for_removal:
                            marked_for_removal = True
                            
                if marked_for_removal:
                    
                    # Check if glycan peak is also bad in other samples to fill gap... must be good in at least half the samples
                    if fill_gaps[0]:
                        if len(sample_groups) == 0:
                            indexes = df2["Sample_Number"]
                        else: # If samples are separated into groups
                            # Identify the sample group
                            sample_group = sample_groups[df2['File_Name'][i_i]]
                            indexes = []
                            for index, _sample_group in enumerate(sample_groups.items()):
                                sample = _sample_group[0]
                                group = _sample_group[1]
                                if group == sample_group:
                                    indexes.append(index)
                                    
                        remove = True
                        good_count = 0
                        for sample_index in indexes:
                            if sample_index == i_i:
                                continue
                            good_in_sample = False
                            for glycan_adduct_index, glycan_adduct in enumerate(df1[sample_index]["Adduct"]):
                                if df1[sample_index]["Glycan"][glycan_adduct_index] == df1[i_i]["Glycan"][j_j] and glycan_adduct == j:
                                    
                                    temp_rt_fill_gap = df1[sample_index]["RT"][glycan_adduct_index]
                                    temp_auc_fill_gap = df1[sample_index]["AUC"][glycan_adduct_index]
                                    temp_ppm_fill_gap = df1[sample_index]["PPM"][glycan_adduct_index]
                                    temp_sn_fill_gap = df1[sample_index]["S/N"][glycan_adduct_index]
                                    temp_fit_fill_gap = df1[sample_index]["Iso_Fitting_Score"][glycan_adduct_index]
                                    temp_curve_fill_gap = df1[sample_index]["Curve_Fitting_Score"][glycan_adduct_index]
                                    
                                    for score_index, rt in enumerate(temp_rt_fill_gap):
                                        if abs(rt-temp_rt[k_k]) < fill_gaps[2]:
                                            peak_good = True
                                            if float(temp_sn_fill_gap[score_index]) < sn and peak_good:
                                                peak_good = False
                                            if float(temp_fit_fill_gap[score_index]) < iso_fit_score and peak_good:
                                                peak_good = False
                                            if float(temp_curve_fill_gap[score_index]) < curve_fit_score and peak_good:
                                                peak_good = False
                                            if type(max_ppm) == float:
                                                if abs(float(temp_ppm_fill_gap[score_index])) > max_ppm and peak_good:
                                                    peak_good = False
                                            else:
                                                if (float(temp_ppm_fill_gap[score_index]) < max_ppm[0] or float(temp_ppm_fill_gap[score_index]) > max_ppm[1]) and peak_good:
                                                    peak_good = False
                                            if peak_good:
                                                good_count += 1
                                                good_in_sample = True
                                                if good_count >= round(len(indexes) * (fill_gaps[1]/100)):
                                                    remove = False
                                                break
                                if not remove or good_in_sample:
                                    break
                            if not remove:
                                break
                        if not remove:
                            continue
                            
                    to_remove.append(k_k)
                    to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                    to_remove_adduct.append(j)
                    
            remove_sub_id[i_i].append(to_remove)
            remove_glycan[i_i].append(to_remove_glycan)
            remove_adduct[i_i].append(to_remove_adduct)
            
    # Here the specific glycans+adduct combos are removed
    for i_i, i in enumerate(df2["Sample_Number"]):
        for j_j, j in enumerate(df1[i_i]["Adduct"]): 
            temp_rt = df1[i_i]["RT"][j_j]
            temp_auc = df1[i_i]["AUC"][j_j]
            temp_ppm = df1[i_i]["PPM"][j_j]
            temp_sn = df1[i_i]["S/N"][j_j]
            temp_fit = df1[i_i]["Iso_Fitting_Score"][j_j]
            temp_curve = df1[i_i]["Curve_Fitting_Score"][j_j]
            
            to_remove = remove_sub_id[i_i][j_j]
            to_remove_glycan = remove_glycan[i_i][j_j]
            to_remove_adduct = remove_adduct[i_i][j_j]
            
            if len(to_remove) != 0:
                to_remove.reverse()
                to_remove_glycan.reverse()
                to_remove_adduct.reverse()
                for k_k in to_remove:
                    del temp_rt[k_k]
                    del temp_auc[k_k]
                    del temp_ppm[k_k]
                    del temp_sn[k_k]
                    del temp_fit[k_k]
                    del temp_curve[k_k]
                    if analyze_ms2:
                        if "Detected_Fragments" not in list(df1[i_i].keys()):
                            print("\nThe data you are trying to reanalyze doesn't\ncontain MS2 data. Set 'analyze_ms2' to 'no' and\ntry again.\n")
                            if os.isatty(0):
                                input("\nPress Enter to exit.")
                                os._exit(1)
                            else:
                                print("Close the window or press CTRL+C to exit.")
                                try:
                                    while True:
                                        time.sleep(3600)
                                except KeyboardInterrupt:
                                    os._exit(1)
                        if len(temp_rt) == 0:
                            df1[i_i]["Detected_Fragments"][j_j] = ""
                            for k in range(len(fragments_dataframes[i_i]["Glycan"])-1, -1, -1):
                                if (fragments_dataframes[i_i]["Glycan"][k] == to_remove_glycan[k_k] 
                                    and fragments_dataframes[i_i]["Adduct"][k] == to_remove_adduct[k_k]
                                    and not unrestricted_fragments):
                                    del fragments_dataframes[i_i]["Glycan"][k]
                                    del fragments_dataframes[i_i]["Adduct"][k]
                                    del fragments_dataframes[i_i]["Fragment"][k]
                                    del fragments_dataframes[i_i]["Fragment_mz"][k]
                                    del fragments_dataframes[i_i]["Fragment_Intensity"][k]
                                    del fragments_dataframes[i_i]["RT"][k]
                                    del fragments_dataframes[i_i]["Precursor_mz"][k]
                                    del fragments_dataframes[i_i]["% TIC explained"][k]
                                    
            # Here specific glycan+adduct combos are selected for removal if they don't meet the minimum peak intensity threshold
            to_remove = []
            to_remove_glycan = []
            to_remove_adduct = []  
            for k_k, k in enumerate(temp_sn): 
                if max(temp_auc) == 0.0:
                    to_remove.append(k_k)
                    to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                    to_remove_adduct.append(j)
                    continue
                if float(temp_auc[k_k]/max(temp_auc)) <= percentage_auc:
                    to_remove.append(k_k)
                    to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                    to_remove_adduct.append(j)
                    continue
                    
            # Here they are removed
            if len(to_remove) != 0:
                to_remove.reverse()
                to_remove_glycan.reverse()
                to_remove_adduct.reverse()
                for k_k in to_remove:
                    del temp_rt[k_k]
                    del temp_auc[k_k]
                    del temp_ppm[k_k]
                    del temp_sn[k_k]
                    del temp_fit[k_k]
                    del temp_curve[k_k]
                    if analyze_ms2:
                        if len(temp_rt) == 0:
                            df1[i_i]["Detected_Fragments"][j_j] = ""
                            for k in range(len(fragments_dataframes[i_i]["Glycan"])-1, -1, -1):
                                if (fragments_dataframes[i_i]["Glycan"][k] == to_remove_glycan[k_k] 
                                    and fragments_dataframes[i_i]["Adduct"][k] == to_remove_adduct[k_k]
                                    and not unrestricted_fragments):
                                    del fragments_dataframes[i_i]["Glycan"][k]
                                    del fragments_dataframes[i_i]["Adduct"][k]
                                    del fragments_dataframes[i_i]["Fragment"][k]
                                    del fragments_dataframes[i_i]["Fragment_mz"][k]
                                    del fragments_dataframes[i_i]["Fragment_Intensity"][k]
                                    del fragments_dataframes[i_i]["RT"][k]
                                    del fragments_dataframes[i_i]["Precursor_mz"][k]
                                    del fragments_dataframes[i_i]["% TIC explained"][k]
                                    
            # Filtered data is added back to the dataframe
            df1[i_i]["RT"][j_j] = str(temp_rt)[1:-1]
            df1[i_i]["AUC"][j_j] = str(temp_auc)[1:-1]
            df1[i_i]["PPM"][j_j] = str(temp_ppm)[1:-1]
            df1[i_i]["S/N"][j_j] = str(temp_sn)[1:-1]
            df1[i_i]["Iso_Fitting_Score"][j_j] = str(temp_fit)[1:-1]
            df1[i_i]["Curve_Fitting_Score"][j_j] = str(temp_curve)[1:-1]
    
    # Empty glycans+adduct combos (with no remaining peaks) are removed
    to_remove = []
    to_remove_glycan = []
    to_remove_adduct = []
    for i_i, i in enumerate(df1[0]["Adduct"]):
        glycan_good = False
        for j_j, j in enumerate(df2["Sample_Number"]):
            check = df1[j_j]["RT"][i_i]
            if len(check) != 0:
                glycan_good = True
                break
        if not glycan_good:
            to_remove.append(i_i)
            to_remove_glycan.append(df1[j_j]["Glycan"][i_i])
            to_remove_adduct.append(i)
    if len(to_remove) != 0:
        to_remove.reverse()
        to_remove_glycan.reverse()
        to_remove_adduct.reverse()
        for j_j, j in enumerate(df2["Sample_Number"]):
            for i_index, i_i in enumerate(to_remove):
                del df1[j_j]["Glycan"][i_i]
                del df1[j_j]["Adduct"][i_i]
                del df1[j_j]["mz"][i_i]
                del df1[j_j]["RT"][i_i]
                del df1[j_j]["AUC"][i_i]
                del df1[j_j]["PPM"][i_i]
                del df1[j_j]["S/N"][i_i]
                del df1[j_j]["Iso_Fitting_Score"][i_i]
                del df1[j_j]["Curve_Fitting_Score"][i_i]
                if analyze_ms2:
                    del df1[j_j]["Detected_Fragments"][i_i]
                    for k in range(len(fragments_dataframes[j_j]["Glycan"])-1, -1, -1):
                        if (fragments_dataframes[j_j]["Glycan"][k] == to_remove_glycan[i_index] 
                            and fragments_dataframes[j_j]["Adduct"][k] == to_remove_adduct[i_index]
                            and not unrestricted_fragments):
                            del fragments_dataframes[j_j]["Glycan"][k]
                            del fragments_dataframes[j_j]["Adduct"][k]
                            del fragments_dataframes[j_j]["Fragment"][k]
                            del fragments_dataframes[j_j]["Fragment_mz"][k]
                            del fragments_dataframes[j_j]["Fragment_Intensity"][k]
                            del fragments_dataframes[j_j]["RT"][k]
                            del fragments_dataframes[j_j]["Precursor_mz"][k]
                            del fragments_dataframes[j_j]["% TIC explained"][k] #QCs cutoff end
                            
    # Creation of the reorganized dataframe, where each peak has an individual line
    for i_i, i in enumerate(df1):
        for j_j, j in enumerate(df1[i_i]["Adduct"]):
            for k_k, k in enumerate(df1[i_i]["RT"][j_j].split(", ")):
                if k != "":
                    df1_refactor[i_i]["Glycan"].append(df1[i_i]["Glycan"][j_j])
                    df1_refactor[i_i]["Adduct"].append(df1[i_i]["Adduct"][j_j])
                    df1_refactor[i_i]["mz"].append(df1[i_i]["mz"][j_j])
                    df1_refactor[i_i]["RT"].append(float(k))
                    df1_refactor[i_i]["AUC"].append(float(df1[i_i]["AUC"][j_j].split(", ")[k_k]))
                    df1_refactor[i_i]["PPM"].append(float(df1[i_i]["PPM"][j_j].split(", ")[k_k]))
                    df1_refactor[i_i]["S/N"].append(float(df1[i_i]["S/N"][j_j].split(", ")[k_k]))
                    df1_refactor[i_i]["Iso_Fitting_Score"].append(float(df1[i_i]["Iso_Fitting_Score"][j_j].split(", ")[k_k]))
                    df1_refactor[i_i]["Curve_Fitting_Score"].append(float(df1[i_i]["Curve_Fitting_Score"][j_j].split(", ")[k_k])) 
    
    # Filter glycans that are not found in x number of samples
    sample_groups_for_filter = {}
    for index, sample in enumerate(df2['File_Name']):
        group = sample_groups.get(sample, 'ungrouped')
        if group not in sample_groups_for_filter:
            sample_groups_for_filter[group] = [index]
        else:
            sample_groups_for_filter[group].append(index)
        
    samples_per_glycan = {}
    for group in sample_groups_for_filter:
        temp_samples_per_glycan = {}
        for i_i in sample_groups_for_filter[group]:
            checked_glycans = []
            for j_j, j in enumerate(df1_refactor[i_i]["Glycan"]):
                if j not in temp_samples_per_glycan.keys():
                    temp_samples_per_glycan[j] = 1
                    checked_glycans.append(j)
                elif j in temp_samples_per_glycan.keys() and j not in checked_glycans:
                    temp_samples_per_glycan[j] += 1
                    checked_glycans.append(j)
        samples_per_glycan[group] = temp_samples_per_glycan
    
    for group in sample_groups_for_filter:
        for i_i in sample_groups_for_filter[group]:
            to_remove = []  
            to_remove_glycan = []          
            for j_j, j in enumerate(df1_refactor[i_i]["Glycan"]):
                if samples_per_glycan[group][j] < round(len(sample_groups_for_filter[group])*(min_samples/100)):
                    to_remove.append(j_j)
                    to_remove_glycan.append(j)
            if len(to_remove) != 0:
                to_remove.reverse()
                to_remove_glycan.reverse()
                for j_j, j in enumerate(to_remove):
                    for k in df1_refactor[i_i]:
                        if k == 'Detected_Fragments':
                            continue
                        del df1_refactor[i_i][k][j]
                    if analyze_ms2:
                        for k in range(len(fragments_dataframes[i_i]["Glycan"])-1, -1, -1):
                            if fragments_dataframes[i_i]["Glycan"][k] == to_remove_glycan[j_j]:
                                del fragments_dataframes[i_i]["Glycan"][k]
                                del fragments_dataframes[i_i]["Adduct"][k]
                                del fragments_dataframes[i_i]["Fragment"][k]
                                del fragments_dataframes[i_i]["Fragment_mz"][k]
                                del fragments_dataframes[i_i]["Fragment_Intensity"][k]
                                del fragments_dataframes[i_i]["RT"][k]
                                del fragments_dataframes[i_i]["Precursor_mz"][k]
                                del fragments_dataframes[i_i]["% TIC explained"][k]
        
    return df1_refactor, fragments_dataframes, noise_levels
        
def make_filtered_ms2_refactor(df1_refactor,
                               fragments_dataframes,
                               reporter_ions,
                               unrestricted_fragments,
                               rt_tolerance_frag):
    '''Refactors and filters the fragments_dataframes. Filters by the reporter ions, if there are any set; Filters out peaks outside a peak retention time, if unrestricted_fragments is not used; and calculates the %TIC of a the annotated MS2 spectra.
    
    Parameters
    ----------
    df1_refactor : list
        A list containing dictionaries for each sample, which contains the organized data for that sample.
        
    fragments_dataframes : list
        Dataframe containing fragments information.
    
    reporter_ions : list
        A list containing the reporter ions you wish to filter your data with.

    unrestricted_fragments : boolean
        Whether or not fragments are restricted to the precursor putative composition.

    rt_tolerance_frag : float
        Same as rt_tolerance, but applies to identify whether the precursor retention time matches the fragmentation and ms1 identification.
        
    Returns
    -------
    df1_refactor : list
        A list containing dictionaries for each sample, which contains the organized data for that sample.
        
    fragments_dataframes : list
        Dataframe containing fragments information.
        
    fragments_refactor_dataframes : list
        A dataframe containing reorganized MS2 dataframes.
    '''
    if len(reporter_ions) != 0: #reporter_ions filtering
        for i_i, i in enumerate(fragments_dataframes):
            to_remove = []
            temp_list_fragments = []
            temp_list_mz = []
            current_checking = ""
            for j_j, j in enumerate(fragments_dataframes[i_i]["Glycan"]):
                to_check = j+"_"+fragments_dataframes[i_i]["Adduct"][j_j]+"_"+str(fragments_dataframes[i_i]["RT"][j_j])
                if to_check != current_checking:
                    found = False
                    current_checking = to_check
                    temp_list_fragments = []
                    temp_list_mz = []
                    for k in range(j_j, len(fragments_dataframes[i_i]["Glycan"])):
                        to_check_2 = fragments_dataframes[i_i]["Glycan"][k]+"_"+fragments_dataframes[i_i]["Adduct"][k]+"_"+str(fragments_dataframes[i_i]["RT"][k])
                        if to_check == to_check_2:
                            temp_list_fragments.append(fragments_dataframes[i_i]["Fragment"][k])
                            temp_list_mz.append(fragments_dataframes[i_i]["Fragment_mz"][k])
                        else:
                            break
                    for k_k, k in enumerate(reporter_ions):
                        try:
                            current_reporter = float(k)
                        except:
                            current_reporter = k
                        if type(current_reporter) == float:
                            for l in temp_list_mz:
                                if abs(current_reporter - l) <= 0.1:
                                    found = True
                                    break
                            if found:
                                break
                        if type(current_reporter) == str:
                            for l in temp_list_fragments:
                                if current_reporter == l.split("+")[0].split("-")[0].split("_")[0]:
                                    found = True
                                    break
                            if found:
                                break
                if not found and current_checking == to_check:
                    to_remove.append(j_j)
            if len(to_remove) != 0:
                for k_k in sorted(to_remove, reverse = True):
                    del fragments_dataframes[i_i]["Glycan"][k_k]
                    del fragments_dataframes[i_i]["Adduct"][k_k]
                    del fragments_dataframes[i_i]["Fragment"][k_k]
                    del fragments_dataframes[i_i]["Fragment_mz"][k_k]
                    del fragments_dataframes[i_i]["Fragment_Intensity"][k_k]
                    del fragments_dataframes[i_i]["RT"][k_k]
                    del fragments_dataframes[i_i]["Precursor_mz"][k_k]
                    del fragments_dataframes[i_i]["% TIC explained"][k_k] #end of reporter ions filtering
                    
    to_keep = [] #filters by retention time, if unrestricted fragments is not used
    for j_j, j in enumerate(df1_refactor):
        to_keep.append([])
        for k_k, k in enumerate(df1_refactor[j_j]["RT"]):
            if k_k == 0:
                df1_refactor[j_j]["%_TIC_explained"] = []
            df1_refactor[j_j]["%_TIC_explained"].append(None)
            found = False
            for l_l, l in enumerate(fragments_dataframes[j_j]["Glycan"]):
                if l == df1_refactor[j_j]["Glycan"][k_k] and fragments_dataframes[j_j]["Adduct"][l_l] == df1_refactor[j_j]["Adduct"][k_k]:
                    if abs(fragments_dataframes[j_j]["RT"][l_l] - k) <= rt_tolerance_frag:
                        found = True
                        to_keep[j_j].append(l_l)
            if found:
                df1_refactor[j_j]["Detected_Fragments"].append("Yes")
            else:
                df1_refactor[j_j]["Detected_Fragments"].append("No")
                
    if not unrestricted_fragments:  #Filters fragments with RT outside the detected peaks range
        for i_i, i in enumerate(to_keep):
            for k_k in range(len(fragments_dataframes[i_i]["Glycan"])-1, -1, -1):
                if k_k not in to_keep[i_i]:
                    del fragments_dataframes[i_i]["Glycan"][k_k]
                    del fragments_dataframes[i_i]["Adduct"][k_k]
                    del fragments_dataframes[i_i]["Fragment"][k_k]
                    del fragments_dataframes[i_i]["Fragment_mz"][k_k]
                    del fragments_dataframes[i_i]["Fragment_Intensity"][k_k]
                    del fragments_dataframes[i_i]["RT"][k_k]
                    del fragments_dataframes[i_i]["Precursor_mz"][k_k]
                    del fragments_dataframes[i_i]["% TIC explained"][k_k]
                    
    for i_i, i in enumerate(fragments_dataframes): #% TIC explained calculation
        fragments_int_sum = 0
        current_checking = ""
        for j_j, j in enumerate(fragments_dataframes[i_i]["Glycan"]):
            to_check = j+"_"+fragments_dataframes[i_i]["Adduct"][j_j]+"_"+str(fragments_dataframes[i_i]["RT"][j_j])
            if to_check != current_checking:
                current_checking = to_check
                fragments_int_sum = 0
                for k in range(j_j, len(fragments_dataframes[i_i]["Glycan"])):
                    to_check_2 = fragments_dataframes[i_i]["Glycan"][k]+"_"+fragments_dataframes[i_i]["Adduct"][k]+"_"+str(fragments_dataframes[i_i]["RT"][k])
                    if to_check == to_check_2:
                        fragments_int_sum += fragments_dataframes[i_i]["Fragment_Intensity"][k]
                    else:
                        break
            if current_checking == to_check and fragments_dataframes[i_i]["% TIC explained"] != 0:
                fragments_dataframes[i_i]["% TIC explained"][j_j] = float("%.2f" % round((fragments_int_sum/fragments_dataframes[i_i]["% TIC explained"][j_j])*100, 2)) #end of annotated_peaks ratio calculation
                
    for i_i, i in enumerate(df1_refactor): #start of ms2 score calculation (at the moment its just % TIC explained)
        for j_j, j in enumerate(i['Glycan']):
            if i['Detected_Fragments'][j_j] == 'Yes':
                glycan = j
                adduct = i['Adduct'][j_j]
                rt = i['RT'][j_j]
                list_tics_explained = []
                for k_k, k in enumerate(fragments_dataframes):
                    for l_l, l in enumerate(k['Glycan']):
                        if glycan == l and adduct == k['Adduct'][l_l] and abs(rt-k['RT'][l_l]) <= rt_tolerance_frag:
                            list_tics_explained.append(k['% TIC explained'][l_l])
                if len(list_tics_explained) > 0:
                    i['%_TIC_explained'][j_j] = max(list_tics_explained)
            else:
                i['%_TIC_explained'][j_j] = 0.0
                
    fragments_refactor_dataframes = [] #here it re-structures the fragments data to an updated format
    for i_i, i in enumerate(fragments_dataframes): #moving through samples
        fragments_refactor_dataframes.append({})
        current_scan = ""
        glycan_number = -1
        for j_j, j in enumerate(i['Glycan']):
            scan_name = j+"_"+i['Adduct'][j_j]+"_"+str(i['RT'][j_j])
            if scan_name != current_scan:
                glycan_number+=1
                current_scan = scan_name
                fragments_refactor_dataframes[i_i]['Glycan_'+str(glycan_number)+':'] = [j, 'RT_'+str(glycan_number)+':', i['RT'][j_j], 'Fragment:']
                fragments_refactor_dataframes[i_i]['Adduct_'+str(glycan_number)+':'] = [i['Adduct'][j_j], '% TIC assigned_'+str(glycan_number)+':', i['% TIC explained'][j_j], 'm/z:']
                fragments_refactor_dataframes[i_i]['m/z_'+str(glycan_number)+':'] = [i['Precursor_mz'][j_j], None, None, 'Intensity:']
                fragments_refactor_dataframes[i_i]['Glycan_'+str(glycan_number)+':'].append(i['Fragment'][j_j])
                fragments_refactor_dataframes[i_i]['Adduct_'+str(glycan_number)+':'].append(i['Fragment_mz'][j_j])
                fragments_refactor_dataframes[i_i]['m/z_'+str(glycan_number)+':'].append(i['Fragment_Intensity'][j_j])
            else:
                fragments_refactor_dataframes[i_i]['Glycan_'+str(glycan_number)+':'].append(i['Fragment'][j_j])
                fragments_refactor_dataframes[i_i]['Adduct_'+str(glycan_number)+':'].append(i['Fragment_mz'][j_j])
                fragments_refactor_dataframes[i_i]['m/z_'+str(glycan_number)+':'].append(i['Fragment_Intensity'][j_j])
    for i in fragments_refactor_dataframes: #makes all lists in the dataframe equal size so it can be ported to excel
        for j in i:
            while len(i[j]) < 1000:
                i[j].append(None)
    return df1_refactor, fragments_dataframes, fragments_refactor_dataframes
    
def make_total_dataframes(df1_refactor,
                          rt_tolerance):
    '''Makes and reorganizes the total_dataframes from df1_refactor data.
    
    Parameters
    ----------
    df1_refactor : list
        A list of dictionaries containing organized and filtered data.
        
    rt_tolerance : float
        Retention time tolerance used to consider two different peaks RT as distinct peaks.
        
    Returns
    -------
    total_dataframes : list
        A list of dictionaries containing the total AUC dataframes for each sample.
    
    '''
    total_dataframes = [] #total glycans AUC dataframe
    for i_i, i in enumerate(df1_refactor):
        ambiguities_solved = []
        total_dataframes.append({"Glycan": [], "RT": [], "AUC": []})
        current_glycan = ""
        RTs = []
        AUCs = []
        for j_j in range(len(i["Glycan"])):
            j = i["Glycan"][j_j]
            found = False
            if j == current_glycan:
                for k_k, k in enumerate(RTs):
                    if abs(i["RT"][j_j] - k) <= rt_tolerance:
                        RTs[k_k] = (k+i["RT"][j_j])/2
                        AUCs[k_k] = AUCs[k_k]+i["AUC"][j_j]
                        found = True
                        break
                if not found:
                    RTs.append(i["RT"][j_j])
                    AUCs.append(i["AUC"][j_j])
            if j != current_glycan:
                if j_j != 0:
                    if current_glycan not in ambiguities_solved:
                        glycan_name = current_glycan
                        if i["Ambiguity"][current_glycan_index] != "No":
                            ambiguities = i["Ambiguity"][current_glycan_index].split(", ")
                            for k in ambiguities:
                                current_ambiguity = k.split("_")[0]
                                if current_ambiguity not in ambiguities_solved:
                                    ambiguities_solved.append(current_ambiguity)
                                    glycan_name+=f"/{current_ambiguity}"
                        for k_k, k in enumerate(RTs):
                                total_dataframes[i_i]["Glycan"].append(glycan_name)
                                total_dataframes[i_i]["RT"].append(k)
                                total_dataframes[i_i]["AUC"].append(AUCs[k_k])
                    RTs = []
                    AUCs = []
                current_glycan = j
                current_glycan_index = j_j
                RTs.append(i["RT"][j_j])
                AUCs.append(i["AUC"][j_j])
            if j_j == len(i["Glycan"])-1:
                if current_glycan not in ambiguities_solved:
                    glycan_name = current_glycan
                    if i["Ambiguity"][current_glycan_index] != "No":
                        ambiguities = i["Ambiguity"][current_glycan_index].split(", ")
                        for k in ambiguities:
                            current_ambiguity = k.split("_")[0]
                            if current_ambiguity not in ambiguities_solved:
                                ambiguities_solved.append(current_ambiguity)
                                glycan_name+=f"/{current_ambiguity}"
                    for k_k, k in enumerate(RTs):
                        total_dataframes[i_i]["Glycan"].append(glycan_name)
                        total_dataframes[i_i]["RT"].append(k)
                        total_dataframes[i_i]["AUC"].append(AUCs[k_k])
                    RTs = []
                    AUCs = [] #total glycans AUC dataframe
    
    arranged_total_dataframes = []
    for i_i, i in enumerate(total_dataframes):
        arranged_total_dataframes.append({'Glycan' : [], 'RT' : [], 'AUC' : []})
        list_of_glycans = []
        glycans_dict_rt = {}
        glycans_dict_AUC = {}
        current_glycan = ''
        for j_j, j in enumerate(i['Glycan']):
            if j != current_glycan:
                current_glycan = j
                list_of_glycans.append(j)
                glycans_dict_rt[j] = []
                glycans_dict_AUC[j] = []
                glycans_dict_rt[j].append(i['RT'][j_j])
                glycans_dict_AUC[j].append(i['AUC'][j_j])
            else:
                glycans_dict_rt[j].append(i['RT'][j_j])
                glycans_dict_AUC[j].append(i['AUC'][j_j])
        if "Internal Standard" in list_of_glycans:
            highest = glycans_dict_AUC["Internal Standard"].index(max(glycans_dict_AUC["Internal Standard"]))
            to_remove = []
            for j_j, j in enumerate(glycans_dict_AUC["Internal Standard"]):
                if j_j != highest:
                    to_remove.append(j_j)
            for j in sorted(to_remove, reverse = True):
                del glycans_dict_rt["Internal Standard"][j]
                del glycans_dict_AUC["Internal Standard"][j]
        list_of_glycans = sorted(list_of_glycans)
        for j_j, j in enumerate(list_of_glycans):
            for k in range(len(glycans_dict_rt[j])):
                arranged_total_dataframes[i_i]['Glycan'].append(j)
            zipped = zip(glycans_dict_rt[j], glycans_dict_AUC[j])
            zipped = sorted(zipped)
            current_RTs, current_AUCs = zip(*zipped)
            arranged_total_dataframes[i_i]['RT'] += list(current_RTs)
            arranged_total_dataframes[i_i]['AUC'] += list(current_AUCs)
    total_dataframes = arranged_total_dataframes
    return total_dataframes
        
def determine_nglycan_class(total_dataframes,
                            compositions_dataframes):
    '''Determines the class of N-Glycans based on the number of H and N in the glycans composition.
    
    Parameters
    ----------
    total_dataframes : list
        A list of dictionaries containing the total AUC dataframes for each sample.
    
    compositions_dataframes : list
        A list containing dictionaries for each sample with combined peak AUC for each composition identified.
        
    Returns
    -------
    glycan_class : dictionary
        A dictionary containing the class of each glycan.
        
    proportion_classes : dictionaries
        A dictionary containing the proportion of each class per sample.
    '''
    glycan_class = {}
    for i_i, i in enumerate(total_dataframes):
        for j_j, j in enumerate(i['Glycan']):
            if j != 'Internal Standard' and j not in glycan_class.keys():
                j_list = j.split("/")
                temp_class = ''
                for k in j_list:
                    comp = General_Functions.form_to_comp_glycans(k)
                    if len(temp_class) > 0:
                        temp_class += '/'
                    if 'N' in comp.keys() and 'H' in comp.keys() and comp['N'] == 2 and comp['H'] <= 3:
                        temp_class += 'Paucimannose'
                        continue
                    if 'N' in comp.keys() and 'H' in comp.keys() and comp['H'] > comp['N']+1 and comp['N'] > 2:
                        temp_class += 'Hybrid'
                        continue
                    if 'N' in comp.keys() and 'H' in comp.keys() and comp['N'] == 2 and comp['H'] > 3:
                        temp_class += 'High-Mannose'
                        continue
                    else:
                        temp_class += 'Complex'
                glycan_class[j] = temp_class
            if j == 'Internal Standard':
                glycan_class[j] = j
    for i_i, i in enumerate(total_dataframes):
        total_dataframes[i_i]['Class'] = []
        for j_j, j in enumerate(i['Glycan']):
            i['Class'].append(glycan_class[j])
    for i_i, i in enumerate(compositions_dataframes):
        compositions_dataframes[i_i]['Class'] = []
        for j_j, j in enumerate(i['Glycan']):
            i['Class'].append(glycan_class[j])
    proportion_classes = {'Paucimannose' : [], 'Hybrid' : [], 'High-Mannose' : [], 'Complex' : []}        
    for i_i, i in enumerate(compositions_dataframes):
        total_sample = sum(i['AUC'])
        if total_sample == 0:
            total_sample = inf
        total_pauci = 0
        total_hybrid = 0
        total_oligo = 0
        total_complex = 0
        for j_j, j in enumerate(i['Class']):
            if 'Paucimannose' in j:
                total_pauci+= i['AUC'][j_j]
            if 'Hybrid' in j:
                total_hybrid+= i['AUC'][j_j]
            if 'High-Mannose' in j:
                total_oligo+= i['AUC'][j_j]
            if 'Complex' in j:
                total_complex+= i['AUC'][j_j]
        proportion_classes['Paucimannose'].append(float("%.2f" % round((total_pauci/total_sample)*100, 2)))
        proportion_classes['Hybrid'].append(float("%.2f" % round((total_hybrid/total_sample)*100, 2)))
        proportion_classes['High-Mannose'].append(float("%.2f" % round((total_oligo/total_sample)*100, 2)))
        proportion_classes['Complex'].append(float("%.2f" % round((total_complex/total_sample)*100, 2)))
        
    return glycan_class, proportion_classes
    
def calculate_fucosylation_sialylation(compositions_dataframes):
    '''Determines which glycans are fucosylated, sialylated or both and calculates their proportion by sample.
    
    Parameters
    ----------
    compositions_dataframes : list
        A list containing dictionaries for each sample with combined peak AUC for each composition identified.
        
    Returns
    -------
    glycans_fucsia : dictionary
        A dictionary containing whether each glycan is fucosylated, sialylated or both.
        
    proportion_fucsia : dictionaries
        A dictionary containing the proportion of fucosylated, sialylated and fucosylated+sialylated glycans per sample.
    '''
    glycans_fucsia = {}
    for i_i, i in enumerate(compositions_dataframes):
        for j_j, j in enumerate(i['Glycan']):
            fucsia = []
            if 'S' in j or 'Am' in j or 'E' in j or 'G' in j:
                fucsia.append(True)
            else:
                fucsia.append(False)
            if 'F' in j:
                fucsia.append(True)
            else:
                fucsia.append(False)
            glycans_fucsia[j] = fucsia
            
    proportion_fucsia = {'Fucosylated' : [], 'Sialylated' : [], 'Fuc+Sia' : []}       
    for i_i, i in enumerate(compositions_dataframes):
        total_sample = sum(i['AUC'])
        if total_sample == 0:
            total_sample = inf
        total_fuc = 0
        total_sia = 0
        total_fucsia = 0
        for j_j, j in enumerate(i['Glycan']):
            if glycans_fucsia[j][0] and glycans_fucsia[j][1]:
                total_fucsia+= i['AUC'][j_j]
            elif glycans_fucsia[j][0] and not glycans_fucsia[j][1]:
                total_sia+= i['AUC'][j_j]
            elif not glycans_fucsia[j][0] and glycans_fucsia[j][1]:
                total_fuc+= i['AUC'][j_j]
        proportion_fucsia['Fucosylated'].append(float("%.2f" % round((total_fuc/total_sample)*100, 2)))
        proportion_fucsia['Sialylated'].append(float("%.2f" % round((total_sia/total_sample)*100, 2)))
        proportion_fucsia['Fuc+Sia'].append(float("%.2f" % round((total_fucsia/total_sample)*100, 2)))
        
    return glycans_fucsia, proportion_fucsia
    
def create_metaboanalyst_files(plot_metaboanalyst,
                               df2,
                               total_dataframes,
                               all_glycans_list,
                               compositions,
                               compositions_dataframes,
                               save_path,
                               begin_time,
                               rt_tolerance,
                               from_GUI = False,
                               metab_groups = {},
                               fill_gaps_noise_level = False,
                               noise_levels = [],
                               glycans_mz = {}):
    '''Creates the metaboanalyst-compatible .csv files.
    
    Parameters
    ----------
    plot_metaboanalyst : list
        A list containing a boolean for whether or not to plot metaboanalyst data and a list of the groups.
    
    df2 : list
        A dictionary containing information about analyzed files.
        
    total_dataframes : list
        A list of dictionaries containing the total AUC dataframes for each sample.
        
    all_glycans_list : list
        A list containing all glycans found in all samples.
    
    compositions : boolean
        If set to True, also outputs the compositions analysis.
    
    compositions_dataframes : list
        A list containing dictionaries for each sample with combined peak AUC for each composition identified.
        
    save_path : string
        A string containing the path to the working directory of the script.
        
    begin_time : string
        Time at which the analysis began.
    
    rt_tolerance : float
        Tolerance of retention time (in minutes) at which an MS2 feature can be attributed to a specific retention time peak and also for peaks in different samples to be regarded as the same peak (and thus be compared with each other).
    
    from_GUI : boolean
        Whether or not the execution of this function came from GUI
                     
    metab_groups : list
        List of metaboanalyst groups to assign the samples to.
        
    Returns
    -------
    nothing
        Just creates the .csv files.
    '''
    if len(noise_levels) == 0:
        fill_gaps_noise_level = False
    local_noises_dict = {}
    with open(os.path.join(save_path, begin_time+"_glycan_abundance_table.csv"), "w") as f:
        # Make samples line
        samples_line = ["Sample"]
        for i_i, i in enumerate(df2["File_Name"]):
            samples_line.append(i)
        
        # Make groups line
        groups_line = ["Group"]
        if from_GUI:
            for i in df2['File_Name']:
                found = False
                for group, samples in metab_groups.items():
                    if i in samples:
                        found = True
                        groups_line.append(group)
                        break
                if not found:
                    groups_line.append("Ungrouped")
        else:
            if len(plot_metaboanalyst[1]) > 0:
                sample_group_dict = {}
                with open(plot_metaboanalyst[1], "r") as file:
                    for index, line in enumerate(file):
                        if index == 0:
                            continue
                        sample, group = line.strip().split(",")
                        sample_group_dict[sample] = group
                for i in df2["File_Name"]:
                    groups_line.append(sample_group_dict.get(i, "Ungrouped"))
            else:
                for i in df2["File_Name"]:
                    groups_line.append("Ungrouped")
        
        # Check if internal standard is present and pick the biggest area found within them for each sample
        found_int_std = False
        for i in total_dataframes:
            if "Internal Standard" in i["Glycan"]:
                found_int_std = True
                break
        if found_int_std:
            # Create the normalized data file if IS is present
            with open(os.path.join(save_path, begin_time+"_glycan_abundance_table_normalized.csv"), "a") as g:
                g.write(",".join(samples_line)+"\n")
                g.write(",".join(groups_line)+"\n")
                g.close()
            is_areas = []
            for i in total_dataframes:
                if "Internal Standard" in i["Glycan"]:
                    temp_areas = []
                    for j_j, j in enumerate(i["Glycan"]):
                        if j == "Internal Standard":
                            temp_areas.append(i["AUC"][j_j])
                    is_areas.append(max(temp_areas))
                else:
                    is_areas.append(0.0)
        
        # Write to peak-separated file
        f.write(",".join(samples_line)+"\n")
        f.write(",".join(groups_line)+"\n")
        
        for glycan_rt in all_glycans_list:
            glycan_line = []
            glycan_line_IS = []
            target_glycan, target_rt = glycan_rt.split("_")
            glycan_line_IS.append(glycan_rt)
            glycan_line.append(glycan_rt)
            for sample_index, sample in enumerate(total_dataframes): #moving through samples
                found = False
                temp_AUC = 0
                for glycan_index, glycan in enumerate(sample["Glycan"]):
                    if glycan == "Internal Standard":
                        continue
                    if glycan == target_glycan and abs(sample["RT"][glycan_index] - float(target_rt)) <= rt_tolerance:
                        found = True
                        if "Internal Standard" in sample["Glycan"]:
                            if is_areas[sample_index] > 0.0:
                                temp_AUC_IS = sample["AUC"][glycan_index]/is_areas[sample_index]
                            else:
                                temp_AUC_IS = 0.0
                            temp_AUC+= sample["AUC"][glycan_index]
                        else:
                            temp_AUC += sample["AUC"][glycan_index]
                if found:
                    if "Internal Standard" in sample["Glycan"]:
                        glycan_line_IS.append(str(temp_AUC_IS))
                    glycan_line.append(str(temp_AUC))
                    continue
                    
                if not found:
                    if fill_gaps_noise_level:
                        local_noise = []
                        
                        single_glycan = target_glycan.split("/")[0]
                        mz_values = glycans_mz[single_glycan]
                                
                        for noise_specs in noise_levels[1][sample_index]:
                            noises = []
                            for mz_value in mz_values:
                                noises.append(General_Functions.local_noise_calc(noise_specs, float(mz_value), noise_levels[0][sample_index]))
                            local_noise.append(sum(noises)/len(noises))
                        local_noise = sum(local_noise)/len(local_noise)
                        
                        if "Internal Standard" in sample["Glycan"]:
                            if is_areas[sample_index] > 0.0:
                                glycan_line_IS.append(str(local_noise/is_areas[sample_index]))
                        else:
                            glycan_line_IS.append("0.0")
                        glycan_line.append(str(local_noise))
                        if single_glycan in local_noises_dict.keys():
                            local_noises_dict[single_glycan].append(local_noise)
                        else:
                            local_noises_dict[single_glycan] = [local_noise]
                    else:
                        glycan_line_IS.append("")
                        glycan_line.append("")
                    continue
            if found_int_std:
                with open(os.path.join(save_path, begin_time+"_glycan_abundance_table_normalized.csv"), "a") as g:
                    g.write(",".join(glycan_line_IS)+"\n")
                    g.close()
            f.write(",".join(glycan_line)+"\n")
        f.close()
    
    # Make compositional metaboanalyst
    if compositions:
        total_glycans_compositions = []
        with open(os.path.join(save_path, begin_time+"_glycan_abundance_table_compositions.csv"), "w") as f:
            found_int_std = False
            for sample in compositions_dataframes:
                if "Internal Standard" in sample["Glycan"]:
                    found_int_std = True
                    break
            if found_int_std:
                with open(os.path.join(save_path, begin_time+"_glycan_abundance_table_compositions_normalized.csv"), "w") as g:
                    g.write(",".join(samples_line)+"\n")
                    g.write(",".join(groups_line)+"\n")
                    g.close()
            f.write(",".join(samples_line)+"\n")
            f.write(",".join(groups_line)+"\n")
            for sample_index, sample in enumerate(compositions_dataframes):
                for glycan_index, glycan in enumerate(sample['Glycan']):
                    if glycan not in total_glycans_compositions and glycan != 'Internal Standard':
                        total_glycans_compositions.append(glycan)
                        
            for glycan_index, glycan in enumerate(sorted(total_glycans_compositions)):
                glycan_line = [glycan]
                glycan_line_IS = [glycan]
                for sample_index, sample in enumerate(compositions_dataframes):
                    if glycan in sample['Glycan']:
                        glycan_line.append(str(sample['AUC'][sample['Glycan'].index(glycan)]))
                        if 'Internal Standard' in sample['Glycan']:
                            glycan_line_IS.append(str(sample['AUC'][sample['Glycan'].index(glycan)]/sample['AUC'][sample['Glycan'].index('Internal Standard')]))
                        else:
                            glycan_line_IS.append("0.0")
                    else:
                        if fill_gaps_noise_level:
                            single_glycan = glycan.split("/")[0]
                            glycan_line.append(str(sum(local_noises_dict[single_glycan])/len(local_noises_dict[single_glycan])))
                            if 'Internal Standard' in sample['Glycan']:
                                glycan_line_IS.append(str((sum(local_noises_dict[single_glycan])/len(local_noises_dict[single_glycan]))/sample['AUC'][sample['Glycan'].index('Internal Standard')]))
                            else:
                                glycan_line_IS.append("0.0")
                        else:
                            glycan_line.append('')
                            glycan_line_IS.append('')
                f.write(",".join(glycan_line)+"\n")
                if found_int_std:
                    with open(os.path.join(save_path, begin_time+"_glycan_abundance_table_compositions_normalized.csv"), "a") as g:
                        g.write(",".join(glycan_line_IS)+"\n")
                        g.close()
            f.close()
        
def output_filtered_data(args_dict):
    '''This function filters and converts raw results data into human readable
    excel files.
    
    Parameters
    ----------
    curve_fit_score : float
        The minimum curve fitting score to consider when outputting results.
        
    iso_fit_score : float
        The minimum isotopic fitting score to consider when outputting results.
        
    sn : int
        The minimum signal-to-noise ration to consider when outputting results.
        
    max_ppm : int
        The maximum amount of PPM difference to consider when outputting results.
        
    percentage_auc : float
        Percentage of the highest peak AUC that another peak AUC must have in a same chromatogram in order to be saved.
        
    reanalysis : tuple
        Contains two indexes: The first one indicates whether this is only a reanalysis execution and the second one indicates whether or not to produce plot data excel files. The second one is there because the plot data is more heavy and doesn't change on reanalysis, so should only be set to True if you lost original plot data.
        
    gg_file : string
        String containing the path to the .gg analysis file.
        
    save_path : string
        A string containing the path to the working directory of the script.
    
    analyze_ms2 : tuple
        A tuple with two indexes: The first one indicates whether to analyze ms2 data and the  second one indicates whether ms2 data should be forced to fit glycans composition.
        
    unrestricted_fragments : boolean
        Whether or not fragments are restricted to the precursor putative composition.
        
    reporter_ions : list
        A list containing the reporter ions you wish to filter your data with.
        
    plot_metaboanalyst : tuple
        A tuple with two indexes: The first one indicates whether or not to output a file to be used in metaboanalyst and the second one indicates the groups for the samples to be separated in.
        
    compositions : boolean
        If set to True, also outputs the compositions analysis.
        
    align_chromatograms : boolean
        Whether or not to align results and chromatograms.

    forced : string
        Indicates whether the function should force strict conditions based on the
        biological knowledge of glycans in order to avoid possible false positives when
        analysing N-glycans, O-glycans or GAGs.

    rt_tolerance : float
        Tolerance of retention time (in minutes) at which an MS2 feature can be attributed to a specific retention time peak and also for peaks in different samples to be regarded as the same peak (and thus be compared with each other).
    
    rt_tolerance_frag : float
        Same as rt_tolerance, but applies to identify whether the precursor retention time matches the fragmentation and ms1 identification.
        
    output_isotopic_fittings : boolean
        Whether or not to output curve-fitting and isotopic-fitting data for data check, if desired.
        
    output_plot_data : boolean
        Whether or not to output COMPLETE chromatogram plotting data (always outputs regular data of found glycans).
        
    multithreaded : boolean
        Whether the execution is multithreaded or not.
        
    number_cores : int
        Number of cores used if execution is multithreaded.
        
    temp_time : str
        Time at which execution of GG started.
        
    from_GUI : boolean
        Whether or not the execution of this function came from GUI
                         
    metab_groups : list
        List of metaboanalyst groups to assign the samples to.
        
    Uses
    ----
    datetime.datetime.now : Time object
        Returns the current date and time.
        
    dill.dump : None
        Pickle the current state of __main__ or another module to a file.
    
    dill.load : Module object
        Update the selected module (default is __main__) with the state saved at filename.
        
    pathlib.Path : Path object
        A subclass of PurePath, this class represents concrete paths of the system’s path flavour
        
    os.unlink : None
        Remove (delete) the file path. 
        
    pandas.DataFrame : Dataframe object
        Two-dimensional, size-mutable, potentially heterogeneous tabular data.
        
    pandas.ExcelWriter : ExcelWriter object
        Class for writing DataFrame objects into excel sheets.
        
    Returns
    -------
    nothing
        Creates excel files of processed data.
    '''
    # Grab all the args
    curve_fit_score = args_dict.get('curve fitting score threshold', 0.0)
    iso_fit_score = args_dict.get('isotopic fitting score threshold', 0.8)
    sn = args_dict.get('signal to noise ratio threshold', 3)
    max_ppm = args_dict.get('maximum ppm threshold', [-10, 10])
    percentage_auc = args_dict.get('minimum percentage auc threshold', 0.1)
    reanalysis = args_dict.get('reanalysis', False)
    gg_file = args_dict.get('reanalysis_path', None)
    save_path = args_dict.get('save path', None)
    analyze_ms2 = args_dict.get('analyze ms2', False)
    unrestricted_fragments = args_dict.get('unrestricted fragments', False)
    reporter_ions = args_dict.get('reporter ions', [])
    plot_metaboanalyst = args_dict.get('output abundance table', False)
    compositions = args_dict.get('output composition-separated data', True)
    align_chromatograms = args_dict.get('align chromatograms', False)
    forced = args_dict.get('glycan class', None)
    rt_tolerance = args_dict.get('retention time tolerance', 0.2)
    rt_tolerance_frag = args_dict.get('fragmentation retention time tolerance', 1)
    output_isotopic_fittings = args_dict.get('output isotopic fitting data', False)
    output_plot_data = args_dict.get('output chromatogram plotting data', False)
    multithreaded = args_dict.get('multithreaded', True)
    number_cores = args_dict.get('number of cpu cores', 99)
    temp_time = args_dict.get('analysis done time', datetime.datetime.now())
    min_samples = args_dict.get('minimum percentage of samples', 0)
    temp_folder = args_dict.get('temporary folder', None)
    from_GUI = args_dict.get('from GUI', False)
    metab_groups = args_dict.get('sample groups', {})
    fill_gaps = args_dict.get('fill data gaps', (False, 50, 0.2, False))
    
    # Preparation for arranging the data
    date = datetime.datetime.now() #gets date information
    begin_time = str(date)[2:4]+str(date)[5:7]+str(date)[8:10]+"_"+str(date)[11:13]+str(date)[14:16]+str(date)[17:19] #arranges the date information for the filename
    
    # Checks whether it's a reanalysis or not
    if reanalysis:
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        if from_GUI:
            print(time_formatted+"Saving results files...")
        else:
            print(time_formatted+"Reanalyzing raw data...")
        try:
            General_Functions.open_gg(gg_file, temp_folder)
        except:
            print("\nAnalysis file not found. If you're reanalyzing\nan existing analysis result file, check if the\npath in 'analysis_file' is correct, otherwise\nchange 'reanalysis' to 'no' and try again.\n")
            return
    else: #if its not, it opens the already loaded raw data files
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        print(time_formatted+"Analyzing raw data...")
        
    with open(os.path.join(temp_folder, 'results'), 'rb') as f:
        file = dill.load(f)
        df1 = file[0]
        df2 = file[1]
        if len(file) > 3: #checks whether MS2 data is present
            analyze_ms2 = True
            fragments_dataframes = file[2]
            if reanalysis and ".".join(version.split('.')[:2]) != ".".join(file[3].split('.')[:2]): #checks if version is compatible by checking the second version number
                print("Raw data files version incompatible with\ncurrent version (Current version: "+version+";\nRaw data version: "+file[3]+")")
                return
        else:
            fragments_dataframes = []
            if reanalysis and ".".join(version.split('.')[:2]) != ".".join(file[2].split('.')[:2]):
                print("Raw data files version incompatible with\ncurrent version (Current version: "+version+";\nRaw data version: "+file[2]+")")
                return
        f.close()
    
    if 'noise_levels' in os.listdir(temp_folder):
        with open(os.path.join(temp_folder, 'noise_levels'), 'rb') as f:
            noise_levels = dill.load(f)
    else:
        noise_levels = []
    #reorganizes and filters the raw data based on the quality thresholds
    if from_GUI:
        df1_refactor, fragments_dataframes, noise_levels = make_df1_refactor(df1, df2, curve_fit_score, iso_fit_score, sn, max_ppm, percentage_auc, analyze_ms2, unrestricted_fragments, min_samples, fragments_dataframes, fill_gaps, metab_groups, noise_levels)
    else:
        if len(plot_metaboanalyst[1]) != 0:
            df1_refactor, fragments_dataframes, noise_levels = make_df1_refactor(df1, df2, curve_fit_score, iso_fit_score, sn, max_ppm, percentage_auc, analyze_ms2, unrestricted_fragments, min_samples, fragments_dataframes, fill_gaps, plot_metaboanalyst[1], noise_levels)
        else:
            df1_refactor, fragments_dataframes, noise_levels = make_df1_refactor(df1, df2, curve_fit_score, iso_fit_score, sn, max_ppm, percentage_auc, analyze_ms2, unrestricted_fragments, min_samples, fragments_dataframes, fill_gaps, {}, noise_levels)
    
    #filters ms2 data by reporter ions, calculates %TIC of MS2 spectra and reorganizes MS2 data
    if analyze_ms2:
        df1_refactor, fragments_dataframes, fragments_refactor_dataframes = make_filtered_ms2_refactor(df1_refactor, fragments_dataframes, reporter_ions, unrestricted_fragments, rt_tolerance_frag)
    
    #checks the ambiguities
    for i_i, i in enumerate(df1_refactor):
        i['Ambiguity'] = []
        for j in i['Glycan']:
            i['Ambiguity'].append([])
        for j_j, j in enumerate(i['Glycan']):
            glycan_j = j+'_'+i['Adduct'][j_j]
            for k_k, k in enumerate(i['Glycan'][j_j+1:]):
                k_k = j_j+k_k+1
                glycan_k = k+'_'+i['Adduct'][k_k]
                if j != k and i['mz'][j_j] == i['mz'][k_k]:
                    i['Ambiguity'][j_j].append(i['Glycan'][k_k]+'_'+i['Adduct'][k_k])
                    i['Ambiguity'][k_k].append(i['Glycan'][j_j]+'_'+i['Adduct'][j_j])
        for j_j, j in enumerate(i['Ambiguity']):
            if len(j) > 0:
                i['Ambiguity'][j_j] = ', '.join(j)
            else:
                i['Ambiguity'][j_j] = 'No'
    
    #counts ambiguities without counting twice
    ambiguity_count = []
    for i_i, i in enumerate(df1_refactor):
        noted_ambiguities = []
        ambiguity_count.append(0)
        for j_j, j in enumerate(i['Glycan']):
            if j not in noted_ambiguities:
                if i['Ambiguity'][j_j] != 'No':
                    noted_ambiguities.append(j)
                    ambiguity_count[-1] += 1
                    for k in i['Ambiguity'][j_j].split(", "):
                        noted_ambiguities.append(k.split("_")[0])
    
    #makes dataframes containing the combined adducts AUC for each peak
    total_dataframes = make_total_dataframes(df1_refactor, rt_tolerance)
    
    compositions_dataframes = [] #compositions_dataframes
    for i_i, i in enumerate(total_dataframes):
        compositions_dataframes.append({'Glycan' : [], 'AUC' : []})
        glycans_lib = {}
        for j_j, j in enumerate(i['Glycan']):
            if j not in glycans_lib.keys():
                glycans_lib[j] = i['AUC'][j_j]
            else:
                glycans_lib[j] += i['AUC'][j_j]
        for j_j, j in enumerate(glycans_lib):
            compositions_dataframes[i_i]['Glycan'].append(j)
            compositions_dataframes[i_i]['AUC'].append(glycans_lib[j])
    
    if forced == 'n_glycans': #if N-Glycans, determines its class
        glycan_class, proportion_classes = determine_nglycan_class(total_dataframes, compositions_dataframes)
    
    #hook for alignment tool. it'll use the total_dataframes (total_glycans) 
    if align_chromatograms:
        if len(total_dataframes) > 1:
            good_alignment = False
            for i in total_dataframes:
                if len(i['Glycan']) > 0:
                    good_alignment = True
                    break
            if not good_alignment:
                print("No good glycans to align the data with.")
            else:
                aligned_total_glycans = align_assignments(total_dataframes, 'total_glycans', multithreaded, number_cores, rt_tol = rt_tolerance)
                total_dataframes = aligned_total_glycans[0]
                df2["Average Delta t"] = []
                for i_i, i in enumerate(aligned_total_glycans[1]):
                    temp_list = []
                    for k_k, k in enumerate(i):
                        temp_list.append(i[k][0])
                    if len(temp_list) != 0:
                        average_delta = float("%.2f" % round(sum(temp_list)/len(temp_list), 2))
                        df2["Average Delta t"].append(average_delta)
                    else:
                        if i_i == aligned_total_glycans[2]:
                            df2["Average Delta t"].append("Reference Sample")
                        else:
                            df2["Average Delta t"].append("No peaks to align")
    #hook for alignment tool. it'll use the total_dataframes (total_glycans)         
                
    glycans_count = [] #glycans composition counts
    for i_i, i in enumerate(df1_refactor):
        current_glycan = ""
        glycans_count.append(0)
        for j_j, j in enumerate(i["Glycan"]):
            if j == 'Internal Standard':
                continue
            if j != current_glycan:
                glycans_count[i_i]+= 1
                current_glycan = j
    df2["MS1_Glycans_Compositions"] = glycans_count
    if analyze_ms2:
        glycans_count = []
        for i_i, i in enumerate(fragments_dataframes):
            current_glycan = ""
            glycans_count.append(0)
            for j_j, j in enumerate(i["Glycan"]):
                if j != current_glycan:
                    glycans_count[i_i]+= 1
                    current_glycan = j
        df2["MS2_Glycans_Compositions"] = glycans_count #end of glycans counts
        
    df2["Ambiguities"] = ambiguity_count
    
    glycans_fucsia, proportion_fucsia = calculate_fucosylation_sialylation(compositions_dataframes)
    
    df2["Fucosylated %"] = proportion_fucsia["Fucosylated"]
    df2["Sialylated %"] = proportion_fucsia["Sialylated"]
    df2["Fuc+Sia %"] = proportion_fucsia["Fuc+Sia"]
    
    if forced == 'n_glycans':
        df2["Paucimannose %"] = proportion_classes["Paucimannose"]
        df2["Hybrid %"] = proportion_classes["Hybrid"]
        df2["High-Mannose %"] = proportion_classes["High-Mannose"]
        df2["Complex %"] = proportion_classes["Complex"]
        
    all_glycans_list = [] #here it makes a list of ALL the glycans found for use in other parts of the data arrangement workflow
    for sample in total_dataframes:
        for rt_index, rt in enumerate(sample["RT"]):
            if sample["Glycan"][rt_index] == "Internal Standard":
                continue
            found = False
            if len(all_glycans_list) > 0:
                for glycan_rt in all_glycans_list:
                    target_glycan, target_rt = glycan_rt.split("_")
                    if sample["Glycan"][rt_index] == target_glycan and abs(rt-float(target_rt)) <= rt_tolerance:
                        found = True
            if found:
                continue
            all_glycans_list.append(f"{sample['Glycan'][rt_index]}_{rt:.2f}")
                
    if plot_metaboanalyst[0]: #start of metaboanalyst plot
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        print(time_formatted+"Creating glycan abundance table...", end="", flush=True)
        
        glycans_mzs = {}
        if fill_gaps[3]:
            for sample in df1_refactor:
                for glycan, mz in (zip(sample['Glycan'], sample['mz'])):
                    if glycan not in glycans_mzs:
                        glycans_mzs[glycan] = [mz]
                    else:
                        if mz not in glycans_mzs[glycan]:
                            glycans_mzs[glycan].append(mz)
        
        create_metaboanalyst_files(plot_metaboanalyst, df2, total_dataframes, all_glycans_list, compositions, compositions_dataframes, save_path, begin_time, rt_tolerance, from_GUI, metab_groups, fill_gaps[3], noise_levels, glycans_mzs)
        
        print("Done!") #end of metaboanalyst plot
    
    total_glycans_compositions = []  #start to build meta_dataframe
    for i_i, i in enumerate(compositions_dataframes): #get all possible compositions
        for j_j, j in enumerate(i['Glycan']):
            if j not in total_glycans_compositions and j != 'Internal Standard':
                total_glycans_compositions.append(j)
    adducts = []
    for i_i, i in enumerate(df1_refactor): #get all possible adducts
        for j_j, j in enumerate(i['Adduct']):
            if j not in adducts and i['Glycan'][j_j] != 'Internal Standard':
                adducts.append(j)
                
    meta_dataframe = {'No.' : [], 'Composition' : []}
    if forced == 'n_glycans':
        meta_dataframe['Class'] = []
    for i_i, i in enumerate(sorted(adducts)):
        meta_dataframe['[M+'+i+']'] = []
        meta_dataframe['Avg PPM Error - '+i] = []
        meta_dataframe['In Samples - '+i] = []
    meta_dataframe['No. Samples'] = []
    for i_i, i in enumerate(sorted(total_glycans_compositions)):
        raw_i = i
        i = i.split("/")[0]
        total_replicates = []
        if forced == 'n_glycans':
            meta_dataframe['Class'].append(glycan_class[raw_i])
        meta_dataframe['No.'].append(i_i+1)
        meta_dataframe['Composition'].append(raw_i)
        for j_j, j in enumerate(sorted(adducts)):
            meta_dataframe['[M+'+j+']'].append("-")
            current_adduct_PPM_Error = []
            replicates_present = []
            for k_k, k in enumerate(df1_refactor):
                found_replicate = False
                ppm_error_data = []
                for l_l, l in enumerate(k['Glycan']):
                    if l == i and k['Adduct'][l_l] == j:
                        if meta_dataframe['[M+'+j+']'][i_i] == "-":
                            meta_dataframe['[M+'+j+']'][i_i] = k['mz'][l_l]
                        found_replicate = True
                        ppm_error_data.append(k['PPM'][l_l])
                if found_replicate:
                    replicates_present.append(k_k)
                    if k_k not in total_replicates:
                        total_replicates.append(k_k)
                if len(ppm_error_data) != 0:
                    current_adduct_PPM_Error.append(sum(ppm_error_data)/len(ppm_error_data))
            if len(current_adduct_PPM_Error) != 0:
                meta_dataframe['Avg PPM Error - '+j].append(float("%.3f" % round(sum(current_adduct_PPM_Error)/len(current_adduct_PPM_Error), 3)))
            else:
                meta_dataframe['Avg PPM Error - '+j].append("-")
            if len(replicates_present) != 0:
                meta_dataframe['In Samples - '+j].append(str(replicates_present)[1:-1])
            else:
                meta_dataframe['In Samples - '+j].append("-")
        meta_dataframe['No. Samples'].append(len(total_replicates))  #end of meta_dataframe building
    
    #start of excel data printing
    df2 = DataFrame(df2) 
    meta_df = DataFrame(meta_dataframe)
    if type(max_ppm) == float:
        ppm_title_label = str(max_ppm)
    else:
        ppm_title_label = str(max_ppm[0])+"-"+str(max_ppm[1])
    with ExcelWriter(os.path.join(save_path, begin_time+'_Results_'+ppm_title_label+'_'+str(iso_fit_score)+'_'+str(curve_fit_score)+'_'+str(sn)+'.xlsx')) as writer:
        dfs = [df2, meta_df]
        sheets_names = ['Index references', 'Detected Glycans']
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        print(time_formatted+"Creating results file...", end="", flush=True)
        df2.to_excel(writer, sheet_name="Index references", index = False)
        meta_df.to_excel(writer, sheet_name="Detected Glycans", index = False)
        for i_i, i in enumerate(df1_refactor):
            result_df = DataFrame(i)
            result_df.to_excel(writer, sheet_name="Sample_"+str(i_i), index = False)
            total_aucs_df = DataFrame(total_dataframes[i_i])
            total_aucs_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Total_AUCs", index = False)
            dfs.append(result_df)
            dfs.append(total_aucs_df)
            sheets_names.append("Sample_"+str(i_i))
            sheets_names.append("Sample_"+str(i_i)+"_Total_AUCs")
            if compositions:
                compositions_df = DataFrame(compositions_dataframes[i_i])
                compositions_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Compositions_AUCs", index = False)
                dfs.append(compositions_df)
                sheets_names.append("Sample_"+str(i_i)+"_Compositions_AUCs")
            if analyze_ms2:
                if len(fragments_dataframes[i_i]["Glycan"]) > 0:
                    fragments_df = DataFrame(fragments_refactor_dataframes[i_i])
                    fragments_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Fragments", index = False)
                    dfs.append(fragments_df)
                    sheets_names.append("Sample_"+str(i_i)+"_Fragments")
        for index, sheet in enumerate(sheets_names):
            General_Functions.autofit_columns_excel(dfs[index], writer.sheets[sheet])
    del df1
    del result_df
    del total_dataframes
    del total_aucs_df
    del dfs
    del sheets_names
    if compositions:
        del compositions_dataframes
        del compositions_df
    
    if analyze_ms2:
        if len(fragments_dataframes[i_i]["Glycan"]) > 0:
            del fragments_df
        del fragments_refactor_dataframes
        del fragments_dataframes
    print("Done!")
         
    samples_aligned = False
    if align_chromatograms:
        if len(df2['Sample_Number']) > 1 and good_alignment: #aligns the chromatograms, may take some time (around 1 minute per sample, depending on run length)
                
            time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
            print(time_formatted+"Aligning chromatograms...", end='', flush=True)
            align_assignments(df2, 'chromatograms', multithreaded, number_cores, temp_folder, gg_file, aligned_total_glycans[1], None, iso_fit_score, curve_fit_score, max_ppm, sn)
            samples_aligned = True
            print("Done!")
        
    if not from_GUI:
        # Select found glycans EICs
        found_eic_processed_dataframes = []
        for i_i, i in enumerate(df1_refactor):
            found_eic_processed_dataframes.append({})
            
            eic_name = 'RTs'
            if samples_aligned:
                with open(os.path.join(temp_folder, f"{i_i}_aligned_{eic_name}_{iso_fit_score}_{curve_fit_score}_{max_ppm}_{sn}"), "rb") as f:
                    found_eic_processed_dataframes[i_i]['RTs_'+str(i_i)] = dill.load(f)
                    f.close()
            else:
                found_eic_processed_dataframes[i_i]['RTs_'+str(i_i)] = General_Functions.access_chromatogram(i_i, f"{i_i}_{eic_name}", temp_folder, gg_file)
            
            for j_j, j in enumerate(i['Glycan']):
                query = j+"+"+i['Adduct'][j_j]+" - "+str(i['mz'][j_j])
                try:
                    found_eic_processed_dataframes[i_i][query] = General_Functions.access_chromatogram(i_i, f"{i_i}_smoothed_{query}", temp_folder, gg_file)
                except:
                    pass
        
        # Combines adducts of found EICs (deconvolution)
        print(time_formatted+"Creating data plotting files...", end='', flush=True)
        found_eic_processed_dataframes_simplified = []
        found_eic_processed_dataframes_copy = copy.deepcopy(found_eic_processed_dataframes)
        for i_i, i in enumerate(found_eic_processed_dataframes_copy):
            current_glycan = ""
            found_eic_processed_dataframes_simplified.append({})
            for j_j, j in enumerate(i):
                working_glycan = j.split("+")[0].split("_")[0].split("-")[0]
                if j_j == 0:
                    found_eic_processed_dataframes_simplified[i_i][j] = i[j]
                    continue
                elif working_glycan != current_glycan:
                    current_glycan = working_glycan
                    found_eic_processed_dataframes_simplified[i_i][working_glycan] = i[j]
                else:
                    for k_k, k in enumerate(i[j]):
                        found_eic_processed_dataframes_simplified[i_i][working_glycan][k_k] += k
        del found_eic_processed_dataframes_copy
        
        # Print found EICs to excel files
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        with ExcelWriter(os.path.join(save_path, begin_time+'_Found_Glycans_EICs.xlsx')) as writer:
            df2.to_excel(writer, sheet_name="Index references", index = False)
            General_Functions.autofit_columns_excel(df2, writer.sheets["Index references"])
            for i_i, i in enumerate(found_eic_processed_dataframes_simplified):
                found_eic_processed_dataframes_simplified_df = DataFrame(found_eic_processed_dataframes_simplified[i_i])
                found_eic_processed_dataframes_simplified_df.to_excel(writer, sheet_name="Processed_Sample_"+str(i_i), index = False)
        del found_eic_processed_dataframes_simplified
        del found_eic_processed_dataframes_simplified_df
        
        # Plot isotopic fittings and curve fittings on demand
        if output_isotopic_fittings:
            with open(os.path.join(temp_folder, 'isotopic_fittings'), 'rb') as f: #start of isotopic fits output
                isotopic_fits_dataframes = dill.load(f)
                f.close()
                
            isotopic_fits_dataframes_arranged = []
            biggest_len = 0
            results = []
            if multithreaded:
                if number_cores == 'all':
                    cpu_count = (os.cpu_count())-2
                    if cpu_count <= 0:
                        cpu_count = 1
                else:
                    number_cores = int(number_cores)
                    if number_cores > (os.cpu_count())-2:
                        cpu_count = (os.cpu_count())-2
                        if cpu_count <= 0:
                            cpu_count = 1
                    else:
                        cpu_count = number_cores
            else:
                cpu_count = 1
                
            with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_count if cpu_count < 60 else 60) as executor:
                for i_i, i in enumerate(isotopic_fits_dataframes): #sample
                    result = executor.submit(arrange_iso_outputs, i_i, i, isotopic_fits_dataframes)
                    results.append(result)
                    isotopic_fits_dataframes_arranged.append(None)
                for index, i in enumerate(results):
                    current_result = i.result()
                    isotopic_fits_dataframes_arranged[current_result[2]] = current_result[0]
                    if current_result[1] > biggest_len:
                        biggest_len = current_result[1]
                    results[index] = None
            del isotopic_fits_dataframes
            
            with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_count if cpu_count < 60 else 60) as executor:
                for i_i, i in enumerate(isotopic_fits_dataframes_arranged):
                    executor.submit(write_iso_to_excel, save_path, begin_time, i_i, i, isotopic_fits_dataframes_arranged, biggest_len)
            del isotopic_fits_dataframes_arranged
            
            with open(os.path.join(temp_folder, 'curve_fittings'), 'rb') as f:
                curve_fitting_dataframes = dill.load(f)
                f.close()
            biggest_len = 0
            for i in curve_fitting_dataframes: #finds out the biggest len
                for j in i:
                    if len(i[j]) > biggest_len:
                        biggest_len = len(i[j])
            for i in curve_fitting_dataframes: #elongates the smaller dataframes so that they are all the same size
                for j in i:
                    if len(i[j]) < biggest_len:
                        while len(i[j]) < biggest_len:
                            i[j].append(None)
                            
            with ExcelWriter(os.path.join(save_path, begin_time+'_curve_fitting_Plot_Data.xlsx')) as writer:
                df2.to_excel(writer, sheet_name="Index references", index = False)
                General_Functions.autofit_columns_excel(df2, writer.sheets["Index references"])
                for i_i, i in enumerate(curve_fitting_dataframes):
                    if len(curve_fitting_dataframes[i_i]) > 16384:
                        for j in range(int(len(curve_fitting_dataframes[i_i])/16384)+1):
                            if j == 0:
                                curve_df = DataFrame(dict(islice(curve_fitting_dataframes[i_i].items(), 16384)))
                                curve_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Curve_Fits_0", index = False)
                            else:
                                if len(dict(islice(curve_fitting_dataframes[i_i].items(), j*16384, len(curve_fitting_dataframes[i_i])))) <= 16384:
                                    curve_df = DataFrame(dict(islice(curve_fitting_dataframes[i_i].items(), j*16384, len(curve_fitting_dataframes[i_i]))))
                                    curve_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Curve_Fits_"+str(j), index = False)
                                else:
                                    curve_df = DataFrame(dict(islice(curve_fitting_dataframes[i_i].items(), j*16384, (j+1)*16384)))
                                    curve_df.to_excel(writer, sheet_name="Sample_"+str(i_i)+"_Curve_Fits_"+str(j), index = False)
                    else:
                        curve_df = DataFrame(curve_fitting_dataframes[i_i])
                        curve_df.to_excel(writer, sheet_name="Sample_"+str(i_i), index = False)
                        
            del curve_fitting_dataframes
            del curve_df        
        
        if not output_plot_data:
            print("Done!")
        
    # Prints EIC of all glycans
    if output_plot_data:
        if from_GUI:
            print(time_formatted+"Creating data plotting files...", end='', flush=True)
        with open(os.path.join(temp_folder, "eics_list"), "rb") as f:
            eics = dill.load(f)
            f.close()
    
        with ExcelWriter(os.path.join(save_path, begin_time+'_processed_EIC_Plot_Data.xlsx')) as writer: #smoothed eic, now changed to processed to avoid TMI
            df2.to_excel(writer, sheet_name="Index references", index = False)
            General_Functions.autofit_columns_excel(df2, writer.sheets["Index references"])
            
            for i_i in eics:
                smoothed_eic_dataframes = {}
                
                eic_name = "RTs"
                if samples_aligned:
                    rts_name = f"{i_i}_aligned_{eic_name}_{iso_fit_score}_{curve_fit_score}_{max_ppm}_{sn}"
                else:
                    rts_name = f"{i_i}_{eic_name}"
                    
                smoothed_eic_dataframes[f"RTs_{i_i}"] = General_Functions.access_chromatogram(i_i, rts_name, temp_folder, gg_file)
                
                for j_j, j in enumerate(eics[i_i]):
                    if j_j == 0:
                        continue
                    file_name = str(i_i)+"_smoothed_"+j
                    
                    smoothed_eic_dataframes[j] = General_Functions.access_chromatogram(i_i, file_name, temp_folder, gg_file)
                    
                smoothed_eic_df = DataFrame(smoothed_eic_dataframes)
                smoothed_eic_df.to_excel(writer, sheet_name="Sample_"+str(i_i), index = False)
        del smoothed_eic_dataframes
        del smoothed_eic_df
        
        with ExcelWriter(os.path.join(save_path, begin_time+'_raw_EIC_Plot_Data.xlsx')) as writer:
            df2.to_excel(writer, sheet_name="Index references", index = False)
            General_Functions.autofit_columns_excel(df2, writer.sheets["Index references"])
                
            for i_i in eics:
                raw_eic_dataframes = {}
                rts_name = str(i_i)+"_RTs"
                
                raw_eic_dataframes[f"RTs_{i_i}"] = General_Functions.access_chromatogram(i_i, rts_name, temp_folder, gg_file)
                
                for j_j, j in enumerate(eics[i_i]):
                    if j_j == 0:
                        continue
                    file_name = str(i_i)+"_raw_"+j
                    
                    raw_eic_dataframes[j] = General_Functions.access_chromatogram(i_i, file_name, temp_folder, gg_file)
                        
                raw_eic_df = DataFrame(raw_eic_dataframes)
                raw_eic_df.to_excel(writer, sheet_name="Sample_"+str(i_i), index = False)
        del raw_eic_dataframes
        del raw_eic_df
        del eics
        
        print("Done!")
    
def arrange_iso_outputs(i_i,
                        i,
                        isotopic_fits_dataframes):
    '''Function to organize the isotopic fits output of a single sample. To be used in multithreading for faster execution.
    
    Parameters
    ----------
    i : str
        Sample data.
    
    i_i : int
        Sample index.
        
    isotopic_fits_dataframes : list
        A list containing dictionaries with information of all the isotopic fits calculated for all the samples.
        
    Returns
    -------
    temp_fits_dataframes : dict
        Dictionary containing all the isotopic fittings organized.
        
    biggest_len : int
        The biggest len of a given column in the temp_fits_dataframes.
        
    i_i : int
        Sample index.
    '''
    try:
        biggest_len = 0
        temp_fits_dataframes = {}
        for j_j, j in enumerate(i): #glycan
            temp_fits_dataframes[j] = {}
            for k_k, k in enumerate(i[j]): #peaks of glycan
                temp_fits_dataframes[j]['RT_'+str(k_k)+':'] = []
                temp_fits_dataframes[j]['Score_'+str(k_k)+':'] = []
                temp_fits_dataframes[j]['fit_'+str(k_k)] = []
                temp_fits_dataframes[j]['RT_'+str(k_k)+':'].append(k)
                temp_fits_dataframes[j]['Score_'+str(k_k)+':'].append(isotopic_fits_dataframes[i_i][j][k][3])
                temp_fits_dataframes[j]['fit_'+str(k_k)].append(None)
                temp_fits_dataframes[j]['RT_'+str(k_k)+':'].append('mz:')
                temp_fits_dataframes[j]['Score_'+str(k_k)+':'].append('Ideal:')
                temp_fits_dataframes[j]['fit_'+str(k_k)].append('Actual:')
                temp_fits_dataframes[j]['RT_'+str(k_k)+':'] = temp_fits_dataframes[j]['RT_'+str(k_k)+':']+isotopic_fits_dataframes[i_i][j][k][0]
                temp_fits_dataframes[j]['Score_'+str(k_k)+':'] = temp_fits_dataframes[j]['Score_'+str(k_k)+':']+isotopic_fits_dataframes[i_i][j][k][1]
                temp_fits_dataframes[j]['fit_'+str(k_k)] = temp_fits_dataframes[j]['fit_'+str(k_k)]+isotopic_fits_dataframes[i_i][j][k][2]
                if len(temp_fits_dataframes[j]['RT_'+str(k_k)+':']) > biggest_len:
                    biggest_len = len(temp_fits_dataframes[j]['RT_'+str(k_k)+':'])
        return temp_fits_dataframes, biggest_len, i_i
    except KeyboardInterrupt:
        print("\n\n----------Execution cancelled by user.----------\n", flush=True)
        raise SystemExit(1)
    
def write_iso_to_excel(save_path,
                       begin_time,
                       i_i,
                       i,
                       isotopic_fits_dataframes_arranged,
                       biggest_len):
    '''Writes the organized isotopic fittings to an excel file.
    
    Parameters
    ----------
    save_path : string
        A string containing the path to the working directory of the script.
    
    begin_time : string
        Time at which the file output started, to save to file name.
        
    i_i : int
        Sample index.
        
    i : list
        Sample data.
        
    isotopic_fits_dataframes_arranged : dict
        The dictionary containing the organized isotopic fittings data.
        
    biggest_len : int
        The biggest length of all the columns in the isotopic_fits_dataframes_arranged.
        
    Uses
    ----
    pandas.DataFrame : Dataframe object
        Two-dimensional, size-mutable, potentially heterogeneous tabular data.
        
    pandas.ExcelWriter : ExcelWriter object
        Class for writing DataFrame objects into excel sheets.
    
    Returns
    -------
    nothing
        Creates excel files with the data.
    '''
    try:
        with ExcelWriter(os.path.join(save_path, begin_time+'_Isotopic_Fits_Sample_'+str(i_i)+'.xlsx')) as writer:
            for j_j, j in enumerate(i): #navigating glycans
                for k_k, k in enumerate(list(i[j].keys())):
                    while len(i[j][k]) < biggest_len:
                        i[j][k].append(None)
                isotopic_fits_df = DataFrame(isotopic_fits_dataframes_arranged[i_i][j])
                isotopic_fits_df.to_excel(writer, sheet_name=j, index = False)
        del isotopic_fits_df
    except KeyboardInterrupt:
        print("\n\n----------Execution cancelled by user.----------\n", flush=True)
        raise SystemExit(1)
        
def combine_raw_data_per_sample(sample_no, temp_folder):
    '''
    '''
    files_list = os.listdir(temp_folder)
    with zipfile.ZipFile(os.path.join(temp_folder, f"{sample_no}_eics"), 'w', compression=zipfile.ZIP_DEFLATED) as zipf:
        for file in files_list:
            if file.split("_")[0] == f"{sample_no}":
                zipf.write(os.path.join(temp_folder, file), arcname=file)
                os.remove(os.path.join(temp_folder, file))
        zipf.close()

def arrange_raw_data(args_dict):
    '''Arrange the raw results data into pickled files to be processed by output_filtered_data.

    Parameters
    ----------
    analyzed_data : tuple
        Tuple containing multiple informations about the glycans analysis, as outputted by
        analyze_files or analyze_ms2.
        
    samples_names : list
        List of samples names extracted from file names.
        
    analyze_ms2 : tuple
        A tuple with two indexes: The first one indicates whether to analyze ms2 data and the
        second one indicates whether ms2 data should be forced to fit glycans composition.
    
    save_path : string
        A string containing the path to the working directory of the script.
        
    parameters : list
        A list of parameters for saving as metadata into the .gg file.
        
    from_GUI : boolean
        Whether or not this function is being executed by the GUI.
        
    file_name : string
        Allows to use a custom file name for the .gg file.
        
    Uses
    ----
    dill.dump : None
        Pickle the current state of __main__ or another module to a file.
    
    dill.load : Module object
        Update the selected module (default is __main__) with the state saved at filename.
        
    Returns
    -------
    nothing
        Creates raw_data files.
    '''
    # Grab args
    analyzed_data = args_dict.get('analyzed data', None)
    samples_names = args_dict.get('sample names', [])
    analyze_ms2 = args_dict.get('ms2 analyzed', False)
    save_path = args_dict.get('save path', None)
    parameters = args_dict.get('analysis parameters', {})
    adduct_combos = args_dict.get('adduct combos', [])
    library = args_dict.get('library', None)
    temp_folder = args_dict.get('temporary folder', None)
    file_name = args_dict.get('results file name', None)
    from_GUI = args_dict.get('from GUI', False)
    erase_files = args_dict.get('erase files', True)
    noise = args_dict.get('calculated noise levels', [])
    
    date = datetime.datetime.now()
    begin_time = str(date)[2:4]+str(date)[5:7]+str(date)[8:10]+"_"+str(date)[11:13]+str(date)[14:16]+str(date)[17:19]
    date, time = begin_time.split("_")
    time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
    print(time_formatted+'Arranging raw data...', end='', flush = True)
    df1 = []
    df2 = {"Sample_Number" : [], "File_Name" : [], "Average_Noise_Level" : []}
    eics_list = {}
    curve_fitting_dataframes = []
    isotopic_fits_dataframes = []
    if analyze_ms2:
        fragments_dataframes = []
    for i_i, i in enumerate(samples_names):
        isotopic_fits_dataframes.append({})
        
        temp_eic_rt = []
        for j in analyzed_data[1][i_i]:
            temp_eic_rt.append(j)
        
        # Name of EIC        
        eic_name = 'RTs'
        eics_list[i_i] = [f'RTs_{i_i}']
        
        # Create the retention time list for the samples
        with open(os.path.join(temp_folder, f"{i_i}_{eic_name}"), 'wb') as f:
            dill.dump(temp_eic_rt, f)
            f.close()
            
        curve_fitting_dataframes.append({})
        df2["Sample_Number"].append(i_i)
        df2["File_Name"].append(i)
        df2["Average_Noise_Level"].append(float("%.1f" % round(analyzed_data[2][i_i],1)))
        if analyze_ms2:
            df1.append({"Glycan" : [], "Adduct" : [], "mz" : [], "RT" : [], "AUC" : [], "PPM" : [], "S/N" : [], "Iso_Fitting_Score" : [], "Curve_Fitting_Score" : [], "Detected_Fragments" : []})
            fragments_dataframes.append({"Glycan" : [], "Adduct" : [], "Precursor_mz" : [], "Fragment" : [], "Fragment_mz" : [], "Fragment_Intensity" : [], "RT" : [], "% TIC explained" : []})
        else:
            df1.append({"Glycan" : [], "Adduct" : [], "mz" : [], "RT" : [], "AUC" : [], "PPM" : [], "S/N" : [], "Iso_Fitting_Score" : [], "Curve_Fitting_Score" : []})
            
    for i_i, i in enumerate(library): #i = glycan (key)
    
        # Load data from file
        with open(os.path.join(temp_folder, i), 'rb') as f:
            glycan = dill.load(f)
            f.close()
            
        for j_j, j in enumerate(glycan['Adducts_mz_data']): #j = adduct (key)
            for k_k, k in enumerate(glycan['Adducts_mz_data'][j]): #k = sample number (key)
                isotopic_fits_dataframes[k_k][i+'_'+j] = glycan['Adducts_mz_data'][j][k][4]
                
                # Determine names of EICs
                eic_name = str(i)+'+'+str(j)+' - '+str(float("%.4f" % round(glycan['Adducts_mz'][j], 4)))
                eics_list[k_k].append(eic_name)
                
                # Raw EIC
                temp_eic_int = []
                for l in glycan['Adducts_mz_data'][j][k][3]:
                    temp_eic_int.append(int(l))
                
                # Create the Raw EIC files for the samples and glycans
                with open(os.path.join(temp_folder, f"{k_k}_raw_{eic_name}"), 'wb') as f:
                    dill.dump(temp_eic_int, f)
                    f.close()
                
                temp_eic_int = []
                for l in glycan['Adducts_mz_data'][j][k][0]:
                    temp_eic_int.append(int(l))
                
                # Create the Filtered EIC files for the samples and glycans
                # with open(os.path.join(temp_folder, f"{k_k}_eic_{eic_name}"), 'wb') as f:
                    # dill.dump(temp_eic_int, f)
                    # f.close()
                    
                temp_eic_int = []
                for l in glycan['Adducts_mz_data'][j][k][2]:
                    temp_eic_int.append(int(l))
                    
                # Create the Smoothed EIC files for the samples and glycans
                with open(os.path.join(temp_folder, f"{k_k}_smoothed_{eic_name}"), 'wb') as f:
                    dill.dump(temp_eic_int, f)
                    f.close()
                    
            found = False
            for k_k, k in enumerate(glycan['Adducts_mz_data'][j]):
                if len(glycan['Adducts_mz_data'][j][k][1]) != 0:
                    found = True
            if not found:
                continue
            for k_k, k in enumerate(glycan['Adducts_mz_data'][j]): #k = sample (key)
                df1[k_k]["Glycan"].append(i)
                df1[k_k]["Adduct"].append(j)
                df1[k_k]["mz"].append(float("%.4f" % round(glycan['Adducts_mz'][j], 4)))
                temp_rts = []
                temp_aucs = []
                temp_ppm = []
                temp_s_n = []
                temp_iso_score = []
                temp_curve_score = []
                temp_curve_data_total = []
                for l_l, l in enumerate(glycan['Adducts_mz_data'][j][k][1]):
                    temp_rts.append(float("%.4f" % round(l['rt'], 4)))
                    temp_aucs.append(float("%.2f" % round(l['AUC'], 2)))
                    temp_ppm.append(float("%.2f" % round(l['Average_PPM'][0], 2)))
                    temp_s_n.append(float("%.1f" % round(l['Signal-to-Noise'], 1)))
                    if isnan(l['Iso_Fit_Score']):
                        temp_iso_score.append(0.0)
                    else:
                        temp_iso_score.append(float("%.4f" % round(l['Iso_Fit_Score'], 4)))
                    if isnan(l['Curve_Fit_Score'][0]):
                        temp_curve_score.append(0.0)
                    else:
                        temp_curve_score.append(float("%.4f" % round(l['Curve_Fit_Score'][0], 4)))
                    temp_curve_data_rt = []
                    temp_curve_data_actual = []
                    temp_curve_data_ideal = []
                    for m_m, m in enumerate(l['Curve_Fit_Score'][1]):
                        temp_curve_data_rt.append(m)
                        temp_curve_data_actual.append(l['Curve_Fit_Score'][2][m_m])
                        temp_curve_data_ideal.append(l['Curve_Fit_Score'][3][m_m])
                    temp_curve_data_total.append((temp_curve_data_rt, temp_curve_data_actual, temp_curve_data_ideal))
                if analyze_ms2:
                
                    # Load MS2 data
                    with open(os.path.join(temp_folder, 'frag_data_'+i), 'rb') as f:
                        glycan_fragments = dill.load(f)
                        f.close()
                    temp_fragments = glycan_fragments[j][k_k]
                    
                if len(temp_rts) == 0:
                    df1[k_k]["RT"].append([0.0])
                    df1[k_k]["AUC"].append([0.0])
                    df1[k_k]["PPM"].append([0.0])
                    df1[k_k]["S/N"].append([0.0])
                    df1[k_k]["Iso_Fitting_Score"].append([0.0])
                    df1[k_k]["Curve_Fitting_Score"].append([0.0])
                    if analyze_ms2:
                        df1[k_k]["Detected_Fragments"].append('Glycan+Adduct not found in sample')
                        temp_fragments = []
                else:
                    df1[k_k]["RT"].append(temp_rts)
                    df1[k_k]["AUC"].append(temp_aucs)
                    df1[k_k]["PPM"].append(temp_ppm)
                    df1[k_k]["S/N"].append(temp_s_n)
                    df1[k_k]["Iso_Fitting_Score"].append(temp_iso_score)
                    df1[k_k]["Curve_Fitting_Score"].append(temp_curve_score)
                    if analyze_ms2:
                        if len(temp_fragments) != 0:
                            for m in temp_fragments:
                                fragments_dataframes[k_k]["Glycan"].append(m[0])
                                fragments_dataframes[k_k]["Adduct"].append(m[1])
                                fragments_dataframes[k_k]["Fragment"].append(m[2])
                                fragments_dataframes[k_k]["Fragment_mz"].append(float("%.4f" % round(m[3], 4)))
                                fragments_dataframes[k_k]["Fragment_Intensity"].append(float("%.2f" % round(m[4], 2)))
                                fragments_dataframes[k_k]["RT"].append(float("%.4f" % round(m[5],4)))
                                fragments_dataframes[k_k]["Precursor_mz"].append(float("%.4f" % round(m[6], 4)))
                                fragments_dataframes[k_k]["% TIC explained"].append(float(m[7]))
                            df1[k_k]["Detected_Fragments"].append('Yes')
                        else:
                            df1[k_k]["Detected_Fragments"].append('No')
                    for m_m, m in enumerate(temp_rts):
                        temp_array = []
                        for n in temp_curve_data_total[m_m][0]:
                            temp_array.append(float("%.4f" % round(n,4)))
                        curve_fitting_dataframes[k_k][str(i)+"+"+str(j)+"_"+str(m)+"_RTs"] = temp_array
                        temp_array = []
                        for n in temp_curve_data_total[m_m][1]:
                            temp_array.append(int(n))
                        curve_fitting_dataframes[k_k][str(i)+"+"+str(j)+"_"+str(m)+"_Found_ints"] = temp_array
                        temp_array = []
                        for n in temp_curve_data_total[m_m][2]:
                            temp_array.append(int(n))
                        curve_fitting_dataframes[k_k][str(i)+"+"+str(j)+"_"+str(m)+"_Ideal_ints"] = temp_array
        try:
            if erase_files:
                os.remove(os.path.join(temp_folder, i))
                if analyze_ms2:
                    os.remove(os.path.join(temp_folder, 'frag_data_'+i))
        except:
            pass
            
    for i_i, i in enumerate(samples_names):
        combine_raw_data_per_sample(i_i, temp_folder)
                        
    with open(os.path.join(temp_folder, 'results'), 'wb') as f:
        if analyze_ms2:
            dill.dump([df1, df2, fragments_dataframes, version], f)
            del df1
            del df2
            del fragments_dataframes
        else:
            dill.dump([df1, df2, version], f)
            del df1
            del df2
        f.close()
        
    with open(os.path.join(temp_folder, 'eics_list'), 'wb') as f:
        dill.dump(eics_list, f)
        f.close()
        
    with open(os.path.join(temp_folder, 'curve_fittings'), 'wb') as f:
        dill.dump(curve_fitting_dataframes, f)
        del curve_fitting_dataframes
        f.close()
    with open(os.path.join(temp_folder, 'isotopic_fittings'), 'wb') as f:
        dill.dump(isotopic_fits_dataframes, f)
        del isotopic_fits_dataframes
        f.close()
    with open(os.path.join(temp_folder, 'noise_levels'), 'wb') as f:
        dill.dump(noise, f)
        del noise
        f.close
    if analyze_ms2 and len(analyzed_data[3]) > 0:
        with open(os.path.join(temp_folder, 'fragments_library'), 'wb') as f:
            dill.dump(analyzed_data[3], f)
            f.close
        with open(os.path.join(temp_folder, 'spectra_score'), 'wb') as f:
            dill.dump(analyzed_data[4], f)
            f.close
    with open(os.path.join(temp_folder, 'metadata'), 'wb') as f:
        parameters['analysis time'] = begin_time
        parameters['adduct combos'] = adduct_combos 
        parameters['analysis settings']['samples path'] = str(parameters['analysis settings']['samples path'])
        parameters['analysis settings']['save path'] = str(parameters['analysis settings']['save path'])
        dill.dump(parameters, f)
        f.close()
        
    if file_name == None:
        gg_name = begin_time+"_Analysis"
    else:
        if len(file_name.split('.')) > 1:
            file_name = file_name.split('.')[0]
        if "<date>" in file_name or "<time>" in file_name:
            temp_gg_name = []
            temp_gg_name_first = file_name.split('>')
            for i in temp_gg_name_first:
                temp_gg_name += i.split('<')
            gg_name = ''
            for word in temp_gg_name:
                if word == 'date':
                    gg_name += str(date)
                elif word == 'time':
                    gg_name += str(time)
                else:
                    gg_name += word
        else:
            gg_name = file_name
        counter = 0
        while True:
            if counter == 0 and os.path.isfile(os.path.join(save_path, gg_name+'.gg')):
                counter+=1
                continue
            elif counter != 0 and os.path.isfile(os.path.join(save_path, gg_name+'('+str(counter)+').gg')):
                counter+=1
                continue
            else:
                if counter != 0:
                    gg_name = gg_name+'('+str(counter)+')'
                else:
                    gg_name = gg_name
                break
                
    General_Functions.make_gg(temp_folder, save_path, gg_name, True)
    
    # Clean up temp dir in case of safeguard
    if not erase_files:
        for i_i, i in enumerate(samples_names):
            os.remove(os.path.join(temp_folder, f"{i_i}_eics"))
        os.remove(os.path.join(temp_folder, 'results'))
        os.remove(os.path.join(temp_folder, 'eics_list'))
        os.remove(os.path.join(temp_folder, 'curve_fittings'))
        os.remove(os.path.join(temp_folder, 'isotopic_fittings'))
        os.remove(os.path.join(temp_folder, 'noise_levels'))
        os.remove(os.path.join(temp_folder, 'metadata'))
        
    print("Done!")
    if from_GUI:
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        print(time_formatted+"File name is '"+gg_name+".gg'.")
    return begin_time

def print_sep(): ##Complete
    '''Prints a separator consisting of 48 '-' character.
    
    Parameters
    ----------
    none
    
    Uses
    ----
    nothing
    
    Returns
    -------
    nothing
        Just prints a fixed length separator.
    '''
    print('------------------------------------------------')
    
def pre_processing(data,
                   ms1_index,
                   ret_time_interval,
                   custom_noise,
                   data_id,
                   from_GUI = False):
    '''Calculates the noise level of samples and creates dummy empty arrays for use down the pipeline.
    
    Parameters
    ----------
    data : list
        A list containing the sample files.
        
    ms1_index : dict
        A dictionary containing the indexes of MS1 spectra for each file.
        
    ret_time_interval : list
        A list with the beggining and end of retention times to analyze.
        
    custom_noise : list
        If custom noise levels are inputted, uses it.
        
    data_id : int
        The ID of one file to be analyzed.
        
    Uses
    ----
    General_Functions.rt_noise_level_parameters_set : float, tuple
        Receives 2 combined arrays containing the x and y information of a spectrum
        and calculate parameters for dynamic noise calculation down the pipeline.
        
    Returns
    -------
    tuple
        A tuple containing the dummy arrays as well as the calculated noise for each sample.
    '''
    try:
        zeroes_arrays= [0.0]*len(ms1_index[data_id])
        inf_arrays = [inf]*len(ms1_index[data_id])
        rt_array_report = [None]*len(ms1_index[data_id])
        threads_arrays = []
        ms1_id = []
        temp_noise = [1.0]*len(ms1_index[data_id])
        temp_avg_noise = [1.0]*len(ms1_index[data_id])
        list_for_avg = []        
        for j_j, j in enumerate(ms1_index[data_id]):
            rt_array_report[j_j] = data[j]['retentionTime']
            mz_ints = [data[j]['m/z array'], data[j]['intensity array']]
            
            # If custom noise is used...
            if custom_noise[0]:
                temp_noise[j_j] = custom_noise[1][data_id]
                temp_avg_noise[j_j] = custom_noise[1][data_id]
            
            # Else it calculates the noise level
            elif data[j]['retentionTime'] >= ret_time_interval[0] and data[j]['retentionTime'] <= ret_time_interval[1]:
                if len(data[j]['intensity array']) == 0:
                    temp_noise[j_j] = (1.0, 0.0, 0.0)
                    temp_avg_noise[j_j] = 1.0
                else:
                    threads_arrays.append(j)
                    ms1_id.append(j_j)
                    temp_noise_segments = General_Functions.rt_noise_level_parameters_set(mz_ints, "segments")
                    temp_noise[j_j] = temp_noise_segments
                    temp_noise_whole = General_Functions.rt_noise_level_parameters_set(mz_ints, "whole")
                    temp_avg_noise[j_j] = temp_noise_whole
                    if temp_noise_whole != 1.0:
                        list_for_avg.append(temp_noise_whole)
            
            # Outside the RT range, dummy noise levels
            else:
                temp_noise[j_j] = (1.0, 0.0, 0.0)
                temp_avg_noise[j_j] = 1.0
                
        if len(list_for_avg) != 0:
            avg = np.mean(list_for_avg)      
            for i_i, i in enumerate(temp_avg_noise):
                if i == 1.0:
                    temp_avg_noise[i_i] = avg
        
        acquisition_interval = rt_array_report[len(rt_array_report)//2]-rt_array_report[(len(rt_array_report)//2)-1]
        sampling_rate = round((1/60)/acquisition_interval)
        
        return zeroes_arrays, inf_arrays, threads_arrays, ms1_id, rt_array_report, temp_noise, temp_avg_noise, data_id, sampling_rate
        
    except KeyboardInterrupt:
        if not from_GUI:
            print("\n\n----------Execution cancelled by user.----------\n", flush=True)
        raise SystemExit(1)
    
def analyze_files(args_dict):
    '''Integrates all the file-accessing associated functions in this script to go
    through the files data, draw and process eic of hypothetical glycans, does 
    peak-picking and calculates AUC of the peaks.

    Parameters
    ----------
    library : dict
        A glycans library, as generated by Library_Tools.full_glycans_library.
        
    lib_size : int
        The length of the library.
        
    data : list
        A list with each index containing a generator object of the sample file
        to be parsed.
        
    ms1_index : dict
        A dictionary containing the ms1 indexes of each sample file.
        
    tolerance : tuple
        First index contains the unit of the tolerance and the second one is the value of 
        that unit.
        
    ret_time_interval : tuple
        A tuple where the first index contains the beggining time of the retention time
        interval you wish to analyze and the second contains the end time.
        
    min_isotops : int
        The minimum amount of isotopologues required to consider an RT mz peak valid.
    
    min_ppp : tuple
        A tuple where the first index contains a boolean indicating whether or not to
        consider this parameter and the second one containing the minimum amount of
        data points per peak. If min_ppp[0] is set to False, calculates this automatically.
        
    max_charges : int
        The maximum amount of charges the queried mz should have.
        
    custom_noise : tuple
        A tuple containing two indexes: The first one indicates whether or not to use a custom
        noise set by the user and the second one is a list of noise levels, with each index
        indicating the noise level of each set file.
        
    close_peaks : tuple
        A tuple where the first index contains a boolean indicating whether or not to
        consider this parameter and the second one contains the amount of peaks to save.
        If close_peaks[0] is set to True, selects only the most intense peak and its 
        close_peaks[1] surrounding peaks.
        
    multithreaded : boolean
        Whether or not to use multiple threads.
        
    number_cores : string or int
        Number of cores to be used.

    Uses
    ----
    analyze_glycan : tuple
        Returns a tuple containing a series of informations for a given glycan.

    Returns
    -------
    analyzed_data : dict
        A dictionary similar to the one generated by the full glycans library generating
        function, with the added information of each adducts' peaks' retention time and
        AUC.
        
    rt_array_report : dict
        A dictionary with lists containing the retention times for each sample.
        
    noise : dict
        A dictionary containing the noise level for each sample.
    '''
    # Grab the args
    library = args_dict.get('library', None)
    lib_size = args_dict.get('library size', None)
    data = args_dict.get('raw data', None)
    ms1_index = args_dict.get('ms1 index', None)
    tolerance = args_dict.get('mass tolerance', ['mz', 0.01])
    ret_time_interval = args_dict.get('retention time interval', [0, 999])
    min_isotops = args_dict.get('min isotopologue peaks', 2)
    min_ppp = args_dict.get('min points per peak', [False, 0])
    adduct_combos = args_dict.get('adduct combos', None)
    max_charges = args_dict.get('maximum charges', 3)
    custom_noise = args_dict.get('custom noise level', [False, []])
    close_peaks = args_dict.get('only x most intense peaks', [False, 3])
    multithreaded = args_dict.get('multithreaded', True)
    number_cores = args_dict.get('number of cpu cores', 'all')
    begin_time = args_dict.get('analysis start time', None)
    temp_folder = args_dict.get('temporary folder', None)
    from_GUI = args_dict.get('from GUI', False)
    
    adduct_combos_dict = {General_Functions.comp_to_formula(adduct[0]): {'charges': adduct[1], 'comp': adduct[0]} for adduct in adduct_combos}
                  
    mean_sel_peaks = 0.0
    noise = {}
    noise_avg = {}
    rt_array_report = {}
    zeroes_arrays = [[]]*len(data)
    inf_arrays = [[]]*len(data)
    threads_arrays = [[]]*len(data)
    ms1_id = [[]]*len(data)
    sampling_rates = [1]*len(data)
    
    if not custom_noise[0]:
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        print(time_formatted+'Calculating noise levels...', end='', flush = True)
    
    results = []
    if multithreaded:
        if number_cores == 'all':
            cpu_count = (os.cpu_count())-2
            if cpu_count <= 0:
                cpu_count = 1
        else:
            number_cores = int(number_cores)
            if number_cores > (os.cpu_count())-2:
                cpu_count = (os.cpu_count())-2
                if cpu_count <= 0:
                    cpu_count = 1
            else:
                cpu_count = number_cores
    else:
        cpu_count = 1
    
    with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_count if cpu_count < 60 else 60) as executor:
        for i_i, i in enumerate(data):
            result = executor.submit(pre_processing,
                                     i,
                                     ms1_index,
                                     ret_time_interval,
                                     custom_noise,
                                     i_i,
                                     from_GUI)
            results.append(result)
            
        for index, i in enumerate(results):
            result_data = i.result()
            zeroes_arrays[result_data[7]] = result_data[0]
            inf_arrays[result_data[7]] = result_data[1]
            threads_arrays[result_data[7]] = result_data[2]
            ms1_id[result_data[7]] = result_data[3]
            rt_array_report[result_data[7]] = result_data[4]
            noise[result_data[7]] = result_data[5]
            noise_avg[result_data[7]] = percentile(result_data[6], 66.8)
            sampling_rates[result_data[7]] = result_data[8]
            results[index] = None
    
    ambiguities = {}
    for i_i, i in enumerate(library):
        for k_k, k in enumerate(library):
            if k_k >= i_i:
                break
            if library[k]['Neutral_Mass'] == library[i]['Neutral_Mass']:
                if i not in ambiguities.keys():
                    ambiguities[i] = [k]
                else:
                    ambiguities[i].append(k)
        
    print('Done!')
    time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
    print(time_formatted+"Pre-processing done in "+str(datetime.datetime.now() - begin_time).split(".")[0]+"!")
    print_sep()
    print(time_formatted+"Starting MS1 tracing...")
    begin_time = datetime.datetime.now()
    
    is_result = None
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_count if cpu_count < 60 else 60) as executor:
        for i_i, i in enumerate(library):
            if i in ambiguities.keys(): #skips ambiguities
                results.append(i)
                continue
            if i == 'Internal Standard':
                is_result = executor.submit(analyze_glycan, 
                                            library,
                                            lib_size,
                                            data,
                                            ms1_index,
                                            tolerance,
                                            ret_time_interval,
                                            min_isotops,
                                            min_ppp,
                                            adduct_combos_dict,
                                            max_charges,
                                            noise,
                                            noise_avg,
                                            close_peaks,
                                            zeroes_arrays,
                                            inf_arrays,
                                            threads_arrays,
                                            rt_array_report,
                                            ms1_id,
                                            i,
                                            i_i,
                                            sampling_rates,
                                            from_GUI)
            else:
                result = executor.submit(analyze_glycan, 
                                         library,
                                         lib_size,
                                         data,
                                         ms1_index,
                                         tolerance,
                                         ret_time_interval,
                                         min_isotops,
                                         min_ppp,
                                         adduct_combos_dict,
                                         max_charges,
                                         noise,
                                         noise_avg,
                                         close_peaks,
                                         zeroes_arrays,
                                         inf_arrays,
                                         threads_arrays,
                                         rt_array_report,
                                         ms1_id,
                                         i,
                                         i_i,
                                         sampling_rates,
                                         from_GUI)
                results.append(result)
        for index, i in enumerate(results):
            if type(i) == str: #ambiguity
                time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
                print(time_formatted+'Traced glycan '+i+': '+str(index+1)+'/'+str(lib_size))
            else:
                result_data = i.result()
                time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
                print(time_formatted+'Traced glycan '+str(result_data[1])+': '+str(index+1)+'/'+str(lib_size))
                
                # Pickling all the data into separate files
                with open(os.path.join(temp_folder, result_data[1]), 'wb') as f:
                    dill.dump(result_data[0], f)
                    f.close()
                
            results[index] = None
    
    if len(ambiguities) > 0:
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        print(time_formatted+'Sorting ambiguities...')
        for i in ambiguities: #sorts ambiguities
            shutil.copy(os.path.join(temp_folder, ambiguities[i][0]), os.path.join(temp_folder, i))
            with open(os.path.join(temp_folder, i), 'rb') as f:
                glycan = dill.load(f)
                f.close()
            glycan['Monos_Composition'] = General_Functions.form_to_comp_glycans(i)
            with open(os.path.join(temp_folder, i), 'wb') as f:
                dill.dump(glycan, f)
                f.close()
    
    if is_result != None:
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
        print(time_formatted+'Traced Internal Standard: '+str(lib_size)+'/'+str(lib_size))
        with open(os.path.join(temp_folder, 'Internal Standard'), 'wb') as f:
            dill.dump(is_result.result()[0], f)
            f.close()
        del is_result
        
    time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
    print(time_formatted+'MS1 tracing done in '+str(datetime.datetime.now() - begin_time).split(".")[0]+'!')
    return None, rt_array_report, noise_avg, noise
    
def analyze_glycan(library,
                  lib_size,
                  data,
                  ms1_index,
                  tolerance,
                  ret_time_interval,
                  min_isotops,
                  min_ppp,
                  adduct_combos_dict,
                  max_charges,
                  noise,
                  noise_avg,
                  close_peaks,
                  zeroes_arrays,
                  inf_arrays,
                  threads_arrays,
                  rt_arrays,
                  ms1_id,
                  i,
                  i_i,
                  sampling_rates,
                  from_GUI = False):
    '''Analyzes one single glycan. Core function of analyze_files.
    
    Parameters
    ----------
    library : dict
        A glycans library, as generated by Library_Tools.full_glycans_library.
        
    lib_size : int
        The length of the library.
        
    data : list
        A list with each index containing a generator object of the sample file
        to be parsed.
        
    ms1_index : dict
        A dictionary containing the ms1 indexes of each sample file.
        
    tolerance : tuple
        First index contains the unit of the tolerance and the second one is the value of 
        that unit.
        
    ret_time_interval : tuple
        A tuple where the first index contains the beggining time of the retention time
        interval you wish to analyze and the second contains the end time.
        
    min_isotops : int
        The minimum amount of isotopologues required to consider an RT mz peak valid.
    
    min_ppp : tuple
        A tuple where the first index contains a boolean indicating whether or not to
        consider this parameter and the second one containing the minimum amount of
        data points per peak. If min_ppp[0] is set to False, calculates this automatically.
        
    max_charges : int
        The maximum amount of charges the queried mz should have.
        
    noise : dict
        A dictionary containing the parameters for noise calculation for each sample.
        
    noise_avg : dict
        A dictionary containing the average noise of each sample.
        
    close_peaks : tuple
        A tuple where the first index contains a boolean indicating whether or not to
        consider this parameter and the second one contains the amount of peaks to save.
        If close_peaks[0] is set to True, selects only the most intense peak and its 
        close_peaks[1] surrounding peaks.
        
    zeroes_arrays : list
        List of correctly sized arrays for each sample, containing only zeroes.
        
    inf_arrays : list
        List of correctly sized arrays for each sample, containing only infinities.
        
    threads_arrays : list
        List of the IDs of the spectra of the chromatogram that will be analyzed.
        
    rt_arrays : dict
        A dictionary containing the retention time arrays for each file.
    
    ms1_id : list
        List of the MS1 spectra that will be analyzed, synchronized to the threads_array.
        
    i : string
        The glycan to be analyzed.
        
    i_i : int
        The ID of the glycan to be analyzed.
        
    total_glycans : int
        The total amount of glycans to be analyzed, to be used for progress reporting.
        
    Uses
    ----
    File_Accessing.eic_from_glycan : tuple
        Generates a very processed EIC for each adduct of each glycan of each sample.
        Removes non-monoisotopic peaks, check charges and calculates multiple quality 
        scoring data.
        
    File_Accessing.eic_smoothing : list
        Smoothes the EIC using the Savitsky Golay algorithm. Smoothed EIC may be
        used for peak-picking and curve-fitting scoring afterwards.

    File_Accessing.peaks_from_eic() : list
        Does multi peak-picking in a given smoothed EIC.

    File_Accessing.average_ppm_calc : tuple
        Calculates the arithmetic mean of the PPM differences of a given peak.
        
    File_Accessing.iso_fit_score_calc : float
        Calculates the mean isotopic fitting score of a given peak.
        
    File_Accessing.peak_curve_fit : tuple
        Calculates the fitting between the actual peak and an ideal peak based
        on a calculated gaussian bell curve.
        
    Returns
    -------
    glycan_data : dict
        A modified dictionary entry of the library, including a lot of information of the glycan obtained by the program.
        
    i : string
        The glycan analyzed.
    '''
    try:
        glycan_data = library[i]
        glycan_data['Adducts_mz_data'] = {}
        temp_eic = File_Accessing.eic_from_glycan(data,
                                                  i,
                                                  glycan_data,
                                                  ms1_index,
                                                  ret_time_interval,
                                                  tolerance,
                                                  min_isotops,
                                                  noise,
                                                  noise_avg,
                                                  adduct_combos_dict,
                                                  max_charges,
                                                  zeroes_arrays,
                                                  inf_arrays,
                                                  threads_arrays,
                                                  rt_arrays,
                                                  ms1_id,
                                                  sampling_rates)
        for j in temp_eic[0]: #moving through adducts
            glycan_data['Adducts_mz_data'][j] = {}
            for k in temp_eic[0][j]: #moving through samples
                temp_eic_smoothed = File_Accessing.eic_smoothing(temp_eic[0][j][k])
                glycan_data['Adducts_mz_data'][j][k] = []
                glycan_data['Adducts_mz_data'][j][k].append(temp_eic[0][j][k][1]) #processed chromatogram
                glycan_data['Adducts_mz_data'][j][k].append([]) #placeholder for inserting data about the glycan and adduct
                glycan_data['Adducts_mz_data'][j][k].append(temp_eic_smoothed[1]) #smoothed chromatogram
                glycan_data['Adducts_mz_data'][j][k].append(temp_eic[4][j][k][1]) #raw chromatogram
                glycan_data['Adducts_mz_data'][j][k].append(temp_eic[5][j][k]) #isotopic fits data
                if max(temp_eic[0][j][k][1]) < noise_avg[k] and i != "Internal Standard":
                    continue
                temp_peaks = File_Accessing.peaks_from_eic(temp_eic[0][j][k],
                                                           temp_eic_smoothed,
                                                           temp_eic[4][j][k],
                                                           ret_time_interval,
                                                           min_ppp,
                                                           close_peaks,
                                                           i)
                if len(temp_peaks) == 0:
                    continue
                temp_peaks_auc = File_Accessing.peaks_auc_from_eic(temp_eic[0][j][k],
                                                                   ms1_index[k],
                                                                   temp_peaks)
                for l_l, l in enumerate(temp_peaks):
                    if temp_peaks_auc[l_l] >= noise_avg[k]:
                        l['Curve_Fit_Score'] = File_Accessing.peak_curve_fit(temp_eic_smoothed, l)
                        l['Average_PPM'] = File_Accessing.average_ppm_calc(temp_eic[1][j][k], (tolerance[0], tolerance[1], glycan_data['Adducts_mz'][j]), l, l['Curve_Fit_Score'][3])
                        l['Iso_Fit_Score'] = File_Accessing.iso_fit_score_calc(temp_eic[2][j][k], l, l['Curve_Fit_Score'][3])
                        l['AUC'] = temp_peaks_auc[l_l]
                        l['Signal-to-Noise'] = l['int']/(General_Functions.local_noise_calc(noise[k][l['id']], glycan_data['Adducts_mz'][j], noise_avg[k]))
                        glycan_data['Adducts_mz_data'][j][k][1].append(l)
        return glycan_data, i
        
    except KeyboardInterrupt:
        if not from_GUI:
            print("\n\n----------Execution cancelled by user.----------\n", flush=True)
        raise SystemExit(1)
    
def analyze_ms2(args_dict):
    '''Analyzes the MS2 data in the sample files, outputting the found matches.
    
    Parameters
    ----------
    ms2_index : dict
        A dictionary containing the ms2 indexes of each sample file.
        
    data : list
        A list with each index containing a generator object of the sample file
        to be parsed.
        
    analyzed_data : tuple
        Data outputted by the analyze_data function.
        
    ret_time_interval : tuple
        A tuple where the first index contains the beggining time of the retention time
        interval you wish to analyze and the second contains the end time.
        
    tolerance : tuple
        First index contains the unit of the tolerance and the second one is the value of 
        that unit.
        
    min_max_mono : tuple
        Minimum and maximum amount of monosaccharides for the hypotethical glycans in the
        library. ie. (5, 20).

    min_max_hex : tuple
        Minimum and maximum amount of hexoses for the hypotethical glycans in the library.
        ie. (5, 20).

    min_max_hexnac : tuple
        Minimum and maximum amount of N-Acetyl hexosamines for the hypotethical glycans
        in the library. ie. (5, 20).

    min_max_sia : tuple
        Minimum and maximum amount of sialic acids for the hypotethical glycans in the
        library. ie. (5, 20).

    min_max_fuc : tuple
        Minimum and maximum amount of deoxyhexoses for the hypotethical glycans in the
        library. ie. (5, 20).

    min_max_ac : tuple
        Minimum and maximum amount of N-Acetyl Neuraminic acids for the hypotethical
        glycans in the library. ie. (5, 20).

    min_max_gc : tuple
        Minimum and maximum amount of N-Glicolyl Neuraminic acids for the hypotethical
        glycans in the library. ie. (5, 20).
        
    min_max_hn : tuple
        Minimum and maximum amount of Hexosamines for the hypotethical glycans in the library. ie. (5, 20).
        
    min_max_ua : tuple
        Minimum and maximum amount of Uronic Acids for the hypotethical glycans in the library. ie. (5, 20).
        
    max_charges : int
        The maximum amount of charges to calculate per glycan.
        
    tag_mass : float
        The tag's added mass to the glycans, if the glycans are tagged.
        Default = 0 (No Tag).

    forced : string
        Indicates whether the function should force strict conditions based on the
        biological knowledge of glycans in order to avoid possible false positives when
        analysing N-glycans, O-glycans or GAGs.

    permethylated : boolean
        Whether or not the sample was permethylated.
        
    reduced : boolean
        Whether or not the sample was reduced.
        
    lactonized_ethyl_esterified : boolean
        Whether the glycans were submitted to lactonization/ethyl-esterification
        derivatization, which differentiates the mass of alpha2-3 and alpha2-6 -bound 
        sialic acids.
        
    filter_output : boolean
        Whether or not to force the output to fit glycans compositions.
        
    unrestricted_fragments : boolean
        Whether or not should take any fragment found
    
    rt_tolerance : float
        Tolerance of retention time (in minutes) at which an MS2 feature can be attributed to a specific retention time peak and also for peaks in different samples to be regarded as the same peak (and thus be compared with each other).
        
    multithreaded : boolean
        Whether or not to use multiple threads.
        
    number_cores : string or int
        Number of cores to be used.
        
    Uses
    ----
    Library_Tools.fragments_library : list
        Generates a list of combinatorial analysis of monosaccharides from the minimum
        amount of monosaccharides to the maximum amount of monosaccharides set, then uses 
        the generated library and increments it with calculations with a series of information,
        clumps the duplicated mzs together as one query and then the script can use this
        library as a fragments library for MS2 analysis.
        
    General_Functions.form_to_charge : int
        Converts adducts formula into raw charge.
        
    analyze_glycan_ms2 : tuple
        Returns a series of information on the MS2 data analyzed.
        
    Returns
    -------
    analyzed_data : tuple
        Data outputted by the analyze_data function.
        
    fragments_data : dict
        Dictionary containing the fragments data.
    '''
    # Grab the args
    ms2_index = args_dict.get('ms2 index', None) 
    data = args_dict.get('raw data', None)
    analyzed_data = args_dict.get('analyzed data', None) 
    rt_interval = args_dict.get('retention time interval', [0, 999])
    tolerance = args_dict.get('mass tolerance', ['mz', 0.01])
    min_max_monos = args_dict.get('min/max monosaccharides', [0,0])
    min_max_hex = args_dict.get('min/max hexoses', [0,0])
    min_max_hexnac = args_dict.get('min/max hexnac', [0,0])
    min_max_xyl = args_dict.get('min/max xyloses', [0,0])
    min_max_sia = args_dict.get('min/max sialic acids', [0,0])
    min_max_fuc = args_dict.get('min/max fucoses', [0,0])
    min_max_ac = args_dict.get('min/max acetyl sialic acids', [0,0])
    min_max_gc = args_dict.get('min/max glycolyl sialic acids', [0,0])
    min_max_hn = args_dict.get('min/max hexosamines', [0,0])
    min_max_ua = args_dict.get('min/max uronic acids', [0,0])
    adduct_combos = args_dict.get('adduct combos', None)
    max_charges = args_dict.get('maximum charges', 3)
    tag_mass = args_dict.get('reducing end tag', 0.0)
    forced = args_dict.get('glycan class', None)
    permethylated = args_dict.get('permethylated', False)
    reduced = args_dict.get('reducing end reduced', False)
    lactonized_ethyl_esterified = args_dict.get('lactonized/ethyl-esterified', False)
    filter_output = args_dict.get('assign fragments compatible with precursor', True)
    unrestricted_fragments = args_dict.get('annotate all glycans in library', False)
    rt_tolerance = args_dict.get('fragmentation retention time tolerance', 0.2)
    multithreaded = args_dict.get('multithreaded', True)
    number_cores = args_dict.get('number of cpu cores', 'all')
    library = args_dict.get('library', None)
    temp_folder = args_dict.get('temporary folder', None)
    from_GUI = args_dict.get('from GUI', False)
    custom_monos = args_dict.get('custom monosaccharides', [])
    
    adduct_combos_dict = {General_Functions.comp_to_formula(adduct[0]): {'charges': adduct[1], 'comp': adduct[0]} for adduct in adduct_combos}
                
    begin_time = datetime.datetime.now()
    
    # Check if any spectra file has MS2 data
    no_ms2 = True
    for i in ms2_index:
        if len(ms2_index[i]) != 0:
            no_ms2 = False
            break
            
    # If no spectra file has MS2 data, create dummy data and finish execution
    if no_ms2:
        print('No MS2 data to analyze...')
        dummy_fragment_data = {}
        for i in library:
            dummy_fragment_data[i] = {}
            for j in library[i]['Adducts_mz']:
                dummy_fragment_data[i][j] = {}
                for k_k, k in enumerate(data):
                    dummy_fragment_data[i][j][k_k] = []
            with open(os.path.join(temp_folder, 'frag_data_'+i), 'wb') as f:
                dill.dump(dummy_fragment_data[i], f)
                f.close()
            dummy_fragment_data[i] = None
        return library, analyzed_data[1], analyzed_data[2], {}, {}
        
    # Otherwise, analyze the MS2 data
    time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
    print(time_formatted+'Analyzing MS2 data...')
    
    # Create the multithreading safe dictionaries to store the fragments
    fragments_dict = multiprocessing.Manager().dict()
    fragments_per_glycan = multiprocessing.Manager().dict()
    
    # Calculate the number of usable cores
    if multithreaded:
        if number_cores == 'all':
            cpu_count = (os.cpu_count())-2
            if cpu_count <= 0:
                cpu_count = 1
        else:
            number_cores = int(number_cores)
            if number_cores > (os.cpu_count())-2:
                cpu_count = (os.cpu_count())-2
                if cpu_count <= 0:
                    cpu_count = 1
            else:
                cpu_count = number_cores
    else:
        cpu_count = 1
    
    # Resolve the reducing end tag
    tag = General_Functions.determine_tag_comp(tag_mass)
    
    # Warn about fragment library being built
    begin_time_frag_library_building = datetime.datetime.now()
    time_formatted = str(begin_time_frag_library_building).split(" ")[-1].split(".")[0]+" - "
    print(f"{time_formatted}Building fragments library...")
    
    # Start fragment library building
    with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_count if cpu_count < 60 else 60) as executor:
        for index, glycan in enumerate(library):
            executor.submit(Library_Tools.calculate_glycan_fragments,
                            glycan,
                            library[glycan]['Monos_Composition'],
                            adduct_combos_dict,
                            tolerance,
                            tag,
                            permethylated,
                            reduced,
                            lactonized_ethyl_esterified,
                            forced,
                            fragments_dict,
                            fragments_per_glycan,
                            custom_monos)
            
    # Index the fragments, arranged by mz, to allow for binary search
    indexed_fragments = {}
    for fragment, frag_data in fragments_dict.items():
        for adduct, mz in frag_data['Adducts_mz'].items():
            indexed_fragments[mz] = indexed_fragments.get(mz, []) + [f"{fragment}_{adduct}"]
    indexed_fragments = dict(sorted(indexed_fragments.items()))
    indexed_fragments_list = list(indexed_fragments.keys())
                                  
    # Notify that library building is complete, how long it took and its size
    time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
    print(f"{time_formatted}Fragments library built in {str(datetime.datetime.now()-begin_time_frag_library_building).split('.')[0]}.")
    print(f"{time_formatted}Fragments library length: {len(fragments_dict)}")
    
    # Prepare for MS2 spectra scanning
    fragments_data = {}
    print(time_formatted+'Scanning MS2 spectra...')
    scan_begin_time = datetime.datetime.now()
    
    # Create dictionary to store whether a spectrum was annotated or not
    all_samples_analyzed_spectra = {}
    
    # Create dictionary for fragments ranking
    fragments_ranking = {}
    
    # Scan the MS2 spectra
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_count if cpu_count < 60 else 60) as executor:
        for i_i, i in enumerate(library): #goes through each glycan found in analysis
            with open(os.path.join(temp_folder, i), 'rb') as f:
                glycan = dill.load(f)
                f.close()
            result = executor.submit(analyze_glycan_ms2,
                                     ms2_index,
                                     fragments_dict,
                                     fragments_per_glycan,
                                     indexed_fragments,
                                     indexed_fragments_list,
                                     data, 
                                     glycan,
                                     adduct_combos_dict,
                                     lactonized_ethyl_esterified,
                                     rt_interval,
                                     tolerance,
                                     filter_output,
                                     unrestricted_fragments,
                                     rt_tolerance,
                                     i_i,
                                     i,
                                     from_GUI)
            results.append(result)
            
        for index, i in enumerate(results):

            result_data = i.result()
            
            for sample, spectra in result_data[2].items():
                if sample not in all_samples_analyzed_spectra.keys():
                    all_samples_analyzed_spectra[sample] = {}
                for spectrum, analyzed in spectra.items():
                    if spectrum not in all_samples_analyzed_spectra[sample].keys():
                        all_samples_analyzed_spectra[sample][spectrum] = False
                    if analyzed:
                        all_samples_analyzed_spectra[sample][spectrum] = True
            
            for fragment, number in result_data[3].items():
                fragments_ranking[fragment] = fragments_ranking.get(fragment, 0) + number
                            
            time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
            print(time_formatted+'Analyzed glycan '+str(result_data[1])+': '+str(index+1)+'/'+str(len(library)))
            with open(os.path.join(temp_folder, 'frag_data_'+result_data[1]), 'wb') as f:
                dill.dump(result_data[0], f)
                f.close()
            results[index] = None
    
    sorted_fragments_ranking = dict(sorted(fragments_ranking.items(), key=lambda item: item[1], reverse=True))
    top_10_fragments = list(sorted_fragments_ranking.keys())[:10]
    
    spectra_to_check = {}
    for file, spectra in all_samples_analyzed_spectra.items():
        spectra_index_to_check = [spectrum for spectrum, analyzed in spectra.items() if not analyzed]
        spectra_to_check[file] = spectra_index_to_check
    
    results = []
    spectra_scores = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_count if cpu_count < 60 else 60) as executor:
        for file, spectra in spectra_to_check.items():
            result = executor.submit(score_remaining_spectra,
                                     file,
                                     spectra,
                                     data,
                                     top_10_fragments,
                                     fragments_dict,
                                     tolerance)
            results.append(result)
        
        for index, i in enumerate(results):
            result_data = i.result()
            
            spectra_scores[result_data[0]] = result_data[1]
            
            results[index] = None
        
    time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
    print(time_formatted+'MS2 analysis done in '+str(datetime.datetime.now() - begin_time).split(".")[0]+'!')
    return library, analyzed_data[1], analyzed_data[2], dict(fragments_dict), spectra_scores
                                 
def analyze_glycan_ms2(ms2_index,
                       fragments,
                       fragments_per_glycan,
                       indexed_fragments,
                       indexed_fragments_list,
                       data, 
                       glycan_data,
                       adduct_combos_dict,
                       lactonized_ethyl_esterified,
                       rt_interval,
                       tolerance,
                       filter_output,
                       unrestricted_fragments,
                       rt_tolerance,
                       i_i,
                       i,
                       from_GUI = False):
    '''Core function of analyze_ms2. Analyze a single glycan.
    
    Parameters
    ----------
    ms2_index : dict
        A dictionary containing the ms2 indexes of each sample file.
        
    data : list
        A list with each index containing a generator object of the sample file
        to be parsed.
        
    analyzed_data : tuple
        Data outputted by the analyze_data function.
        
    lactonized_ethyl_esterified : boolean
        Whether the glycans were submitted to lactonization/ethyl-esterification
        derivatization, which differentiates the mass of alpha2-3 and alpha2-6 -bound 
        sialic acids.
        
    rt_interval : list
        A list with the beggining and end of retention times to analyze.
        
    tolerance : tuple
        First index contains the unit of the tolerance and the second one is the value of 
        that unit.
        
    filter_output : boolean
        Whether or not to force the output to fit glycans compositions.
        
    unrestricted_fragments : boolean
        Whether or not should take any fragment found
    
    rt_tolerance : float
        Tolerance of retention time (in minutes) at which an MS2 feature can be attributed to a specific retention time peak and also for peaks in different samples to be regarded as the same peak (and thus be compared with each other).
        
    i : string
        The glycan to be analyzed.
        
    i_i : int
        The ID of the glycan to be analyzed.
    
    Uses
    ----
    nothing
    
    Returns
    -------
    tuple
        A series of information on the MS2 data analyzed.
    '''
    try:
        # Create Analyzed Spectra dictionary
        analyzed_spectra = {}
        
        # Create fragments ranking
        fragments_ranking = {}
        
        # Superscripts for better annotation
        superscripts = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹', '+': '⁺', '-': '⁻', '=': '⁼', '(': '⁽', ')': '⁾', 'n': 'ⁿ', 'i': 'ⁱ'}
        
        # Temporary dictionary to save the data on
        glycan_fragments_data = {}
        
        # Go through each adduct of the glycan
        for adduct in glycan_data['Adducts_mz_data']:
            
            # Calculate the charge of the adduct
            adduct_name = adduct_combos_dict[adduct]['comp']
            adduct_charge = adduct_combos_dict[adduct]['charges']
            glycan_fragments_data[adduct] = {}
            
            # Go through each file
            for file_index, file in enumerate(data):
        
                glycan_fragments_data[adduct][file_index] = []
                
                # Just to be sure, check again whether the file has MS2 data
                if len(ms2_index[file_index]) == 0:
                    continue
                    
                # Checks whether the glycan+adduct were found in MS1 and whether unrestricted_fragments is on
                if len(glycan_data['Adducts_mz_data'][adduct][file_index][1]) == 0 and not unrestricted_fragments:
                    continue
                
                # Add sample to analyzed spectra
                if file_index not in analyzed_spectra.keys():
                    analyzed_spectra[file_index] = {}
                    
                # Go through each MS2 spectrum within the spectra file
                for spectrum in ms2_index[file_index]:
                    
                    # If the retention time of that spectrum is outside the chosen analysis range, ignore it
                    if file[spectrum]['retentionTime'] < rt_interval[0] or file[spectrum]['retentionTime'] > rt_interval[1]:
                        continue
                        
                    # Add spectrum to analyzed spectra
                    if spectrum not in analyzed_spectra[file_index].keys():
                        analyzed_spectra[file_index][spectrum] = False
                        
                    # If there are no peaks in this MS2 spectrum, skip it
                    if len(file[spectrum]['intensity array']) == 0:
                        continue
                        
                    # If not unrestricted_fragments and there are no identified MS1 peaks in which this spectrum falls within boundaries of, skip this spectrum
                    if not unrestricted_fragments:
                        
                        peak_intervals = [[peak['peak_interval'][0], peak['peak_interval'][1]] for peak in glycan_data['Adducts_mz_data'][adduct][file_index][1]]
                        
                        found = False
                        for peak_interval in peak_intervals:
                            if file[spectrum]['retentionTime'] > peak_interval[0] and file[spectrum]['retentionTime'] < peak_interval[1]:
                                found = True
                                break
                                
                        if not found:
                            continue
                            
                    # Check whether the precursor mz matches the adduct mz up to the fourth isotopic peak within five times the regular tolerance
                    found_matching_mz = False
                    for isotopic_peak_index, isotopic_peak_mass in enumerate(glycan_data['Isotopic_Distribution_Masses']): 
                        
                        if isotopic_peak_index > 4:
                            break
                            
                        target_mz = (isotopic_peak_mass+(General_Functions.h_mass*adduct_charge))/abs(adduct_charge)
                        
                        tolerance_calculated = General_Functions.tolerance_calc(tolerance[0], tolerance[1], target_mz)*5
                        
                        if 'isolation window lower offset' in file[spectrum] and 'isolation window upper offset' in file[spectrum]:
                            lower_boundary = file[spectrum]['precursorMz'][0]['precursorMz'] - file[spectrum]['isolation window lower offset'] - tolerance_calculated
                            upper_boundary = file[spectrum]['precursorMz'][0]['precursorMz'] + file[spectrum]['isolation window upper offset'] + tolerance_calculated
                            
                            if lower_boundary <= target_mz <= upper_boundary:
                                found_matching_mz = True
                                break
                        else:
                            if abs((file[spectrum]['precursorMz'][0]['precursorMz']) - target_mz) <= tolerance_calculated:
                                found_matching_mz = True
                                break
                    
                    # If the MS2 spectrum percursor mz matches the mz of the adduct, analyze the spectrum
                    if found_matching_mz:
                        
                        # The array summed intensities information to calculate % TIC assigned
                        total_array_intensity = sum(file[spectrum]['intensity array'])
                        
                        # The last peak information
                        former_peak_mz = 0
                        former_peak_intensity = 0
                        former_peak_identified_mz = 0
                        
                        # The maximum intensity of the array and intensity cutoff base for array trimming
                        max_int = max(file[spectrum]['intensity array'])
                        intensity_cutoff = 0
                        
                        # Trim extremely large arrays
                        max_size = 10000
                        if len(file[spectrum]['intensity array']) > max_size:
                            
                            while np.sum(file[spectrum]['intensity array'] > intensity_cutoff) > max_size:
                                intensity_cutoff += max_int*0.01
                                
                            mask = file[spectrum]['intensity array'] > intensity_cutoff
                            mz_array = file[spectrum]['m/z array'][mask]
                            int_array = file[spectrum]['intensity array'][mask]
                            
                        else:
                            mz_array = file[spectrum]['m/z array']
                            int_array = file[spectrum]['intensity array']
                            
                        # Check the spectrum peak by peak
                        for mz_peak_index, mz_peak in enumerate(mz_array):
                            
                            # Moving intensity threshold to ignore minuscule peaks that are between isotopologues
                            if int_array[mz_peak_index] < former_peak_intensity*0.05:
                                continue
                                
                            # Check whether this mz peak is part of the isotopic envelope of the last checked peak. Checks for possibly singly, doubly or triply charged
                            if (
                                abs(mz_peak-(former_peak_mz+General_Functions.h_mass)) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak) 
                                or abs(mz_peak-(former_peak_mz+(General_Functions.h_mass/2))) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak) 
                                or abs(mz_peak-(former_peak_mz+(General_Functions.h_mass/3))) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak)
                                ): 
                                # And here it checks if it is part of the isotopic envelope of the last IDENTIFIED peak. If it is, it reduces its intensity from the total array intensity so to better represent the % TIC annotated
                                if (
                                    abs(mz_peak-(former_peak_identified_mz+General_Functions.h_mass)) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak) 
                                    or abs(mz_peak-(former_peak_identified_mz+(General_Functions.h_mass/2))) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak) 
                                    or abs(mz_peak-(former_peak_identified_mz+(General_Functions.h_mass/3))) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak)
                                    ):
                                    former_peak_identified_mz = mz_peak
                                    total_array_intensity -= int_array[mz_peak_index]
                                    
                                former_peak_mz = mz_peak
                                continue
                                
                            former_peak_mz = mz_peak
                            former_peak_intensity = int_array[mz_peak_index]
                            
                            # Checks whether this peak matches the mz of a fragment
                            fragment_id = General_Functions.binary_search_with_tolerance(indexed_fragments_list, mz_peak, 0, len(indexed_fragments_list)-1, General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak))
                            
                            # If no matching fragment is found, continue to the next mz peak
                            if fragment_id == -1:
                                continue
                                
                            # Identify the fragment mz
                            identified_fragment_mz = indexed_fragments_list[fragment_id]
                            
                            # Otherwise, fetch the list of fragments that match the mz and, if the user wants to limit the output to the composition of precursor, filter it
                            if not filter_output:
                                possible_fragments = indexed_fragments[identified_fragment_mz]
                            else:
                                possible_fragments = [fragment for fragment in indexed_fragments[identified_fragment_mz] if (fragment.split("_")[0] in fragments_per_glycan[i] and set(adduct_combos_dict[fragment.split("_")[1]]['comp'].keys()).issubset(set(adduct_name.keys())))]
                            
                            # If after filtering, no possible fragments are found, move on to the next peak
                            if len(possible_fragments) == 0:
                                continue
                                
                            # Mark this as a successful identification and save the mz as the last identified
                            former_peak_identified_mz = mz_peak 
                            
                            # Format the fragments name to look prettier
                            fragment_name_list = []
                            for fragment in possible_fragments:
                                fragment_name, fragment_adduct = fragment.split("_")
                                
                                fragments_ranking[fragment_name] = fragments_ranking.get(fragment_name, 0) + 1
                                
                                adduct_comp_frag = adduct_combos_dict[fragment_adduct]['comp']
                                adduct_charge_frag = adduct_combos_dict[fragment_adduct]['charges']
                                
                                adduct_str = ""
                                for o in adduct_comp_frag:
                                    polarity = '+' if adduct_comp_frag[o] > 0 else ''
                                    adduct_str += f"{polarity}{adduct_comp_frag[o]}{o}"
                                    
                                superscript_polarity = superscripts['+'] if adduct_charge_frag > 0 else superscripts['-']
                                
                                fragment_name_list.append(f"{fragment_name}[M{adduct_str}]{superscript_polarity}{superscripts[str(abs(adduct_charge_frag))]}")
                                
                            fragment_name = "/".join(fragment_name_list)
                            
                            glycan_fragments_data[adduct][file_index].append([i, adduct, fragment_name, mz_peak, int_array[mz_peak_index], file[spectrum]['retentionTime'], file[spectrum]['precursorMz'][0]['precursorMz'], total_array_intensity])
                            
                        # Update the total intensities after it's done with the array
                        for fragment in glycan_fragments_data[adduct][file_index]:
                            if fragment[5] == file[spectrum]['retentionTime']:
                                fragment[7] = total_array_intensity
                        
                        analyzed_spectra[file_index][spectrum] = True
                       
        return glycan_fragments_data, i, analyzed_spectra, fragments_ranking
        
    except KeyboardInterrupt:
        if not from_GUI:
            print("\n\n----------Execution cancelled by user.----------\n", flush=True)
        raise SystemExit(1)

def score_remaining_spectra(file,
                            spectra_to_check,
                            data,
                            top_10_fragments,
                            fragments_dict,
                            tolerance):
    '''
    '''
    spectra_score = {}
    
    # Go through each spectrum
    for spectrum_number in spectra_to_check:
        # Extract the spectrum data
        mz_array = data[file][spectrum_number]['m/z array']
        int_array = data[file][spectrum_number]['intensity array']
        
        # Fragment match counter for this spectrum
        match_counter = 0
        
        # Cycle through the top 10 fragments
        for fragment_name in top_10_fragments:
            # Acquire fragment information
            fragment_data = fragments_dict[fragment_name]
            
            # Cycle through the fragments adduct mz
            for adduct_name, mz in fragment_data['Adducts_mz'].items():
                # Binary search the fragment within the mz array
                fragment_id = General_Functions.binary_search_with_tolerance(mz_array, mz, 0, len(mz_array)-1, General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz), int_arr = int_array)
                
                # Check if it was found and increment the score, move to the next fragment
                if fragment_id != -1:
                    match_counter += 1
                    break
        
        spectra_score[spectrum_number] = match_counter
        
    return file, spectra_score