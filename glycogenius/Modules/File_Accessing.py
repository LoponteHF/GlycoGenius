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

from . import General_Functions
from pyteomics import mzxml, mzml, mass, auxiliary
from itertools import combinations_with_replacement
from scipy.sparse.linalg import splu
from scipy import sparse
from statistics import mean
from re import split
from math import inf, atan, pi, exp, sqrt
import concurrent.futures
import pathlib
import importlib
import numpy
import sys
import datetime
import copy
import os

##---------------------------------------------------------------------------------------
##File accessing-associated functions (these are functions that deal with the data in
##some way. They make vast use of general functions and require a library to work).

class make_mzxml(object):
    '''A wrapper that takes the output of pyteomics mzML parser and converts it to
    the mzXML pyteomics parser standard to be used within the script. Allows for full
    support of mzML.
    
    Parameters
    ----------
    path_to_mzml : string
        A string containing the path to the mzML file.
        
    Uses
    ----
    pyteomics.mzml.MzML : generator
        Indexes the mzML file into a generator, allowing you to parse the file for
        analysis.
        
    Returns
    -------
    dict
        If class is queried for a single index or iterated, returns the dictionary of the
        converted output (mzML -> mzXML)
        
    data : list
        If class is queried for slices, returns a list of dictionaries containing
        the converted output (mzML -> mzXML)
    '''
    def __init__(self,it):
        self.data = mzml.MzML(it)
    def __iter__(self):
        return self.make_mzxml_iterator(self.data)
    def __getitem__(self,index):
        if type(index) == int:
            pre_data = self.data[index]
            data_to_return = None
            if pre_data['ms level'] == 2:
                if float(self.data[-1]['scanList']['scan'][0]['scan start time']) > 300: #300 scan time should allow for the correct evaluation of scan time being in seconds or minutes for every run that lasts between 5 minutes and 5 hours
                    data_to_return = {'num': pre_data['id'].split('=')[-1], 'retentionTime': float(pre_data['scanList']['scan'][0]['scan start time'])/60, 'msLevel': pre_data['ms level'], 'm/z array': pre_data['m/z array'], 'intensity array': pre_data['intensity array']}
                else:
                    data_to_return = {'num': pre_data['id'].split('=')[-1], 'retentionTime': float(pre_data['scanList']['scan'][0]['scan start time']), 'msLevel': pre_data['ms level'], 'm/z array': pre_data['m/z array'], 'intensity array': pre_data['intensity array']}
                        
                # Some MzML have isolation window, some don't
                if 'isolationWindow' in pre_data['precursorList']['precursor'][0] and 'isolation window target m/z' in pre_data['precursorList']['precursor'][0]['isolationWindow'] and 'isolation window lower offset' in pre_data['precursorList']['precursor'][0]['isolationWindow'] and 'isolation window upper offset' in pre_data['precursorList']['precursor'][0]['isolationWindow']:
                    data_to_return['precursorMz'] = [{'precursorMz': pre_data['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']}]
                    data_to_return['isolation window lower offset'] = pre_data['precursorList']['precursor'][0]['isolationWindow']['isolation window lower offset']
                    data_to_return['isolation window upper offset'] = pre_data['precursorList']['precursor'][0]['isolationWindow']['isolation window upper offset']
                else:
                    data_to_return['precursorMz'] = [{'precursorMz': pre_data['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']}]
            else:
                if float(self.data[-1]['scanList']['scan'][0]['scan start time']) > 300:
                    data_to_return = {'num': pre_data['id'].split('=')[-1], 'retentionTime': float(pre_data['scanList']['scan'][0]['scan start time'])/60, 'msLevel': pre_data['ms level'], 'm/z array': pre_data['m/z array'], 'intensity array': pre_data['intensity array']}
                else:
                    data_to_return = {'num': pre_data['id'].split('=')[-1], 'retentionTime': float(pre_data['scanList']['scan'][0]['scan start time']), 'msLevel': pre_data['ms level'], 'm/z array': pre_data['m/z array'], 'intensity array': pre_data['intensity array']}
                    
            return data_to_return
        else:
            first_index = str(index).split('(')[1].split(', ')[0]
            if first_index != 'None':
                first_index = int(first_index)
            else:
                first_index = 0
            last_index = str(index).split('(')[1].split(', ')[1]
            if last_index != 'None':
                last_index = int(last_index)
            else:
                last_index = len(self.data)
            data = []
            for index in range(first_index, last_index):
                if self.data[index]['ms level'] == 2:
                    if float(self.data[-1]['scanList']['scan'][0]['scan start time']) > 300:
                        data.append({'num': self.data[index]['id'].split('=')[-1], 'retentionTime': float(self.data[index]['scanList']['scan'][0]['scan start time'])/60, 'msLevel': self.data[index]['ms level'], 'm/z array': self.data[index]['m/z array'], 'intensity array': self.data[index]['intensity array']})
                    else:
                        data.append({'num': self.data[index]['id'].split('=')[-1], 'retentionTime': float(self.data[index]['scanList']['scan'][0]['scan start time']), 'msLevel': self.data[index]['ms level'], 'm/z array': self.data[index]['m/z array'], 'intensity array': self.data[index]['intensity array']})
                        
                    # Some MzML have isolation window, some don't
                    if 'isolationWindow' in pre_data['precursorList']['precursor'][0] and 'isolation window target m/z' in pre_data['precursorList']['precursor'][0]['isolationWindow'] and 'isolation window lower offset' in pre_data['precursorList']['precursor'][0]['isolationWindow'] and 'isolation window upper offset' in pre_data['precursorList']['precursor'][0]['isolationWindow']:
                        data[-1]['precursorMz'] = [{'precursorMz': self.data[index]['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']}]
                        data[-1]['isolation window lower offset'] = self.data[index]['precursorList']['precursor'][0]['isolationWindow']['isolation window lower offset']
                        data[-1]['isolation window upper offset'] = self.data[index]['precursorList']['precursor'][0]['isolationWindow']['isolation window upper offset']
                    else:
                        data[-1]['precursorMz'] = [{'precursorMz': self.data[index]['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']}]
                else:
                    if float(self.data[-1]['scanList']['scan'][0]['scan start time']) > 300:
                        data.append({'num': self.data[index]['id'].split('=')[-1], 'retentionTime': float(self.data[index]['scanList']['scan'][0]['scan start time'])/60, 'msLevel': self.data[index]['ms level'], 'm/z array': self.data[index]['m/z array'], 'intensity array': self.data[index]['intensity array']})
                    else:
                        data.append({'num': self.data[index]['id'].split('=')[-1], 'retentionTime': float(self.data[index]['scanList']['scan'][0]['scan start time']), 'msLevel': self.data[index]['ms level'], 'm/z array': self.data[index]['m/z array'], 'intensity array': self.data[index]['intensity array']})
            return data
            
    class make_mzxml_iterator:
        def __init__(self, data):
            self.data = data
            self.index = 0

        def __iter__(self):
            return self

        def __next__(self):
            if self.index < len(self.data):
                pre_data = self.data[self.index]
                self.index += 1
                if pre_data['ms level'] == 2:
                    if float(self.data[-1]['scanList']['scan'][0]['scan start time']) > 300: #300 scan time should allow for the correct evaluation of scan time being in seconds or minutes for every run that lasts between 5 minutes and 5 hours
                        return {'num': pre_data['id'].split('=')[-1], 'retentionTime': float(pre_data['scanList']['scan'][0]['scan start time'])/60, 'msLevel': pre_data['ms level'], 'm/z array': pre_data['m/z array'], 'intensity array': pre_data['intensity array'], 'precursorMz': [{'precursorMz': pre_data['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']}]}
                    else:
                        return {'num': pre_data['id'].split('=')[-1], 'retentionTime': float(pre_data['scanList']['scan'][0]['scan start time']), 'msLevel': pre_data['ms level'], 'm/z array': pre_data['m/z array'], 'intensity array': pre_data['intensity array'], 'precursorMz': [{'precursorMz': pre_data['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']}]}
                else:
                    if float(self.data[-1]['scanList']['scan'][0]['scan start time']) > 300:
                        return {'num': pre_data['id'].split('=')[-1], 'retentionTime': float(pre_data['scanList']['scan'][0]['scan start time'])/60, 'msLevel': pre_data['ms level'], 'm/z array': pre_data['m/z array'], 'intensity array': pre_data['intensity array']}
                    else:
                        return {'num': pre_data['id'].split('=')[-1], 'retentionTime': float(pre_data['scanList']['scan'][0]['scan start time']), 'msLevel': pre_data['ms level'], 'm/z array': pre_data['m/z array'], 'intensity array': pre_data['intensity array']}
            else:
                raise StopIteration
       
def eic_from_glycan(files,
                    glycan,
                    glycan_info,
                    ms1_indexes,
                    rt_interval,
                    tolerance,
                    min_isotops,
                    noise,
                    avg_noise,
                    max_charges,
                    zeroes_arrays,
                    inf_arrays,
                    threads_arrays,
                    rt_arrays,
                    ms1_id,
                    sampling_rates):
    '''Generates a very processed EIC for each adduct of each glycan of each sample.
    Removes non-monoisotopic peaks, check charges and calculates multiple quality scoring
    data.

    Parameters
    ----------
    files : list
        A list of generators obtained from the function pyteomics.mzxml.MzXML().
        
    glycan : str
        The name of the target glycan.

    glycan_info : dict
        A dictionary containing multiple glycan info, including the needed 'Adducts_mz'
        key required for this function. Generated by full_glycans_library().

    ms1_indexes : dict
        A dictionary where each key is the ID of a file and each value is a list
        containing the indexes of all MS1 scans in the file.
        
    rt_interval : tuple
        A tuple where the first index contains the beggining time of the retention time
        interval you wish to analyze and the second contains the end time.
        
    tolerance : tuple
        First index contains the unit of the tolerance and the second one is the value of 
        that unit.
        
    min_isotops : int
        The minimum amount of isotopologues required to consider an RT mz peak valid.
        
    noise : list
        A list containing the calculated noise level of each sample.
    
    avg_noise : float
        Fallback average noise in case local noise can't be calculated.
        
    max_charges : int
        The maximum amount of charges the queried mz should have.
        
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
        
    rewind : int
        The maximum number of spectras it can go back to for retroactively checking after finding a peak.

    Uses
    ----
    General_Functions.form_to_charge : int
        Converts adducts formula into raw charge.
        
    pyteomics.mass.calculate_mass(*args, **kwargs) : float
        Calculates the monoisotopic mass of a polypeptide defined by a sequence string,
        parsed sequence, chemical formula or Composition object.
        
    General_Functions.form_to_comp(string) : dict
        Separates a molecular formula or monosaccharides composition of glycans into a
        dictionary with each atom/monosaccharide as a key and its amount as value.
        
    General_Functions.h_mass : float
        The mass of Hydrogen-1.
        
    General_Functions.calculate_ppm_diff : float
        Calculates the PPM difference between a mz and a target mz.
        
    numpy.corrcoef : matrix
        Return Pearson product-moment correlation coefficients.

    Returns
    -------
    data : dict
        A dictionary with keys for each adduct combo of a given glycan, which value is a
        dict with keys for each file id, at which value is a list of lists, with the
        first one being the rt array and the second one the processed int array.
        
    ppm_info : dict
        A dictionary with keys for each adduct combo of a given glycan, which value is a
        dict with keys for each file id, at which each value corresponds to a given 
        calculated ppm difference for the retention time at the same index.
        
    iso_fitting_quality : dict
        A dictionary with keys for each adduct combo of a given glycan, which value is a
        dict with keys for each file id, at which each value corresponds to a given 
        calculated isotopic fitting score for the retention time at the same index.
        
    verbose_info : list
        Legacy list used to save verbose info to. Unused and in need of CLEANUP.
        
    raw_data : dict
        A dictionary with keys for each adduct combo of a given glycan, which value is a
        dict with keys for each file id, at which value is a list of lists, with the
        first one being the rt array and the second one the raw int array.
        
    isotopic_fits : dict
        A dictionary containing all the isotopic fitting information for reporting and 
        checking data.
    '''
    global buffer, buffer_good
    data = {}
    ppm_info = {}
    iso_fitting_quality = {}
    isotopic_fits = {}
    verbose_info = []
    raw_data = {}
    for i in glycan_info['Adducts_mz']:
        
        adduct_comp, adduct_charge = General_Functions.fix_adduct_determine_charge(i)
        
        adduct_mass = mass.calculate_mass(composition=adduct_comp)
        ppm_info[i] = {}
        iso_fitting_quality[i] = {}
        data[i] = {}
        raw_data[i] = {}
        isotopic_fits[i] = {}
        for j_j, j in enumerate(files):
            ppm_info[i][j_j] = copy.deepcopy(inf_arrays[j_j])
            iso_fitting_quality[i][j_j] = copy.deepcopy(zeroes_arrays[j_j])
            data[i][j_j] = [copy.deepcopy(rt_arrays[j_j]), copy.deepcopy(zeroes_arrays[j_j])]
            raw_data[i][j_j] = [copy.deepcopy(rt_arrays[j_j]), copy.deepcopy(zeroes_arrays[j_j])]
            isotopic_fits[i][j_j] = {}
            thread_numbers = threads_arrays[j_j]
            
            #checked possibility of parallelization here, too much overhead (over 10 more time to run, even if using 1 core)
            buffer = []
            buffer_good = 0
            for k_k, k in enumerate(thread_numbers):
                analyze_mz_array(j[k]['m/z array'],
                                 j[k]['intensity array'],
                                 glycan_info,
                                 tolerance,
                                 min_isotops,
                                 noise,
                                 avg_noise,
                                 max_charges,
                                 ppm_info,
                                 iso_fitting_quality,
                                 data,
                                 raw_data,
                                 isotopic_fits,
                                 i,
                                 j_j,
                                 k,
                                 j[k]['retentionTime'],
                                 ms1_id[j_j][k_k],
                                 ms1_id[-1],
                                 adduct_mass,
                                 adduct_charge,
                                 sampling_rates)
                if len(buffer) > max([4*sampling_rates[j_j], 4]):
                    rewind = False
                    found_count = 0
                    for l in range(-1, -max([(4*sampling_rates[j_j]), 4])-2, -1):
                        if buffer[l] != None:
                            found_count += 1
                        else:
                            if found_count >= 2:
                                rewind = True
                                break
                            else:
                                break
                    if rewind:
                        for l_l in range(-found_count-1, -len(buffer)-1, -1):
                            analyze_mz_array(j[thread_numbers[k_k+l_l+1]]['m/z array'],
                                             j[thread_numbers[k_k+l_l+1]]['intensity array'],
                                             glycan_info,
                                             tolerance,
                                             min_isotops,
                                             noise,
                                             avg_noise,
                                             max_charges,
                                             ppm_info,
                                             iso_fitting_quality,
                                             data,
                                             raw_data,
                                             isotopic_fits,
                                             i,
                                             j_j,
                                             thread_numbers[k_k+l_l+1],
                                             j[thread_numbers[k_k+l_l+1]]['retentionTime'],
                                             ms1_id[j_j][k_k+l_l+1],
                                             ms1_id[-1],
                                             adduct_mass,
                                             adduct_charge,
                                             sampling_rates,
                                             retest = True,
                                             retest_no = l_l)
                            if buffer[l_l] == None:
                                break
    return data, ppm_info, iso_fitting_quality, verbose_info, raw_data, isotopic_fits

    
def analyze_mz_array(sliced_mz,
                     sliced_int,
                     glycan_info,
                     tolerance,
                     min_isotops,
                     noise,
                     avg_noise,
                     max_charges,
                     ppm_info,
                     iso_fitting_quality,
                     data,
                     raw_data,
                     isotopic_fits,
                     glycan_id,
                     file_id,
                     thread_id,
                     ret_time,
                     ms1_id,
                     last_ms1_id,
                     adduct_mass,
                     adduct_charge,
                     sampling_rates,
                     filtered = True,
                     retest = False,
                     retest_no = 0):
    '''The core function of eic_from_glycan. Analyzes a single spectrum and outputs
    relevant information, such as the PPM-error of the monoisotopic peak, isotopic fitting
    information, etc.
    
    Parameters
    ----------
    sliced_mz : list
        The mz array of the spectrum.
        
    sliced_int : list
        The intensity array of the spectrum, synchronized with sliced_mz.
        
    glycan_info : dict
        A dictionary containing multiple glycan info, including the needed 'Adducts_mz'
        key required for this function. Generated by full_glycans_library().
        
    tolerance : tuple
        First index contains the unit of the tolerance and the second one is the value of 
        that unit.
        
    min_isotops : int
        The minimum amount of isotopologues required to consider an RT mz peak valid.
        
    noise : list
        A list containing the calculated noise level of each sample.
    
    avg_noise : float
        Fallback average noise in case local noise can't be calculated.
        
    max_charges : int
        The maximum amount of charges the queried mz should have.
        
    data : dict
        A dictionary with keys for each adduct combo of a given glycan, which value is a
        dict with keys for each file id, at which value is a list of lists, with the
        first one being the rt array and the second one the processed int array.
        
    ppm_info : dict
        A dictionary with keys for each adduct combo of a given glycan, which value is a
        dict with keys for each file id, at which each value corresponds to a given 
        calculated ppm difference for the retention time at the same index.
        
    iso_fitting_quality : dict
        A dictionary with keys for each adduct combo of a given glycan, which value is a
        dict with keys for each file id, at which each value corresponds to a given 
        calculated isotopic fitting score for the retention time at the same index.
        
    raw_data : dict
        A dictionary with keys for each adduct combo of a given glycan, which value is a
        dict with keys for each file id, at which value is a list of lists, with the
        first one being the rt array and the second one the raw int array.
        
    isotopic_fits : dict
        A dictionary containing all the isotopic fitting information for reporting and 
        checking data.
        
    glycan_id : str
        The adduct of the glycan ie. H1, Na1, etc.
        
    file_id : int
        The number of the sample file.
        
    thread_id : int
        The id of the spectrum to be analyzed.
        
    ret_time : float
        The retention time correspondent to the thread_id.
        
    ms1_id : int
        The ID of the MS1 spectra (important when you also have MS2 in file, otherwise
        matches thread_id).
        
    adduct_mass : float
        The mass of the adduct (ie. 1.0074 for H1 adduct)
        
    adduct_charge : int
        The charge of the adduct.
        
    filtered : boolean
        Whether the output will be filtered or not. Non-filtered output allows for MUCH FASTER tracing. Might get used in the future.
        
    retest : boolean
        Whether this is a retest check or not. For use with rewind on eic_from_glycan.
        
    retest_no : int
        Buffer ID of the retest.
    
    Uses
    ----
    General_Functions.h_mass : float
        The mass of a hydrogen atom.
        
    General_Functions.tolerance_calc : float
        An accurate way to convert 'ppm' mass accuracy into 'pw' (peak width, aka. mz
        tolerance).
    
    General_Functions.local_noise_calc : float
        Uses the noise_specs produced by rt_noise_level_parameters_set to 
        calculate the local noise levels. If any of the parameters are considered
        abnormal (ie. absurdly high or no peaks on the array) it defaults to avg_noise.
        
    General_Functions.calculate_ppm_diff : float
        Calculates the PPM difference between a mz and a target mz.
        
    numpy.linalg.norm : float or ndarray
        Matrix or vector normalization.
    
    numpy.dot : ndarray
        Dot product of two arrays.
        
    Returns
    -------
    nothing
        Edits the target dictionaries/lists directly.
    '''
    global buffer, buffer_good
    target_mz = glycan_info['Adducts_mz'][glycan_id]
    local_noise = General_Functions.local_noise_calc(noise[file_id][ms1_id], target_mz, avg_noise[file_id])
    sliced_mz_length = len(sliced_mz)-1
    # print(f"Analyzing {ret_time}... retest? {retest}")
    if target_mz > sliced_mz[-1]:
        mz_id = -1
        # print(f"Target mz {target_mz} outside mz range")
    else:
        mz_id = General_Functions.binary_search_with_tolerance(sliced_mz, target_mz, 0, sliced_mz_length, General_Functions.tolerance_calc(tolerance[0], tolerance[1], target_mz), sliced_int)
        # if mz_id == -1:
            # print(f"Target not found in this retention time")
    
    # print(f"Target found... Local noise: {local_noise}, Threshold for detection: {local_noise*0.5}, Intensity: {sliced_int[mz_id]}")
    if mz_id != -1 and sliced_int[mz_id] >= local_noise*0.5:
        found_mz = sliced_mz[mz_id]
        charge_range = range(1, max(4, abs(adduct_charge)*2))
        intensity = sliced_int[mz_id] #variable to sum the total deisotopotized intensity
        mono_int = sliced_int[mz_id] #variable to sum the total intensity of the monoisotopic peak
        raw_data[glycan_id][file_id][1][ms1_id] = mono_int #unfiltered EIC
        ppm_error = General_Functions.calculate_ppm_diff(sliced_mz[mz_id], target_mz)
        
        if filtered == True:
            iso_actual = [1]
            bad = False #here starts quality checks
            margin = 0.2 #0.2 = 20%, 0.4 = 40% - margin for checking, higher the margin, more strict is the checking and less glycans will probably be found
            
            temp_id = General_Functions.binary_search_with_tolerance(sliced_mz, found_mz+(General_Functions.h_mass/abs(adduct_charge)), mz_id, sliced_mz_length, General_Functions.tolerance_calc(tolerance[0], tolerance[1], found_mz+(General_Functions.h_mass/abs(adduct_charge))), sliced_int) #check if second isotopic is actually present or not
            if temp_id == -1:
                # print(f"Second isotopic peak not found")
                bad = True
                
            if not bad:
                # print(f"Checking if target is monoisotopic...")
                for i in charge_range: #check if it's monoisotopic and correct charge
                    # print(f"Testing charge {i}")
                    if (len(buffer) <= round(max([2*sampling_rates[file_id], 2])) or (len(buffer) > round(max([2*sampling_rates[file_id], 2])) and buffer[round(-max([2*sampling_rates[file_id], 2]))-1] == None)) and not retest: #only check if it's monoisotopic if at least one of the last 3 RT got nothing...
                        temp_id = General_Functions.binary_search_with_tolerance(sliced_mz, found_mz-(General_Functions.h_mass/i), 0, mz_id, General_Functions.tolerance_calc(tolerance[0], tolerance[1], found_mz-(General_Functions.h_mass/i)), sliced_int) #check monoisotopic
                        if temp_id != -1 and sliced_int[temp_id] > 0:
                            expected_value = (sliced_mz[temp_id]*i*0.0006)+0.1401 #based on linear regression of the relationship between masses and the second isotopic peak relative intensity of the average of different organic macromolecules
                            # print(f"Found possible monoisotopic peak: Expected value greater than: {expected_value*(1+(margin*2))}, theoretical mass of monoisotopic: {sliced_mz[temp_id]*i}, charges: {i}, second isotopic actual: {mono_int/sliced_int[temp_id]}")
                            if (mono_int/sliced_int[temp_id] < expected_value*(1+(margin*2))):
                                # print(f"Target not monoisotopic")
                                bad = True
                                break
                            # print(f"Putative monoisotopic peak discarded")
                        # else:
                            # print(f"No candidate monoisotopic found for this charge")
                    # else:
                        # print(f"Check not necessary")
                    
                    # print(f"Monoisotopic check passed. Checking if it's correctly charged isotopic envelope...")
                    if i == 1 or i == abs(adduct_charge) or (i == 2 and abs(adduct_charge) == 4) or (i == 3 and abs(adduct_charge) == 6): #ignores charge 1 due to the fact that any charge distribution will find a hit on that one
                        # print(f"Check not necessary for this charge")
                        continue
                    temp_id = General_Functions.binary_search_with_tolerance(sliced_mz, found_mz+(General_Functions.h_mass/i), mz_id, sliced_mz_length, General_Functions.tolerance_calc(tolerance[0], tolerance[1], found_mz+(General_Functions.h_mass/i)), sliced_int) #check for correct charge
                    if temp_id != -1:
                        expected_value = (target_mz*i*0.0006)+0.1401 #based on linear regression of the relationship between masses and the second isotopic peak relative intensity of the average of different organic macromolecules
                        # print(f"Found possible second peak: Expected value smaller than: {expected_value*(1-margin)}, theoretical mass of monoisotopic: {target_mz*i}, charges: {i}, second isotopic actual: {sliced_int[temp_id]/mono_int}")
                        if (sliced_int[temp_id]/mono_int > expected_value*(1-margin)):
                            # print(f"Target not correctly charged")
                            bad = True
                            break
                        # print("Putative second peak of isotopic envelope discarded")
                    # else:
                        # print(f"No candidate second peak found for this charge")
            if not bad:
                # print("Checks passed! Checking isotopic envelope")
                isos_found = 0
                mz_isos = []
                for i_i, i in enumerate(glycan_info['Isotopic_Distribution_Masses']): #check isotopic peaks and add to the intensity
                    if i_i == 0: #ignores monoisotopic this time around
                        continue
                    # print(f"Looking for isotopic peak no. {i_i+1}, mz {found_mz+(i_i*(General_Functions.h_mass/abs(adduct_charge)))}")
                    temp_id = General_Functions.binary_search_with_tolerance(sliced_mz, found_mz+(i_i*(General_Functions.h_mass/abs(adduct_charge))), mz_id, sliced_mz_length, General_Functions.tolerance_calc(tolerance[0], tolerance[1], found_mz+(i_i*(General_Functions.h_mass/abs(adduct_charge))))*2, sliced_int, mz_isos)
                    if temp_id != -1 and sliced_int[temp_id] > 0:
                        # print(f"Found! Intensity {sliced_int[temp_id]/mono_int}")
                        isos_found += 1
                        mz_isos.append(sliced_mz[temp_id])
                        iso_actual.append(sliced_int[temp_id]/mono_int)
                        if sliced_int[temp_id] > mono_int*glycan_info['Isotopic_Distribution'][i_i]:
                            intensity += mono_int*glycan_info['Isotopic_Distribution'][i_i]
                        if sliced_int[temp_id] <= mono_int*glycan_info['Isotopic_Distribution'][i_i]:
                            intensity += sliced_int[temp_id]
                    else:
                        # print(f"Not found...")
                        if isos_found == 0: #a compound needs at least 2 identifiable peaks (monoisotopic + 1 from isotopic envelope)
                            bad = True
                            # print("Not enough isotopic envelope peaks found, discarded this RT.")
                            break
                        else:
                            break
            
            if not bad and (iso_actual[1] < 0.2 or iso_actual[1] > 5): #this should avoid situations where it's obvious that it's picking the wrong charge because the second peak is almost invisible compared to the third and first, which when z=2 means that it's very likely actually a singly charge compound, for example
                # print(f"Last tests on isotopic envelope...")
                if len(iso_actual) > 2:
                    if iso_actual[2] > iso_actual[1]*10 or iso_actual[1] < glycan_info['Isotopic_Distribution'][1]*0.5 or iso_actual[1] > glycan_info['Isotopic_Distribution'][1]*2:
                        # print(f"Failed last tests on isotopic envelope: Third isotopic peak much bigger than second or second isotopic peak too big or too small... Second isotopic peak intensity: {iso_actual[1]}, Expected: {glycan_info['Isotopic_Distribution'][1]}")
                        bad = True
                else:
                    if iso_actual[1] < glycan_info['Isotopic_Distribution'][1]*0.5 or iso_actual[1] > glycan_info['Isotopic_Distribution'][1]*2: #smallest glycan should have the second iso_actual somewhere around 0.5, so a cutoff lower than that is fine
                        # print(f"Failed last tests on isotopic envelope: Second isotopic peak too big or too small... Second isotopic peak intensity: {iso_actual[1]}, Expected: {glycan_info['Isotopic_Distribution'][1]}")
                        bad = True
            
            if bad:
                # print(f"Failed. RT set to zero.\n")
                if not retest:
                    buffer.append(None)
            else:
                # print(f"Passed! RT intensity saved to buffer!\n")
                iso_target = glycan_info['Isotopic_Distribution'][:isos_found+1]
                
                ratios = []
                weights = []
                for i_i, i in enumerate(iso_target):
                    if i_i == 0:
                        continue
                    intensities = [i, iso_actual[i_i]]
                    ratio = min(intensities)/max(intensities)
                    
                    #scales the score in a sigmoid, with steepness determine by k_value
                    k_value = 10
                    corrected_ratio = 1 / (1 + numpy.exp(-k_value * (ratio - 0.5)))
                    
                    ratios.append(corrected_ratio)
                    weights.append(1/(exp(1.25*i_i)))
                
                iso_quali = numpy.average(ratios, weights = weights)
            
                #reduces score if fewer isotopic peaks are found: punishing for only 1 peaks, normal score from 2 and over (besides the monoisotopic)
                if len(iso_actual) == 2:
                    iso_quali = (iso_quali*0.8)
                
                mz_isos = [found_mz]+mz_isos
                
                if not retest:
                    buffer.append(([glycan_id, file_id, ms1_id, float("%.4f" % round(ret_time, 4))], [ppm_error, iso_quali, intensity, [mz_isos, iso_target, iso_actual, iso_quali]]))
                else:
                    buffer[retest_no] = ([glycan_id, file_id, ms1_id, float("%.4f" % round(ret_time, 4))], [ppm_error, iso_quali, intensity, [mz_isos, iso_target, iso_actual, iso_quali]])
                
    else:
        if filtered == True:
            if not retest:
                buffer.append(None)
        else:
            info = ([glycan_id, file_id, ms1_id, float("%.4f" % round(ret_time, 4))], [inf, 1.0, 0.0, [[], [], [], 1.0]])
            isotopic_fits[info[0][0]][info[0][1]][info[0][3]] = info[1][3]
    
    # print(f"Buffer before clean-up: {buffer}\n")
    
    #dynamical clean-up of buffer
    min_in_a_row = round(max([4*sampling_rates[file_id], 4]))
    buffer_size = round(max([10*sampling_rates[file_id], 10]))
    if filtered:
        in_a_row = 0
        for i_i, i in enumerate(buffer):
            if i != None:
                in_a_row += 1
            else:
                in_a_row = 0
                for k_k, k in enumerate(buffer): #this clears the buffer up to the point where it was checked, only once "None" is found
                    if k_k >= i_i:
                        break
                    buffer[k_k] = None
        if in_a_row >= min_in_a_row:
            for i in buffer:
                if i != None:
                    ppm_info[i[0][0]][i[0][1]][i[0][2]] = i[1][0]
                    iso_fitting_quality[i[0][0]][i[0][1]][i[0][2]] = i[1][1]
                    data[i[0][0]][i[0][1]][1][i[0][2]] = i[1][2]
                    isotopic_fits[i[0][0]][i[0][1]][i[0][3]] = i[1][3]
        if len(buffer) >= buffer_size: #this means that the buffer will be worked on once it gets to the buffer_size or over it or the end of the MS1 array is reached
            buffer.pop(0)
            
    elif not filtered:
        if not retest:
            info = ([glycan_id, file_id, ms1_id, float("%.4f" % round(ret_time, 4))], [ppm_error, 1.0, mono_int, [[], [], [], 1.0]])
            ppm_info[info[0][0]][info[0][1]][info[0][2]] = info[1][0]
            iso_fitting_quality[info[0][0]][info[0][1]][info[0][2]] = info[1][1]
            data[info[0][0]][info[0][1]][1][info[0][2]] = info[1][2]
            isotopic_fits[info[0][0]][info[0][1]][info[0][3]] = info[1][3]
        
    # print(f"Buffer after clean-up: {buffer}")
    
def eic_smoothing(y, lmbd = 100, d = 2):
    '''Implementation of the Whittaker smoothing algorithm,
    based on the work by Eilers [1].

    [1] P. H. C. Eilers, "A perfect smoother", Anal. Chem. 2003, (75), 3631-3636
    
    The larger 'lmbd', the smoother the data.
    For smoothing of a complete data series, sampled at equal intervals

    This implementation uses sparse matrices enabling high-speed processing
    of large input vectors
    
    Parameters
    ----------    
    y : list
        Two synchronized lists, one with retention times and second with intensities.
        Second list works as vector containing raw data to be smoothed by this function.
    lmbd : int
        Parameter for the smoothing algorithm (roughness penalty).
    d : int
        Order of the smoothing.

    Returns
    -------
    y[0] : list
        Retention time array
    z : list
        Smoothed intensity array.
    '''
    datapoints_per_min = 1/(y[0][y[1].index(max(y[1]))]-y[0][y[1].index(max(y[1]))-1])
    lmbd = exp(datapoints_per_min/20)
    array = numpy.array(y[1])
    m = len(array)
    E = sparse.eye(m, format='csc')
    D = General_Functions.speyediff(m, d, format='csc')
    coefmat = E + lmbd * D.conj().T.dot(D)
    z = splu(coefmat).solve(array)
    for i_i, i in enumerate(z):
        if i < 0:
            z[i_i] = 0.0
    return y[0], list(z)
    
def peak_curve_fit(rt_int, 
                   peak):
    '''Calculates the fitting between the actual peak and an ideal peak based
    on a calculated gaussian bell curve.
    
    Parameters
    ----------
    rt_int : list
        A list containing two lists: the first one has all the retention times, the
        second one has all the intensities.
        
    peak : dict
        A dictionary containing all sorts of identified peaks information.
        
    Uses
    ----
    General_Functions.normpdf : float
        Calculates the intensity of a gaussian bell curve at the x-axis point x
        with the set parameters.
        
    numpy.corrcoef : matrix
        Return Pearson product-moment correlation coefficients.
        
    Returns
    -------
    fits_list : tuple
        A tuple containing the R_sq of the curve fitting and plotting information
        of the actual and ideal peak curves.
    '''
    x = rt_int[0][peak['peak_interval_id'][0]:peak['peak_interval_id'][1]+1]
    y = rt_int[1][peak['peak_interval_id'][0]:peak['peak_interval_id'][1]+1]
    baseline_correction = min(y[0], y[-1])
    interval = x[len(x)//2]-x[(len(x)//2)-1]
    fits_list = []
    for m in range(1, 11):
        for j in range(len(x)):
            before = []
            after = []
            y_adds = []
            for k in range(j):
                before.append(x[0]-k*interval)
                after.append(x[-1]+k*interval)
                y_adds.append(0.0)
            before = sorted(before)
            temp_x = before+x+after
            temp_y = y_adds+y+y_adds
            max_amp = max(temp_y)
            maximums = []
            for index, i in enumerate(temp_y):
                if i > max_amp*0.8:
                    maximums.append(index)
            max_amp_id = round(sum(maximums)/len(maximums))
            y_gaussian = []
            for i in temp_x:
                y_gaussian.append(General_Functions.normpdf(i, temp_x[max_amp_id], (temp_x[-1]-temp_x[0])/m))
            y_gaussian = [i-min(y_gaussian) for i in y_gaussian]
            scaler = (max_amp/y_gaussian[max_amp_id])
            y_gaussian_scaled = []
            for j in y_gaussian:
                y_gaussian_scaled.append(j*scaler)
            if len(temp_y[len(y_adds):len(temp_y)-len(y_adds)]) <= 2:
                temp_relation = []
                for l in range(len(temp_y[len(y_adds):len(temp_y)-len(y_adds)])):
                    if temp_y[len(y_adds):len(temp_y)-len(y_adds)][l] >= y_gaussian_scaled[len(y_adds):len(temp_y)-len(y_adds)][l]:
                        temp_relation.append(temp_y[len(y_adds):len(temp_y)-len(y_adds)][l]/y_gaussian_scaled[len(y_adds):len(temp_y)-len(y_adds)][l])
                    else:
                        temp_relation.append(y_gaussian_scaled[len(y_adds):len(temp_y)-len(y_adds)][l]/temp_y[len(y_adds):len(temp_y)-len(y_adds)][l])
                R_sq = mean(temp_relation)
            else:
                corr_matrix = numpy.corrcoef(temp_y[len(y_adds):len(temp_y)-len(y_adds)], y_gaussian_scaled[len(y_adds):len(temp_y)-len(y_adds)])
                corr = corr_matrix[0,1]
                R_sq = corr**2
            fits_list.append((R_sq, temp_x[len(y_adds):len(temp_y)-len(y_adds)], temp_y[len(y_adds):len(temp_y)-len(y_adds)], y_gaussian_scaled[len(y_adds):len(temp_y)-len(y_adds)]))
    fits_list = sorted(fits_list, reverse=True)
    return fits_list[0]
    
def iso_fit_score_calc(iso_fits,
                       peak,
                       weights_list):
    '''Calculates the mean isotopic fitting score of a given peak.
    
    Parameters
    ----------
    iso_fits : list
        A list of isotopic fitting scores calculated from eic_from_glycan.
        
    peak : dict
        A dictionary containing all sorts of identified peaks information.
        
    Uses
    ----
    numpy.average : float
        Calculates the weighted average of an array based on another array of weights.
    
    weights_list : list
        A list containing the weights used to calculate the scores. Comes from the curve fitting.
        
    Returns
    -------
    iso_fit_score : float
        Average of the isotopic fits score calculated for the interval of
        the peak of the glycan, weighted by the gaussian fit of it.
    '''
    new_weights = [i**2 for i in weights_list]
    iso_fit_score = numpy.average(iso_fits[peak['peak_interval_id'][0]:peak['peak_interval_id'][1]+1], weights = new_weights)
    return iso_fit_score
    
def average_ppm_calc(ppm_array,
                     tolerance,
                     peak,
                     weights_list):
    '''Calculates the arithmetic mean of the PPM differences of a given peak.
    
    Parameters
    ----------
    ppm_array : list
        A list containing ppm differences.
        
    tolerance : tuple
        First index contains the unit of the tolerance and the second one is the value of 
        that unit and the third index is the mz of the glycan.
        
    peak : dict
        A dictionary containing all sorts of identified peaks information.
    
    weights_list : list
        A list containing the weights used to calculate the scores. Comes from the curve fitting.

    Uses
    ----
    General_Functions.calculate_ppm_diff : float
        Calculates the PPM difference between a mz and a target mz.
        
    numpy.average : float
        Calculates the weighted average of an array based on another array of weights.
        
    Returns
    -------
    tuple
        A tuple containing the mean ppm differences (averaged by the gaussian fit of the peak)
        and the number of missing points.
        The number of missing points in the peak has their ppm difference set to the
        tolerance of the analysis.
    '''
    new_weights = [i**2 for i in weights_list]
    ppms = []
    ppm_default = General_Functions.calculate_ppm_diff(tolerance[2]-General_Functions.tolerance_calc(tolerance[0], tolerance[1], tolerance[2]), tolerance[2])
    missing_points = 0
    for i in ppm_array[peak['peak_interval_id'][0]:peak['peak_interval_id'][1]+1]:
        if i != inf:
            ppms.append(i)
        else:
            ppms.append(ppm_default)
            missing_points+= 1
    regular_mean = numpy.average(ppms, weights = new_weights)
    return regular_mean, missing_points

def peaks_from_eic(rt_int, 
                   rt_int_smoothed,
                   raw_rt_int,
                   rt_interval,
                   min_ppp,
                   close_peaks,
                   glycan):
    '''Does multi peak-picking in a given smoothed EIC.
    
    Parameters
    ----------
    rt_int : list
        A list containing two lists: the first one has all the retention times, the
        second one has all the FILTERD BUT UNSMOOTHED intensities.
        
    rt_int_smoothed : list
        A list containing two lists: the first one has all the retention times, the
        second one has all the SMOOTHED intensities.
        
    raw_rt_int : list
        A list containing two lists: the first one has all the retention times, the
        second one has all the UNFILTERED intensities. Used to calculate signal-to-noise ratio.
        
    rt_interval : tuple
        A tuple where the first index contains the beggining time of the retention time
        interval you wish to analyze and the second contains the end time.
        
    min_ppp : tuple
        A tuple where the first index contains a boolean indicating whether or not to
        consider this parameter and the second one containing the minimum amount of
        data points per peak. If min_ppp[0] is set to False, calculates this automatically.
        
    close_peaks : tuple
        A tuple where the first index contains a boolean indicating whether or not to
        consider this parameter and the second one contains the amount of peaks to save.
        If close_peaks[0] is set to True, selects only the most intense peak and its 
        close_peaks[1] surrounding peaks.
        
    glycan : string
        Glycan name to identify when dealing with the Internal Standard.
        
    Returns
    -------
    peaks : list
        A list of dictionaries with each one containing numerous peaks information.
    '''
    peaks = []
    peaks_ranges = []
    datapoints_per_time = int((0.2/(rt_int[0][rt_int[1].index(max(rt_int[1]))]-rt_int[0][rt_int[1].index(max(rt_int[1]))-1]))*(rt_int[0][-1]/60))
    maximums_index = []
    min_relative_int_peak = 0.001
    threshold = min_relative_int_peak*max(rt_int_smoothed[1])
    for i_i, i in enumerate(rt_int_smoothed[1]):
        if (rt_int[0][i_i] > rt_interval[1] or rt_int[0][i_i] == rt_int[0][-2]):
            break    
        if rt_int[0][i_i] >= rt_interval[0]:
            if i > threshold and rt_int_smoothed[1][i_i-1] <= i and rt_int_smoothed[1][i_i+1] <= i:
                if len(maximums_index) == 0:
                    maximums_index.append(i_i)
                elif len(maximums_index) > 0 and i_i-maximums_index[-1] >= datapoints_per_time:
                    maximums_index.append(i_i)
    
    former_peak_limit = 0
    for i_i, i in enumerate(maximums_index):
        #print('starting id: ', i, 'former peak limit: ', former_peak_limit)
        peak_limits = []
        temp_min = inf
        temp_min_rt_index = 0
        for k_k in range(i, 0, -1):
            #print('before peak: ', k_k, rt_int[0][k_k], rt_int_smoothed[1][k_k], temp_min)
            if rt_int_smoothed[1][k_k] < temp_min:
                temp_min = rt_int_smoothed[1][k_k]
                temp_min_rt_index = k_k
            if rt_int_smoothed[1][k_k] < rt_int_smoothed[1][i]*min_relative_int_peak:
                #print('breaking due to low intensity')
                break
            if rt_int[0][k_k] < rt_interval[0] or k_k == former_peak_limit:
                #print('breaking due to stumbling upon last peak')
                break
        if temp_min != inf:
            peak_limits.append(temp_min_rt_index)
                    
        temp_min = inf
        temp_min_rt_index = 0
        next_peak_limit = sorted(maximums_index)[i_i+1] if i_i != len(maximums_index)-1 else len(rt_int[0])-1
        
        for k_k in range(i, len(rt_int[0])):
            #print('after peak: ', k_k, rt_int[0][k_k], rt_int_smoothed[1][k_k], temp_min)
            if rt_int_smoothed[1][k_k] < temp_min:
                temp_min = rt_int_smoothed[1][k_k]
                temp_min_rt_index = k_k
            if rt_int_smoothed[1][k_k] < rt_int_smoothed[1][i]*min_relative_int_peak:
                #print('breaking due to low intensity')
                break
            if rt_int[0][k_k] > rt_interval[1] or k_k == next_peak_limit:
                #print('breaking due to stumbling upon next peak')
                break
        if temp_min != inf:
            peak_limits.append(temp_min_rt_index)
            former_peak_limit = temp_min_rt_index
        #print(peak_limits)
        if len(peak_limits) == 2:
            peaks_ranges = peaks_ranges + peak_limits
            
    #print(peaks_ranges)
    
    removal = []        
    for i_i, i in enumerate(peaks_ranges):
        if i_i > 1 and i != peaks_ranges[-1]:
            #print(peaks_ranges[i_i-1], peaks_ranges[i_i])
            if peaks_ranges[i_i] == peaks_ranges[i_i-1]:
                #print("True!")
                if (rt_int_smoothed[1][i] > max(rt_int_smoothed[1][peaks_ranges[i_i-2]:peaks_ranges[i_i+1]])*0.8 
                    or abs(max(rt_int_smoothed[1][peaks_ranges[i_i-2]:peaks_ranges[i_i-1]])-rt_int_smoothed[1][peaks_ranges[i_i]]) < rt_int_smoothed[1][peaks_ranges[i_i]]*0.1 
                    or abs(max(rt_int_smoothed[1][peaks_ranges[i_i]:peaks_ranges[i_i+1]])-rt_int_smoothed[1][peaks_ranges[i_i]]) < rt_int_smoothed[1][peaks_ranges[i_i]]*0.1):
                    
                    removal.append(i_i-1)
                    removal.append(i_i)
    for i in sorted(removal, reverse=True):
        del peaks_ranges[i]
        
    #print(peaks_ranges)
    
    for i_i, i in enumerate(peaks_ranges):
        if i_i % 2 == 0:
            peak_limits = [i, peaks_ranges[i_i+1]]
            if peak_limits[0] == peak_limits[1] or (min_ppp[0] and peak_limits[1] - peak_limits[0] < min_ppp[1]):
                continue
            temp_peak_width = (rt_int[0][peak_limits[1]]-rt_int[0][peak_limits[0]])
            peaks.append({'id': i, 'rt': rt_int[0][rt_int_smoothed[1].index(max(rt_int_smoothed[1][peak_limits[0]:peak_limits[1]+1]))], 'int': max(raw_rt_int[1][peak_limits[0]:peak_limits[1]+1]), 'peak_width': temp_peak_width, 'peak_interval': (rt_int[0][peak_limits[0]], rt_int[0][peak_limits[1]]), 'peak_interval_id': (peak_limits[0], peak_limits[1])})
    
    #print(peaks)
        
    if close_peaks[0] or glycan == "Internal Standard":
        peaks = sorted(peaks, key=lambda x: x['int'], reverse = True)
        for i_i, i in enumerate(peaks):
            peaks[i_i]['proximity'] = abs(i['rt']-peaks[0]['rt'])
        peaks = sorted(peaks, key=lambda x: x['proximity'])
        if glycan == "Internal Standard":
            close_peaks = (True, 1)
        return sorted(peaks[:close_peaks[1]], key=lambda x: x['rt'])
    return sorted(peaks, key=lambda x: x['rt'])

def peaks_auc_from_eic(rt_int,
                       ms1_index,
                       peaks):
    '''Calculates the area under curve of the picked peaks based on the EIC given.
    Overlapped portions of peaks are calculated as fractional proportion, based on peaks
    max intensity.

    Parameters
    ----------
    rt_int : list
        A list containing two synchronized lists, the first one containing retention
        times and the second one containing intensity.
        
    ms1_indexes : dict
        A dictionary where each key is the ID of a file and each value is a list
        containing the indexes of all MS1 scans in the file.

    peaks : list
        A list of dictionaries with each one containing numerous peaks information.

    Returns
    -------
    auc : list
        A list of AUCs, with each index containing a float of the AUC of a peak.
    '''
    auc = []
    for i in peaks:
        temp_auc = 0
        for j in range(i['peak_interval_id'][0], i['peak_interval_id'][1]):
            temp_auc+=rt_int[1][j]
        auc.append(temp_auc)
    return auc
