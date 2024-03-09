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
from math import inf, atan, pi, exp
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
    object : string
        A string containing the path to the mzML file.
        
    Uses
    ----
    pyteomics.mzml.MzML : generator
        Indexes the mzML file into a generator, allowing you to parse the file for
        analysis.
        
    Returns
    -------
    standard mzML output
        If class is iterated, returns regular mzml.MzML output.
        UNFORTUNATELY COULDN'T CONVERT THIS FUNCTIONALITY
        CIRCUNVENTED USING THE INDEX FUNCTION TO ACQUIRE SPECIFIC SPECTRA INFORMATION
        AND DIRECTLY QUERYING EACH SPECTRUM
    
    dict
        If class is queried for a single index, returns the dictionary of the
        converted output (mzML -> mzXML)
        
    data : list
        If class is queried for slices, returns a list of dictionaries containing
        the converted output (mzML -> mzXML)
    '''
    def __init__(self,it):
        self.it = mzml.MzML(it)
    def __iter__(self):
        return self.it
    def __getitem__(self,index):
        if type(index) == int:
            pre_data = self.it[index]
            if pre_data['ms level'] == 2:
                if float(self.it[-1]['scanList']['scan'][0]['scan start time']) > 210: #210 scan time should allow for the correct evaluation of scan time being in seconds or minutes for every run that lasts between 3.5 minutes and 3.5 hours
                    return {'num': pre_data['id'].split('=')[-1], 'retentionTime': float(pre_data['scanList']['scan'][0]['scan start time'])/60, 'msLevel': pre_data['ms level'], 'm/z array': pre_data['m/z array'], 'intensity array': pre_data['intensity array'], 'precursorMz': [{'precursorMz': pre_data['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']}]}
                else:
                    return {'num': pre_data['id'].split('=')[-1], 'retentionTime': float(pre_data['scanList']['scan'][0]['scan start time']), 'msLevel': pre_data['ms level'], 'm/z array': pre_data['m/z array'], 'intensity array': pre_data['intensity array'], 'precursorMz': [{'precursorMz': pre_data['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']}]}
            else:
                if float(self.it[-1]['scanList']['scan'][0]['scan start time']) > 210:
                    return {'num': pre_data['id'].split('=')[-1], 'retentionTime': float(pre_data['scanList']['scan'][0]['scan start time'])/60, 'msLevel': pre_data['ms level'], 'm/z array': pre_data['m/z array'], 'intensity array': pre_data['intensity array']}
                else:
                    return {'num': pre_data['id'].split('=')[-1], 'retentionTime': float(pre_data['scanList']['scan'][0]['scan start time']), 'msLevel': pre_data['ms level'], 'm/z array': pre_data['m/z array'], 'intensity array': pre_data['intensity array']}
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
                last_index = len(self.it)
            data = []
            for index in range(first_index, last_index):
                if self.it[index]['ms level'] == 2:
                    if float(self.it[-1]['scanList']['scan'][0]['scan start time']) > 210:
                        data.append({'num': self.it[index]['id'].split('=')[-1], 'retentionTime': float(self.it[index]['scanList']['scan'][0]['scan start time'])/60, 'msLevel': self.it[index]['ms level'], 'm/z array': self.it[index]['m/z array'], 'intensity array': self.it[index]['intensity array'], 'precursorMz': [{'precursorMz': self.it[index]['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']}]})
                    else:
                        data.append({'num': self.it[index]['id'].split('=')[-1], 'retentionTime': float(self.it[index]['scanList']['scan'][0]['scan start time']), 'msLevel': self.it[index]['ms level'], 'm/z array': self.it[index]['m/z array'], 'intensity array': self.it[index]['intensity array'], 'precursorMz': [{'precursorMz': self.it[index]['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']}]})
                else:
                    if float(self.it[-1]['scanList']['scan'][0]['scan start time']) > 210:
                        data.append({'num': self.it[index]['id'].split('=')[-1], 'retentionTime': float(self.it[index]['scanList']['scan'][0]['scan start time'])/60, 'msLevel': self.it[index]['ms level'], 'm/z array': self.it[index]['m/z array'], 'intensity array': self.it[index]['intensity array']})
                    else:
                        data.append({'num': self.it[index]['id'].split('=')[-1], 'retentionTime': float(self.it[index]['scanList']['scan'][0]['scan start time']), 'msLevel': self.it[index]['ms level'], 'm/z array': self.it[index]['m/z array'], 'intensity array': self.it[index]['intensity array']})
            return data
       
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
                    ms1_id):
    '''Generates a very processed EIC for each adduct of each glycan of each sample.
    Removes non-monoisotopic peaks, check charges and calculates multiple quality scoring
    data.

    Parameters
    ----------
    files : list
        A list of generators obtained from the function pyteomics.mzxml.MzXML().

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
        
    raw_data : dict
        A dictionary with keys for each adduct combo of a given glycan, which value is a
        dict with keys for each file id, at which value is a list of lists, with the
        first one being the rt array and the second one the raw int array.
        
    isotopic_fits : dict
        A dictionary containing all the isotopic fitting information for reporting and 
        checking data.
    '''
    data = {}
    ppm_info = {}
    iso_fitting_quality = {}
    isotopic_fits = {}
    verbose_info = []
    raw_data = {}
    for i in glycan_info['Adducts_mz']:
        adduct_mass = mass.calculate_mass(composition=General_Functions.form_to_comp(i))
        adduct_charge = General_Functions.form_to_charge(i)
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
            
            #this is a possible point of parallelization
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
                                 adduct_mass,
                                 adduct_charge)
                
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
                     adduct_mass,
                     adduct_charge):
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
    intensity = 0.0 #variable to sum the total deisotopotized intensity
    mono_int = 0.0 #variable to sum the total intensity of the monoisotopic peak
    iso_distro = 1 #starting isotopic distribution peak, always 1
    charge_range = range(1, abs(max_charges)*2)
    target_mz = glycan_info['Adducts_mz'][glycan_id]
    second_peak = (glycan_info['Isotopic_Distribution_Masses'][1]+adduct_mass)/adduct_charge
    sec_peak_rel_int = glycan_info['Isotopic_Distribution'][1]
    last_peak = (glycan_info['Isotopic_Distribution_Masses'][-1]+adduct_mass)/adduct_charge
    no_iso_peaks = len(glycan_info['Isotopic_Distribution_Masses'])
    before_target = glycan_info['Adducts_mz'][glycan_id]-General_Functions.h_mass-General_Functions.tolerance_calc(tolerance[0], tolerance[1], glycan_info['Adducts_mz'][glycan_id])
    mono_ppm = []
    mz_isos = []
    iso_actual = []
    iso_target = []
    bad_peaks_before_target = []
    max_int = max(sliced_int)
    nearby_id = 0
    checked_bad_peaks_before_target = False
    iso_found = False
    current_iso_peak1 = []
    current_iso_peak2 = []
    current_mz = []
    not_good = False
    found = False
    for l_l, l in enumerate(sliced_mz):
        local_noise = General_Functions.local_noise_calc(noise[file_id][ms1_id], l, avg_noise[file_id])
        current_tolerance = General_Functions.tolerance_calc(tolerance[0], tolerance[1], l)
        if iso_found and l > ((glycan_info['Isotopic_Distribution_Masses'][iso_distro]+adduct_mass)/adduct_charge) + current_tolerance:
            if len(current_iso_peak1) > 5: #here data is profile (NOT RECOMMENDED, but it will work.... very slowly....)
                iso_actual.append(max(current_iso_peak1))
                iso_target.append(max(current_iso_peak2))
            else:
                iso_actual.append(sum(current_iso_peak1))
                iso_target.append(sum(current_iso_peak2))
            mz_isos.append(sum(current_mz)/len(current_mz))
            current_iso_peak1 = []
            current_iso_peak2 = []
            current_mz = []
            iso_distro += 1
            iso_found = False
        if not_good: #Here are checks for quick skips
            break
        if max_int < avg_noise[file_id]*0.5: #if there's no peak with intensity higher than 50% of the average noise of sample, don't even check it
            break
        if l_l == len(sliced_mz)-1 and not found:
            not_good = True
            break
        if target_mz > sliced_mz[-1]:
            not_good = True
            break
        if l < before_target: #This is the last check for quick skips
            continue
        if sliced_int[l_l] < local_noise*0.5: #This ignores peaks that are far below the noise level (less than 50% of the calculated noise level)
            continue
        if l > second_peak + current_tolerance and iso_distro == 1: #This is the first check for quick skips that's dependent on mz array acquired data but independent of noise
            not_good = True
            break
        if l > (target_mz + current_tolerance) and not found:
            not_good = True
            break
        if l > last_peak + current_tolerance and found:
            if iso_distro < min_isotops:
                not_good = True
            break
        if l > target_mz + current_tolerance and l < target_mz + General_Functions.h_mass and found:
            for m in charge_range:
                if m == 1:
                    continue
                expected_value = (l*m*0.0006)+0.1401 #based on linear regression of the relationship between masses and the second isotopic peak relative intensity of the average of different organic macromolecules
                if m != adduct_charge and (sliced_int[l_l] > mono_int*expected_value*0.6)  and abs(l-(target_mz+(General_Functions.h_mass/m))) <= current_tolerance:
                    not_good = True
                    break
            if not_good:
                break
        if l > target_mz - General_Functions.h_mass - current_tolerance and l < target_mz - current_tolerance:
            nearby_id = l_l
            for m in charge_range:
                if abs(l-(target_mz-(General_Functions.h_mass/m))) <= current_tolerance:
                    bad_peaks_before_target.append((sliced_int[l_l], m))
                    break
        if l > target_mz + current_tolerance and abs(l-((glycan_info['Isotopic_Distribution_Masses'][iso_distro]+adduct_mass)/adduct_charge)) <= current_tolerance:
            if sliced_int[l_l] > mono_int*glycan_info['Isotopic_Distribution'][iso_distro]:
                intensity += mono_int*glycan_info['Isotopic_Distribution'][iso_distro]
            if sliced_int[l_l] <= mono_int*glycan_info['Isotopic_Distribution'][iso_distro]:
                intensity += sliced_int[l_l]
            current_mz.append(l)
            current_iso_peak1.append(sliced_int[l_l]/mono_int)
            current_iso_peak2.append(glycan_info['Isotopic_Distribution'][iso_distro])
            iso_found = True
            continue
        if not checked_bad_peaks_before_target and l > target_mz + current_tolerance:
            checked_bad_peaks_before_target = True
            for m in bad_peaks_before_target:
                expected_value = 1/((l*m[1]*0.0006)+0.1401) #based on linear regression of the relationship between masses and the second isotopic peak relative intensity of the average of different organic macromolecules
                if (m[0] > mono_int*expected_value*0.6):
                    not_good = True
                    break
            if not_good:
                break
        if sliced_int[l_l] < local_noise: #Everything from here is dependent on noise level (currently: monoisotopic peak detection only)
            continue
        if l >= target_mz - current_tolerance and abs(l-target_mz) <= current_tolerance:
            mono_ppm.append(General_Functions.calculate_ppm_diff(l, target_mz))
            intensity += sliced_int[l_l]
            mono_int += sliced_int[l_l]
            found = True
            continue
            
    for l_l in range(nearby_id, len(sliced_mz)):
        if sliced_mz[l_l] > target_mz+General_Functions.tolerance_calc(tolerance[0], tolerance[1], target_mz) or l_l == len(sliced_mz)-1:
            raw_data[glycan_id][file_id][1][ms1_id] = 0.0
            break
        if abs(sliced_mz[l_l] - target_mz) <= General_Functions.tolerance_calc(tolerance[0], tolerance[1], target_mz):
            raw_data[glycan_id][file_id][1][ms1_id] = sliced_int[l_l]
            break
            
    if len(iso_actual) > 0:
        dotp = []
        number = range(1, len(mz_isos)+1)
        starting_points_actual = [0, 1]
        starting_points_theoretical = [0, 1]
        for m_m, m in enumerate(number):
            vector_actual = [m-starting_points_actual[0], iso_actual[m_m]-starting_points_actual[1]]
            vector_target = [m-starting_points_theoretical[0], iso_target[m_m]-starting_points_theoretical[1]]
            normalized_actual = vector_actual/numpy.linalg.norm(vector_actual)
            normalized_target = vector_target/numpy.linalg.norm(vector_target)
            starting_points_actual = [m, iso_actual[m_m]]
            starting_points_theoretical = [m, iso_target[m_m]]
            dotproduct = numpy.dot(normalized_actual, normalized_target)
            dotp.append(dotproduct)
        iso_quali = mean(dotp)
        
        #reduces score if fewer isotopic peaks are found: very punishing for only 2 peaks, much less punishing for three, normal score from 4 and over
        if len(iso_actual) == 1:
            iso_quali = (iso_quali+(iso_quali*0.70))/2
        if len(iso_actual) == 2:
            iso_quali = (iso_quali+(iso_quali*0.95))/2
            
    if len(iso_actual) == 0:
        iso_quali = 0.0
        
    if len(mono_ppm) > 0:
        ppm_error = mean(mono_ppm)
    else:
        ppm_error = inf
        
    ppm_info[glycan_id][file_id][ms1_id] = ppm_error
    iso_fitting_quality[glycan_id][file_id][ms1_id] = iso_quali
    data[glycan_id][file_id][1][ms1_id] = intensity
    isotopic_fits[glycan_id][file_id][ret_time] = [[glycan_info['Adducts_mz'][glycan_id]]+mz_isos, [1]+iso_target, [1]+iso_actual, iso_quali]
    
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
    interval = x[len(x)//2]-x[(len(x)//2)-1]
    fits_list = []
    for j in range(int(0.2/interval)):
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
        max_amp_id = temp_y.index(max(temp_y))
        counts_max = 0
        for i in temp_y:
            if i == temp_y[max_amp_id]:
                counts_max += 1
        max_amp_id += int(counts_max/2)
        y_gaussian = []
        for i in temp_x:
            y_gaussian.append(General_Functions.normpdf(i, temp_x[max_amp_id], (temp_x[-1]-temp_x[0])/6))
        scaler = (temp_y[max_amp_id]/y_gaussian[max_amp_id])
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
                       peak):
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
        
    Returns
    -------
    iso_fit_score : float
        Average of the isotopic fits score calculated for the interval of
        the peak of the glycan, weighted by the gaussian fit of it.
    '''
    weights_list = []
    for i in range(peak['peak_interval_id'][1]-peak['peak_interval_id'][0]+1):
        weights_list.append(General_Functions.normpdf(i, (peak['peak_interval_id'][1]-peak['peak_interval_id'][0]+1)/2, (peak['peak_interval_id'][1]-peak['peak_interval_id'][0])/6))
    scaler = (1/max(weights_list))
    weights_list_scaled = []
    for i in weights_list:
        weights_list_scaled.append(i*scaler)
    iso_fit_score = numpy.average(iso_fits[peak['peak_interval_id'][0]:peak['peak_interval_id'][1]+1], weights = weights_list_scaled)
    return iso_fit_score
    
def average_ppm_calc(ppm_array,
                     tolerance,
                     peak):
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
    ppms = []
    ppm_default = General_Functions.calculate_ppm_diff(tolerance[2]-General_Functions.tolerance_calc(tolerance[0], tolerance[1], tolerance[2]), tolerance[2])
    missing_points = 0
    weights_list = []
    for i in range(peak['peak_interval_id'][1]-peak['peak_interval_id'][0]+1):
        weights_list.append(General_Functions.normpdf(i, (peak['peak_interval_id'][1]-peak['peak_interval_id'][0]+1)/2, (peak['peak_interval_id'][1]-peak['peak_interval_id'][0])/6))
    scaler = (1/max(weights_list))
    weights_list_scaled = []
    for i in weights_list:
        weights_list_scaled.append(i*scaler)
    for i in ppm_array[peak['peak_interval_id'][0]:peak['peak_interval_id'][1]+1]:
        if i != inf:
            ppms.append(i)
        else:
            ppms.append(ppm_default)
            missing_points+= 1
    regular_mean = numpy.average(ppms, weights = weights_list_scaled)
    return regular_mean, missing_points

def peaks_from_eic(rt_int, 
                   rt_int_smoothed,
                   rt_interval,
                   min_ppp,
                   close_peaks,
                   glycan):
    '''Does multi peak-picking in a given smoothed EIC.
    
    Parameters
    ----------
    rt_int : list
        A list containing two lists: the first one has all the retention times, the
        second one has all the intensities.
        
    rt_int_smoothed : list
        A list containing two lists: the first one has all the retention times, the
        second one has all the SMOOTHED intensities.
        
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
    temp_start = 0
    temp_max = 0
    temp_max_id_iu = 0
    going_up = False
    going_down = False
    datapoints_per_time = int(0.2/(rt_int[0][rt_int[1].index(max(rt_int[1]))]-rt_int[0][rt_int[1].index(max(rt_int[1]))-1]))
    slope_threshold = 5*datapoints_per_time #this number indicates the minimum angle to consider a peak
    for i_i, i in enumerate(rt_int_smoothed[1]):
        if (rt_int[0][i_i] >= rt_interval[1] or rt_int[0][i_i] == rt_int[0][-2]):
            break    
        if rt_int[0][i_i] >= rt_interval[0]:
            if i > 0 and rt_int_smoothed[1][i_i+1] > 0:
                arc_tan_slope = atan(((rt_int_smoothed[1][i_i+1]/i)-1)/(1/datapoints_per_time))*(180/pi)
            elif i > 0 and rt_int_smoothed[1][i_i+1] == 0:
                arc_tan_slope = atan(((1/(i+1))-1)/(1/datapoints_per_time))*(180/pi)
            elif i == 0 and rt_int_smoothed[1][i_i+1] > 0:
                arc_tan_slope = atan((rt_int_smoothed[1][i_i+1])/(1/datapoints_per_time))*(180/pi)
            elif i == 0 and rt_int_smoothed[1][i_i+1] == 0:
                arc_tan_slope = 0.0
            # print(rt_int[0][i_i], i, rt_int_smoothed[1][i_i+1], arc_tan_slope, going_up, going_down, datapoints_per_time, slope_threshold)
            if (going_up or going_down) and rt_int[1][i_i] > temp_max:
                temp_max = rt_int[1][i_i]
                temp_max_id_iu = i_i
            if (going_up and (arc_tan_slope < 0 or rt_int_smoothed[1][i_i] == 0)):
                going_up = False
                going_down = True
            if (going_down and (arc_tan_slope >= 0 or rt_int_smoothed[1][i_i] == 0.0)):
                min_id = i_i
                if rt_int_smoothed[1][i_i+1] < i:
                    min_id = i_i+1
                temp_peak_width = (rt_int[0][min_id]-rt_int[0][temp_start])
                going_down = False
                good = False
                if min_ppp[0]:
                    if i_i-temp_start >= min_ppp[1] or glycan == "Internal Standard":
                        good = True
                else:
                    if min_id-temp_start >= datapoints_per_time or glycan == "Internal Standard":
                        good = True
                if good:
                    peaks.append({'id': temp_max_id_iu, 'rt': rt_int[0][temp_max_id_iu], 'int': temp_max, 'peak_width': temp_peak_width, 'peak_interval': (rt_int[0][temp_start], rt_int[0][min_id]), 'peak_interval_id': (temp_start, min_id)})
                    temp_max = 0
                    temp_max_id_iu = 0
                if not good:
                    temp_max = 0
                    temp_max_id_iu = 0
            if (i_i != 0 and arc_tan_slope >= slope_threshold and not going_up and not going_down):
                temp_start = i_i-1
                going_up = True
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
