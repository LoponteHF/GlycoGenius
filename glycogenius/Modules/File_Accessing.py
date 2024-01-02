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

import pathlib
import importlib
import_path = str(pathlib.Path(__file__).parent.resolve())
for i_i, i in enumerate(import_path):
    if i == "\\":
        import_path = import_path[:i_i]+"/"+import_path[i_i+1:]
        
#Absolute import of General_Functions
spec1 = importlib.util.spec_from_file_location("General_Functions", import_path+'/General_Functions.py')
General_Functions = importlib.util.module_from_spec(spec1)
spec1.loader.exec_module(General_Functions)

from pyteomics import mzxml, mzml, mass, auxiliary
from itertools import combinations_with_replacement
from scipy.signal import savgol_filter
from statistics import mean
from re import split
from math import inf
import numpy
import sys
import datetime

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
                    fast_iso,
                    verbose = False):
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
        
    verbose : boolean
        Allows to print to file debugging information of each EIC traced so that you
        may find out why some specific pattern or behavior is ocurring.

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
        A list containing the conclusion of the processing for each retention time of each
        chromatogram processed, for each file.
        
    raw_data : dict
        A dictionary with keys for each adduct combo of a given glycan, which value is a
        dict with keys for each file id, at which value is a list of lists, with the
        first one being the rt array and the second one the raw int array.
    '''
    data = {}
    ppm_info = {}
    iso_fitting_quality = {}
    verbose_info = []
    raw_data = {}
    for i in glycan_info['Adducts_mz']:
        if verbose:
            print('Adduct: '+str(i)+" mz: "+str(glycan_info['Adducts_mz'][i]))
            verbose_info.append('Adduct: '+str(i)+" mz: "+str(glycan_info['Adducts_mz'][i]))
        adduct_mass = mass.calculate_mass(composition=General_Functions.form_to_comp(i))
        adduct_charge = General_Functions.form_to_charge(i)
        ppm_info[i] = {}
        iso_fitting_quality[i] = {}
        data[i] = {}
        raw_data[i] = {}
        for j_j, j in enumerate(files):
            if verbose:
                print("--Drawing EIC for Sample: "+str(j_j))
                verbose_info.append("--Sample: "+str(j_j))
            ppm_info[i][j_j] = []
            iso_fitting_quality[i][j_j] = []
            data[i][j_j] = [[], []]
            raw_data[i][j_j] = [[], []]
            for k_k, k in enumerate(ms1_indexes[j_j]):
                iso_fitting_quality[i][j_j].append(0.0)
                ppm_info[i][j_j].append(inf)
                not_good = False
                found = False
                rt = j[k]['retentionTime']
                data[i][j_j][0].append(rt)
                raw_data[i][j_j][0].append(rt)
                if (j[k]['retentionTime'] < rt_interval[0] or j[k]['retentionTime'] > rt_interval[1]):
                    data[i][j_j][1].append(0.0)
                    raw_data[i][j_j][1].append(0.0)
                    continue
                else:   
                    if verbose:
                        verbose_info.append("----RT: "+str(rt))
                    l_l = 0
                    last_skip = 2
                    intensity = 0.0
                    mono_int = 0.0
                    iso_distro = 1
                    sliced_mz = j[k]['m/z array']
                    sliced_int = j[k]['intensity array']
                    m_range = range(1, abs(max_charges)+1)
                    target_mz = glycan_info['Adducts_mz'][i]
                    second_peak = (glycan_info['Isotopic_Distribution_Masses'][1]+adduct_mass)/adduct_charge
                    sec_peak_rel_int = glycan_info['Isotopic_Distribution'][1]
                    last_peak = (glycan_info['Isotopic_Distribution_Masses'][-1]+adduct_mass)/adduct_charge
                    no_iso_peaks = len(glycan_info['Isotopic_Distribution_Masses'])
                    before_target = glycan_info['Adducts_mz'][i]-General_Functions.h_mass-General_Functions.tolerance_calc(tolerance[0], tolerance[1], glycan_info['Adducts_mz'][i])
                    mono_ppm = []
                    iso_actual = []
                    iso_target = []
                    bad_peaks_before_target = []
                    nearby_id = 0
                    for l_l, l in enumerate(sliced_mz):
                        if not_good: #Here are checks for quick skips
                            break
                        if l_l == len(sliced_mz)-1 and not found:
                            if verbose:
                                verbose_info.append("------m/z "+str(l)+", int "+str(sliced_int[l_l]))
                                verbose_info.append("--------Couldn't find peak and array ends shortly afterwards.")
                            not_good = True
                            break
                        if target_mz > sliced_mz[-1]:
                            if verbose:
                                verbose_info.append("------m/z "+str(l)+", int "+str(sliced_int[l_l]))
                                verbose_info.append("--------Target mz outside acquired mz interval.")
                            not_good = True
                            break
                        if l < before_target: #This is the last check for quick skips
                            continue
                        if l > second_peak + General_Functions.tolerance_calc(tolerance[0], tolerance[1], l) and iso_distro == 1: #This is the first check for quick skips that's dependent on mz array acquired data but independent of noise
                            if verbose:
                                verbose_info.append("------m/z "+str(l)+", int "+str(sliced_int[l_l]))
                                verbose_info.append("--------Couldn't check correct charge.")
                            not_good = True
                            break
                        if iso_distro == no_iso_peaks:
                            if verbose:
                                verbose_info.append("------m/z "+str(l)+", int "+str(sliced_int[l_l]))
                                verbose_info.append("--------Reached the end of isotopic distribution without errors.")
                            break
                        if l > (target_mz + General_Functions.tolerance_calc(tolerance[0], tolerance[1], l)) and not found:
                            if verbose:
                                verbose_info.append("------m/z "+str(l)+", int "+str(sliced_int[l_l]))
                                verbose_info.append("--------Couldn't find monoisotopic peak.")
                            not_good = True
                            break
                        if l > last_peak + General_Functions.tolerance_calc(tolerance[0], tolerance[1], l) and found:
                            if verbose:
                                verbose_info.append("------m/z "+str(l)+", int "+str(sliced_int[l_l]))
                                verbose_info.append("--------Checked all available isotopologue peaks with no error.")
                                verbose_info.append("--------Isotopologue peaks found: "+str(iso_distro))
                            if iso_distro < min_isotops:
                                not_good = True
                            break
                        if l > target_mz + General_Functions.tolerance_calc(tolerance[0], tolerance[1], l) and l < target_mz + General_Functions.h_mass and found:
                            for m in m_range:
                                if abs(l-(target_mz+(General_Functions.h_mass/m))) <= General_Functions.tolerance_calc(tolerance[0], tolerance[1], l) and m != adduct_charge and sliced_int[l_l] < mono_int*(sec_peak_rel_int*1.2) and sliced_int[l_l] > mono_int*(sec_peak_rel_int*0.8):
                                    if verbose:
                                        verbose_info.append("------m/z "+str(l)+", int "+str(sliced_int[l_l]))
                                        verbose_info.append("--------Incorrect charge assigned.")
                                    not_good = True
                                    break
                        if l > target_mz - General_Functions.h_mass - General_Functions.tolerance_calc(tolerance[0], tolerance[1], l) and l < target_mz - General_Functions.tolerance_calc(tolerance[0], tolerance[1], l):
                            nearby_id = l_l
                            for m in m_range:
                                if abs(l-(target_mz-(General_Functions.h_mass/m))) <= General_Functions.tolerance_calc(tolerance[0], tolerance[1], l):
                                    if verbose:
                                        verbose_info.append("------m/z "+str(l)+", int "+str(sliced_int[l_l]))
                                        verbose_info.append("--------Not monoisotopic.")
                                    bad_peaks_before_target.append(sliced_int[l_l])
                                    break
                        if l > target_mz + General_Functions.tolerance_calc(tolerance[0], tolerance[1], l) and abs(l-((glycan_info['Isotopic_Distribution_Masses'][iso_distro]+adduct_mass)/adduct_charge)) <= General_Functions.tolerance_calc(tolerance[0], tolerance[1], l):
                            if sliced_int[l_l] > mono_int*glycan_info['Isotopic_Distribution'][iso_distro]:
                                intensity += mono_int*glycan_info['Isotopic_Distribution'][iso_distro]
                            if sliced_int[l_l] <= mono_int*glycan_info['Isotopic_Distribution'][iso_distro]:
                                intensity += sliced_int[l_l]
                            if iso_distro <= min_isotops - 1:
                                iso_actual.append(sliced_int[l_l])
                                iso_target.append(mono_int*glycan_info['Isotopic_Distribution'][iso_distro])
                            iso_distro += 1
                            continue
                        if sliced_int[l_l] < General_Functions.local_noise_calc(noise[j_j][k_k], l, avg_noise[j_j]): #Everything from here is dependent on noise level
                            continue
                        if l >= target_mz - General_Functions.tolerance_calc(tolerance[0], tolerance[1], l) and abs(l-target_mz) <= General_Functions.tolerance_calc(tolerance[0], tolerance[1], l):
                            mono_ppm.append(General_Functions.calculate_ppm_diff(l, target_mz))
                            intensity += sliced_int[l_l]
                            mono_int += sliced_int[l_l]
                            for m in bad_peaks_before_target:
                                if m > mono_int*(1/sec_peak_rel_int*0.8):
                                    not_good = True
                                    break
                            found = True
                            continue
                    for l_l in range(nearby_id, len(sliced_mz)):
                        if sliced_mz[l_l] > target_mz+General_Functions.tolerance_calc(tolerance[0], tolerance[1], sliced_mz[l_l]) or l_l == len(sliced_mz)-1:
                            raw_data[i][j_j][1].append(0.0)
                            break
                        if abs(sliced_mz[l_l] - target_mz) <= General_Functions.tolerance_calc(tolerance[0], tolerance[1], sliced_mz[l_l]):
                            raw_data[i][j_j][1].append(sliced_int[l_l])
                            break
                if not_good:
                    ppm_info[i][j_j][-1] = inf
                    iso_fitting_quality[i][j_j][-1] = 0.0
                    data[i][j_j][1].append(0.0)
                    continue
                else:
                    ppm_info[i][j_j][-1] = mean(mono_ppm)
                    if len(iso_actual) > 0:
                        weights = [0.5, 0.35, 0.15, 0.1, 0.01, 0.001, 0.0001, 0.001, 0.0001, 0.00001]
                        temp_relation = []
                        for l in range(len(iso_actual)):
                            if iso_actual[l] >= iso_target[l]:
                                temp_relation.append(((iso_target[l]/iso_actual[l])+1)/2)
                            else:
                                temp_relation.append(iso_actual[l]/iso_target[l])
                        R_sq = numpy.average(temp_relation, weights = weights[:len(temp_relation)])
                    if len(iso_actual) == 0:
                        R_sq = 0.0
                    iso_fitting_quality[i][j_j][-1] = R_sq
                    data[i][j_j][1].append(intensity)
    return data, ppm_info, iso_fitting_quality, verbose_info, raw_data
    
def eic_smoothing(rt_int):
    '''Smoothes the EIC using the Savitsky Golay algorithm. Smoothed EIC may be
    used for peak-picking and curve-fitting scoring afterwards.
    
    Parameters
    ----------
    rt_int : list
        A list containing two lists: the first one has all the retention times, the
        second one has all the intensities.
        
    Uses
    ----
    scipy.savgol_filter : array
        Apply a Savitzky-Golay filter to an array.
        
    Returns
    -------
    rt_int[0] : list
        A list of each retention time.
        
    filtered_ints : list
        A list containing the processed intensities.
    '''
    points = 21
    polynomial_degree = 3
    while polynomial_degree >= points:
        points+=1
    filtered_ints = list(savgol_filter(rt_int[1], points, polynomial_degree))
    for i_i, i in enumerate(filtered_ints):
        if i < 0:
            filtered_ints[i_i] = 0.0
    return rt_int[0], filtered_ints
    
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
    interval = x[-1]-x[-2]
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
    statistics.mean : float
        Return the sample arithmetic mean of data which can be a sequence or iterable.
        
    Returns
    -------
    float
       The arithmetic mean of the isotopic fitting scores for the given peak.
    '''
    missing_points = 0
    calculated_iso_points = []
    max_missing_points = (peak['peak_interval_id'][1]-peak['peak_interval_id'][0])/2
    regular_mean = mean(iso_fits[peak['peak_interval_id'][0]:peak['peak_interval_id'][1]+1])
    iter_range = range(peak['peak_interval_id'][0], peak['peak_interval_id'][1]+1)
    for i in iter_range:
        if iso_fits[i] == 0.0:
            if missing_points == 0:
                missing_points+= 1
                continue
            if missing_points > 0 and missing_points <= max_missing_points:
                calculated_iso_points.append(regular_mean-(missing_points*(regular_mean/max_missing_points)))
                missing_points+= 1
                continue
            if missing_points > max_missing_points:
                return 0.0
        calculated_iso_points.append(iso_fits[i])
    if len(calculated_iso_points) > 0:
        return mean(calculated_iso_points)
    else:
        return 0.0
    
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
        
    Returns
    -------
    tuple
        A tuple containing the mean ppm differences and the number of missing points.
        The number of missing points in the peak has their ppm difference set to the
        tolerance of the analysis.
    '''
    ppms = []
    ppm_default = General_Functions.calculate_ppm_diff(tolerance[2]-General_Functions.tolerance_calc(tolerance[0], tolerance[1], tolerance[2]), tolerance[2])
    missing_points = 0
    for i in ppm_array[peak['peak_interval_id'][0]:peak['peak_interval_id'][1]+1]:
        if i != inf:
            ppms.append(i)
        else:
            ppms.append(ppm_default)
            missing_points+= 1
    return mean(ppms), missing_points

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
    counter = 0
    max_counter = int(0.1/(rt_int[0][-1]-rt_int[0][-2]))
    for i_i, i in enumerate(rt_int_smoothed[1]):
#        print("RT:", "%.2f" % round(rt_int[0][i_i], 2),"INT:", "%.2f" % round(rt_int[1][i_i], 2), "TEMP_MAX:", temp_max,"TEMP_MAX_ID:", temp_max_id_iu, going_up, going_down, counter)
        if (rt_int[0][i_i] >= rt_interval[1] or rt_int[0][i_i] == rt_int[0][-1]):
            break
        if rt_int[0][i_i] >= rt_interval[0]:
            if (going_up or going_down) and rt_int[1][i_i] > temp_max:
                temp_max = rt_int[1][i_i]
                temp_max_id_iu = i_i
            if (going_up and (i < rt_int_smoothed[1][i_i-1] or rt_int_smoothed[1][i_i] == 0)):
                counter+=1
                if counter <= max_counter:
                    continue
                elif counter > max_counter:
                    counter = 0
                    going_up = False
                    going_down = True
            if (going_down and (i > rt_int_smoothed[1][i_i-1] or rt_int_smoothed[1][i_i] == 0.0)):
                temp_peak_width = (rt_int[0][i_i]-rt_int[0][temp_start])
                going_down = False
                good = False
                if min_ppp[0]:
                    if i_i-temp_start >= min_ppp[1] or glycan == "Internal Standard":
                        good = True
                else:
                    if i_i-temp_start >= int(0.2/(rt_int[0][-1]-rt_int[0][-2])) or glycan == "Internal Standard":
                        good = True
                if good:
                    peaks.append({'id': temp_max_id_iu, 'rt': rt_int[0][temp_max_id_iu], 'int': temp_max, 'peak_width': temp_peak_width, 'peak_interval': (rt_int[0][temp_start], rt_int[0][i_i]), 'peak_interval_id': (temp_start, i_i)})
                    temp_max = 0
                    temp_max_id_iu = 0
                if not good:
                    temp_max = 0
                    temp_max_id_iu = 0
            if (i_i != 0 and i > rt_int_smoothed[1][i_i-1] and not going_up and not going_down):
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
