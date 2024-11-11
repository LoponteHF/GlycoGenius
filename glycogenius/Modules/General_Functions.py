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

from pyteomics import mzxml, mzml, mass, auxiliary
from itertools import combinations_with_replacement
from numpy import percentile, arange, zeros, array, polyfit, std, where
from re import split
from math import inf, atan, acos, exp, pi
from statistics import stdev, mean
from scipy.stats import linregress
from scipy.sparse.linalg import splu
from scipy import sparse
import dill
import numpy
import sys
import datetime
import xlsxwriter
import zipfile
import pathlib
import tempfile
import shutil
import os
import copy

##---------------------------------------------------------------------------------------
##Hard-coded permanent information

monosaccharides = {
    "H": ("Hexose", "C6O6H12", {"C": 6, "O": 5, "N": 0, "H": 10}, "H"),
    "N": ("N-Acetyl Hexosamine", "C8O6NH15", {"C": 8, "O": 5, "N": 1, "H": 13}, "N"),
    "X": ("Xylose", "C5O5H10", {"C": 5, "O": 4, "N": 0, "H": 8}, "X"),
    "S": ("Acetyl Neuraminic Acid", "C11O9NH19", {"C": 11, "O": 8, "N": 1, "H": 17}, "S"),
    "Am": ("Lactonized Acetyl Neuraminic Acid alpha2,3 bound", "C11O8N2H20", {"C": 11, "O": 7, "N": 2, "H": 18}, "L"),
    "E": ("Ethyl-Esterified Acetyl Neuraminic Acid alpha2,6 bound", "C13O9NH23", {"C": 13, "O": 8, "N": 1, "H": 21}, "E"),
    "F": ("Deoxyhexose", "C6O5H12", {"C": 6, "O": 4, "N": 0, "H": 10}, "F"),
    "G": ("Glycolyl Neuraminic Acid", "C11O10NH19", {"C": 11, "O": 9, "N": 1, "H": 17}, "G"),
    "AmG": ("Lactonized Glycolyl Neuraminic Acid alpha2,3 bound", "C11O8N2H20", {"C": 11, "O": 8, "N": 2, "H": 18}, "A"),
    "EG": ("Ethyl-Esterified Glycolyl Neuraminic Acid alpha2,6 bound", "C13O9NH23", {"C": 13, "O": 9, "N": 1, "H": 21}, "R"),
    "HN": ("Hexosamine", "C6O5NH13", {"C": 6, "O": 4, "N": 1, "H": 11}, "M"),
    "UA": ("Uronic Acid", "C6O7H10", {"C": 6, "O": 6, "N": 0, "H": 8}, "U")
    }
'''A hardcoded dictionary containing each single letter code for monosaccharides as key
and a tuple containing the full monosaccharide name, its full molecular formula and its
residue composition in dict form as value.
'''

default_composition = {"H": 0, "N": 0, "X": 0, "S": 0, "Am": 0, "E": 0, "F": 0, "G": 0, "AmG": 0, "EG": 0, "HN": 0, "UA": 0}

h_mass = mass.calculate_mass(composition={'H' : 1})
'''The mass of an hydrogen-1 atom. Pre-calculated here to avoid calculating too many
times during a run.
'''

##---------------------------------------------------------------------------------------
##General functions (these functions use only external libraries, such as itertools and
##pyteomics).

def binary_search_with_tolerance(arr, target, low, high, tolerance, int_arr = [], black_list = []):
    '''A function to quickly find a target in an array by splitting the array recurssively in two and looking for the mid point, finding out if value is bigger or smaller than  target, then splitting again. It also checks if found target in target array is within a tolerance, and picks the most intense one within the tolerance if intensity array is available, else picks the closest one to target.
    
    Parameters
    ----------
    arr : list/np.array
        Target array to search for target.
        
    target : float
        Float to find in target array.
        
    low : int
        Index of the first element in target array.
        
    high : int
        Index of the last element in target array.
        
    tolerance : float
        Tolerance to check for target in target array.
        
    int_arr : list/np.array
        List of intensities, synchronized with arr, to help choose the best target within tolerance.
        
    Uses
    ----
    numpy.argmax : int
        Outputs the index of the highest value in a given array.
        
    numpy.argmin : int
        Outputs the index of the value with smallest difference to a given target, in an array.
        
    Returns
    -------
    selected_id : index
        The index of the selected target.
    '''
    # Base case: if the range is invalid, the target is not in the array
    if low > high:
        return -1  # Target not found
    
    # Find the middle index
    mid = (low + high) // 2
    
    # Check if the target is within the tolerance range of the middle element
    if abs(arr[mid] - target) <= tolerance:
    
        range_width = 5
        range_search = [mid, mid+1]
        for i in range(mid, mid-range_width, -1):
            if i == 0 or arr[i] < target-tolerance:
                break
            range_search[0] = i
        for i in range(mid+1, mid+range_width+1):
            range_search[1] = i
            if i > high or arr[i] > target+tolerance:
                break
                
        if len(int_arr) != 0:
            array_slice = int_arr[range_search[0]:range_search[1]]
        else:
            array_slice = arr[range_search[0]:range_search[1]]
            
        if len(array_slice) == 0:
            return -1

        if len(int_arr) != 0:
            relative_id = numpy.argmax(array_slice)
        else:
            relative_id = (numpy.abs(array_slice - target).argmin())
        selected_id = range_search[0]+relative_id
            
        # This avoids picking the same peak twice
        forbidden_ids = []
        while arr[selected_id] in black_list:
            forbidden_ids.append(relative_id)
            if len(int_arr) != 0:
                array_slice[relative_id] = 0
            else:
                array_slice[relative_id] += 1000
            if len(int_arr) != 0:
                relative_id = numpy.argmax(array_slice)
            else:
                relative_id = (numpy.abs(array_slice - target).argmin())
            if relative_id in forbidden_ids:
                return -1
            selected_id = range_search[0]+relative_id
            
        return selected_id
    elif arr[mid] < target:
        # If target is greater, ignore the left half
        return binary_search_with_tolerance(arr, target, mid + 1, high, tolerance, int_arr, black_list)
    else:
        # If target is smaller, ignore the right half
        return binary_search_with_tolerance(arr, target, low, mid - 1, tolerance, int_arr, black_list)

def linear_regression(x, y, th = 2.5):
    '''Traces a linear regression of supplied 2d data points and returns the slope,
    y-intercept and the indices of the outliers outside the determined threshold.
       
    Parameters
    ----------
    x : list
        A list containing the x elements of the datapoints.
        
    y : list
        A list containing the y elements of the datapoints.
        
    th : float
        The threshold for calculating outliers.
        
    Uses
    ----
    numpy.array : ndarray
        Transforms a list into a numpy array.
        
    scipy.stats.linregress : tuple
        Takes the x and y arrays as inputs and returns a series of results in a
        tuple, such as the slope and the y-intercept of the fit of the datapoints
        to a linear equation.
        
    Returns
    -------
    m : float
        Slope of the fitted linear equation curve.
        
    b : float
        y-intercept of the fitted linear equation curve.
        
    outlier_indices : list
        An index list of the residuals outside the given threshold of the linear equation.
    '''
    if len(x) != len(y): # Ensure x and y have the same length
        raise ValueError("Input arrays x and y must have the same length.")
    if len(x) == 0:
        return 0, 0, []
    if len(x) == 1:
        return 0, y[0], []
    x = array(x) # Convert input data to numpy arrays
    y = array(y) # Convert input data to numpy arrays
    m, b, r_value, p_value, std_err = linregress(x, y) # Calculate the slope (m) and y-intercept (b) using numpy's polyfit function
    predicted_y = m * x + b
    residuals = y - predicted_y
    threshold = th * std(residuals)
    outlier_indices = where(abs(residuals) > threshold)[0]
    # Final equation is used as (y = mx + b)
    return m, b, outlier_indices
    
def cleanup_by_average(x, y, th = 2.5):
    '''Used to find the outliers of a given fit average.
    
    Parameters
    ----------
    x : list
        A list containing the x elements of the datapoints.
        
    y : list
        A list containing the y elements of the datapoints.
        
    th : float
        The threshold for calculating outliers.
    
    Uses
    ----
    numpy.array : ndarray
        Transforms a list into a numpy array.
        
    Returns
    -------
    mean : float
        The mean of the fit.
        
    residuals : float
        Values of residuals.
        
    outlier_indices : list
        An index list of the residuals outside the given threshold.
    '''
    if len(x) != len(y): # Ensure x and y have the same length
        raise ValueError("Input arrays x and y must have the same length.")
    if len(x) == 0:
        return 0, 0, []
    if len(x) == 1:
        return 0, y[0], []
    x = array(x) # Convert input data to numpy arrays
    y = array(y) # Convert input data to numpy arrays
    predicted_y = mean(y)
    residuals = y - predicted_y
    threshold = th * std(residuals)
    outlier_indices = where(abs(residuals) > threshold)[0]
    # Final equation is used as (y = mx + b)
    return mean, residuals, outlier_indices
    
def make_gg(temp_dir, save_dir, filename, ignore_files = False):
    '''Creates a zipped file with extension .gg containing raw_data files.
    
    Parameters
    ----------
    temp_dir : string
        A string containing the path to the temp_dir where raw_data files are.
        
    save_dir : string
        A string containing the path to save the .gg file.
        
    filename : string
        The name of the .gg file.
        
    Uses
    ----
    zipfile.ZipFile
        Creates a zipfile object to add files to.
    
    zipfile.write
        Adds files to the zipfile object.
        
    Returns
    -------
    nothing
        Creates filename.gg file, containing files in temp_dir, at the save_dir.
    '''
    files_list = os.listdir(temp_dir)
    with zipfile.ZipFile(os.path.join(save_dir, filename+".gg"), 'w', compression=zipfile.ZIP_DEFLATED) as zipf:
        for file in files_list:
            if ignore_files:
                if file != "metadata" and file != "results" and "_" not in file:
                    continue
            zipf.write(os.path.join(temp_dir, file), arcname=file)
        zipf.close()
            
def open_gg(gg_file, temp_path, file = 'all'):
    '''Unzipts the .gg file to temp_path.
    
    Parameters
    ----------
    gg_file : string
        Path to the .gg file.
        
    temp_path : string
        String containing the path to the folder where the files will be extracted.
        
    Uses
    ----
    zipfile.ZipFile
        Creates a zipfile object to add files to.
        
    zipfile.extractall
        Extract all files in .gg file to the target directory.
        
    Returns
    -------
    nothing
        Extracts raw_data files from gg_file to temp_path.
    '''
    with zipfile.ZipFile(gg_file, 'r') as zip_ref:
        if file == 'all':
            zip_ref.extractall(temp_path)
        else:
            zip_ref.extract(file, temp_path)
        zip_ref.close()
        
def access_chromatogram(file_number, chromatogram_name, temp_folder, reanalysis_path):
    '''
    '''
    if chromatogram_name not in os.listdir(temp_folder):
        if f'{file_number}_eics' not in os.listdir(temp_folder):
            open_gg(reanalysis_path, temp_folder, f'{file_number}_eics')
        open_gg(os.path.join(temp_folder, f'{file_number}_eics'), temp_folder, chromatogram_name)
    with open(os.path.join(temp_folder, chromatogram_name), 'rb') as f:
        chromatogram = dill.load(f)
        f.close()
    return chromatogram

def calculate_ppm_diff(mz, target):
    '''Calculates the PPM difference between a mz and a target mz.
    
    Parameters
    ----------
    mz : float
        A float of the mz you want to check for the difference.
        
    target : float
        A float of the target mz you're comparing your mz to.
        
    Returns
    -------
    float
        A float of the PPM difference between the mz and target mz.
    '''
    return ((target-mz)/target)*(10**6)
    
def tolerance_calc(unit,
                   value, 
                   mz = 1000.0):
    '''An accurate way to convert 'ppm' mass accuracy into 'pw' (peak width, aka. mz
    tolerance).

    Parameters
    ----------
    unit : string
        Can be either "ppm" (particles-per-million) or "pw" (peak width [tolerance]).

    value : float
        Float value of the tolerance, based on the unit inputted.
        
    mz : float
        mz at which to calculate the PPM difference.

    Returns
    -------
    tolerance : float
        If unit == "ppm", converts value into "pw", if unit == "pw", outputs value as is.
        ie. 10 ppm gets converted to 0.01 pw tolerance.
    '''
    if unit == "ppm":
        return (mz-(-((value*mz)/10**6)+mz))
    elif unit == "mz":
        return value
    else:
        return("Unit for tolerance not 'ppm' or 'mz'.")
        
def autofit_columns_excel(df, worksheet):
    '''Autofits the column width in a excel worksheet based on a dataframe used to make it.
    
    Parameters
    ----------
    df : Pandas Dataframe
        Dataframe containing the data used to make the worksheet with Pandas.
        
    worksheet : Worksheet XLSXwriter object
        The worksheet to autofit the columns.
        
    Uses
    ----
    xlsxwriter
        Library to write and edit excel files and objects.
        
    Returns
    -------
    nothing
        Directly edits the worksheet.
    '''
    for idx, col in enumerate(df):  # loop through all columns
        series = df[col]
        max_len = max((
            series.astype(str).map(len).max(),  # len of largest item
            len(str(series.name))  # len of column name/header
            )) + 1  # adding a little extra space
        worksheet.set_column(idx, idx, max_len)
        
def speyediff(N, d, format='csc'):
    '''Construct a d-th order sparse difference matrix based on an initial 
    N x N identity matrix. Obtained from https://github.com/mhvwerts/whittaker-eilers-smoother, 
    applied as is. Due credits given.
    
    Parameters
    ----------
    N : int
        Length of vector containing raw data
    
    d : int
        Order of smoothing
        
    format : string
        To be used by scipy.sparse.diags
        
    Uses
    ----
    numpy.zeros : ndarray
        Return a new array of given shape and type, filled with zeros.
    
    numpy.arange : ndarray
        Return evenly spaced values within a given interval.
    
    scipy.sparse.diags :  ndarrays
        Construct a sparse matrix from diagonals.
    
    Returns
    -------
    spmat : ndarray
        Final matrix (N-d) x N
    '''
    
    assert not (d < 0), "d must be non negative"
    shape     = (N-d, N)
    diagonals = zeros(2*d + 1)
    diagonals[d] = 1.
    for i in range(d):
        diff = diagonals[:-1] - diagonals[1:]
        diagonals = diff
    offsets = arange(d+1)
    spmat = sparse.diags(diagonals, offsets, shape, format=format)
    return spmat        
    
def rt_noise_level_parameters_set(mz_int, mode):
    '''Receives 2 combined arrays containing the x and y information of a spectrum
    and calculate parameters for dynamic noise calculation down the pipeline.
    
    Parameters
    ----------
    mz_int : list
        A list containing two synchronized lists: the first one contains the
        mzs and the second one the intensity.
        
    mode : string
        Whether you're analyzing segments or the whole mz_int.
        
    Uses
    ----
    numpy.std : ndarray
        Compute the standard deviation along the specified axis.
        
    Returns
    -------
    float
        Returns a noise threshold if the operation runs on the whole array.
        
    tuple
        If analyzing segments, gives the noise threshold of the first quarter and
        the last quarter, as well as the last mz in the spectrum.
    '''
    if mode == "segments":
        first_quarter_end = int(len(mz_int[1])/4)
        last_quarter_begin = first_quarter_end*3
        int_list_first_quarter = mz_int[1][:first_quarter_end]
        int_list_last_quarter = mz_int[1][last_quarter_begin:]
        segments_list = [sorted(int_list_first_quarter), sorted(int_list_last_quarter)]
        
    if mode == "whole":
        segments_list = [sorted(list(mz_int[1]))]
    
    noise = []
    for j_j, j in enumerate(segments_list):
        if len(j) == 0:
            noise.append(1.0)
            continue
        intensity_std = numpy.std(j)
        noise_threshold = 2.0 * intensity_std
        if (min(j) != 0 and noise_threshold > min(j)*5) or noise_threshold > max(j)*0.5: #this means that the data is denoised already, so it picks really high intensity as possible noise
            # if mode == "whole": print("picked minimum", min(j), max(j), len(j), noise_threshold) 
            if min(j) != 0:
                noise.append(min(j))
            else:
                noise.append(1.0)
        else:
            # if mode == "whole": print("picked 2 std", min(j), max(j), len(j), noise_threshold) 
            noise.append(noise_threshold)
    if len(noise) == 1:
        return noise[0]
    else:
        return noise[0], noise[1], mz_int[0][-1]
    
def local_noise_calc(noise_specs, x, avg_noise):
    '''Uses the noise_specs produced by rt_noise_level_parameters_set to 
    calculate the local noise levels. If any of the parameters are considered
    abnormal (ie. absurdly high or no peaks on the array) it defaults to avg_noise.
    
    Parameters
    ----------
    noise_specs : tuple
        A tuple containing scalars, as produced by rt_noise_level_parameters_set
        
    x : float
        The mz at which you want to calculate the noise.
        
    avg_noise : float
        The average noise of the file, used for fallback.
        
    Returns
    -------
    float
        A float containing the local or average noise level.
    '''
    if noise_specs[2] == 0.0 or noise_specs[1] > noise_specs[0]*5 or noise_specs[0] > noise_specs[1]*5:
        return avg_noise
    return noise_specs[0] + (((noise_specs[1]-noise_specs[0])/noise_specs[2])*x)
    
def normpdf(x, mean, sd):
    '''Calculates the intensity of a gaussian bell curve at the x-axis point x
    with the set parameters.
    
    Parameters
    ----------
    x : float
        Given point in the x axis of the gaussian bell curve.
        
    mean : float
        Determines the maximum intensity of the mean point of the gaussian
        bell curve.
        
    sd : float
        Determines the standard deviation of the gaussian.
        
    Uses
    ----
    math.pi : float
        The mathematical constant pi = 3.141592…, to available precision.
        
    math.exp : float
        Return e raised to the power x, where e = 2.718281… is the base of 
        natural logarithms.
        
    Returns
    -------
    float
        A float containing the intensity of the gaussian bell curve at the 
        given x point.
    '''
    var = float(sd)**2
    denom = (2*pi*var)**.5
    num = exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom

def form_to_comp(string, form_type = 'glycan'):
    '''Separates a molecular formula or monosaccharides formula of glycans into a
    dictionary with each atom/monosaccharide as a key and its amount as value.

    Parameters
    ----------
    string : str
        A string in the form of C6O6N0H12 or H5N4S1F1G1.

    Returns
    -------
    counts : dict
        A dictionary with keys containing the monosaccharides/atoms letter(s) and values
        containing the amounts of each. ie. {"H": 5, "N": 4, "S": 1, "F": 1, "G": 1}.
    '''
    string = string.split("+")[0] #avoids getting phosphorylation and sulfation symbols
    counts = {}
    split_str = split('(\\d+)', string)
    negative = False
    for i_i, i in enumerate(split_str):
        if i != '' and i[-1] == '-':
            split_str[i_i] = i[:-1]
            negative = True
        if i_i%2 != 0 and i != '' and negative:
            split_str[i_i] = '-'+i
    if len(split_str)%2 != 0:
        split_str.append('1')
    for i in range(len(split_str)-1):
        if i%2 == 0:
            counts[split_str[i]] = int(split_str[i+1])
    if '' in counts:
        del counts['']
    return counts

def form_to_charge(string):
    '''Converts adducts formula into raw charge.

    Parameters
    ----------
    string : str
        A string containing adducts formula.

    Uses
    ----
    form_to_comp() : dict
        Separates a molecular formula or monosaccharides formula of glycans into a
        dictionary with each atom/monosaccharide as a key and its amount as value

    Returns
    -------
    charge : int
        The raw charge of the adduct.
    '''
    comp = form_to_comp(string)
    charge = 0
    for i in comp:
        if i != '':
            charge+=comp[i]
    return charge

def glycan_to_atoms(glycan_composition, permethylated):
    '''Calculates the amounts of atoms based on glycan monosaccharides.

    Parameters
    ----------
    glycan_composition : dict
        Accepts as input the glycan monosaccharides formula in the format of {"H": 5,
        "N": 4, "S": 1, "F": 1, "G": 1}.
        
    permethylated : boolean
        Whether the glycans are permethylated or not.

    Uses
    ----
    monosaccharides : dict
        A hardcoded dictionary containing each single letter code for monosaccharides as
        key and a tuple containing the full monosaccharide name, its full molecular
        formula and its residue composition in dict form.

    Returns
    -------
    atoms : dict
        Returns a dictionary with atoms as keys and amounts as values. ie. {"C": 6,
        "O": 6, "N": 0, "H": 12}.
    '''
    atoms = {"C": 0, "O": 0, "N": 0, "H": 0}
    monosaccharides_local = copy.deepcopy(monosaccharides)
    if permethylated:
        for i in monosaccharides_local:
            if i == 'H' or i == 'N':
                monosaccharides_local[i][2]['C'] += 3
                monosaccharides_local[i][2]['H'] += 6
            if i == 'F':
                monosaccharides_local[i][2]['C'] += 2
                monosaccharides_local[i][2]['H'] += 4
            if i == 'S' or i == 'G':
                monosaccharides_local[i][2]['C'] += 5
                monosaccharides_local[i][2]['H'] += 10
    for i in glycan_composition:
        if i == "T":
            continue
        for j in atoms:
            atoms[j] += monosaccharides_local[i][2][j]*glycan_composition[i]
    return atoms

def count_seq_letters(string):
    '''If you make anything with itertools for combinatorial analysis, it will produce a
    string that's not very human readable. This converts it into a human readable form.

    Parameters
    ----------
    string : str
        A string of atoms or glycans in the form o CCCCCCOOOONH or HHHHHNNNNFSG.

    Returns
    -------
    friendly_letters : dict
        A dictionary containing the count for each letter in the string. ie. CCCCOONH
        returns {"C": 4, "O": 2, "N": 1, "H": 1}.
    '''
    friendly_letters = {}
    current_letter = ""
    for i in string:
        if i != current_letter:
            current_letter = i
            count = string.count(i)
            if i == "L":
                friendly_letters['Am'] = count
            if i == "A":
                friendly_letters['AmG'] = count
            if i == "R":
                friendly_letters['EG'] = count
            else:
                friendly_letters[i] = count
    return friendly_letters

def sum_atoms(*compositions):
    '''Sums the atoms of two compositions.

    Parameters
    ----------
    compositions : dict
        Dictionaries containing the atomic compositions of the items to be summed.
        ie. {"C": 4, "O": 2, "N": 1, "H": 1}.

    Returns
    -------
    summed_comp : dict
        Dictionary containing the sum of each atom of the compositions.
    '''
    summed_comp = {"C": 0, "O": 0, "N": 0, "H": 0}
    for i in compositions:
        for j in i:
            if j not in list(summed_comp.keys()):
                summed_comp[j] = i[j]
            else:
                summed_comp[j]+=i[j]
    return summed_comp

def sum_monos(*compositions):
    '''Sums the monosaccharides of two glycan compositions. 'T' stands for TAG, used
    in fragments library calculation.

    Parameters
    ----------
    compositions : dict
        Dictionaries containing the monosaccharrides compositions of the glycans to be
        summed. ie. {"H": 5, "N": 4, "S": 1, "F": 1, "G": 1}.

    Returns
    -------
    summed_comp : dict
        Dictionary containing the sum of each monosaccharides of the compositions.
    '''
    summed_comp = {"H": 0, "N": 0, "X": 0, "S": 0, "Am": 0, "E": 0, "F": 0, "G": 0, "AmG": 0, "EG": 0, "T": 0, "HN": 0, "UA": 0}
    for i in compositions:
        for j in i:
            summed_comp[j]+=i[j]
    return summed_comp

def comp_to_formula(composition):
    '''Transforms a composition dictionary into string formula.

    Parameters
    ----------
    composition : dict
        Dictionary containing the composition of the molecule or glycan.

    Returns
    -------
    formula : string
        Formula of the atomic or monosaccharides composition in string form.
    '''
    formula = ''
    for i in composition:
        if composition[i] != 0:
            formula+=i+str(composition[i])
    return formula

def calculate_comp_from_mass(tag_mass):
    '''Calculates the composition of a molecule based on its mass. Intended to use with
    small tags added to the glycans.

    Parameters
    ----------
    tag_mass : float
        The monoisotopic molecular weight of the tag's added mass.

    Uses
    ----
    pyteomics.mass.calculate_mass(*args, **kwargs) : float
        Calculates the monoisotopic mass of a polypeptide defined by a sequence string,
        parsed sequence, chemical formula or Composition object.

    itertools.combinations_with_replacement : generator
        Return r length subsequences of elements from the input iterable allowing
        individual elements to be repeated more than once.

    Returns
    -------
    closest : tuple
        Returns the proposed composition of the molecule with the calculated mass closest
        to the tag's added mass in dictionary form and the calculated hypothetical mass
        of the molecule as float.
    '''
    elements = "CONH"
    closest = ({}, float('inf'))
    
    # Estimate minimum and maximum number of carbon atoms
    min_carbons = int(((0.05 * tag_mass) + 0.5) * 0.7)
    max_carbons = int(((0.05 * tag_mass) + 0.5) * 1.3)

    max_nitrogens = max(3, int(tag_mass/200))  # Maximum number of nitrogen atoms allowed
    min_atoms = max(min_carbons, int(tag_mass/10))  # Min total atoms considered per combination
    max_atoms = int(tag_mass/5)  # Max total atoms considered per combination
    
    # print(f"Carbon range: {min_carbons} to {max_carbons}, Max Nitrogens: {max_nitrogens}, Min atoms: {min_atoms}, Max atoms: {max_atoms}")

    # Generate combinations within the carbon range
    for c_count in range(min_carbons, max_carbons + 1):
        for n_count in range(0, max_nitrogens + 1):
            for total_atoms in range(min_atoms, max_atoms + 1):
                non_carbon_nitrogen_count = total_atoms - c_count
                if non_carbon_nitrogen_count < n_count:
                    continue
                
                # Generate combinations of 'O' and 'H', with fixed 'N' count
                combos = combinations_with_replacement('OH', non_carbon_nitrogen_count - n_count)
                
                for combo in combos:
                    seq_readable = count_seq_letters('C' * c_count + 'N' * n_count + ''.join(combo))
                    test_tag_mass = mass.calculate_mass(composition=seq_readable)
                    
                    if abs(test_tag_mass - tag_mass) < abs(closest[1] - tag_mass):
                        closest = (seq_readable, test_tag_mass)

                    # Early stopping if we find a very close match
                    if abs(test_tag_mass - tag_mass) < 1e-3:
                        return closest

    return closest
    
def calculate_isotopic_pattern(glycan_atoms,
                               fast=True,
                               high_res=False):
    '''Calculates up to 5 isotopic pattern peaks relative abundance in relation with the
    monoisotopic one.

    Parameters
    ----------
    glycan_atoms : dict
        A dictionary containing glycan atomic composition in the form of {"C": 4, "O": 2,
        "N": 1, "H": 1}.
    
    fast : boolean
        If True, only calculates the isotopic pattern based on isotopes of carbon and
        nitrogen, thus very inaccurate, but enough for most uses. If fast = False, it
        will take a significantly higher amount of time to produce results (from 10x to
        1000x longer, depending on number of atoms).
        Default = True.
    
    high_res : boolean
        If True, doesn't clump similar masses together. Useful for analysis of data
        obtained on FT MS instruments.

    Uses
    ----
    pyteomics.mass.isotopologues(*args, **kwargs) : iterator
        Iterate over possible isotopic states of a molecule. The molecule can be defined
        by formula, sequence, parsed sequence, or composition.
        
    pyteomics.mass.calculate_mass(*args, **kwargs) : float
        Calculates the monoisotopic mass of a polypeptide defined by a sequence string,
        parsed sequence, chemical formula or Composition object.
    
    Returns
    -------
    relative_isotop_pattern : list
        A list with the isotopic pattern, with the first element being the abundance of
        the monoisotopic peak (1.0, or 100%) and the following ones are the isotopologues
        relative abundance in relation to the monoisotopic peak (around 1 Da apart).
        
    relative_isotop_mass : list
        A list of isotopologue masses synchronized with relative_isotop_pattern list.
    '''
    if fast:
        isotopologue = mass.isotopologues(glycan_atoms, report_abundance = True,
                                          elements_with_isotopes = ["C", "P"],
                                          overall_threshold = 1e-4)
    else:
        isotopologue = mass.isotopologues(glycan_atoms, report_abundance = True,
                                          overall_threshold = 1e-4)
    isotop_arranged = []
    relative_isotop_pattern = []
    relative_isotop_mass = []
    try:
        for i in isotopologue:
            isotop_arranged.append({'mz' : mass.calculate_mass(i[0]), 'Ab' : i[1]})
    except Exception:
        pass
    isotop_arranged = sorted(isotop_arranged, key=lambda x: x['mz'])
    for i_i, i in enumerate(isotop_arranged):
        relative_isotop_pattern.append(i['Ab']/isotop_arranged[0]['Ab'])
        relative_isotop_mass.append(i['mz'])
    if not high_res and not fast:
        relative_isotop_pattern_low_res = []
        relative_isotop_mass_low_res = []
        for i_i, i in enumerate(relative_isotop_mass):
            if i_i == 0:
                relative_isotop_pattern_low_res.append(relative_isotop_pattern[i_i])
                relative_isotop_mass_low_res.append(i)
            else:
                if abs(i-relative_isotop_mass[i_i-1]) < h_mass/2:
                    relative_isotop_mass_low_res[-1] = (relative_isotop_mass_low_res[-1]+i)/2
                    relative_isotop_pattern_low_res[-1]+= relative_isotop_pattern[i_i]
                else:
                    relative_isotop_mass_low_res.append(i)
                    relative_isotop_pattern_low_res.append(relative_isotop_pattern[i_i])
        return relative_isotop_pattern_low_res, relative_isotop_mass_low_res
    return relative_isotop_pattern, relative_isotop_mass

def gen_adducts_combo(adducts,
                      exclusions=[],
                      max_charge=3):
    '''Generates a list of dictionaries with compositions of adducts combinations, based
    on parameters set.

    Parameters
    ----------
    adducts : dict
        A dictionary with each key containing a single atom adduct and its value
        containing the maximum amount of such adduct.
        
    exclusions : list
        A list containing undesired adducts combinations to be excluded from the function
        result.

    max_charge : int
        The maximum amount of charges for the adducts.

    Uses
    ----
    itertools.combinations_with_replacement : generator
        Return r length subsequences of elements from the input iterable allowing
        individual elements to be repeated more than once.

    Returns
    -------
    adducts_combo_dict : list
        A list of dictionaries containing the composition of each adducts combination.
    '''
    adducts_list = []
    adducts_combo = []
    adducts_combo_dict = []
    for i in adducts:
        adducts_list.append(i)
    for i in range(1, abs(max_charge)+1):
        for j in combinations_with_replacement(adducts_list, i):
            adducts_combo.append(j)
    for i in adducts_combo:
        temp_dict = {}
        for j in i:
            if j not in temp_dict:
                if max_charge > 0:
                    temp_dict[j] = 1
                else:
                    temp_dict[j] = -1
            else:
                if max_charge > 0:
                    temp_dict[j]+= 1
                else:
                    temp_dict[j]-= 1
        adducts_combo_dict.append(temp_dict)
    to_remove = []
    for i in adducts_combo_dict:
        charges = 0
        for atom in i:
            if atom == "Cl":
                if j[atom] > 0:
                    charges -= j[atom]
                else:
                    charges += j[atom]
            if atom == 'Fe':
                charges += 2*i[atom]
            else:
                charges += i[atom]
        if max_charge > 0 and charges > max_charge:
            to_remove.append(i)
            continue
        elif max_charge < 0 and charges < max_charge:
            to_remove.append(i)
            continue
        if i in exclusions:
            to_remove.append(i)
            continue
        for j in i:
            if abs(i[j]) > abs(adducts[j]):
                to_remove.append(i)
                break
    for i in to_remove:
        adducts_combo_dict.remove(i)
    return adducts_combo_dict
