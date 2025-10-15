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
from . import Config_Handler
from . import CLI
import time
import os
import datetime
import tempfile
import copy
import shutil

#-----------------------------------------------------------------------------

def main(args={}):
    try:
        # Check whether is running interactively or not
        if len(args) == 0:
            from_GUI = False
            if not os.isatty(0):
                args = Config_Handler.config_handler()
            else:
                args = CLI.CLI()
        else:
            from_GUI = True
        
        # Set the start time for the analysis
        begin_time = datetime.datetime.now()
        begin_time_string = str(begin_time)[2:4]+str(begin_time)[5:7]+str(begin_time)[8:10]+"_"+str(begin_time)[11:13]+str(begin_time)[14:16]+str(begin_time)[17:19]
        
        # Create a temporary folder for files created during the run
        temp_folder = os.path.join(tempfile.gettempdir(), "gg_"+begin_time_string)
        os.makedirs(temp_folder, exist_ok=True)
        
        # Add the temporary folder to the args of output_filtered_data
        args["output_filtered_data_args"]['temporary folder'] = temp_folder
        
        # If only reanalysis is set, only reanalyze, else do the complete workflow
        if args["output_filtered_data_args"]['reanalysis']:
            Execution_Functions.output_filtered_data(args["output_filtered_data_args"])
        
        else:
            # Add the temporary folder to the imp_exp_gen_library args
            args['imp_exp_gen_library_args']['temporary folder'] = temp_folder
            
            # Generate the library
            library, adduct_combos = Execution_Functions.imp_exp_gen_library(args['imp_exp_gen_library_args'])
            
            # Print a separator
            Execution_Functions.print_sep()
            
            # Save the beginning time for pre processing and output it
            pre_processing_begin_time = datetime.datetime.now()
            time_formatted = str(pre_processing_begin_time).split(" ")[-1].split(".")[0]+" - "
            print(time_formatted+"Starting pre-processing...")
            print(time_formatted+"Loading files...", end = "", flush = True)
            
            # Load the data from the raw spectra files
            data = Execution_Functions.list_of_data(args['list_of_data_args'])
            print("Done!")
            
            # Print the updated time and start indexing the spectra
            time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
            print(time_formatted+"Indexing spectra...", end = "", flush = True)
            
            # Add the loaded data to the spectra indexing args
            args['index_spectra_from_file_ms1_args']['raw data'] = data
            
            # Index the ms1 spectra
            ms1_index = Execution_Functions.index_spectra_from_file(args['index_spectra_from_file_ms1_args'])
            
            # If set to analyze ms2 spectra, reload the data and index ms2 spectra
            if args["output_filtered_data_args"]['analyze ms2']:
                data = Execution_Functions.list_of_data(args['list_of_data_args'])
                args['index_spectra_from_file_ms2_args']['raw data'] = data
                ms2_index = Execution_Functions.index_spectra_from_file(args['index_spectra_from_file_ms2_args'])
                
            print("Done!")
            
            # Check the size of the library
            lib_size = len(library)
            
            # Add the necessary args to the analyze_files function args
            args['analyze_files_args']['library'] = library
            args['analyze_files_args']['library size'] = lib_size
            args['analyze_files_args']['raw data'] = data
            args['analyze_files_args']['ms1 index'] = ms1_index
            args['analyze_files_args']['analysis start time'] = pre_processing_begin_time
            args['analyze_files_args']['temporary folder'] = temp_folder
            args['analyze_files_args']['adduct combos'] = adduct_combos
            
            # Analyze the data
            analyzed_data = Execution_Functions.analyze_files(args['analyze_files_args'])
            
            # Take the calculated noise levels from the analyzed data
            noise = [analyzed_data[2], analyzed_data[3]]
            
            # If set to analyze ms2, start the workflow
            if args["output_filtered_data_args"]['analyze ms2']:
                
                # Safeguard MS1 save in case of crashes during MS2 analysis
                if datetime.datetime.now()-begin_time > datetime.timedelta(hours=2):
                    
                    # Print a separator
                    Execution_Functions.print_sep()
                    
                    # Print the current time and start saving the .gg file
                    time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "   
                    print(time_formatted+"Saving MS1 .gg file as safeguard\nin case of crashes.")
                    
                    # Make a deep copy of the arrange_raw_data args
                    temp_args_raw_data = copy.deepcopy(args['arrange_raw_data_args'])
                    
                    # Add the required args to it
                    temp_args_raw_data['analyzed data'] = analyzed_data
                    temp_args_raw_data['ms2 analyzed'] = False
                    temp_args_raw_data['analysis parameters']['analysis settings']['analyze ms2'] = [False, False, False]
                    temp_args_raw_data['adduct combos'] = adduct_combos
                    temp_args_raw_data['library'] = library
                    temp_args_raw_data['temporary folder'] = temp_folder
                    temp_args_raw_data['results file name'] = '<date>_<time>_MS1_Analysis_safeguard'
                    temp_args_raw_data['from GUI'] = False
                    temp_args_raw_data['erase files'] = False
                    temp_args_raw_data['calculated noise levels'] = noise
                    
                    # Save the .gg file
                    Execution_Functions.arrange_raw_data(temp_args_raw_data)
                    
                # Print a separator
                Execution_Functions.print_sep()
                
                # Add the required args to analyze ms2
                args['analyze_ms2_args']['ms2 index'] = ms2_index
                args['analyze_ms2_args']['raw data'] = data
                args['analyze_ms2_args']['analyzed data'] = analyzed_data
                args['analyze_ms2_args']['library'] = library
                args['analyze_ms2_args']['adduct combos'] = adduct_combos
                args['analyze_ms2_args']['temporary folder'] = temp_folder
                
                # Start ms2 analysis
                analyzed_data = Execution_Functions.analyze_ms2(args['analyze_ms2_args'])
                
            # Print a separator
            Execution_Functions.print_sep()
            
            # Add the required args to arrange the raw data and save .gg file
            args['arrange_raw_data_args']['analyzed data'] = analyzed_data
            args['arrange_raw_data_args']['adduct combos'] = adduct_combos
            args['arrange_raw_data_args']['library'] = library
            args['arrange_raw_data_args']['temporary folder'] = temp_folder
            args['arrange_raw_data_args']['calculated noise levels'] = noise
            
            # Arrange the raw data and generate .gg file
            temp_time = Execution_Functions.arrange_raw_data(args['arrange_raw_data_args'])
            
            # Print a separator
            Execution_Functions.print_sep()
            
            # Add the .gg file save time to the output filtered data args
            args['output_filtered_data_args']['analysis done time'] = temp_time
            if not from_GUI:
                Execution_Functions.output_filtered_data(args['output_filtered_data_args'])
                
        # Get the current time and output the conclusion, with elapsed time
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "   
        print(time_formatted+'Finished! Time elapsed: '+str(datetime.datetime.now() - begin_time).split(".")[0])
        
        # Remove the temporary folder
        shutil.rmtree(temp_folder)
        
        # If running interactively, press enter to exit, else just close
        if os.isatty(0):
            input('\nPress Enter to exit.')
        else:
            print("Close the window or press CTRL+C to exit.")
            try:
                while True:
                    time.sleep(3600)
            except KeyboardInterrupt:
                os._exit(1)
    except KeyboardInterrupt:
        print("\n\n----------Execution cancelled by user.----------\n", flush=True)
        raise SystemExit(1)