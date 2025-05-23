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

def main(args=[]):
    try:
        if len(args) == 0:
            if not os.isatty(0):
                args = Config_Handler.config_handler()
            else:
                args = CLI.CLI()
            
        begin_time = datetime.datetime.now()
        begin_time_string = str(begin_time)[2:4]+str(begin_time)[5:7]+str(begin_time)[8:10]+"_"+str(begin_time)[11:13]+str(begin_time)[14:16]+str(begin_time)[17:19]
        temp_folder = os.path.join(tempfile.gettempdir(), "gg_"+begin_time_string)
        os.makedirs(temp_folder, exist_ok=True)
        
        args[0][23] = temp_folder

        if args[9]:
            Execution_Functions.output_filtered_data(*args[0])

        else:
            args[1][30] = temp_folder
            library = Execution_Functions.imp_exp_gen_library(*args[1])
            Execution_Functions.print_sep()
            pre_processing_begin_time = datetime.datetime.now()
            time_formatted = str(pre_processing_begin_time).split(" ")[-1].split(".")[0]+" - "
            print(time_formatted+"Starting pre-processing...")
            print(time_formatted+"Loading files...", end = "", flush = True)
            data = Execution_Functions.list_of_data(*args[2])
            print("Done!")
            time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "
            print(time_formatted+"Indexing spectra...", end = "", flush = True)
            args[3][0] = data
            ms1_index = Execution_Functions.index_spectra_from_file(*args[3])
            if args[10]:
                data = Execution_Functions.list_of_data(*args[2])
                args[4][0] = data
                ms2_index = Execution_Functions.index_spectra_from_file(*args[4])
            print("Done!")
            lib_size = len(library)
            args[5][0] = library
            args[5][1] = lib_size
            args[5][2] = data
            args[5][3] = ms1_index
            args[5][13] = pre_processing_begin_time
            args[5][14] = temp_folder
            analyzed_data = Execution_Functions.analyze_files(*args[5])
            noise = [analyzed_data[2], analyzed_data[3]]
            if args[10]:
                # Safeguard MS1 save in case of crashes during MS2 analysis
                if datetime.datetime.now()-begin_time > datetime.timedelta(hours=2):
                    Execution_Functions.print_sep()
                    time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "   
                    print(time_formatted+"Saving MS1 .gg file as safeguard\nin case of crashes.")
                    temp_args_raw_data = copy.deepcopy(args[7])
                    temp_args_raw_data[0] = analyzed_data
                    temp_args_raw_data[2] = False
                    temp_args_raw_data[4][1][2] = [False, False, False]
                    temp_args_raw_data[5] = library
                    temp_args_raw_data[6] = temp_folder
                    temp_args_raw_data[7] = '<date>_<time>_MS1_Analysis_safeguard'
                    temp_args_raw_data[8] = False
                    temp_args_raw_data[9] = False
                    temp_args_raw_data[10] = noise
                    Execution_Functions.arrange_raw_data(*temp_args_raw_data)
                Execution_Functions.print_sep()
                args[6][0] = ms2_index
                args[6][1] = data
                args[6][2] = analyzed_data
                args[6][26] = library
                args[6][27] = temp_folder
                analyzed_data = Execution_Functions.analyze_ms2(*args[6])
            Execution_Functions.print_sep()
            args[7][0] = analyzed_data
            args[7][5] = library
            args[7][6] = temp_folder
            args[7][10] = noise
            temp_time = Execution_Functions.arrange_raw_data(*args[7])
            Execution_Functions.print_sep()
            args[0][21] = temp_time
            if len(args) == 11:
                Execution_Functions.output_filtered_data(*args[0])
                
        time_formatted = str(datetime.datetime.now()).split(" ")[-1].split(".")[0]+" - "   
        print(time_formatted+'Finished! Time elapsed: '+str(datetime.datetime.now() - begin_time).split(".")[0])
        shutil.rmtree(temp_folder)
        
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