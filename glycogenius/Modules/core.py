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

#-----------------------------------------------------------------------------

def main():
    try:
        if not os.isatty(0):
            args = Config_Handler.config_handler()
        else:
            args = CLI.CLI()
            
        begin_time = datetime.datetime.now()

        if args[9]:
            Execution_Functions.output_filtered_data(*args[0])

        else:
            library = Execution_Functions.imp_exp_gen_library(*args[1])
            print('Library length: '+str(len(library)))
            Execution_Functions.print_sep()
            print("Starting pre-processing...")
            print("Preparing files for analysis...", end = "", flush = True)
            data = Execution_Functions.list_of_data(*args[2])
            print("Done!")
            print("Indexing spectra...", end = "", flush = True)
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
            analyzed_data = Execution_Functions.analyze_files(*args[5])
            if args[10]:
                Execution_Functions.print_sep()
                args[6][0] = ms2_index
                args[6][1] = data
                args[6][2] = analyzed_data
                analyzed_data = Execution_Functions.analyze_ms2(*args[6])
            Execution_Functions.print_sep()
            args[7][0] = analyzed_data
            temp_time = Execution_Functions.arrange_raw_data(*args[7])
            Execution_Functions.print_sep()
            args[0].append(temp_time)
            Execution_Functions.output_filtered_data(*args[0])
                                                     
        if os.isatty(0):
            input('Execution complete. Time elapsed: '+str(datetime.datetime.now() - begin_time)+'\nPress Enter to exit.')
        else:
            print('Execution complete. Time elapsed: '+str(datetime.datetime.now() - begin_time))
            print("Close the window or press CTRL+C to exit.")
            try:
                while True:
                    time.sleep(3600)
            except KeyboardInterrupt:
                os._exit(1)
    except KeyboardInterrupt:
        print("\n\n----------Execution cancelled by user.----------\n", flush=True)
        raise SystemExit(1)