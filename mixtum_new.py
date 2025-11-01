#    Mixtum: the geometry of admixture in population genetics.
#    Copyright (C) 2025  Jose Maria Castelo Ares
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

import argparse, sys
from pathlib import Path

from gui.core import Core



# Raise error if file not found
def check_file_path(file_path):
    if not file_path.is_file():
        raise FileNotFoundError(f'File {file_path} does not exist')



# Raise error if dir not found
def check_dir_path(dir_path):
    if not dir_path.is_dir():
        raise FileNotFoundError(f'Directory {dir_path} does not exist')



class Helper():
    def __init__(self, core):
        # Core
        self.core = core

        self.input_files_messages = { 'geno': '', 'ind': '', 'snp': '', 'pops': '' }

    def print_input_files_progress(self, key, message):
        self.input_files_messages[key] = message

    def print_computation_progress(self, *args):
        if len(args) == 1 and isinstance(args[0], int):
            print(f'{100 * args[0] / len(self.core.selected_pops):.2f}%')
        elif len(args) == 3 and isinstance(args[0], str) and isinstance(args[1], str) and isinstance(args[2], int):
            print(args[1])

    def run(self):
        print(f'Mixtum v{self.core.version}\n')
        print('Parsing and checking input files...')

        self.check_input_files()
        self.parse_pops_file()

        for key, message in self.input_files_messages.items():
            print(f'{key} file: {message}')

        print('Parsing and checking finished.\n')

        self.compute_frequencies()

        self.compute_results()

    def check_input_files(self):
        self.core.geno_table_shape(self.print_input_files_progress)
        self.core.parse_ind_file(self.print_input_files_progress)
        self.core.parse_snp_file(self.print_input_files_progress)

        valid = True
        if not self.core.check_geno_file():
            print('Error in .geno file: not all rows have the same number of columns.')
            valid = False
        if not self.core.check_ind_and_geno():
            print('Error: Number of populations in .ind file is not equal to number of columns in .geno file.')
            valid = False
        if not self.core.check_snp_and_geno():
            print('Error: Number of alleles in .snp file is not equal to number of rows in .geno file.')
            valid = False

        if not valid:
            sys.exit(1)

        print('Parsed input files seem to have a valid structure.')

    def parse_pops_file(self):
        self.core.parse_selected_populations(self.print_input_files_progress)
        missing_pops = self.core.check_parsed_pops()
        if len(missing_pops) > 0:
            print(f'Error: The following populations are missing from .ind file and were deselected: {','.join(missing_pops)}')

    def compute_frequencies(self):
        self.core.parallel_compute_populations_frequencies(self.print_computation_progress)

    def compute_results(self):
        print('\nComputing...')

        self.core.init_admixture_model()

        n = 1
        ntotal = 9

        print(f"({n}/{ntotal}) Mixing coefficient pre-jl")
        self.core.mixing_coefficient_pre_jl()
        n += 1

        print(f"({n}/{ntotal}) Admixture angle pre-jl")
        self.core.admixture_angle_pre_jl()
        n += 1

        print(f"({n}/{ntotal}) f3")
        self.core.f3()
        n += 1

        print(f"({n}/{ntotal}) f4'")
        self.core.f4_prime()
        n += 1

        print(f"({n}/{ntotal}) Alpha'")
        self.core.alpha_prime()
        n += 1

        print(f"({n}/{ntotal}) f4-standard")
        self.core.f4_std()
        n += 1

        print(f"({n}/{ntotal}) Alpha-standard")
        self.core.alpha_standard()
        n += 1

        print(f"({n}/{ntotal}) Admixture angle post-jl")
        self.core.admixture_angle_post_jl()
        n += 1

        print(f"({n}/{ntotal}) f4-ratio")
        self.core.f4_ratio()

        print('\nResults:')
        print(self.core.admixture_data())

    def plot(self):
        self.core.plot()


if __name__ == '__main__':
    core = Core()

    parser = argparse.ArgumentParser(description = f'Mixtum v{core.version}: The geometry of admixture in population genetics')
    parser.add_argument('--geno', type = str, required = True, help = 'Path of .geno file')
    parser.add_argument('--ind', type = str, required = True, help = 'Path of .ind file')
    parser.add_argument('--snp', type = str, required = True, help = 'Path of .snp file')
    parser.add_argument('--pops', type = str, required = True, help = 'Path of selected populations file')
    parser.add_argument('--outdir', type = str, required = True, help = 'Path of output dir')
    parser.add_argument('--nprocs', type = int, default = 1, help = 'Number of parallel computation processes (default 1)')
    parser.add_argument('--plot', action=argparse.BooleanOptionalAction, help='Plot fits and histogram')

    args = parser.parse_args()

    num_procs = args.nprocs

    geno_file_path = Path(args.geno)
    ind_file_path = Path(args.ind)
    snp_file_path = Path(args.snp)
    pops_file_path = Path(args.pops)
    out_dir_path = Path(args.outdir)

    check_file_path(geno_file_path)
    check_file_path(ind_file_path)
    check_file_path(snp_file_path)
    check_file_path(pops_file_path)
    check_dir_path(out_dir_path)

    core.set_num_procs(num_procs)

    core.set_geno_file_path(args.geno)
    core.set_ind_file_path(args.ind)
    core.set_snp_file_path(args.snp)
    core.set_pops_file_path(args.pops)

    helper = Helper(core)

    helper.run()

    if args.plot:
        helper.plot()