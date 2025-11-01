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

from PySide6.QtCore import QCoreApplication, QObject, Signal, Slot, QThreadPool

from gui.core import Core
from gui.worker import Worker
from gui.log_system import LogSystem



# Raise error if file not found
def check_file_path(file_path):
    if not file_path.is_file():
        raise FileNotFoundError(f'File {file_path} does not exist')



# Raise error if dir not found
def check_dir_path(dir_path):
    if not dir_path.is_dir():
        raise FileNotFoundError(f'Directory {dir_path} does not exist')



class Helper(QObject):
    def __init__(self, core, app):
        QObject.__init__(self)

        # Core
        self.core = core

        # App
        self.app = app

        # Check error
        self.core.geno_file_error.connect(self.geno_check_failed)
        self.core.ind_file_error.connect(self.ind_check_failed)
        self.core.snp_file_error.connect(self.snp_check_failed)
        self.core.parsed_pops_error.connect(self.pops_check_failed)

        # Thread pool
        self.thread_pool = QThreadPool()
        self.worker_finished = {'geno': False, 'ind': False, 'snp': False}

        # Log system
        self.log_text = ''
        self.log = LogSystem(['main', 'geno', 'ind', 'snp', 'pops', 'check', 'progress', 'timing', 'percent'])
        self.log.set_entry('main', 'Checking input files...')
        self.log.append_entry('geno', '')
        self.log.append_entry('ind', '')
        self.log.append_entry('snp', '')
        self.log.append_entry('pops', '')
        self.log.append_entry('timing', '')
        self.log.append_entry('percent', '')
        self.log.append_entry('check', '')

        self.log.changed.connect(self.set_log_text)

        self.core.file_path_set.connect(self.log_file_path)

    def print_log_text(self):
        print(self.log_text)

    @Slot(str)
    def set_log_text(self, text):
        self.log_text = text

    @Slot(str, str)
    def log_file_path(self, key, file_path):
        self.log.set_entry(key, file_path, 0)

    @Slot(str, str)
    def log_progress(self, key, message):
        self.log.set_entry(key, message, 0)

    @Slot(str)
    def checking_finished(self, worker_name):
        self.worker_finished[worker_name] = True

        if all([status for name, status in self.worker_finished.items()]):
            if self.core.check_input_files():
                self.log.append_entry('check', 'Parsed input files seem to have a valid structure.')
                self.parse_pops_file()
            else:
                self.print_log_text()
                self.app.quit()

    @Slot()
    def geno_check_failed(self):
        self.log.set_entry('main', 'Checking error!')
        self.log.append_entry('check', 'Error in .geno file: not all rows have the same number of columns.')

    @Slot(int, int)
    def ind_check_failed(self, num_pops, num_cols):
        self.log.set_entry('main', 'Checking error!')
        self.log.append_entry('check', f'Error: Number of populations ({num_pops}) in .ind file is not equal to number of columns ({num_cols}) in .geno file.')

    @Slot(int, int)
    def snp_check_failed(self, num_alleles, num_rows):
        self.log.set_entry('main', 'Checking error!')
        self.log.append_entry('check', f'Error: Number of alleles ({num_alleles}) in .snp file is not equal to number of rows ({num_rows}) in .geno file.')

    def check_input_files(self):
        self.print_log_text()

        self.worker_finished = {'geno': False, 'ind': False, 'snp': False}

        geno_worker = Worker('geno', self.core.geno_table_shape)
        geno_worker.signals.progress[str, str].connect(self.log_progress)
        geno_worker.signals.finished.connect(self.checking_finished)

        ind_worker = Worker('ind', self.core.parse_ind_file)
        ind_worker.signals.progress[str, str].connect(self.log_progress)
        ind_worker.signals.finished.connect(self.checking_finished)

        snp_worker = Worker('snp', self.core.parse_snp_file)
        snp_worker.signals.progress[str, str].connect(self.log_progress)
        snp_worker.signals.finished.connect(self.checking_finished)

        self.thread_pool.start(geno_worker)
        self.thread_pool.start(ind_worker)
        self.thread_pool.start(snp_worker)

    @Slot(str)
    def parsing_finished(self, worker_name):
        self.core.check_parsed_pops()

        self.log.set_entry('main', 'Parsing and checking finished.')
        self.print_log_text()

        self.compute_frequencies()

    @Slot(list)
    def pops_check_failed(self, missing_pops):
        self.log.set_entry('pops', f'Error: The following populations are missing from .ind file and were unselected: {','.join(missing_pops)}')

    def parse_pops_file(self):
        pops_worker = Worker('pops', self.core.parse_selected_populations)
        pops_worker.signals.progress[str, str].connect(self.log_progress)
        pops_worker.signals.finished.connect(self.parsing_finished)

        self.thread_pool.start(pops_worker)

    @Slot(str, str, int)
    def log_computation_progress(self, key, message, line_num):
        # self.log.set_text_from_entry(key, message, line_num)
        # self.print_log_text()
        if len(message) > 0:
            print(message, flush=True)

    @Slot(int)
    def computation_progress(self, index):
        if index > 0:
            percent = f'{100 * index / len(self.core.selected_pops):.2f}%'
            # self.log.set_text_from_entry('percent', percent, 0)
            # self.print_log_text()
            print(percent, flush=True)

    @Slot(bool)
    def computation_result(self, status):
        if not status:
            self.app.quit()

    @Slot(str)
    def computation_finished(self, worker_name):
        self.compute_results()

    def compute_frequencies(self):
        self.log.clear_entry('main')
        self.log.clear_entry('geno')
        self.log.clear_entry('ind')
        self.log.clear_entry('snp')
        self.log.clear_entry('pops')
        self.log.clear_entry('check')

        print('')

        worker = Worker('freqs', self.core.parallel_compute_populations_frequencies)
        worker.signals.progress[int].connect(self.computation_progress)
        worker.signals.progress[str, str, int].connect(self.log_computation_progress)
        worker.signals.result.connect(self.computation_result)
        worker.signals.finished.connect(self.computation_finished)

        self.thread_pool.start(worker)

    def compute_results(self):
        print('\nComputing:')

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

        self.app.quit()



def main():
    parser = argparse.ArgumentParser(description = 'Mixtum: The geometry of admixture in population genetics')
    parser.add_argument('--geno', type = str, required = True, help = 'Path of .geno file')
    parser.add_argument('--ind', type = str, required = True, help = 'Path of .ind file')
    parser.add_argument('--snp', type = str, required = True, help = 'Path of .snp file')
    parser.add_argument('--pops', type = str, required = True, help = 'Path of selected populations file')
    parser.add_argument('--outdir', type = str, required = True, help = 'Path of output dir')
    parser.add_argument('--nprocs', type = int, default = 1, help = 'Number of parallel computation processes')

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

    app = QCoreApplication()

    core = Core()

    core.set_num_procs(num_procs)

    helper = Helper(core, app)

    core.set_geno_file_path(args.geno)
    core.set_ind_file_path(args.ind)
    core.set_snp_file_path(args.snp)
    core.set_pops_file_path(args.pops)

    helper.check_input_files()

    sys.exit(app.exec())

if __name__ == '__main__':
    main()