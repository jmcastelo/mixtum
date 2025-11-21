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

from pathlib import Path

from PySide6.QtCore import QObject, Signal, Slot

class InputFilesChecker(QObject):
    # Signals
    file_path_set = Signal(str, str)

    input_file_paths_state = Signal(bool)
    pops_file_path_state = Signal(bool)

    geno_file_error = Signal()
    ind_file_error = Signal(int, int)
    snp_file_error = Signal(int, int)

    parsed_pops_error = Signal(list)

    def __init__(self, core):
        QObject.__init__(self)

        self.core = core

    def check_file_paths(self):
        self.input_file_paths_state.emit(bool(self.core.geno_file_path.is_file() and self.core.ind_file_path.is_file() and self.core.snp_file_path.is_file()))

    @Slot(str)
    def set_geno_file_path(self, file_path):
        self.core.set_geno_file_path(file_path)
        self.file_path_set.emit('geno', file_path)
        self.check_file_paths()

    @Slot(str)
    def set_ind_file_path(self, file_path):
        self.core.set_ind_file_path(file_path)
        self.file_path_set.emit('ind', file_path)
        self.check_file_paths()

    @Slot(str)
    def set_snp_file_path(self, file_path):
        self.core.set_snp_file_path(file_path)
        self.file_path_set.emit('snp', file_path)
        self.check_file_paths()

    @Slot(str)
    def set_pops_file_path(self, file_path):
        self.core.set_pops_file_path(file_path)
        self.file_path_set.emit('pops', file_path)
        self.pops_file_path_state.emit(self.core.pops_file_path.is_file())

    # Check input file consistency
    def check_input_files(self):
        valid = True

        if self.core.is_geno_file_ascii():
            if not self.core.check_geno_file():
                self.geno_file_error.emit()
                valid = False
            if not self.core.check_ind_and_geno():
                self.ind_file_error.emit(self.core.num_ind_rows, self.core.num_geno_cols[0])
                valid = False
            if not self.core.check_snp_and_geno():
                self.snp_file_error.emit(self.core.num_snp_rows, self.core.num_geno_rows)
                valid = False
        else:
            if not self.core.check_ind_and_geno_packed():
                self.ind_file_error.emit(self.core.num_ind_rows, self.core.num_ind)
                valid = False
            if not self.core.check_snp_and_geno_packed():
                self.snp_file_error.emit(self.core.num_snp_rows, self.core.num_snp)
                valid = False

        return valid

    def check_parsed_pops(self):
        missing_pops = self.core.check_parsed_pops()
        if len(missing_pops) > 0:
            self.parsed_pops_error.emit(missing_pops)