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

from gui.log_system import LogSystem
from gui.select_file_widget import SelectFileWidget
from gui.worker import Worker
from gui.about_dialog import AboutDialog
from gui.input_files_checker import InputFilesChecker

from PySide6.QtCore import Qt, Signal, Slot, QThreadPool, QSize
from PySide6.QtWidgets import QWidget, QPushButton, QSizePolicy, QGroupBox, QVBoxLayout, QHBoxLayout



class InputFilesWidget(QWidget):
    snp_file_parsed = Signal()
    ind_file_parsed = Signal()
    parsed_pops_changed = Signal()
    min_snp_cutoff_check_failed = Signal()

    def __init__(self, core):
        QWidget.__init__(self)

        # Core
        self.core = core

        # File checker
        self.input_files_checker = InputFilesChecker(core)
        self.input_files_checker.file_path_set.connect(self.log_file_path)

        # Thread pool
        self.thread_pool = QThreadPool()
        self.worker_finished = {'geno': False, 'ind': False, 'snp': False}

        # Log system
        self.log = LogSystem(['main', 'geno', 'ind', 'snp', 'pops', 'check'])
        self.log.set_entry('main', 'Select input file triad and optionally a selected populations file.')
        self.log.append_entry('geno', '')
        self.log.append_entry('ind', '')
        self.log.append_entry('snp', '')
        self.log.append_entry('pops', '')

        # Stylesheets
        stylesheet_11 = 'color: white; background-color: rgb(128, 45, 0); font-size: 24pt;'
        stylesheet_12 = 'color: white; background-color: rgb(128, 64, 0); font-size: 24pt;'
        stylesheet_21 = 'color: white; background-color: rgb(0, 45, 128); font-size: 24pt;'
        stylesheet_22 = 'color: white; background-color: rgb(0, 64, 128); font-size: 24pt;'
        stylesheet_31 = 'color: white; background-color: rgb(32, 128, 32); font-size: 20pt;'

        # Select file widgets
        self.geno_file_widget = SelectFileWidget('Select .geno file', '(*.geno)', stylesheet_11)
        self.ind_file_widget = SelectFileWidget('Select .ind file', '(*.ind)', stylesheet_11)
        self.snp_file_widget = SelectFileWidget('Select .snp file', '(*.snp)', stylesheet_11)
        self.pops_file_widget = SelectFileWidget('Select populations file', None, stylesheet_21)

        self.geno_file_widget.file_path_selected.connect(self.input_files_checker.set_geno_file_path)
        self.ind_file_widget.file_path_selected.connect(self.input_files_checker.set_ind_file_path)
        self.snp_file_widget.file_path_selected.connect(self.input_files_checker.set_snp_file_path)
        self.pops_file_widget.file_path_selected.connect(self.input_files_checker.set_pops_file_path)

        # Check files button
        self.check_button = QPushButton('Parse and check files')
        self.check_button.setMinimumSize(QSize(500, 100))
        self.check_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.check_button.setStyleSheet(stylesheet_12)
        self.check_button.setEnabled(False)
        self.input_files_checker.input_file_paths_state.connect(self.check_button.setEnabled)
        self.check_button.clicked.connect(self.check_input_files)

        # Check error
        self.input_files_checker.geno_file_error.connect(self.geno_check_failed)
        self.input_files_checker.ind_file_error.connect(self.ind_check_failed)
        self.input_files_checker.snp_file_error.connect(self.snp_check_failed)

        # Parse selected pops file button
        self.parse_pops_button = QPushButton('Load selected populations')
        self.parse_pops_button.setMinimumSize(QSize(500, 100))
        self.parse_pops_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.parse_pops_button.setStyleSheet(stylesheet_22)
        self.parse_pops_button.setEnabled(False)
        self.input_files_checker.pops_file_path_state.connect(self.parse_pops_button.setEnabled)
        self.parse_pops_button.clicked.connect(self.parse_pops_file)

        # Parse error
        self.input_files_checker.parsed_pops_error.connect(self.pops_check_failed)

        # Required files buttons group box
        req_group_box = QGroupBox('Required files')
        req_group_box.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        req_layout = QVBoxLayout()
        req_layout.addWidget(self.geno_file_widget, 0, Qt.AlignmentFlag.AlignCenter)
        req_layout.addWidget(self.ind_file_widget, 0, Qt.AlignmentFlag.AlignCenter)
        req_layout.addWidget(self.snp_file_widget, 0, Qt.AlignmentFlag.AlignCenter)
        req_layout.addWidget(self.check_button, 0, Qt.AlignmentFlag.AlignCenter)
        req_group_box.setLayout(req_layout)

        # Optional files buttons group box
        opt_group_box = QGroupBox('Optional files')
        opt_group_box.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        opt_layout = QVBoxLayout()
        opt_layout.addWidget(self.pops_file_widget, 0, Qt.AlignmentFlag.AlignCenter)
        opt_layout.addWidget(self.parse_pops_button, 0, Qt.AlignmentFlag.AlignCenter)
        opt_group_box.setLayout(opt_layout)

        # About dialog button
        self.about_button = QPushButton('About Mixtum')
        self.about_button.setMinimumSize(QSize(250, 50))
        self.about_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.about_button.setStyleSheet(stylesheet_31)
        self.about_button.clicked.connect(self.show_about_dialog)

        self.about_dialog = AboutDialog(core)
        self.about_dialog.setModal(False)
        self.about_dialog.finished.connect(self.enable_about_button)

        # Child layout
        child_layout = QVBoxLayout()
        child_layout.addWidget(opt_group_box)
        child_layout.addWidget(self.about_button)

        # Layout
        layout = QHBoxLayout(self)
        layout.addWidget(req_group_box)
        layout.addLayout(child_layout)

    @Slot(str, str)
    def log_file_path(self, key, file_path):
        self.log.set_entry(key, file_path, 0)

    @Slot(str, str)
    def log_progress(self, key, message):
        self.log.set_entry(key, message, 1)

    @Slot(str)
    def checking_finished(self, worker_name):
        self.worker_finished[worker_name] = True

        if all([status for name, status in self.worker_finished.items()]):
            if self.input_files_checker.check_input_files():
                self.log.set_entry('main', 'Checking finished.')
                self.log.append_entry('check', 'Parsed input files seem to have a valid structure.')

                self.core.check_parsed_pops()

        if worker_name == 'ind':
            self.ind_file_parsed.emit()
        elif worker_name == 'snp':
            if self.core.check_min_snp_cutoff():
                self.snp_file_parsed.emit()
            else:
                self.min_snp_cutoff_check_failed.emit()

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

    @Slot()
    def check_input_files(self):
        self.log.clear_entry('check')
        self.log.set_entry('main', 'Checking input files...')

        self.worker_finished = {'geno': False, 'ind': False, 'snp': False}

        geno_is_ascii = self.core.is_geno_file_ascii()

        if geno_is_ascii:
            geno_worker = Worker('geno', self.core.geno_table_shape)
            geno_worker.signals.progress[str, str].connect(self.log_progress)
            geno_worker.signals.finished.connect(self.checking_finished)
        elif self.core.read_geno_file_header(self.log_progress):
            self.checking_finished('geno')
        else:
            self.log.set_entry('main', 'Checking error!')
            self.log.append_entry('check', 'Error: unsupported .geno file format.')
            return

        ind_worker = Worker('ind', self.core.parse_ind_file)
        ind_worker.signals.progress[str, str].connect(self.log_progress)
        ind_worker.signals.finished.connect(self.checking_finished)

        snp_worker = Worker('snp', self.core.parse_snp_file)
        snp_worker.signals.progress[str, str].connect(self.log_progress)
        snp_worker.signals.finished.connect(self.checking_finished)

        if geno_is_ascii:
            self.thread_pool.start(geno_worker)
        self.thread_pool.start(ind_worker)
        self.thread_pool.start(snp_worker)

    @Slot(str)
    def parsing_finished(self, worker_name):
        self.core.check_parsed_pops()
        self.log.set_entry('main', 'Parsing finished.')
        self.parsed_pops_changed.emit()

    @Slot(list)
    def pops_check_failed(self, missing_pops):
        self.log.set_entry('pops', f'Error: The following populations are missing from .ind file and were deselected: {','.join(missing_pops)}')
        self.parsed_pops_changed.emit()

    @Slot()
    def parse_pops_file(self):
        self.log.set_entry('main', 'Parsing selected populations file...')

        pops_worker = Worker('pops', self.core.parse_selected_populations)
        pops_worker.signals.progress[int].connect(self.log_progress)
        pops_worker.signals.finished.connect(self.parsing_finished)

        self.thread_pool.start(pops_worker)

    @Slot()
    def show_about_dialog(self):
        self.about_button.setEnabled(False)
        self.about_dialog.show()

    @Slot()
    def enable_about_button(self):
        self.about_button.setEnabled(True)