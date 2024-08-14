from pathlib import Path

from gui.select_file_widget import SelectFileWidget
from gui.worker import Worker

from PySide6.QtCore import Signal, Slot, QThreadPool
from PySide6.QtWidgets import QWidget, QPushButton, QLabel, QVBoxLayout



class InputFilesWidget(QWidget):
    log_text_changed = Signal(str)

    def __init__(self, core):
        QWidget.__init__(self)

        # Core
        self.core = core

        # Thread pool
        self.thread_pool = QThreadPool()
        self.worker_finished = {'geno': False, 'ind': False, 'snp': False}

        # Log text
        self.log_text = ''
        self.log_text_lines = {'main' : 'Please, select input files.', 'geno': '', 'ind': '', 'snp': '', 'pops': '', 'check': []}
        self.set_log_text()

        # Select file widgets
        self.geno_file_widget = SelectFileWidget('Select .geno file', '(*.geno)')
        self.ind_file_widget = SelectFileWidget('Select .ind file', '(*.ind)')
        self.snp_file_widget = SelectFileWidget('Select .snp file', '(*.snp)')
        self.pops_file_widget = SelectFileWidget('Select populations file', None)

        self.geno_file_widget.file_path_selected.connect(self.core.set_geno_file_path)
        self.ind_file_widget.file_path_selected.connect(self.core.set_ind_file_path)
        self.snp_file_widget.file_path_selected.connect(self.core.set_snp_file_path)
        self.pops_file_widget.file_path_selected.connect(self.core.set_pops_file_path)

        # Check files button
        self.check_button = QPushButton('Parse and check')
        self.check_button.setEnabled(False)
        self.core.input_file_paths_state.connect(self.check_button.setEnabled)
        self.check_button.clicked.connect(self.check_input_files)

        # Check error
        self.core.geno_file_error.connect(self.geno_check_failed)
        self.core.ind_file_error.connect(self.ind_check_failed)
        self.core.snp_file_error.connect(self.snp_check_failed)

        # Parse selected pops file button
        self.parse_pops_button = QPushButton('Load selected populations')
        self.parse_pops_button.setEnabled(False)
        self.core.pops_file_path_state.connect(self.parse_pops_button.setEnabled)
        self.parse_pops_button.clicked.connect(self.parse_pops_file)

        # Layout
        layout = QVBoxLayout(self)
        layout.addWidget(self.geno_file_widget)
        layout.addWidget(self.ind_file_widget)
        layout.addWidget(self.snp_file_widget)
        layout.addWidget(self.pops_file_widget)
        layout.addWidget(self.check_button)
        layout.addWidget(self.parse_pops_button)

    def set_log_text(self):
        text_list = [text for key, text in self.log_text_lines.items() if key != 'check' and len(text) > 0] + self.log_text_lines['check']
        self.log_text = '\n'.join(text_list)
        self.log_text_changed.emit(self.log_text)

    @Slot(int)
    def set_geno_log_text(self, num_rows):
        self.log_text_lines['geno'] = f'Number of rows in .geno file: {num_rows}'
        self.set_log_text()

    @Slot(int)
    def set_ind_log_text(self, num_rows):
        self.log_text_lines['ind'] = f'Number of rows in .ind file: {num_rows}'
        self.set_log_text()

    @Slot(int)
    def set_snp_log_text(self, num_rows):
        self.log_text_lines['snp'] = f'Number of rows in .snp file: {num_rows}'
        self.set_log_text()

    @Slot(int)
    def set_pops_log_text(self, num_rows):
        self.log_text_lines['pops'] = f'Number of selected populations: {num_rows}'
        self.set_log_text()

    @Slot(str)
    def checking_finished(self, worker_name):
        self.worker_finished[worker_name] = True
        if all([status for name, status in self.worker_finished.items()]):
            if self.core.check_input_files():
                self.log_text_lines['main'] = 'Check finished.'
                self.log_text_lines['check'].append('Parsed input files seem to have a valid structure.')
                self.set_log_text()

    @Slot()
    def geno_check_failed(self):
        self.log_text_lines['main'] = 'Checking error!'
        self.log_text_lines['check'].append('Error in .geno file: not all rows have the same number of columns.')
        self.set_log_text()

    @Slot(int, int)
    def ind_check_failed(self, num_pops, num_cols):
        self.log_text_lines['main'] = 'Checking error!'
        self.log_text_lines['check'].append(f'Error: Number of populations ({num_pops}) in .ind file is not equal to number of columns ({num_cols}) in .geno file.')
        self.set_log_text()

    @Slot(int, int)
    def snp_check_failed(self, num_alleles, num_rows):
        self.log_text_lines['main'] = 'Checking error!'
        self.log_text_lines['check'].append(f'Error: Number of alleles ({num_alleles}) in .snp file is not equal to number of rows ({num_rows}) in .geno file.')
        self.set_log_text()

    @Slot()
    def check_input_files(self):
        self.log_text_lines['main'] = 'Checking input files...'
        self.log_text_lines['check'] = []
        self.set_log_text()

        self.worker_finished = {'geno': False, 'ind': False, 'snp': False}

        geno_worker = Worker('geno', self.core.geno_table_shape)
        geno_worker.signals.progress[int].connect(self.set_geno_log_text)
        geno_worker.signals.finished.connect(self.checking_finished)

        ind_worker = Worker('ind', self.core.parse_ind_file)
        ind_worker.signals.progress[int].connect(self.set_ind_log_text)
        ind_worker.signals.finished.connect(self.checking_finished)

        snp_worker = Worker('snp', self.core.parse_snp_file)
        snp_worker.signals.progress[int].connect(self.set_snp_log_text)
        snp_worker.signals.finished.connect(self.checking_finished)

        self.thread_pool.start(geno_worker)
        self.thread_pool.start(ind_worker)
        self.thread_pool.start(snp_worker)

    @Slot(str)
    def parsing_finished(self, worker_name):
        self.log_text_lines['main'] = 'Parsing finished.'
        self.set_log_text()

    @Slot()
    def parse_pops_file(self):
        self.log_text_lines['main'] = 'Parsing selected populations file...'
        self.set_log_text()

        pops_worker = Worker('pops', self.core.parse_selected_populations)
        pops_worker.signals.progress[int].connect(self.set_pops_log_text)
        pops_worker.signals.finished.connect(self.parsing_finished)

        self.thread_pool.start(pops_worker)