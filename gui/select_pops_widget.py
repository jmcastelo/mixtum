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
from gui.searchable_table_widget import SearchableTableWidget
from gui.worker import Worker

from pathlib import Path

from PySide6.QtCore import Qt, Signal, Slot, QThreadPool
from PySide6.QtWidgets import QWidget, QTableWidget, QTableWidgetItem, QPushButton, QSizePolicy, QFrame, QSpinBox
from PySide6.QtWidgets import QProgressBar, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QHeaderView, QFileDialog



class SelectPopsWidget(QWidget):
    selected_pops_changed = Signal()
    computation_result = Signal(bool)

    def __init__(self, core):
        QWidget.__init__(self)

        # Core
        self.core = core

        # Thread pool
        self.thread_pool = QThreadPool()

        # Log
        self.log = LogSystem(['main', 'progress', 'timing', 'check'])
        self.log.append_entry('main', 'Select entries from the available populations on the left table and compute their allele frequencies.')
        self.log.append_entry('main', 'A minimum of 7 populations is required (3 for admixture model + 4 auxiliary populations).')
        self.log.append_entry('main', 'To perform bootstrap, a minimum of 11 populations is required (of which 8 are auxiliaries).')
        self.log.append_entry('timing', '')
        self.log.append_entry('check', '')

        # Searchable table containing available populations
        self.search_widget = SearchableTableWidget()

        # Table containing selected populations
        self.selected_table = QTableWidget()
        self.selected_table.setSortingEnabled(True)
        self.selected_table.setColumnCount(1)
        self.selected_table.verticalHeader().setVisible(False)
        self.selected_table.horizontalHeader().setVisible(True)
        self.selected_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.selected_table.setHorizontalHeaderLabels(['Selected populations'])

        # Select populations button
        self.select_button = QPushButton('Select populations')
        self.select_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.select_button.clicked.connect(self.select_populations)

        # Remove populations button
        self.remove_button = QPushButton('Deselect populations')
        self.remove_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.remove_button.clicked.connect(self.remove_populations)

        # Reset populations button
        self.reset_button = QPushButton('Reset populations')
        self.reset_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.reset_button.clicked.connect(self.reset_populations)

        # Number of computing processes
        self.num_procs_spin_box = QSpinBox()
        self.num_procs_spin_box.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.num_procs_spin_box.setMinimum(1)
        self.num_procs_spin_box.setMaximum(9999)
        self.num_procs_spin_box.valueChanged.connect(self.set_num_procs)
        self.num_procs_spin_box.setValue(1)

        # Number of computing processes
        self.snp_cutoff_spin_box = QSpinBox()
        self.snp_cutoff_spin_box.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.snp_cutoff_spin_box.setMinimum(self.core.min_snp_cutoff)
        self.snp_cutoff_spin_box.setMaximum(9999999)
        self.snp_cutoff_spin_box.valueChanged.connect(self.set_snp_cutoff)
        self.snp_cutoff_spin_box.setValue(self.core.min_snp_cutoff)
        self.snp_cutoff_spin_box.setEnabled(False)

        # Compute allele frequencies files button
        self.comp_button = QPushButton('Compute frequencies')
        self.comp_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.comp_button.setEnabled(False)
        self.comp_button.clicked.connect(self.compute_frequencies)

        # Stop computation button
        self.stop_button = QPushButton('Stop computation')
        self.stop_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self.stop_computation)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(1)
        self.progress_bar.setValue(0)

        # Save frequencies button
        self.save_freqs_button = QPushButton('Save frequencies')
        self.save_freqs_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.save_freqs_button.setEnabled(False)
        self.save_freqs_button.clicked.connect(self.save_frequencies)

        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setFrameShadow(QFrame.Shadow.Sunken)

        # Select table buttons layout
        hlayout = QHBoxLayout()
        hlayout.addWidget(self.remove_button, 0, Qt.AlignmentFlag.AlignJustify)
        hlayout.addWidget(self.reset_button, 0, Qt.AlignmentFlag.AlignJustify)

        # Select table layout
        vlayout = QVBoxLayout()
        vlayout.addWidget(self.selected_table)
        vlayout.addLayout(hlayout)

        # Select group box
        select_group = QGroupBox('Selected populations table')
        select_group.setLayout(vlayout)

        # Search table layout
        slayout = QVBoxLayout()
        slayout.addWidget(self.search_widget)
        slayout.addWidget(self.select_button, 0, Qt.AlignmentFlag.AlignHCenter)

        # Search group box
        search_group = QGroupBox('Search and select from available populations')
        search_group.setLayout(slayout)

        # Tables layout
        tlayout = QHBoxLayout()
        tlayout.addWidget(search_group)
        tlayout.addWidget(select_group)

        # Numprocs form layout
        npflayout = QFormLayout()
        npflayout.addRow('Number of processes:', self.num_procs_spin_box)

        # Cutoff form layout
        coflayout = QFormLayout()
        coflayout.addRow('snp cutoff:', self.snp_cutoff_spin_box)

        # Computation layout
        clayout = QHBoxLayout()
        clayout.addLayout(npflayout)
        clayout.addLayout(coflayout)
        clayout.addWidget(self.comp_button)
        clayout.addWidget(self.stop_button)
        clayout.addWidget(self.progress_bar)
        clayout.addWidget(self.save_freqs_button)

        # Layout
        layout = QVBoxLayout(self)
        layout.addLayout(tlayout)
        layout.addWidget(separator)
        layout.addLayout(clayout)

        self.selected_pops_changed.connect(self.reset_controls)

    @Slot()
    def reset_controls(self):
        self.log.clear_entry('main')
        self.log.append_entry('main', 'Select entries from the available populations on the left table and compute their allele frequencies.')
        self.log.append_entry('main', 'A minimum of 7 populations is required (3 for admixture model + 4 auxiliary populations). ')
        self.log.append_entry('main', 'To perform bootstrap, a minimum of 11 populations is required (of which 8 are auxiliaries). ')
        self.log.clear_entry('progress')
        self.log.clear_entry('timing')
        self.log.clear_entry('check')

        self.progress_bar.setValue(0)

    def disable_controls(self):
        self.comp_button.setEnabled(False)
        self.stop_button.setEnabled(False)
        self.save_freqs_button.setEnabled(False)
        self.snp_cutoff_spin_box.setEnabled(False)

    @Slot()
    def init_search_table(self):
        self.search_widget.init_table(self.core.avail_pops)

    @Slot()
    def init_selected_table(self):
        self.populate_selected_table(self.core.parsed_pops)

    @Slot()
    def select_populations(self):
        self.core.append_pops(self.search_widget.selected_items())
        self.populate_selected_table(self.core.selected_pops)
        self.selected_table.scrollToBottom()

    @Slot()
    def remove_populations(self):
        self.core.remove_pops([item.text() for item in self.selected_table.selectedItems()])
        self.populate_selected_table(self.core.selected_pops)

    @Slot()
    def reset_populations(self):
        self.core.reset_pops()
        self.populate_selected_table(self.core.selected_pops)

    @Slot(object)
    def populate_selected_table(self, pops):
        self.selected_table.clearContents()
        self.selected_table.setRowCount(len(pops))

        for index, item in enumerate(pops):
            table_widget_item = QTableWidgetItem(item)
            table_widget_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            self.selected_table.setItem(index, 0, table_widget_item)

        self.comp_button.setEnabled(len(pops) >= 7)
        self.selected_pops_changed.emit()

    @Slot(int)
    def set_num_procs(self, procs):
        self.core.set_num_procs(procs)

    @Slot()
    def set_snp_cutoff_spin_box(self):
        self.snp_cutoff_spin_box.setMaximum(self.core.num_snp)
        self.snp_cutoff_spin_box.setValue(self.core.num_snp)
        self.snp_cutoff_spin_box.setEnabled(True)

    @Slot(int)
    def set_snp_cutoff(self, n):
        self.core.set_snp_cutoff(n)

    @Slot(str, str, int)
    def log_progress(self, key, message, line_num):
        self.log.set_entry(key, message, line_num)

    @Slot(str)
    def computation_finished(self, worker_name):
        if worker_name == 'freqs':
            self.stop_button.setEnabled(False)
            self.save_freqs_button.setEnabled(True)

    @Slot()
    def compute_frequencies(self):
        self.log.clear_entry('main')
        self.log.clear_entry('progress')
        self.log.clear_entry('timing')
        self.log.clear_entry('check')

        self.progress_bar.setMaximum(len(self.core.selected_pops))

        self.stop_button.setEnabled(True)

        worker = Worker('freqs', self.core.parallel_compute_populations_frequencies)
        worker.signals.progress[str, str, int].connect(self.log_progress)
        worker.signals.progress[int].connect(self.progress_bar.setValue)
        worker.signals.result.connect(self.computation_result)
        worker.signals.finished.connect(self.computation_finished)

        self.thread_pool.start(worker)

    @Slot()
    def stop_computation(self):
        self.core.stop_computation()

    @Slot()
    def save_frequencies(self):
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.FileMode.AnyFile)
        dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)

        if dialog.exec():
            file_names = dialog.selectedFiles()
            file_path = file_names[0]
            self.core.save_population_allele_frequencies(Path(file_path))
