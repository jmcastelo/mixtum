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

from gui.core import Core
from gui.log_widget import LogWidget
from gui.input_files_widget import InputFilesWidget
from gui.select_pops_widget import SelectPopsWidget
from gui.mix_model_widget import MixModelWidget
from gui.pca_widget import PCAWidget
from gui.f_statistics_widget import FStatisticsWidget

from PySide6.QtCore import Qt, Slot
from PySide6.QtWidgets import QWidget, QTabWidget, QSplitter, QVBoxLayout



class MainWindow(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle('Mixtum')

        # Core
        self.core = Core()

        # Log widget
        self.log_widget = LogWidget()

        # Input files widget
        self.input_files_widget = InputFilesWidget(self.core)

        # Select pops widget
        self.sel_pops_widget = SelectPopsWidget(self.core)

        # Admixture model widget
        self.mix_model_widget = MixModelWidget(self.core)

        # PCA widget
        self.pca_widget = PCAWidget(self.core)

        # f-statistics widget
        self.f_statistics_widget = FStatisticsWidget(self.core)

        # Connections
        self.input_files_widget.ind_file_parsed.connect(self.sel_pops_widget.init_search_table)
        self.input_files_widget.snp_file_parsed.connect(self.sel_pops_widget.set_snp_cutoff_spin_box)
        self.input_files_widget.min_snp_cutoff_check_failed.connect(self.sel_pops_widget.disable_controls)
        self.input_files_widget.parsed_pops_changed.connect(self.sel_pops_widget.init_selected_table)
        self.sel_pops_widget.computation_result.connect(self.mix_model_widget.init_pop_tables)
        self.sel_pops_widget.computation_result.connect(self.pca_widget.init_sel_pops_table)
        self.sel_pops_widget.computation_result.connect(self.f_statistics_widget.init_pop_tables)

        # Tab widget
        self.tab = QTabWidget()
        self.tab.addTab(self.input_files_widget, 'Input files')
        self.tab.addTab(self.sel_pops_widget, 'Populations')
        self.tab.addTab(self.mix_model_widget, 'Admixture model')
        self.tab.addTab(self.pca_widget, 'PCA')
        self.tab.addTab(self.f_statistics_widget, 'f-statistics')
        self.tab.currentChanged.connect(self.set_log_source)

        self.set_log_source(0)

        # Splitter
        splitter = QSplitter()
        splitter.setOrientation(Qt.Orientation.Vertical)
        splitter.addWidget(self.tab)
        splitter.addWidget(self.log_widget)

        # Layout
        layout = QVBoxLayout(self)
        layout.addWidget(splitter)

    @Slot(int)
    def set_log_source(self, index):
        if index == 0:
            if self.sel_pops_widget.log.is_changed_signal_connected():
                self.sel_pops_widget.log.changed.disconnect(self.log_widget.set_text)
            if self.mix_model_widget.log.is_changed_signal_connected():
                self.mix_model_widget.log.changed.disconnect(self.log_widget.set_text)
            if self.f_statistics_widget.log.is_changed_signal_connected():
                self.f_statistics_widget.log.changed.disconnect(self.log_widget.set_text)
            self.input_files_widget.log.changed.connect(self.log_widget.set_text)
            self.input_files_widget.log.set_text()
        elif index == 1:
            if self.input_files_widget.log.is_changed_signal_connected():
                self.input_files_widget.log.changed.disconnect(self.log_widget.set_text)
            if self.mix_model_widget.log.is_changed_signal_connected():
                self.mix_model_widget.log.changed.disconnect(self.log_widget.set_text)
            if self.f_statistics_widget.log.is_changed_signal_connected():
                self.f_statistics_widget.log.changed.disconnect(self.log_widget.set_text)
            self.sel_pops_widget.log.changed.connect(self.log_widget.set_text)
            self.sel_pops_widget.log.set_text()
        elif index == 2:
            if self.input_files_widget.log.is_changed_signal_connected():
                self.input_files_widget.log.changed.disconnect(self.log_widget.set_text)
            if self.sel_pops_widget.log.is_changed_signal_connected():
                self.sel_pops_widget.log.changed.disconnect(self.log_widget.set_text)
            if self.f_statistics_widget.log.is_changed_signal_connected():
                self.f_statistics_widget.log.changed.disconnect(self.log_widget.set_text)
            self.mix_model_widget.log.changed.connect(self.log_widget.set_text)
            self.mix_model_widget.log.set_text()
        elif index == 4:
            if self.sel_pops_widget.log.is_changed_signal_connected():
                self.sel_pops_widget.log.changed.disconnect(self.log_widget.set_text)
            if self.mix_model_widget.log.is_changed_signal_connected():
                self.mix_model_widget.log.changed.disconnect(self.log_widget.set_text)
            if self.sel_pops_widget.log.is_changed_signal_connected():
                self.sel_pops_widget.log.changed.disconnect(self.log_widget.set_text)
            self.f_statistics_widget.log.changed.connect(self.log_widget.set_text)
            self.f_statistics_widget.log.set_text()

    def closeEvent(self, event):
        self.mix_model_widget.plots_panel = None
        self.input_files_widget.about_dialog.close()
        event.accept()
