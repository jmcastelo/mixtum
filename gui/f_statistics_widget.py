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

from PySide6.QtCore import Qt, Slot
from PySide6.QtWidgets import QWidget, QTableWidget, QTableWidgetItem, QHeaderView
from PySide6.QtWidgets import QPushButton, QSizePolicy, QVBoxLayout, QHBoxLayout



class FStatisticsWidget(QWidget):
    def __init__(self, core):
        QWidget.__init__(self)

        # Core
        self.core = core

        # Log
        self.log = LogSystem(['main', 'f2', 'f3', 'f4'])
        self.log.set_entry('main', 'Choose populations, then compute results.')
        self.log.set_entry('f2', '')
        self.log.set_entry('f3', '')
        self.log.set_entry('f4', '')

        # f2 table widget
        self.f2_table = QTableWidget()
        self.f2_table.setColumnCount(1)
        self.f2_table.setSortingEnabled(True)
        self.f2_table.verticalHeader().setVisible(False)
        self.f2_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.f2_table.setHorizontalHeaderLabels(['f2 Populations'])

        # f2 compute button
        self.f2_compute_button = QPushButton('Compute f2')
        self.f2_compute_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.f2_compute_button.setEnabled(False)

        # f2 layout
        f2_layout = QVBoxLayout()
        f2_layout.addWidget(self.f2_table)
        f2_layout.addWidget(self.f2_compute_button, 0, Qt.AlignmentFlag.AlignCenter)

        # f3 table widget
        self.f3_table = QTableWidget()
        self.f3_table.setColumnCount(1)
        self.f3_table.setSortingEnabled(True)
        self.f3_table.verticalHeader().setVisible(False)
        self.f3_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.f3_table.setHorizontalHeaderLabels(['f3 Populations'])

        # f3 compute button
        self.f3_compute_button = QPushButton('Compute f3')
        self.f3_compute_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.f3_compute_button.setEnabled(False)

        # f3 layout
        f3_layout = QVBoxLayout()
        f3_layout.addWidget(self.f3_table)
        f3_layout.addWidget(self.f3_compute_button, 0, Qt.AlignmentFlag.AlignCenter)

        # f4 table widget
        self.f4_table = QTableWidget()
        self.f4_table.setColumnCount(1)
        self.f4_table.setSortingEnabled(True)
        self.f4_table.verticalHeader().setVisible(False)
        self.f4_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.f4_table.setHorizontalHeaderLabels(['f4 Populations'])

        # f4 compute button
        self.f4_compute_button = QPushButton('Compute f4')
        self.f4_compute_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.f4_compute_button.setEnabled(False)

        # f4 layout
        f4_layout = QVBoxLayout()
        f4_layout.addWidget(self.f4_table)
        f4_layout.addWidget(self.f4_compute_button, 0, Qt.AlignmentFlag.AlignCenter)

        layout = QHBoxLayout(self)
        layout.addLayout(f2_layout)
        layout.addLayout(f3_layout)
        layout.addLayout(f4_layout)

        # Connections
        self.f2_table.itemSelectionChanged.connect(self.f2_pops_changed)
        self.f3_table.itemSelectionChanged.connect(self.f3_pops_changed)
        self.f4_table.itemSelectionChanged.connect(self.f4_pops_changed)
        self.f2_compute_button.clicked.connect(self.compute_f2)
        self.f3_compute_button.clicked.connect(self.compute_f3)
        self.f4_compute_button.clicked.connect(self.compute_f4)

    @Slot()
    def reset_controls(self):
        self.f2_table.setRowCount(0)
        self.f3_table.setRowCount(0)
        self.f4_table.setRowCount(0)

        self.log.clear_entry('main')
        self.log.append_entry('main', 'Choose populations, then compute results.')
        self.log.clear_entry('f2')
        self.log.clear_entry('f3')
        self.log.clear_entry('f4')

    @Slot(bool)
    def init_pop_tables(self, result):
        if not result:
            return

        self.populate_table_widget(self.f2_table)
        self.populate_table_widget(self.f3_table)
        self.populate_table_widget(self.f4_table)

    def populate_table_widget(self, table):
        table.clearContents()
        table.setRowCount(len(self.core.selected_pops))

        for index, pop in enumerate(self.core.selected_pops):
            item = QTableWidgetItem(pop)
            item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            table.setItem(index, 0, item)

    @Slot()
    def f2_pops_changed(self):
        self.f2_compute_button.setEnabled(len(self.f2_table.selectedItems()) == 2)

    @Slot()
    def f3_pops_changed(self):
        self.f3_compute_button.setEnabled(len(self.f3_table.selectedItems()) == 3)

    @Slot()
    def f4_pops_changed(self):
        self.f4_compute_button.setEnabled(len(self.f4_table.selectedItems()) == 4)

    @Slot()
    def compute_f2(self):
        pops = [item.text() for item in  self.f2_table.selectedItems()]
        f2 = self.core.compute_f2(pops)
        self.log.clear_entry('f2')
        self.log.append_entry('f2', f"f2({pops[0]}, {pops[1]}) = {f2:6.4f}")

    @Slot()
    def compute_f3(self):
        pops = [item.text() for item in self.f3_table.selectedItems()]
        self.log.clear_entry('f3')
        for n in range(3):
            rot_pops = pops[n:] + pops[:n]
            f3, angle = self.core.compute_f3(rot_pops)
            self.log.append_entry('f3', f"f3({rot_pops[0]}, {rot_pops[1]}; {rot_pops[2]}) = {f3:6.4f} , angle = {angle:6.2f} deg")

    @Slot()
    def compute_f4(self):
        pops = [item.text() for item in self.f4_table.selectedItems()]
        self.log.clear_entry('f4')
        right_pops = pops[1:]
        for n in range(3):
            rot_pops = [pops[0]] + right_pops[n:] + right_pops[:n]
            f4, angle = self.core.compute_f4(rot_pops)
            self.log.append_entry('f4', f"f4({rot_pops[0]}, {rot_pops[1]}; {rot_pops[2]}, {rot_pops[3]}) = {f4:6.4f} , angle = {angle:6.2f} deg")