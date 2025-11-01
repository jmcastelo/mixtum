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

from gui.plots import Plot

from pathlib import Path

from PySide6.QtCore import Qt, Slot
from PySide6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QTableWidget, QHeaderView, QPushButton, QSizePolicy
from PySide6.QtWidgets import QTableWidgetItem, QTabWidget, QFileDialog

import numpy as np



class PCAWidget(QWidget):
    def __init__(self, core):
        QWidget.__init__(self)

        # Core
        self.core = core

        # Selected populations table widget
        self.sel_pops_table = QTableWidget()
        self.sel_pops_table.setColumnCount(1)
        self.sel_pops_table.setSortingEnabled(True)
        self.sel_pops_table.verticalHeader().setVisible(False)
        self.sel_pops_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.sel_pops_table.setHorizontalHeaderLabels(['Selected populations'])

        # Selected populations PCA table widget
        self.sel_pops_pca_table = QTableWidget()
        self.sel_pops_pca_table.setColumnCount(1)
        self.sel_pops_pca_table.setSortingEnabled(True)
        self.sel_pops_pca_table.verticalHeader().setVisible(False)
        self.sel_pops_pca_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.sel_pops_pca_table.setHorizontalHeaderLabels(['PCA populations'])

        self.pca_items = []
        self.pca_indices = {}
        self.pca_names = []

        # Compute button
        self.compute_button = QPushButton('Compute PCA')
        self.compute_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.compute_button.setEnabled(False)

        # Compute button
        self.save_button = QPushButton('Save PCA data')
        self.save_button.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        self.save_button.setEnabled(False)

        buttons_layout = QHBoxLayout()
        buttons_layout.addWidget(self.compute_button)
        buttons_layout.addWidget(self.save_button)

        # PCA Plots
        self.pca_plot_3d = Plot('PCA', 'PC1', 'PC2', 5, 5, 100, projection='3d', zlabel='PC3', multi_selectable=True)
        self.pca_plot_2d = Plot('PCA', 'PC1', 'PC2', 5, 5, 100, multi_selectable=True)

        # PCA Plots tab widget
        pca_tab_widget = QTabWidget()
        pca_tab_widget.addTab(self.pca_plot_2d, 'PCA: 2D')
        pca_tab_widget.addTab(self.pca_plot_3d, 'PCA: 3D')

        # Controls
        controls_layout = QVBoxLayout()
        controls_layout.addWidget(self.sel_pops_table)
        controls_layout.addLayout(buttons_layout)

        layout = QHBoxLayout(self)
        layout.addLayout(controls_layout, 1)
        layout.addWidget(self.sel_pops_pca_table, 1)
        layout.addWidget(pca_tab_widget, 2)

        # Connections
        self.sel_pops_table.itemSelectionChanged.connect(self.sel_pops_changed)
        self.sel_pops_pca_table.itemSelectionChanged.connect(self.plot_sel_pca_points)
        self.pca_plot_2d.selected_indices_changed.connect(self.select_pops_pca_2d)
        self.pca_plot_3d.selected_indices_changed.connect(self.select_pops_pca_3d)
        self.compute_button.clicked.connect(self.compute_pca)
        self.compute_button.clicked.connect(self.init_sel_pops_pca_table)
        self.save_button.clicked.connect(self.save_pca_data)

    @Slot(bool)
    def init_sel_pops_table(self, result):
        if not result:
            return

        self.sel_pops_table.clearContents()
        self.sel_pops_table.setRowCount(len(self.core.selected_pops))

        for index, pop in enumerate(self.core.selected_pops):
            item = QTableWidgetItem(pop)
            item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            self.sel_pops_table.setItem(index, 0, item)

    @Slot()
    def init_sel_pops_pca_table(self):
        self.sel_pops_pca_table.clearContents()
        self.sel_pops_pca_table.setRowCount(len(self.pca_items))

        self.pca_indices = {}
        self.pca_names.clear()

        for index, item in enumerate(self.pca_items):
            pca_item = QTableWidgetItem(item.text())
            pca_item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            self.sel_pops_pca_table.setItem(index, 0, pca_item)

            # Store indices and poop names to link points and pops
            self.pca_indices[pca_item.text()] = index
            self.pca_names.append(pca_item.text())

    @Slot()
    def sel_pops_changed(self):
        self.compute_button.setEnabled(len(self.sel_pops_table.selectedItems()) > 2)

    @Slot(list)
    def select_pops_pca_2d(self, indices):
        self.sel_pops_pca_table.clearSelection()
        for index in indices:
            pca_items = self.sel_pops_pca_table.findItems(self.pca_names[index], Qt.MatchFlag.MatchExactly)
            pca_items[0].setSelected(True)

        points = np.take(self.core.principal_components, indices, 0)
        self.pca_plot_3d.plot_multiple_selected_points(points, indices)

    @Slot(list)
    def select_pops_pca_3d(self, indices):
        self.sel_pops_pca_table.clearSelection()
        for index in indices:
            pca_items = self.sel_pops_pca_table.findItems(self.pca_names[index], Qt.MatchFlag.MatchExactly)
            pca_items[0].setSelected(True)

        points = np.take(self.core.principal_components, indices, 0)
        self.pca_plot_2d.plot_multiple_selected_points(points, indices)

    @Slot()
    def plot_sel_pca_points(self):
        indices = []
        for item in self.sel_pops_pca_table.selectedItems():
            indices.append(self.pca_indices[item.text()])

        points = np.take(self.core.principal_components, indices, 0)
        self.pca_plot_2d.plot_multiple_selected_points(points, indices)
        self.pca_plot_3d.plot_multiple_selected_points(points, indices)

    @Slot()
    def compute_pca(self):
        self.pca_items = self.sel_pops_table.selectedItems()
        pops = [item.text() for item in self.pca_items]

        self.core.compute_pca(pops)

        xlabel = f"PC1 {self.core.explained_variance[0]:.1f}%"
        ylabel = f"PC2 {self.core.explained_variance[1]:.1f}%"
        zlabel = f"PC3 {self.core.explained_variance[2]:.1f}%"

        self.pca_plot_2d.plot_pca_2d(self.core.principal_components, 'PCA: 2D', xlabel, ylabel)
        self.pca_plot_3d.plot_pca_3d(self.core.principal_components, 'PCA: 3D', xlabel, ylabel, zlabel)

        self.save_button.setEnabled(True)

    @Slot()
    def save_pca_data(self):
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.FileMode.AnyFile)
        dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)

        if dialog.exec():
            file_names = dialog.selectedFiles()
            file_path = file_names[0]
            self.core.save_pca_data(Path(file_path))