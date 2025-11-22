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

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QWidget, QVBoxLayout

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import numpy as np



class MatplotlibCanvas(FigureCanvasQTAgg):
    def __init__(self, parent = None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super(MatplotlibCanvas, self).__init__(self.fig)



class Plot(QWidget):
    # Signals
    selected_index_changed = Signal(int)
    selected_indices_changed = Signal(list)

    def __init__(self, title, xlabel, ylabel, width, height, dpi, show_axes=True, show_toolbar=True, polar=False, projection=None, zlabel=None, selectable=False, multi_selectable=False):
        QWidget.__init__(self)

        self.canvas = MatplotlibCanvas(self, width, height, dpi)

        self.axes = self.canvas.fig.add_subplot(polar=polar, projection=projection)

        self.projection = projection

        self.axes.set_title(title)

        if polar:
            self.axes.set_thetamin(0)
            self.axes.set_thetamax(180)

            self.axes.set_rticks([])

            self.axes.bar(0, 1, np.pi / 2, color=(1, 0.8, 0.4), alpha=0.5, align='edge')
            self.axes.bar(np.pi / 2, 1, np.pi / 2, color=(0.8, 1, 0.4), alpha=0.5, align='edge')
        else:
            self.axes.set_xlabel(xlabel)
            self.axes.set_ylabel(ylabel)

        if not show_axes:
            self.axes.axis('off')

        if zlabel is not None:
            self.axes.set_zlabel(zlabel)

        if selectable:
            self.selectable_plot = None
            self.sel_point_plot = None
            self.conn = self.canvas.mpl_connect('button_press_event', self.select_point)

        if multi_selectable:
            self.multi_selectable_plots = []
            self.sel_point_plot = None
            self.sel_indices = []
            self.conn = self.canvas.mpl_connect('button_press_event', self.multi_select_points)

        self.canvas.fig.tight_layout()

        layout = QVBoxLayout(self)

        if show_toolbar:
            toolbar = NavigationToolbar(self.canvas, self)
            layout.addWidget(toolbar)

        layout.addWidget(self.canvas)

    def clear(self, title, xlabel, ylabel, show_axes=True, polar=False):
        self.axes.clear()
        self.axes.set_title(title)

        if polar:
            self.axes.set_thetamin(0)
            self.axes.set_thetamax(180)

            self.axes.set_rticks([])

            self.axes.bar(0, 1, np.pi / 2, color=(1, 0.8, 0.4), alpha=0.5, align='edge')
            self.axes.bar(np.pi / 2, 1, np.pi / 2, color=(0.8, 1, 0.4), alpha=0.5, align='edge')
        else:
            self.axes.set_xlabel(xlabel)
            self.axes.set_ylabel(ylabel)

        if not show_axes:
            self.axes.axis('off')

        self.canvas.fig.legends = []

        self.canvas.fig.canvas.draw()

    def select_point(self, event):
        if self.selectable_plot is not None and event.inaxes == self.axes:
            cont, ind = self.selectable_plot.contains(event)
            if cont:
                index = ind['ind'][0]
                self.selected_index_changed.emit(index)

    def multi_select_points(self, event):
        if len(self.multi_selectable_plots) > 0 and event.inaxes == self.axes:
            for plt in self.multi_selectable_plots:
                cont, ind = plt.contains(event)
                if cont:
                    index = ind['ind'][0]
                    if index in self.sel_indices:
                        self.sel_indices.remove(index)
                    else:
                        self.sel_indices.append(index)
                    self.selected_indices_changed.emit(self.sel_indices)

    def plot_selected_point(self, x, y):
        if self.sel_point_plot is not None:
            self.sel_point_plot.remove()
        self.sel_point_plot = self.axes.scatter(x, y, c='cyan', s=50)

    def plot_multiple_selected_points(self, points: np.array, indices):
        if self.sel_point_plot is not None:
            self.sel_point_plot.remove()

        if self.projection == '3d':
            self.sel_point_plot = self.axes.scatter(points[:, 0], points[:, 1], points[:, 2], alpha=0.5, c='cyan', s=120)
        else:
            self.sel_point_plot = self.axes.scatter(points[:, 0], points[:, 1], alpha=0.5, c='cyan', s=120)

        self.canvas.fig.canvas.draw()

        self.sel_indices = indices

    def plot_fit(self, x, y, alpha, title, xlabel, ylabel):
        self.axes.clear()

        self.axes.set_title(title)

        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)

        self.axes.axhline(y=0, c='0.6', lw=0.5)
        self.axes.axvline(x=0, c='0.6', lw=0.5)

        self.axes.axline((0, 0), slope=1, color='0.6', lw=0.5)

        self.selectable_plot = self.axes.scatter(x=x, y=y, s=10)
        self.axes.plot(x, alpha * x, c='r', lw=0.5)

        self.canvas.fig.canvas.draw()

    def plot_histogram(self, histogram, title, xlabel, ylabel):
        counts = histogram[0]
        edges = histogram[1]

        self.axes.clear()

        self.axes.set_title(title)

        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)

        self.axes.bar(edges[:-1], counts, width = np.diff(edges), edgecolor = 'black', align = 'edge')

        self.canvas.fig.canvas.draw()

    def plot_bars(self, hybrid, parent1, parent2, alpha):
        self.axes.clear()

        self.axes.axis('off')

        self.axes.set_title(hybrid)

        percent = 100.0 * alpha

        hbar_p1 = self.axes.barh(y = 1, width = percent, left = 0, label = parent1, align = 'center', color = 'cyan')
        hbar_p2 = self.axes.barh(y = 1, width = (100.0 - percent), left = percent, label = parent2, align = 'center', color='yellow')

        self.axes.bar_label(hbar_p1, label_type = 'center', fmt = '%.2f%%')
        self.axes.bar_label(hbar_p2, label_type = 'center', fmt = '%.2f%%')

        self.axes.set_ymargin(0.3)
        self.axes.set_xmargin(0)

        self.canvas.fig.legends = []
        self.canvas.fig.legend(loc = 'outside lower center', bbox_to_anchor = (0.5, -0.02), ncols = 2, fontsize = 'small')

        self.canvas.fig.canvas.draw()

    def plot_angle(self, title, angles):
        self.axes.clear()

        self.axes.set_title(title)

        self.axes.set_thetamin(0)
        self.axes.set_thetamax(180)

        self.axes.set_rmax(1)
        self.axes.set_rticks([])

        self.axes.bar(0, 1, np.pi / 2, color=(1, 0.8, 0.4), alpha=0.5, align='edge')
        self.axes.bar(np.pi / 2, 1, np.pi / 2, color=(0.8, 1, 0.4), alpha=0.5, align='edge')

        self.axes.vlines(angles[0] * np.pi / 180, 0, 1, colors = 'C0', label = f"Pre-JL angle: {angles[0]:.2f} deg")
        self.axes.vlines(angles[1] * np.pi / 180, 0, 1, colors = 'C1', label = f"Post-JL angle: {angles[1]:.2f} deg")

        self.canvas.fig.legends = []
        self.canvas.fig.legend(loc = 'outside lower center', ncols = 1, fontsize = 'medium')

        self.canvas.fig.canvas.draw()

    def plot_pca_3d(self, pcs: np.array, title, xlabel, ylabel, zlabel):
        self.axes.clear()

        self.axes.set_title(title)

        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)
        self.axes.set_zlabel(zlabel)

        xmin = np.min(pcs[:, 0])
        xmax = np.max(pcs[:, 0])
        self.axes.set_xlim(xmin, xmax)

        ymin = np.min(pcs[:, 1])
        ymax = np.max(pcs[:, 1])
        self.axes.set_ylim(ymin, ymax)

        zmin = np.min(pcs[:, 2])
        zmax = np.max(pcs[:, 2])
        self.axes.set_zlim(zmin, zmax)

        self.multi_selectable_plots = [
            self.axes.scatter(pcs[:, 0], pcs[:, 1], pcs[:, 2], alpha=0.5, color='k', s=40),
            self.axes.scatter(pcs[:, 0], pcs[:, 1], color='r', zdir='z', zs=zmin, s=40),
            self.axes.scatter(pcs[:, 0], pcs[:, 2], color='g', zdir='y', zs=ymin, s=40),
            self.axes.scatter(pcs[:, 1], pcs[:, 2], color='b', zdir='x', zs=xmin, s=40)
        ]
        self.sel_indices = []

        self.sel_point_plot = None

        self.canvas.fig.canvas.draw()

    def plot_pca_2d(self, pcs: np.array, title, xlabel, ylabel):
        self.axes.clear()

        self.axes.set_title(title)

        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)

        self.multi_selectable_plots = [
            self.axes.scatter(pcs[:, 0], pcs[:, 1], color='r', s=40)
        ]
        self.sel_indices = []

        self.sel_point_plot = None

        self.canvas.fig.canvas.draw()