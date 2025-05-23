from PySide6.QtWidgets import QWidget, QVBoxLayout

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import numpy as np



class MatplotlibCanvas(FigureCanvasQTAgg):
    def __init__(self, parent = None, width = 5, height = 4, dpi = 100):
        self.fig = Figure(figsize = (width, height), dpi = dpi)
        super(MatplotlibCanvas, self).__init__(self.fig)



class Plot(QWidget):
    def __init__(self, title, xlabel, ylabel, width, height, dpi, show_axes = True, show_toolbar = True, polar = False):
        QWidget.__init__(self)

        self.canvas = MatplotlibCanvas(self, width, height, dpi)

        self.axes = self.canvas.fig.add_subplot(polar = polar)
        
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

        self.canvas.fig.tight_layout()

        layout = QVBoxLayout(self)

        if show_toolbar:
            toolbar = NavigationToolbar(self.canvas, self)
            layout.addWidget(toolbar)

        layout.addWidget(self.canvas)

    def plot_fit(self, x, y, alpha, title, xlabel, ylabel):
        self.axes.clear()

        self.axes.set_title(title)

        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)

        self.axes.axhline(y = 0, color = '0.6', lw = 0.5)
        self.axes.axvline(x = 0, color = '0.6', lw = 0.5)

        self.axes.axline((0, 0), slope = 1, color = '0.6', lw = 0.5)

        self.axes.plot(x, y, '.')
        self.axes.plot(x, alpha * x)

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