from gui.log_system import LogSystem
from gui.searchable_table_widget import SearchableTableWidget
from gui.worker import Worker

from PySide6.QtCore import Qt, Signal, Slot, QThreadPool
from PySide6.QtWidgets import QWidget, QTableWidget, QAbstractScrollArea, QTableWidgetItem, QPushButton, QSizePolicy, QFrame, QLabel, QSpinBox, QProgressBar, QVBoxLayout, QHBoxLayout, QFormLayout, QGridLayout, QGroupBox



class SelectPopsWidget(QWidget):
    computation_result = Signal(bool)

    def __init__(self, core):
        QWidget.__init__(self)

        # Core
        self.core = core

        # Thread pool
        self.thread_pool = QThreadPool()

        # Log
        self.log = LogSystem(['main', 'progress', 'timing', 'check'])
        self.log.set_entry('main', 'Select entries from the available populations on the left table and compute their allele frequencies.')
        self.log.append_entry('timing', '')
        self.log.append_entry('check', '')

        # Searchable table containing available populations
        self.search_widget = SearchableTableWidget()

        # Table containing selected populations
        self.selected_table = QTableWidget()
        self.selected_table.setColumnCount(1)
        self.selected_table.verticalHeader().setVisible(False)
        self.selected_table.horizontalHeader().setVisible(False)
        self.selected_table.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)

        # Select populations button
        self.select_button = QPushButton('Select populations')
        self.select_button.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        self.select_button.clicked.connect(self.select_populations)

        # Remove populations button
        self.remove_button = QPushButton('Deselect populations')
        self.remove_button.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        self.remove_button.clicked.connect(self.remove_populations)

        # Reset populations button
        self.reset_button = QPushButton('Reset populations')
        self.reset_button.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        self.reset_button.clicked.connect(self.reset_populations)

        # Number of computing processes
        self.num_procs_spin_box = QSpinBox()
        self.num_procs_spin_box.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        self.num_procs_spin_box.setMinimum(1)
        self.num_procs_spin_box.setMaximum(9999)
        self.num_procs_spin_box.valueChanged.connect(self.core.set_num_procs)
        self.num_procs_spin_box.setValue(1)

        # Compute allele frequencies files button
        self.comp_button = QPushButton('Compute frequencies')
        self.comp_button.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        self.comp_button.clicked.connect(self.compute_frequencies)

        # Stop computation button
        self.stop_button = QPushButton('Stop computation')
        self.stop_button.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        self.stop_button.clicked.connect(self.core.stop_computation)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(1)
        self.progress_bar.setValue(0)

        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)

        # Select table buttons layout
        hlayout = QHBoxLayout()
        hlayout.addWidget(self.remove_button, 0, Qt.AlignJustify)
        hlayout.addWidget(self.reset_button, 0, Qt.AlignJustify)

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
        slayout.addWidget(self.select_button, 0, Qt.AlignHCenter)

        # Search group box
        search_group = QGroupBox('Search and select from available populations')
        search_group.setLayout(slayout)

        # Tables layout
        tlayout = QHBoxLayout()
        tlayout.addWidget(search_group)
        tlayout.addWidget(select_group)

        # Form layout
        flayout = QFormLayout()
        flayout.addRow('Number of processes:', self.num_procs_spin_box)

        # Computation layout
        glayout = QGridLayout()
        glayout.addLayout(flayout, 0, 0, 1, 2, Qt.AlignLeft)
        glayout.addWidget(self.comp_button, 1, 0, Qt.AlignRight)
        glayout.addWidget(self.stop_button, 1, 1, Qt.AlignRight)
        glayout.addWidget(self.progress_bar, 0, 2, 2, 1)

        # Layout
        layout = QVBoxLayout(self)
        layout.addLayout(tlayout)
        layout.addWidget(separator)
        layout.addLayout(glayout)

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
            table_widget_item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            self.selected_table.setItem(index, 0, table_widget_item)

    @Slot(str, str, int)
    def log_progress(self, key, message, line_num):
        self.log.set_entry(key, message, line_num)

    @Slot()
    def compute_frequencies(self):
        self.log.clear_entry('main')
        self.log.clear_entry('progress')
        self.log.clear_entry('timing')
        self.log.clear_entry('check')

        self.progress_bar.setMaximum(len(self.core.selected_pops))

        worker = Worker('freqs', self.core.parallel_compute_populations_frequencies)
        worker.signals.progress[str, str, int].connect(self.log_progress)
        worker.signals.progress[int].connect(self.progress_bar.setValue)
        worker.signals.result.connect(self.computation_result)

        self.thread_pool.start(worker)
