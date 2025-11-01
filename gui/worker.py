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

import traceback, sys

from PySide6.QtCore import QObject, QRunnable, Signal, Slot



class WorkerSignals(QObject):
    """
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data

    error
        tuple (exctype, value, traceback.format_exc() )

    result
        object data returned from processing, anything

    progress
        int indicating % progress
    """

    finished = Signal(str)
    error = Signal(tuple)
    result = Signal(bool)
    progress = Signal((int,), (float,), (str,), (str, str), (str, str, int))



class Worker(QRunnable):
    """
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function
    """

    def __init__(self, name, fn, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.name = name
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the callback to our kwargs
        # self.kwargs['progress_callback'] = self.signals.progress
        self.kwargs['progress_callback'] = self.progress_callback

    def progress_callback(self, *args):
        if len(args) == 1:
            if isinstance(args[0], int):
                self.signals.progress[int].emit(args[0])
            elif isinstance(args[0], float):
                self.signals.progress[float].emit(args[0])
            elif isinstance(args[0], str):
                self.signals.progress[str].emit(args[0])
        elif len(args) == 2 and isinstance(args[0], str) and isinstance(args[1], str):
            self.signals.progress[str, str].emit(args[0], args[1])
        elif len(args) == 3 and isinstance(args[0], str) and isinstance(args[1], str) and isinstance(args[2], int):
            self.signals.progress[str, str, int].emit(args[0], args[1], args[2])

    @Slot()
    def run(self):
        """
        Initialise the runner function with passed args, kwargs.
        """

        # Retrieve args/kwargs here; and fire processing using them
        try:
            result = self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            # Return the result of the processing
            self.signals.result.emit(result)
        finally:
            # Done
            self.signals.finished.emit(self.name)
