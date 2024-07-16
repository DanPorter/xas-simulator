from nexpy.gui.widgets import NXMessageBox
from nexpy.gui.utils import report_error
from nexusformat.nexus import NeXusError


def show_info_dialog():
    try:
        dialog = INFODialog()
        dialog.show()
    except NeXusError as error:
        report_error("Info dialogue error", error)

class INFODialog(NXMessageBox):
    def __init__(self, title, text):
        super().__init__(title, text)
    def init_from_file(self, title, filename):
        with open(filename,'r') as f:
            super().__init__(title, f.readlines())
