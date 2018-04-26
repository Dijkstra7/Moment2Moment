import cx_Freeze
import sys
import matplotlib
import os

os.environ['TCL_LIBRARY'] = "C:/ProgramData/Anaconda3/pkgs/python-3.6.2-h6679aeb_11/tcl/tcl8.6"
os.environ['TK_LIBRARY'] = "C:/ProgramData/Anaconda3/pkgs/python-3.6.2-h6679aeb_11/tcl/tk8.6"
base = None

if sys.platform == 'win32':
    base = "Win32GUI"

executables = [cx_Freeze.Executable("GraphsGUI.py", base=base, icon="res/m2micon.ico")]

cx_Freeze.setup(
    name = "Moment to moment graph generator",
    options = {"build_exe": {"packages":["matplotlib", "DataHandler"], "include_files":["res/m2micon.ico"]}},
    version = "0.01",
    description = "Moment to moment graph generator",
    executables = executables
    )