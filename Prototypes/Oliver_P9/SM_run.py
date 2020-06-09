import VIS_col as VC
from md_tools import *


if __name__ == "__main__":
    a = load()
    if a =={}:
        shutil.copyfile("topol_orig.top", "topol.top")
        a = VC.VIS_collection("parameters.smg")
    a.parse_parameters()
    a.prepare_endpoints()
    a.create_base_string()
    a.run_method()
