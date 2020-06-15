# StringMethodMaster
## Molecular Dynamics implementing string method with swarms of trajectories
The main project folder is P10_Production it should contain the latest files.

In order to be anle to run this software the following software and packages
need to be installed.
* Gromacs 2020.1
* Python 3.7
* gmax api 0.2.0a1
* NumPy
* SciPy
* Matplotlib
* Plotly
* pandas

Before running the software, 5 files are needed:
* Start and end state configuration files in .gro format
* topology file in .top format
* protein restraint file in .itp format (not necessary if running in vacuum)
* String method parameter file, according to template.

Once all these files are in your working directory and the parameter file is
updated and named "parameters.smg" the program should be runnable using
>pyhton/python3 SM_run.py
