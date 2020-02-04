import gmxapi as gmx
# import numpy as np

## numpy  and networkx needed for initial testrun

trial_input = gmx.read_tpr("topol.tpr")
md = gmx.mdrun(trial_input)
md.run()
