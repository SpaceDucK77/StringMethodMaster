import gmxapi as gmx
# import numpy as np

## create tpr file:
## gmx grompp -f md.mdp -c conf.gro -p topol.top -o topol.tpr -maxwarn 3

## numpy  and networkx needed for initial testrun

trial_input = gmx.read_tpr("topol.tpr")
md = gmx.mdrun(trial_input)
md.run()
