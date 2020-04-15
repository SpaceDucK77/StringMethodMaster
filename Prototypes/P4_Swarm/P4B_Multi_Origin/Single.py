import gmxapi as gmx
import subprocess as sp

def translate(clean_name, out_name):
    sp.run(["gmx", "pdb2gmx", "-water", "spce", "-ff", "oplsaa", "-f", clean_name,\
            "-o", out_name])
    """
    test=gmx.commandline_operation(executable = "pdb2gmx",
                                   arguments = ["-water", "spce"
                                                ,"-ff", "oplsaa"
                                   ]
                                   #,"-f",clean_name
                                   #,"-o",out_name)
                                   ,
                                   #)
                                   input_files = {"-f":clean_name},
                                   output_files = {"-o":out_name})
    test.run()"""

def create_box():
    sp.run(["gmx", "editconf", "-c", "-d", "1.0", "-bt", "cubic", "-f", out_name, \
            "-o", box_name])

def solvate():
    sp.run(["gmx", "solvate", "-cp", box_name, "-o", solv_name, "-p", "topol.top"])
    
def equilibrate():
    # Equilibration 1
    sp.run(["gmx", "grompp", "-f", "nvt.mdp", "-c", "em.gro", "-r", "em.gro", "-p",\
        "topol.top", "-o",  "nvt.tpr"])
    sp.run(["gmx", "mdrun", "-deffnm", "nvt"])

    # Equilibration 2
    sp.run(["gmx", "grompp", "-f", "npt.mdp", "-c", "nvt.gro", "-r", "nvt.gro", "-p",\
            "topol.top", "-o", "npt.tpr"])
    sp.run(["gmx", "mdrun", "-deffnm", "npt"])

    
def e_minimize():
    sp.run(["gmx", "grompp", "-f", "minim.mdp", "-c", is_name, "-p", "topol.top",\
            "-o", "em.tpr"])
    sp.run(["gmx", "mdrun", "-deffnm", "em"])

def ions():
    sp.run(["gmx", "grompp", "-f", ion_name, "-c", solv_name, "-p", "topol.top",\
            "-o", ion_top])
    sp.run(["gmx", "genion", "-s", ion_top, "-o", is_name, "-p","topol.top",\
            "-pname", pname, "-nname", nname, "-neutral"], input="SOL", text=True)

def single_run():
    sp.run(["gmx", "grompp", "-f", "md.mdp", "-c", "npt.gro", "-t", "npt.cpt", "-p", \
            "topol.top", "-o", "md_0_1.tpr"])
    sp.run(["gmx", "mdrun", "-deffnm", "md_0_1"])

def prep_run():
    #translate()
    #solvate()
    #ions()
    e_minimize()
    equilibrate()

def analyse():
    pass

def multi_run(n):
    results = None
    for i in range(n):
        prep_run()
        single_run()
    analyse()
    return results

def main():
    pass
