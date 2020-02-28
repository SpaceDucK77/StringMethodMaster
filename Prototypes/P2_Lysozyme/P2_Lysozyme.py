## Second prototype to implement String method using gromacs API
## Started: 2020-01-29
## Author: Marko Petrovic
##

import gmxapi as gmx
import subprocess as sp

QUITC = "QQ"
QUIT="If you wish to quit the program, you can do so by typing " +\
      QUITC + " at any input\n"

class stop(Exception):
    pass

""" Asks for a file name and returns an existing one. """
def select_file(query):
    file_found =  False
    while not file_found:
        try:
            file_name = input(query)
            if file_name == QUITC:
                raise stop()
            file = open(file_name)
            file_found = True
            file.close()
        except FileNotFoundError:
            print("""\nThere does not appear to be such a file, common causes are:
 - The file was misspelled
 - The file is in a different directory
 - The file name includes a file ending (.pdb, .txt, .docx) and needs to be added to your input\n""")
            print(QUIT)
    return file_name

""" Asks the user to confirm an action. """
def confirm(action):
    no_answer = True
    while no_answer:
        reply = input("Would you like to " + action + "? (yes/no)")
        if reply == QUITC:
            raise stop()
        reply=reply.lower()
        if reply not in ("y", "yes", "n","no"):
            print("\nPlease use either yes or no\n")
            print(QUIT)
            continue
        no_answer =  False
    return reply in ("y", "yes")


""" Creates a new .pdb file from another without crystal water. """
def clean_h2o(file_name, ask_confirmation = True):
    if ask_confirmation and not confirm("clear crystal water from .pdb file"):
        return file_name
    prefix = "cleaned_h2o_"
    reads = open(file_name)
    out = open(prefix + file_name, "w")
    cleaned = 0
    for row in reads:
        if "HOH" not in row:
            out.write(row)
        else:
            cleaned += 1
    if cleaned == 0:
        print("No crystal water found")
    else:
        print(str(cleaned)+" crystal waters removed")
    reads.close()
    out.close()
    return prefix + file_name

def main():
    print('''Welcome to an automated run of the Gromacs Lysozyme tutorial, this program
should work on most systems containable within the same simulation box.
Some files are necessary for running the program:

 - a *.pdb file describing a protein.
 - an "ions.mdp" file for the ion placement configurations
 - an "minim.mdp" file for energy minimization configurations
 - "nvt.mdp" AND "npt.mdp" files for equilibration"
 - "md.mdp" for the md run configurations.\n\n
 ''')
    try:
        #file_name = select_file("Please type the .pdb files name: ")
        file_name = "1aki.pdb"
        clean_name = clean_h2o(file_name, False)
        out_name = "processed_" + ".".join(file_name.split(".")[:-1]) + ".gro"
        box_name = "newbox_" +  ".".join(file_name.split(".")[:-1]) + ".gro"
        solv_name = "solv_" +  ".".join(file_name.split(".")[:-1]) + ".gro"
        is_name = "ion_" + solv_name

        """ Haven't been able to get this to work lately, meant to delete old backups
            might also unintentionally delete valuable work so not very safe. """
        sp.run(["rm", "\#*"])

        # First command line instruction
        pdb2gmx=gmx.commandline_operation(executable = "gmx",
                                          arguments = ["pdb2gmx", "-water", "spce",
                                                       "-ff", "oplsaa"],
                                          input_files =  {"-f": clean_name},
                                          output_files = {"-o" : out_name})
        pdb2gmx.run()
        print("pdb2gmx",pdb2gmx.output.erroroutput.result())
        file_name = "conf.gro"
        box = gmx.commandline_operation(executable ="gmx",
                                        arguments = [ "editconf","-c", "-d", "1.0", "-bt", "cubic"],
                                        input_files = {"-f": pdb2gmx.output.file["-o"]},
                                        output_files= {"-o": box_name})
        box.run()
        print("box",box.output.erroroutput.result())

        """ Until I get a grip on what I'm doing wrong with commandline_operation """
        solvate = gmx.commandline_operation(executable = "gmx",
                                            arguments = ["solvate"],
                                            input_files = {"-cp": box.output.file["-o"],
                                                           "-p": "topol.top"},
                                            output_files = {"-o": solv_name})

        solvate.run()
        print("solvate",solvate.output.erroroutput.result())

        #ion_name = select_file("Please type the name of you ion .mdp files name")
        ion_name = "ions.mdp"
        ion_top =  ".".join(ion_name.split(".")[:-1]) + ".tpr"
        pname = "NA" #get_pname()
        nname = "CL" #get_nname()
        tpr_assemble = gmx.commandline_operation(executable = "gmx",
                                                 arguments = ["grompp"],
                                                 input_files = {"-f": ion_name,
                                                                "-c": solvate.output.file["-o"],
                                                                "-p": "topol.top"},
                                                 output_files = {"-o": ion_top})
        tpr_assemble.run()
        print("tpr_assemble",tpr_assemble.output.erroroutput.result())


        #sp.run(["gmx", "genion", "-s", ion_top, "-o", is_name, "-p","topol.top",\
        #        "-pname", pname, "-nname", nname, "-neutral"], input="SOL", text=True)



        genion = gmx.commandline_operation(executable = "gmx",
                                           arguments =  ["genion", "-pname", pname, "-nname", nname, "-neutral"],
                                           input_files = {"-s": tpr_assemble.output.file["-o"],
                                                          "-p": "topol.top"},
                                           output_files = {"-o": is_name},
                                           stdin = "SOL")

        genion.run()
        print("genion", genion.output.erroroutput.result())
        # Energy minimization
        emtprprep = gmx.commandline_operation(executable = "gmx",
                                              arguments = ["grompp"],
                                              input_files = {"-f": "minim.mdp",
                                                             "-c": genion.output.file["-o"],
                                                             "-p": "topol.top"},
                                              output_files = {"-o": "em.tpr"})
        emtprprep.run()
        print ("emtprprep", emtprprep.output.erroroutput.result())
        #sp.run(["gmx", "grompp", "-f", "minim.mdp", "-c", is_name, "-p", "topol.top",\
        #        "-o", "em.tpr"])
        emtpr = gmx.read_tpr(emtprprep.output.file["-o"])
        em =  gmx.mdrun(emtpr)
        em.run()


        #sp.run(["gmx", "mdrun", "-deffnm", "em"])
        emgro=em.output.trajectory.result()
        emgro=emgro[:emgro.rfind('/') + 1] + "confout.gro"

        # Equilibration
        eq1tpr =  gmx.commandline_operation(executable = "gmx",
                                            arguments = ["grompp"],
                                            input_files = {"-f": "nvt.mdp",
                                                           "-c": emgro,
                                                           "-r": emgro,
                                                           "-p": "topol.top"},
                                            output_files = {"-o": "nvt.tpr"})
        eq1tpr.run()
        print ("eq1tpr", eq1tpr.output.erroroutput.result())
        eq1tp = gmx.read_tpr(eq1tpr.output.file["-o"])
        eq1=gmx.mdrun(eq1tp)
        eq1.run()
        #sp.run(["gmx", "grompp", "-f", "nvt.mdp", "-c", "em.gro", "-r", "em.gro", "-p",\
        #        "topol.top", "-o",  "nvt.tpr"])
        #sp.run(["gmx", "mdrun", "-deffnm", "nvt"])
        eq1gro=eq1.output.trajectory.result()
        eq1gro=eq1gro[:eq1gro.rfind('/') + 1] + "confout.gro"
        

        
        # Equilibration 2
        eq2tpr =  gmx.commandline_operation(executable = "gmx",
                                            arguments = ["grompp"],
                                            input_files = {"-f": "npt.mdp",
                                                           "-c": eq1gro,
                                                           "-r": eq1gro,
                                                           "-p": "topol.top"},
                                            output_files = {"-o": "npt.tpr"})
        eq2tpr.run()
        print ("eq2tpr", eq2tpr.output.erroroutput.result())
        eq2tp = gmx.read_tpr(eq2tpr.output.file["-o"])
        eq2 = gmx.mdrun(eq2tp)
        eq2.run()
        eq2out = eq2.output.trajectory.result()
        eq2out = eq2out[:eq2out.rfind('/') + 1] 





        mdtpr = gmx.commandline_operation(executable = "gmx",
                                          arguments = ["grompp"],
                                          input_files = {"-f": "md.mdp",
                                                         "-c": eq2out + "confout.gro",
                                                         "-t": eq2out + "state.cpt",
                                                         "-p": "topol.top"},
                                          output_files = {"-o": "md_0_1.tpr"})
        mdtpr.run()

        mdtp = gmx.read_tpr(mdtpr.output.file["-o"])
        md = gmx.mdrun(mdtp)
        md.run()

        print("Finished")
                                                         

        # RUN
        """sp.run(["gmx", "grompp", "-f", "md.mdp", "-c", "npt.gro", "-t", "npt.cpt", "-p", \
                "topol.top", "-o", "md_0_1.tpr"])
        sp.run(["gmx", "mdrun", "-deffnm", "md_0_1"])"""
    except stop:
        print ("*** Exiting program per users wishes ***")
        return

if __name__ == "__main__":
    main()
