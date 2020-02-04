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
    try:
        #file_name = select_file("Please type the .pdb files name: ")
        file_name = "1aki.pdb"
        file_name = clean_h2o(file_name,False)
        out_name = "processed_" + ".".join(file_name.split(".")[:-1]) + ".gro"
        box_name = "newbox_" +  ".".join(file_name.split(".")[:-1]) + ".gro"
        solv_name = "solv_" +  ".".join(file_name.split(".")[:-1]) + ".gro"
        is_name = "ion_" + solv_name
        # First command line instruction, currently not working, could not understand error message
        sp.run(["rm", "\#*"])
        sp.run(["gmx", "pdb2gmx", "-water", "spce", "-ff", "oplsaa", "-f", file_name,\
                "-o", out_name])
        """
        test=gmx.commandline_operation(executable = "pdb2gmx",
                                  arguments = ["-water", "spce"
                                               ,"-ff", "oplsaa"
                                               ]
                                               #,"-f",file_name
                                               #,"-o",out_name)
                                       ,
                                       #)
                                       input_files = {"-f":file_name},
                                       output_files = {"-o":out_name})
        test.run()"""
        #file_name = "conf.gro"
        #test2 = gmx.commandline_operation(executable = "editconf",
        #                                  arguments = ["-c", "-d", "1.0", "-bt", "cubic"])
        #test2.run()

        """ Until I get a grip on what Iäm doing wrong with commandline_operation """
        sp.run(["gmx", "editconf", "-c", "-d", "1.0", "-bt", "cubic", "-f", out_name, \
                "-o", box_name])
        sp.run(["gmx", "solvate", "-cp", box_name, "-o", solv_name, "-p", "topol.top"])
        #ion_name = select_file("Please type the name of you ion .mdp files name")
        ion_name = "ions.mdp"
        ion_top =  ".".join(ion_name.split(".")[:-1]) + ".tpr"
        pname = "NA" #get_pname()
        nname = "CL" #get_nname()
        sp.run(["gmx", "grompp", "-f", ion_name, "-c", solv_name, "-p", "topol.top",\
                "-o", ion_top])
        sp.run(["gmx", "genion", "-s", ion_top, "-o", is_name, "-p", "topol.top", "-pname", \
                pname, "-nname", nname, "-neutral"])
        
    except stop:
        print ("*** Exiting program per users wishes ***")
        return

if __name__ == "__main__":
    main()
