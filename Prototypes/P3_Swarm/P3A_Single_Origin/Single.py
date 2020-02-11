import gmxapi as gmx
import subprocess as sp

def solvate():
    pass

def equilibrate():
    pass

def e_minimize():
    pass

def ions():
    pass

def single_run():
    pass

def prep_run():
    solvate()
    ions()
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
    
