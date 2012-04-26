import os

def executeCMD(cmd,verbose=False):
    if verbose: print cmd

    os.system(cmd)


