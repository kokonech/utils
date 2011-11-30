import os
import argparse



def singleCharacter(string):
    value = str(string)
    if len(string) != 1:
        raise argparse.ArgumentTypeError("%s is not a single character")
    return value


def existsDir(folder):
    if not os.path.isdir(folder):
        raise argparse.ArgumentError("%s is not an existing directory"%folder)
    return folder

def existsFile(f):
    if not os.path.isfile(f):
        raise argparse.ArgumentTypeError("%s is not an existing file"%f)
    return f


