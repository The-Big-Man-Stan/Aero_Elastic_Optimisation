import subprocess
import os

def clear_para_file():
        process = subprocess.call(['cp','Optimisation_scripts/rae2822.bmap','rae2822.bmap'])
	process = subprocess.call(['cp','Optimisation_scripts/rae2822.para_primal','rae2822.para_primal'])
	return
clear_para_file()
