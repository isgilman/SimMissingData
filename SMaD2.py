#!/usr/bin/env python

try: import sys
except ImportError:
    print "\n\tError: sys is not installed/loaded"
    sys.exit()
try: import os
except ImportError:
    print "\n\tError: os is not installed/loaded"
    sys.exit()
try: from optparse import OptionParser
except ImportError:
    print "\n\tError: OptionParser (optparse) is not installed/loaded"
    sys.exit()
try: import shutil
except ImportError:
    print "\n\tError: shutil is not installed/loaded"
    sys.exit()
#try: from tqdm import tqdm, tnrange
#except ImportError:
#    print "\n\tError: tqdm is not installed/loaded."
#    print "\n\tTo install, try 'pip install tqdm' or 'conda install -c conda-forge tqdm'"
#    sys.exit()
try: import subprocess
except ImportError:
    print "\n\tError: subprocess is not installed/loaded."
try: import numpy as np
except ImportError:
    print "\n\tError: numpy is not installed/loaded."
try: import pandas as pd
except ImportError:
    print "\n\tError: pandas is not installed/loaded."
try: from Bio import AlignIO, Alphabet
except ImportError:
    print "\n\tError: Bio is not installed/loaded."
try: import re
except ImportError:
    print "\n\tError: re is not installed/loaded."
try: import dendropy
except ImportError:
    print "\n\tError: dendropy is not installed/loaded."
try: import csv
except ImportError:
    print "\n\tError: csv is not installed/loaded."
try: from SMDfx import listabs, Vividict, phylip_to_df, write_missing_data
except ImportError:
    print "\n\tError: SMD is not installed/loaded."
try: from subprocess import Popen, PIPE
except ImportError:
    print "\n\tError: subprocess is not installed/loaded."
from pprint import pprint
try: import SMDfx as smd
except ImportError:
    print "\n\tError: SMDfx is not installed/loaded."
try: from multiprocessing import Process, current_process, log_to_stderr, Pool
except ImportError:
    print "\n\tError: multiprocessing is not installed/loaded."
import SMaDfx2 as sma
try: from multiprocessing.managers import SyncManager, BaseProxy
except ImportError:
    print "\n\tError: multiprocessing is not installed/loaded."
import signal
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simdir = '/mnt/lfs2/dtank/ian/pSims'

phylips = []
for directory in listabs(simdir):
    if os.path.isdir(directory):
        try: 
            for filename in listabs(os.path.join(directory, 'outfiles')):
                if (filename.endswith('.phy')) and ('bN' in filename) and  ('L10000S' in filename): phylips.append(filename)
        except OSError:
            continue

#MedBalTree = '(((((A:2,B:2):2,C:4):4,D:8):4,(((E:2,F:2):2,G:4):4,H:8):4):4,(((I:2,J:2):2,K:4):4,L:8):8);'
#MedUnbTree = '((((((((((((A:1,B:1):1,C:2):1,D:3):1,E:4):1,F:5):1,G:6):1,H:7):1,I:8):1,J:9):1,K:10):1,L:11):12);'
#ShortBalTree = smd.multiply_strints(MedBalTree, 0.5)
#LongBalTree = smd.multiply_strints(MedBalTree, 2.0)
#ShortUnbTree = smd.multiply_strints(MedUnbTree, 0.5)
#LongUnbTree = smd.multiply_strints(MedUnbTree, 2.0)

#treeshapes = smd.Vividict()

#treeshapes['sb'] = ShortBalTree
#treeshapes['mb'] = MedBalTree
#treeshapes['lb'] = LongBalTree
#treeshapes['su'] = ShortUnbTree
#treeshapes['mu'] = MedUnbTree
#treeshapes['lu'] = LongUnbTree

#MedBalTree = '(((((A_0:2,B_0:2):2,C_0:4):4,D_0:8):4,(((E_0:2,F_0:2):2,G_0:4):4,H_0:8):4):4,(((I_0:2,J_0:2):2,K_0:4):4,L_0:8):8);'
MedUnbTree = '((((((((((((A_0:1,B_0:1):1,C_0:2):1,D_0:3):1,E_0:4):1,F_0:5):1,G_0:6):1,H_0:7):1,I_0:8):1,J_0:9):1,K_0:10):1,L_0:11):12);'
#ShortBalTree = smd.multiply_strints(MedBalTree, 0.5)
ShortBalTree = '(((((A_0:0.5,B_0:0.5):0.5,C_0:1.0):0.5,D_0:1.5):0.5,(((E_0:0.5,F_0:0.5):0.5,G_0:1.0):0.5,H_0:1.5):0.5):0.5,(((I_0:0.5,J_0:0.5):0.5,K_0:1.0):0.5,L_0:1.5):1.0);'
#LongBalTree = smd.multiply_strints(MedBalTree, 2.0)
MedBalTree = '(((((A_0:1.0,B_0:1.0):1.0,C_0:2.0):1.0,D_0:3.0):1.0,(((E_0:1.0,F_0:1.0):1.0,G_0:2.0):1.0,H_0:3.0):1.0):1.0,(((I_0:1.0,J_0:1.0):1.0,K_0:2.0):1.0,L_0:3.0):2.0);'
ShortUnbTree = smd.multiply_strints(MedUnbTree, 0.5)
LongBalTree = '(((((A_0:2.0,B_0:2.0):2.0,C_0:4.0):2.0,D_0:6.0):2.0,(((E_0:2.0,F_0:2.0):2.0,G_0:4.0):2.0,H_0:6.0):2.0):2.0,(((I_0:2.0,J_0:2.0):2.0,K_0:4.0):2.0,L_0:6.0):4.0);'
LongUnbTree = smd.multiply_strints(MedUnbTree, 2.0)

treeshapes = smd.Vividict()

treeshapes['sb'] = ShortBalTree
treeshapes['mb'] = MedBalTree
treeshapes['lb'] = LongBalTree
treeshapes['su'] = ShortUnbTree
treeshapes['mu'] = MedUnbTree
treeshapes['lu'] = LongUnbTree
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# initializer for SyncManager
def mgr_init():
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    print 'Initialized manager'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == '__main__':

    processes = []

    # now using SyncManager vs a Manager
    manager = SyncManager()
    # explicitly starting the manager, and telling it to ignore the interrupt signal
    manager.start(mgr_init)
    try:
        shared_array = manager.list()

        for i, sim in enumerate(phylips):
            p = Process(target=sma.astral_sim, args=(sim, i, shared_array))
            p.start()
            processes.append(p)

        try:
            for process in processes:
                process.join()
        except KeyboardInterrupt:
            print "Keyboard interrupt in main"

        for item in shared_array:
            # we still have access to it!  Yay!
            print item
    finally:
      # to be safe -- explicitly shutting down the manager
        manager.shutdown()
