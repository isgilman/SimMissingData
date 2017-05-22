#!/usr/bin/env python

# ## SMD2ASTRAL
# This script will 
# 1. read in full nexus files from the missing data runs
# 2. simulate missing data with functions from SimMissingData.py
# 3. partition in missing data sets into individual loci with G2Genes
# 4. construct gene trees with paup
# 5. combine gene trees into a single file for input into ASTRAL
# 6. run ASTRAL
# 7. record RF distance to true tree with dendropy

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
try: from tqdm import tqdm, tnrange
except ImportError:
    print "\n\tError: tqdm is not installed/loaded."
    print "\n\tTo install, try 'pip install tqdm' or 'conda install -c conda-forge tqdm'"
    sys.exit()
try: import subprocess
except ImportError:
    print "\n\tError: subprocess is not installed/loaded."
    sys.exit()
try: import numpy as np
except ImportError:
    print "\n\tError: numpy is not installed/loaded."
    sys.exit()
try: import pandas as pd
except ImportError:
    print "\n\tError: pandas is not installed/loaded."
    sys.exit()
try: from Bio import AlignIO, Alphabet, SeqIO
except ImportError:
    print "\n\tError: Bio is not installed/loaded."
    sys.exit()
try: import re
except ImportError:
    print "\n\tError: re is not installed/loaded."
    sys.exit()
try: import dendropy
except ImportError:
    print "\n\tError: dendropy is not installed/loaded."
    sys.exit()
try: import csv
except ImportError:
    print "\n\tError: csv is not installed/loaded."
    sys.exit()
try: from SMDfx import listabs, Vividict, phylip_to_df
except ImportError:
    print "\n\tError: SMDfx is not installed/loaded."
    sys.exit()
try: from subprocess import Popen, PIPE
except ImportError:
    print "\n\tError: subprocess is not installed/loaded."
    sys.exit()
try: import SMDfx as smd
except ImportError:
    print "\n\tError: SMDfx is not installed/loaded."
try: from multiprocessing.managers  import BaseProxy
except ImportError:
    print "\n\tError: multiprocessing is not installed/loaded."
#try: import skbio
#except ImportError:
#    print "\n\tError: skbio is not installed/loaded."
#    sys.exit()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MedBalTree = '(((((A_0:2,B_0:2):2,C_0:4):4,D_0:8):4,(((E_0:2,F_0:2):2,G_0:4):4,H_0:8):4):4,(((I_0:2,J_0:2):2,K_0:4):4,L_0:8):8);'
MedUnbTree = '((((((((((((A_0:1,B_0:1):1,C_0:2):1,D_0:3):1,E_0:4):1,F_0:5):1,G_0:6):1,H_0:7):1,I_0:8):1,J_0:9):1,K_0:10):1,L_0:11):12);'
#ShortBalTree = smd.multiply_strints(MedBalTree, 0.5)
ShortBalTree = '(((((A_0:0.5,B_0:0.5):0.5,C_0:1.0):0.5,D_0:1.5):0.5,(((E_0:0.5,F_0:0.5):0.5,G_0:1.0):0.5,H_0:1.5):0.5):0.5,(((I_0:0.5,J_0:0.5):0.5,K_0:1.0):0.5,L_0:1.5):1.0);'
#LongBalTree = smd.multiply_strints(MedBalTree, 2.0)
LongBalTree = '(((((A_0:1.0,B_0:1.0):1.0,C_0:2.0):1.0,D_0:3.0):1.0,(((E_0:1.0,F_0:1.0):1.0,G_0:2.0):1.0,H_0:3.0):1.0):1.0,(((I_0:1.0,J_0:1.0):1.0,K_0:2.0):1.0,L_0:3.0):2.0);'
ShortUnbTree = smd.multiply_strints(MedUnbTree, 0.5)
LongUnbTree = smd.multiply_strints(MedUnbTree, 2.0)

treeshapes = smd.Vividict()

treeshapes['sb'] = ShortBalTree
treeshapes['mb'] = MedBalTree
treeshapes['lb'] = LongBalTree
treeshapes['su'] = ShortUnbTree
treeshapes['mu'] = MedUnbTree
treeshapes['lu'] = LongUnbTree

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_simpath(phylip):
    simname = os.path.basename(phylip).split('.')[0]
    metadata = re.split(r'[LSN]', simname)
    metadata[0] = metadata[0].replace('l', '')
    
    return simname, int(metadata[0]), int(metadata[1]), metadata[2], metadata[3]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def write_Pmissing_data(fulldf, percent, simname, partitions):
    '''Simulates missing data given a DataFrame of sequence data and writes it to file.
    
    Parameters
    ----------
    fulldf : pandas DataFrame containing the full data
    percent : percent missing data to simulate
    simname : simulation name'''

    # Get number of entries
    entries = fulldf.shape[0]*fulldf.shape[1]
    
    # Substitute rand percent with 'N'
    if percent==0: missingdf=fulldf
    else: missingdf = fulldf.where(np.random.uniform(size=fulldf.shape) > float(percent/100.), 'N')
    
    # Split dataframe by partitions
    for i, part in enumerate(partitions):
        locus = '.locus'+str(i).zfill(len(str(len(partitions))))
        
        # Create string for naming files
        p = '.p'+str(percent).zfill(2)+locus
        # Delete missing data from previous replicate if it exists
        if os.path.exists(simname+p+'.nex'):
            os.remove(simname+p+'.nex')
        if os.path.exists(simname+p+'.phy'):
            os.remove(simname+p+'.phy')
        
        partdf = missingdf[range(part[0], part[1]+1)]
 
        # Convert dataframe to sequence strings
        missinglines =[]
        meta = ('%d %d\n') % (partdf.shape[0], partdf.shape[1])
        missinglines.append(meta)

        for i in xrange(len(missingdf)):
            missinglines.append('       '.join([partdf.index[i],''.join(partdf.ix[i])])+'\n')

        # Write strings to phylip file
        with open(simname+p+'.phy', 'w+') as m:
            m.writelines(missinglines)
        #Convert phlyip to nexus for PAUP
        alignment = AlignIO.read(open(simname+p+'.phy'), "phylip", alphabet=Alphabet.generic_dna)
        with open(simname+p+'.nex', "w+") as n:
            n.write(alignment.format("nexus"))
        os.remove(simname+p+'.phy')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def write_paup_block(nexus, blockpath='/Users/IanGilman/simrrls/ml_paup_block.txt'):
    '''Appends a PAUP block to nexus files to run SVDquartets
    
    Parameters
    ----------
    nexus : a nexus file
    percent : percent missing data
    blockpath : path to text file containing block to ammend and append'''
    
    prefix = '.'.join(os.path.basename(nexus).split('.')[:-1])
    
    # Read default PAUP block for funning SVDquartets
    with open(blockpath, 'r') as b:
        defaultblock = b.readlines()
        
    # Edit block lines
    newblock = ['\n']
    for line in defaultblock:
        if 'SaveTrees' in line:
            line = '\tSaveTrees file=%s.tre format=newick brLens=yes;\n' % (prefix)
        newblock.append(line)

    # Append block to nexus file
    with open(nexus, 'a') as n:
        n.writelines(newblock)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flatten = lambda l: [item for sublist in l for item in sublist]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def fixmissing(nexusfile):
    withdata = []
    for record in SeqIO.parse(nexusfile, 'nexus'):
        if len(set(record.seq)) > 1:
            withdata.append(record)
    SeqIO.write(withdata, nexusfile, 'nexus')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def sym_diff(t, true_shape, treedict):
    '''This is a sloppy fix version of the symmetric difference function in Dendropy written by Jeet Sukumaran
    (see his GitHub)
    
    Parameters
    ----------
    t1 : first dendropy tree object
    t2 : second dendropy tree object'''
    
    taxon_namespace = dendropy.TaxonSet()
    
    t = dendropy.Tree.get_from_string(t, "newick", taxon_set=taxon_namespace)
    true_tree = dendropy.Tree.get_from_string(treedict[true_shape], "newick", taxon_set=taxon_namespace)
    t.is_rooted
    true_tree.is_rooted
    t.symmetric_difference(true_tree)
    t.symmetric_difference(true_tree)
    
    return true_tree.symmetric_difference(t)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def astral_sim(phy, process_number, shared_array):
    
    try:
        while True:
            
            ## To combat a Broken Pipe issue (see StackOverflow: broken-pipe-when-using-python-multiprocessing-managers-basemanager-syncmanager)
            ## Cause: BaseProxy class has thread local storage which caches the connection, which is reused for future connections causing 
            ## "broken pipe" errors even on creating a new Manager
            ## Workaround: Delete the cached connection before reconnecting
            try: 
                if address in BaseProxy._address_to_local:
                    del BaseProxy._address_to_local[address][0].connection
            except NameError: None

            shared_array.append(process_number)
            print('Started thread: %d') % (process_number)
            
            # Get values from simulation name
            proc = os.getpid()
            simname, size, nloci, shape, ne = parse_simpath(phy)
            reps=10

            # Load true tree
            #truetree = '"'+treeshapes[shape]+'"'
            #truetree = TreeNode.read(StringIO(eval('u'+truetree)))
       
            print('Simulation: %s\nLocus size: %d\nNumber of loci: %d\nShape: %s\nNe: %s\nPID: %d\n') % (simname, size, nloci, shape, ne, proc)
                
            # Change to simulation directory
            print('CWD: %s') % (os.getcwd())
            os.chdir(os.path.dirname(os.path.dirname(phy)))
            print('Changed CWD to: %s') % (os.getcwd())
            
            # Check for completed simulation
            finaloutput = os.path.join(os.path.join(os.path.join(os.getcwd(), 'outfiles'), 'loci'), simname+'.ASTRAL.RF.csv')
            if os.path.exists(finaloutput):
                print('\tFound ASTRAL RF file for %s. Checking for completion...') % (simname)
                with open(finaloutput, 'r') as f:
                    if not len(f.readlines())==reps+1:
                        print('\t%s ASTRAL RF file not complete. Rerunning...') % (simname)
                    else:
                        print('\t%s ASTRAL RF file complete. Exitting...\n') % (simname)
                        #print('\t%s ASTRAL RF file deleted\n') % (simname)
                        #os.remove(finaloutput)
                        sys.exit()     
                    
            # Read in partitions
            partpath = phy+'.partitions'
            with open(partpath, 'r') as p:
                parts = p.readlines()
            
            partitions = []
            for line in parts:
                split = re.split(r'[=-]', line)[-2:]
                split[0] = int(split[0])
                split[1] = int(split[1].split('\n')[0])
                partitions.append(split)
            
            # Convert phylip to dataframe
            df=phylip_to_df(simname)
            
            print('CWD: %s') % (os.getcwd())
            locidir = os.path.join(os.path.join(os.getcwd(), 'outfiles'),'loci')
            try: os.mkdir(locidir)
            except OSError: None
            os.chdir(locidir)
            print('CWD: %s') % (os.getcwd())
            
            # Initialize dict for RF distances
            treedist = {}
            
            # Loop over rates of missing data
            for rate in tqdm(range(0,85, 10), desc=simname, leave=False) :
                # Initialize dict list
                pkey = 'p'+str(rate).zfill(2)
                treedist[pkey] = []
        
                # Loop over loci nexus files
                for i, rep in enumerate(range(reps)):
                    
                    # Create missing data nexus files for each locus
                    write_Pmissing_data(fulldf=df, percent=rate, simname=simname, partitions=partitions)
        
                    # Write PAUP block to nexus files
                    nexfiles = [filename for filename in os.listdir('.') if pkey in filename and filename.endswith('.nex')]
                    newicks = []
        
                    for nexusfile in nexfiles:
                        
                        fixmissing(nexusfile)

                        write_paup_block(nexusfile, blockpath='/mnt/lfs2/dtank/ian/pSims/ml_paup_block.txt')
        
                        # Run PAUP
                        # Create newick path to write to, remove old trees if need be
                        newick = nexusfile.replace('nex', 'tre')
                        if os.path.exists(newick):
                            os.remove(newick)
                        # Run PAUP
                        FNULL = open(os.devnull, 'w')
                        subprocess.call(['/opt/modules/biology/paup/4a152/bin/paup', '-n', nexusfile], stdout=FNULL, stderr=FNULL)
                        # Read newick file to list and delete
                        with open(newick, 'r') as n:
                            newicks.append(n.readlines())
                        os.remove(newick)
                    
                    # Flatten sublists
                    newicks = flatten(newicks)
        
                    genetreesfile = os.path.join(os.getcwd(), '.'.join([simname, pkey, str(i),'tre']))
                    
                    if os.path.exists(genetreesfile): 
                        os.remove(genetreesfile)
                    with open(genetreesfile, 'w+') as t:
                        #print('\tWRITING NEWICKS')
                        t.writelines(flatten(newicks))
        
                    # Run ASTRAL, capture output tree to 'outtree'
                    p = subprocess.Popen(['java','-Xmx1000M', '-jar', '/mnt/lfs2/dtank/ian/pSims/ASTRAL/Astral/astral.4.10.12.jar', '-i', genetreesfile, '-t', '0'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
                    outtree, err = p.communicate(b"input data that is passed to subprocess' stdin")
                    # Calculate symmetric difference
                    try:
                        #print(simname, outtree)
                        
                        taxon_namespace = dendropy.TaxonSet()
                        #print(taxon_namespace)
                        true_tree = dendropy.Tree.get_from_string(treeshapes[shape], "newick", taxon_set=taxon_namespace)
                        #print(taxon_namespace)
                        #true_tree.print_plot()
                        estimate = dendropy.Tree.get_from_string(outtree, "newick", taxon_set=taxon_namespace) 
                        #print(taxon_namespace)
                        #estimate.print_plot()
                        true_tree.is_rooted
                        estimate.is_rooted
                        true_tree.symmetric_difference(estimate)
                        true_tree.symmetric_difference(estimate)
                        rfd = estimate.symmetric_difference(true_tree)

                        #rfd = sym_diff(t=outtree, true_shape=shape, treedict=treeshapes)
                        #estimate = '"'+outtree.replace('\n','')+'"'
                        #estimate = TreeNode.read(StringIO(eval('u'+estimate)))
                        #prf = truetree.compare_subsets(estimate)
                        #print('\t%s: %f\n') % (simname, rfd)
                        #treedist[pkey].append(prf)
                        treedist[pkey].append(rfd)
                    except IndexError:
                        print('\tERROR OCCURED WHILE TAKING SYMMETRIC DIFFERENCE IN %s\n') % (simname)
                        print(outtree)
                        print(err)
                        sys.exit() 
                    except DataParseError:
                        print('\DataParseError in %s. Here is the outtree:\n\t%s\n') % (simname, outtree)
                    os.remove(nexusfile)
            pd.DataFrame.from_dict(treedist).to_csv(path_or_buf=simname+'.ASTRAL.RF.csv')
            for filename in listabs(locidir):
                if filename.endswith('.nex'): os.remove(filename)
    except KeyboardInterrupt:
        print "Keyboard interrupt in process: ", process_number
    finally:
        print "cleaning up thread", process_number
        
