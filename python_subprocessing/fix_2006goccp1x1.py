# Why this didn't work in the JupyterHub still confuses me...

import sys
# Add common resources folder to path
sys.path.append('/glade/u/home/jonahshaw/Scripts/git_repos/CESM2_analysis')
sys.path.append('/glade/u/home/jonahshaw/Scripts/git_repos/CESM2_analysis/Common/')

from imports import (os,run)

og_dir = '/glade/u/home/jonahshaw/w/obs/CALIPSO/GOCCP/2Ddata/grid_1x1_2006/'
out_dir = '/glade/u/home/jonahshaw/w/obs/CALIPSO/GOCCP/2Ddata/grid_1x1_2006fix/'

commands = [f'module load nco']
for file in os.listdir(og_dir): # iterate over years (2009-2013)
    _prepath = '%s%s' % (og_dir, file)
    _outpath = '%s%s' % (out_dir, file)
    
    comm = f'ncks -C -x -v time {_prepath} {_outpath}'
    commands.append(comm)
    
for com in commands:
    print(com)
    out = run(com, cwd='/glade/u/home/jonahshaw/w/obs/', shell=True)
    print(out)