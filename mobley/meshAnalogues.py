#!/usr/bin/env python

import glob, re, os, shutil,sys


def parse_mesh_out(filename):
    fh = open(filename,'r')
    for line in fh:
        if line.find('probe_radius 1.4') >= 0:
            line_data = line.rstrip().lstrip().split()
    outfile = line_data[-1]
    fh.close()

    fh = open(outfile,'r')
    in_SAS_section = False
    while not in_SAS_section:
        line = fh.readline()
        if line.find('ANALYTICAL SURFACE AREA')>=0:
            in_SAS_section = True
    next_line = fh.readline()
    SAS_data = fh.readline().rstrip().lstrip().split()
    SAS = SAS_data[-1]
    fh.close()
    return SAS
    
def parse_matlab_out(filename):
    out1= 0.0
    fh = open(filename,'r')
    for line in fh:
        if line.find('nlbc out') >= 0:
            print "line = " + line
            line_data = line.lstrip().rstrip().split()
            out1 = line_data[2]
#            out2 = line_data[3]
    fh.close()
    return out1
    

pqr_script = "/Users/jbardhan/repos/molman/mobley_to_pqr.py"
generate_pqr_xyzr_files = "%s --pqr test.pqr --prmtop test.prmtop --crd test.crd --xyzr test.xyzr"%pqr_script

meshmaker_bin = "/Users/jbardhan/repos/fftsvd/meshmaker"
mesh_1_output_file = "mesh_1.out"
mesh_density_1 = "%s test.xyzr test.xyzr test_1.srf 1.4 2.0 1 1 1 1 . > %s "%(meshmaker_bin,mesh_1_output_file)
mesh_density_2 = "%s test.xyzr test.xyzr test_2.srf 1.4 2.0 2 2 1 1 . > mesh_2.out "%meshmaker_bin
mesh_density_4 = "%s test.xyzr test.xyzr test_4.srf 1.4 2.0 4 4 1 1 . > mesh_4.out "%meshmaker_bin
mesh_density_8 = "%s test.xyzr test.xyzr test_8.srf 1.4 2.0 8 8 1 1 . > mesh_8.out "%meshmaker_bin
mesh_density_16 = "%s test.xyzr test.xyzr test_16.srf 1.4 2.0 16 16 1 1 . > mesh_16.out "%meshmaker_bin

matlab_bin = "/Applications/MATLAB_R2011b.app/bin/matlab"
matlab_output_file = "nlbc_out"
run_matlab = "%s -nodesktop -nojvm -r runtest > %s "%(matlab_bin,matlab_output_file)

root_dir = os.getcwd()
file_list = glob.glob("*/*.crd")
file_list = [re.sub("\.crd","",crd_file) for crd_file in file_list]
file_list = [re.sub(".*/","",crd_file) for crd_file in file_list]

counter = 0
log_lines = []
energy_lines = []
for problem in file_list:
    print "problem is " + problem
    
    try:
        os.chdir(problem)
    except:
        print "Could not change directory??"
        
    shutil.copyfile("%s.crd"%problem,"test.crd")
    shutil.copyfile("%s.prmtop"%problem,"test.prmtop")
    shutil.copyfile("../runtest.m","./runtest.m")

    os.system(generate_pqr_xyzr_files)
    os.system(mesh_density_1)
    os.system(mesh_density_2)
    os.system(mesh_density_4)


    os.chdir(root_dir)
    


