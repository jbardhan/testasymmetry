
import os, errno
from pathlib import Path
import numpy as np
import shutil
from shutil import copy2

def twoRegMeshGen (d0,spacing,df):

    distRange = np.arange(d0,df+0.001,spacing)
    curDir = os.getcwd() #in the two-region folder
    dataFilesPath = os.path.join(curDir,'mesh-data')
    try:
        os.makedirs(dataFilesPath)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    for root, dirs, files in os.walk(dataFilesPath):
        for f in files:
            os.unlink(os.path.join(root, f))
        for d in dirs:
            shutil.rmtree(os.path.join(root, d))
    for distance in distRange:
        dirPath = os.path.join(dataFilesPath,"%3.1f" % distance)
        os.makedirs(dirPath)
        os.chdir(dirPath)
        ClDiel = "Cl_diel_%3.1f.xyzr" % (distance)
        ClStern = "Cl_stern_%3.1f.xyzr" % (distance)
        NaDiel = "Na_diel_%3.1f.xyzr" % (distance)
        NaStern = "Na_stern_%3.1f.xyzr" % (distance)
        
        # Generating mesh input files
        
        # Cl mesh
        f = open(ClDiel,"w+")   
        Cl_D_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (-distance/2,0,0,1.260) # Na-Na or Na-Na(-)
        #Cl_D_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (-distance/2,0,0,2.312) # Cl-Cl or Na-Cl
        f.write(Cl_D_xyzr)
        f.close()

        # Na mesh
        f = open(NaDiel,"w+")   
        Na_D_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (distance/2,0,0,1.260) # Na-Cl or Na-Na or Na-Na(-)
        #Na_D_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (distance/2,0,0,2.312) # Cl-Cl
        f.write(Na_D_xyzr)
        f.close()

        # Cl stern mesh
        f = open(ClStern,"w+")   
        Cl_S_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (-distance/2,0,0,3.260) # Na-Na or Na-Na(-)
        #Cl_S_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (-distance/2,0,0,4.312) # Cl-Cl or Na-Cl
        f.write(Cl_S_xyzr)
        f.close()

        # Na stern mesh
        f = open(NaStern,"w+")   
        Na_S_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (distance/2,0,0,3.260) # Na-Cl or Na-Na
        #Na_S_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (distance/2,0,0,4.312) # Cl-Cl 
        f.write(Na_S_xyzr)
        f.close()

    
        # generating vert and face files using msms
        # change the line below to the absolute path of msms in your system to generate meshes

        msmsPath = '/Users/Ali/Downloads/msms_MacOSX_2.6.1/msms.MacOSX.2.6.1'
        Cl_D_msms = "%s%s%s%s" % (msmsPath,' -if ',ClDiel,' -de 8.0 -prob 1.4 -of Cl')
        Na_D_msms = "%s%s%s%s" % (msmsPath,' -if ',NaDiel,' -de 8.0 -prob 1.4 -of Na')
        Cl_S_msms = "%s%s%s%s" % (msmsPath,' -if ',ClStern,' -de 2.0 -prob 1.4 -of Cl_stern')
        Na_S_msms = "%s%s%s%s" % (msmsPath,' -if ',NaStern,' -de 2.0 -prob 1.4 -of Na_stern')
        os.system(Cl_D_msms)
        os.system(Na_D_msms)
        os.system(Cl_S_msms)
        os.system(Na_S_msms)

        
        f=open("Cl.srf","w+")
        f.write("f\nf\n./%s\n./%s\n\n\n0\n" % ('Cl_stern','Cl'))
        f.close()

        f=open("Na.srf","w+")
        f.write("f\nf\n./%s\n./%s\n\n\n0\n" % ('Na_stern','Na'))
        f.close()
        
        
        Path('%s/%s' % (dirPath,'Cl')).touch()
        Path('%s/%s' % (dirPath,'Cl_Stern')).touch()
        Path('%s/%s' % (dirPath,'Na')).touch()
        Path('%s/%s' % (dirPath,'Na_Stern')).touch()

        
        # writing pqr files
            
        # Charge locations
        x_Cl = "-%4.2f" % (distance/2)
        x_Na = "%4.2f" % (distance/2)

        Cl_pqrName = "Cl.pqr"
        Na_pqrName = "Na.pqr"

        f = open(Cl_pqrName,"w+")
        f.write("%s %5d %4s %5s %5d %s %9.6f %9.6f %9.6f %9.6f" % ('ATOM',1,'Cl','TMP',1,x_Cl,0,0,0,1.0))
        #f.write("%s %5d %4s %5s %5d %s %9.6f %9.6f %9.6f %9.6f" % ('ATOM',1,'NA','TMP',1,x_Cl,0,0,-1,1.2979))
        f.close()

        f = open(Na_pqrName,"w+")
        f.write("%s %5d %4s %5s %5d %s %9.6f %9.6f %9.6f %9.6f" % ('ATOM',1,'NA','TMP',1,x_Na,0,0,1,1.0))
        f.close()

        # Removing the first 3 lines from msms output files
        
        listDir = os.listdir(dirPath)
        threeFirstLines = []      
        for item in listDir:

            if (item.endswith(".vert") or item.endswith(".face")):
                with open(item) as f, open("temp.txt", "w") as out:
                    for x in range(3):
                        threeFirstLines.append(next(f))
                    for line in f:
                        out.write(line)
                os.remove(os.path.join(dirPath,item))
                os.rename("temp.txt",item)
        
        os.chdir(dataFilesPath)
