
import os, errno
from pathlib import Path
import numpy as np
import shutil
from shutil import copy2

def meshGen2D2S (distRange):

    #distRange = np.arange(d0,df+0.001,spacing)
    curDir = os.getcwd() #in the two-region folder
    dataFilesPath = os.path.join(curDir,'mesh-2D2S')
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
        dirPath = os.path.join(dataFilesPath,"%4.2f" % distance)
        os.makedirs(dirPath)
        os.chdir(dirPath)
        mol1Diel = "mol1_diel_%3.1f.xyzr" % (distance)
        mol1Stern = "mol1_stern_%3.1f.xyzr" % (distance)
        mol2Diel = "mol2_diel_%3.1f.xyzr" % (distance)
        mol2Stern = "mol2_stern_%3.1f.xyzr" % (distance)
        
        # Generating mesh input files
        
        # Cl mesh
        f = open(mol1Diel,"w+")   
        mol1_D_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (-distance/2,0,0,1.260) # Na-Na or Na-Na(-)
        #Cl_D_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (-distance/2,0,0,2.312) # Cl-Cl or Na-Cl
        f.write(mol1_D_xyzr)
        f.close()

        # Na mesh
        f = open(mol2Diel,"w+")   
        mol2_D_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (distance/2,0,0,1.260) # Na-Cl or Na-Na or Na-Na(-)
        #Na_D_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (distance/2,0,0,2.312) # Cl-Cl
        f.write(mol2_D_xyzr)
        f.close()

        # Cl stern mesh
        f = open(mol1Stern,"w+")   
        mol1_S_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (-distance/2,0,0,3.260) # Na-Na or Na-Na(-)
        #Cl_S_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (-distance/2,0,0,4.312) # Cl-Cl or Na-Cl
        f.write(mol1_S_xyzr)
        f.close()

        # Na stern mesh
        f = open(mol2Stern,"w+")   
        mol2_S_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (distance/2,0,0,3.260) # Na-Cl or Na-Na
        #Na_S_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (distance/2,0,0,4.312) # Cl-Cl 
        f.write(mol2_S_xyzr)
        f.close()

    
        # generating vert and face files using msms
        # change the line below to the absolute path of msms in your system to generate meshes

        msmsPath = '/Users/Ali/Downloads/msms_MacOSX_2.6.1/msms.MacOSX.2.6.1'
        mol1_D_msms = "%s%s%s%s" % (msmsPath,' -if ',mol1Diel,' -de 8.0 -prob 1.4 -of mol1_diel')
        mol2_D_msms = "%s%s%s%s" % (msmsPath,' -if ',mol2Diel,' -de 8.0 -prob 1.4 -of mol2_diel')
        mol1_S_msms = "%s%s%s%s" % (msmsPath,' -if ',mol1Stern,' -de 2.0 -prob 1.4 -of mol1_stern')
        mol2_S_msms = "%s%s%s%s" % (msmsPath,' -if ',mol2Stern,' -de 2.0 -prob 1.4 -of mol2_stern')
        os.system(mol1_D_msms)
        os.system(mol2_D_msms)
        os.system(mol1_S_msms)
        os.system(mol2_S_msms)

        
        f=open("mol1.srf","w+")
        f.write("f\nf\n./%s\n./%s\n\n\n0\n" % ('mol1_stern','mol1_diel'))
        f.close()

        f=open("mol2.srf","w+")
        f.write("f\nf\n./%s\n./%s\n\n\n0\n" % ('mol2_stern','mol2_diel'))
        f.close()
        
        
        Path('%s/%s' % (dirPath,'mol1_diel')).touch()
        Path('%s/%s' % (dirPath,'mol1_Stern')).touch()
        Path('%s/%s' % (dirPath,'mol2_diel')).touch()
        Path('%s/%s' % (dirPath,'mol2_Stern')).touch()

        
        # writing pqr files
            
        # Charge locations
        x_mol1 = "-%4.2f" % (distance/2)
        x_mol2 = "%4.2f" % (distance/2)

        mol1_pqrName = "mol1.pqr"
        mol2_pqrName = "mol2.pqr"

        f = open(mol1_pqrName,"w+")
        f.write("%s %5d %4s %5s %5d %s %9.6f %9.6f %9.6f %9.6f" % ('ATOM',1,'Cl','TMP',1,x_mol1,0,0,-1,1.0))
        #f.write("%s %5d %4s %5s %5d %s %9.6f %9.6f %9.6f %9.6f" % ('ATOM',1,'NA','TMP',1,x_Cl,0,0,-1,1.2979))
        f.close()

        f = open(mol2_pqrName,"w+")
        f.write("%s %5d %4s %5s %5d %s %9.6f %9.6f %9.6f %9.6f" % ('ATOM',1,'NA','TMP',1,x_mol2,0,0,-1,1.0))
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
