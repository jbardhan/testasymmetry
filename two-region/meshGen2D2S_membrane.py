
import os, errno
from pathlib import Path
import pymesh
import numpy as np
import pandas as pd
import shutil
from shutil import copy2

def insert(originalfile, string, newname):
    with open(originalfile,'r') as f:
        with open('newfile.txt','w') as f2: 
            f2.write(string)
            f2.write(f.read())
    os.rename('newfile.txt',newname)

def memMeshGen2D2S (d0,spacing,df,r,h):

    distRange = np.arange(d0,df+0.001,spacing)
    curDir = os.getcwd() #in the two-region folder
    dataFilesPath = os.path.join(curDir,'I-M-2D2S')
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
        molDiel = "mol_diel_%3.1f.xyzr" % (distance)
        molStern = "mol_stern_%3.1f.xyzr" % (distance)
        
        # Generating mesh input files
        
        # Cl mesh
        f = open(molDiel,"w+")   
        mol_D_xyzr = "%10.5f%10.5f%10.5f%10.5f\n" % (0,0,h+2+distance,1.260) # Na-Na or Na-Na(-)
        #Cl_D_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (-distance/2,0,0,2.312) # Cl-Cl or Na-Cl
        f.write(mol_D_xyzr)
        f.close()

        # Cl stern mesh
        f = open(molStern,"w+")   
        mol_S_xyzr = "%10.5f%10.5f%10.5f%10.5f\n" % (0,0,h+2+distance,3.260) # Na-Na or Na-Na(-)
        #Cl_S_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (-distance/2,0,0,4.312) # Cl-Cl or Na-Cl
        f.write(mol_S_xyzr)
        f.close()

        
        # generating vert and face files using msms
        # change the line below to the absolute path of msms in your system to generate meshes

        msmsPath = '/Users/Ali/Downloads/msms_MacOSX_2.6.1/msms.MacOSX.2.6.1'
        mol_D_msms = "%s%s%s%s" % (msmsPath,' -if ',molDiel,' -de 8.0 -prob 1.4 -of mol_diel')
        mol_S_msms = "%s%s%s%s" % (msmsPath,' -if ',molStern,' -de 2.0 -prob 1.4 -of mol_stern')
        os.system(mol_D_msms)
        os.system(mol_S_msms)
        
        
        f=open("mol.srf","w+")
        f.write("f\nf\n./%s\n./%s\n\n\n0\n" % ('mol_stern','mol_diel'))
        f.close()

        
        
        Path('%s/%s' % (dirPath,'mol_diel')).touch()
        Path('%s/%s' % (dirPath,'mol_Stern')).touch()
        
        
        # writing pqr files
            
        # Charge locations
        z_mol = "%4.2f" % (h+2+distance)
        
        mol_pqrName = "mol.pqr"
        
        f = open(mol_pqrName,"w+")
        f.write("%s %5d %4s %5s %5d %9.6f %9.6f %s %9.6f %9.6f" % ('ATOM',1,'NA','TMP',1,0,0,z_mol,1,1.0))
        #f.write("%s %5d %4s %5s %5d %s %9.6f %9.6f %9.6f %9.6f" % ('ATOM',1,'NA','TMP',1,x_Cl,0,0,-1,1.2979))
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

        #Membrane setup
        dielgeo = "diel%d.geo" % (r)
        sterngeo = "stern%d.geo" % (r)
        dielstl = "diel%d.stl" % (r)
        sternstl = "stern%d.stl" % (r)
        listDir = os.listdir(dirPath)

        # Removing old files        
        '''for item in listDir:

            if (item.endswith(".stl") or item.endswith(".geo") or
                item.endswith(".vert") or item.endswith(".face") or
                item.endswith(".msh") or item.endswith(".pqr") or
                os.path.getsize(os.path.join(curDir, item)) == 0 or item.endswith(".info")):
                os.remove(os.path.join(dirPath, item))
                '''

        copy2(os.path.join(curDir,'template.geo'),dirPath)

        # Generating *.geo file
        lc = 0.8 * np.sqrt(r)
        lc_fine = 0.001 * r
        lc_st = 0.9 * np.sqrt(r+2)
        lc_fine_st = 0.001 * (r+2)
        r_diel = r;
        r_stern = r_diel + 2

        # diel.geo
        f = open(dielgeo,"w+")   
        meshCustomizeDiel = "%s%1d%s\n%s%1d%s\n%s%5.3f%s\n%s%5.3f%s\n%s\n" % ('r = ', r_diel, ';', 'h = ', h, ';', 'lc =', lc, ';', 'lc_fine =', lc_fine, ';', 'st_thickness = 0;')
        f.write(meshCustomizeDiel)
        f.close()

        # stern.geo
        f = open(sterngeo,"w+")   
        meshCustomizeStern = "%s%1d%s\n%s%1d%s\n%s%5.3f%s\n%s%5.3f%s\n%s\n" % ('r = ', r_diel + 2, ';', 'h = ', h, ';', 'lc =', lc_st, ';', 'lc_fine =', lc_fine_st, ';', 'st_thickness = 2;')
        f.write(meshCustomizeStern)
        f.close()

        insert('template.geo',meshCustomizeDiel,dielgeo)
        insert('template.geo',meshCustomizeStern,sterngeo)

        # converting *.geo files to *.stl using Gmsh
        gmsh = '/Applications/Gmsh.app/Contents/MacOS/gmsh ' #MacOS
        #gmsh = 'gmsh ' #Linux
        dielMeshCommand = "%s%s%s%s" % (gmsh,dielgeo,' -2 -format stl ',dielstl)
        sternMeshCommand = "%s%s%s%s" % (gmsh,sterngeo,' -2 -format stl ',sternstl)

        os.system(dielMeshCommand)
        os.system(sternMeshCommand)


        # Each folder corresponds to one radius and has geometric information and
        # other information such as charge distribution for each charge location (*.pqr).
        # Geomtry file (*.geo) that is used to enerate each geometry is also included.
        # These folders are setup in a way that it needs minimal change in order to be read
        # by SLIC Matl                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          km,ab code
        
        
        #Dielectric surface generation
        meshDiel = pymesh.load_mesh(dielstl)
        meshDiel, info = pymesh.remove_duplicated_vertices(meshDiel)
        meshDiel, info = pymesh.remove_duplicated_faces(meshDiel)
        elementsDiel = np.copy(meshDiel.elements)
        for i in range(elementsDiel.shape[0]):
            for j in range(elementsDiel.shape[1]):
                elementsDiel[i,j]+=1

        #Stern surface generation        
        meshStern = pymesh.load_mesh(sternstl)
        meshStern, info = pymesh.remove_duplicated_vertices(meshStern)
        meshStern, info = pymesh.remove_duplicated_faces(meshStern)
        elementsStern = np.copy(meshStern.elements)

        #Because element numbers in PyMesh start from zero!
        for i in range(elementsStern.shape[0]):
            for j in range(elementsStern.shape[1]):
                elementsStern[i,j]+=1
            
        dielFileName = "diel%d" % r
        sternFileName = "stern%d" % r

        np.savetxt("%s.vert" % dielFileName,meshDiel.vertices,fmt='%5.3f',delimiter=' ')
        np.savetxt("%s.face" % dielFileName,elementsDiel,fmt='%5.0f',delimiter=' ')
        np.savetxt("%s.vert" % sternFileName,meshStern.vertices,fmt='%5.3f',delimiter=' ')
        np.savetxt("%s.face" % sternFileName,elementsStern,fmt='%5.0f',delimiter=' ')
        
        # Creating mem.srf and two empty files which run 
        f=open("mem.srf","w+")
        f.write("f\nf\n./%s\n./%s\n\n\n0\n" % (sternFileName,dielFileName))
        f.close()
        
        
        mem_Diel = "diel%d" % r
        mem_Stern = "stern%d" % r
        Path('%s/%s' % (dirPath,mem_Diel)).touch()
        Path('%s/%s' % (dirPath,mem_Stern)).touch()

        
        # writing pqr files in each folder
        lambda_z = [25.0,30.0,35.0,40.0,42.5,45.0,46.0,47.0,47.5,48.0,48.5]
        k=0
        #for lambdaFEP in lambda_z:

        f = open(os.path.join(curDir,'po4_new.pqr'))
        lines = f.readlines()
        f.close()
        list =[];
        for id in range(0,len(lines)):
            list.append(lines[id].split())

        df=pd.DataFrame(list,columns=["Type", "Atom_number", "Atom_type", "AA type", "Chain", "X", "Y", "Z", "q", "r"])
        df[["Atom_number", "Chain"]] = df[["Atom_number", "Chain"]].apply(pd.to_numeric)
        df[["X", "Y", "Z", "q", "r"]] = df[["X", "Y", "Z", "q", "r"]].apply(pd.to_numeric)
        df = df[((df.X)**2 + (df.Y)**2   < (r-5)**2)].sort_values("Chain")
        df = pd.concat(g for _, g in df.groupby("Chain") if len(g) > 24)
        df = df.sort_values("Atom_number")
        df['Type'] = df['Type'].map('{0:5s}'.format)
        df['Atom_number'] = df['Atom_number'].map('{0:6d}'.format)
        df['Atom_type'] = df['Atom_type'].map('{0:6s}'.format)
        df['AA type'] = df['AA type'].map('{0:5s}'.format)
        df['Chain'] = df['Chain'].map('{0:5d}'.format)
        df['X'] = df['X'].map('{0:14.6f}'.format)
        df['Y'] = df['Y'].map('{0:12.6f}'.format)
        df['Z'] = df['Z'].map('{0:12.6f}'.format)
        df['q'] = df['q'].map('{0:11.6f}'.format)
        df['r'] = df['r'].map('{0:11.6f}'.format)
        np.savetxt('mem_w_po4.pqr', df.values, fmt='%s')

        for lambdaFEP in range(1,2):
            
            # Charge locations
            x_q = 0;
            y_q = 0;
            
            # Charge moves from the center of the membrane toward the top plane on z-axis
            # z_q = h/2 * (1 + lambdaFEP)
            
            q = 0
            R_q = 1 # Ask Jay!
            #pqrFileName = "test_lambda_%3.1f.pqr" % k
            pqrFileName = "mem.pqr" 
            f = open(pqrFileName,"w+")
            f.write("%s %5d %4s %5s %5d %15.6f %9.6f %9.6f %9.6f %9.6f" % ('ATOM',1,'NA','TMP',1,x_q,y_q,lambdaFEP,q,R_q))
            f.close()
            #k += 1
        memMeshInfo = "%s%d\n%s%d\n%s%d\n%s%d\n%s%d%s" % ('# Elements in dielectric mesh = ',meshDiel.num_elements,
                                               '# Vertices in dielectric mesh = ',meshDiel.num_vertices,                                               '# Elements in stern layer mesh = ',meshStern.num_elements,
                                               '# Vertices in stern layer mesh = ',meshStern.num_vertices,
                                               'Surface area = ',np.pi*(2 * r ** 2 + 2 * r * h),' AA')
        f = open('memMesh.info',"w+")
        f.write(memMeshInfo)
        f.close()
        os.remove(os.path.join(dirPath,'template.geo'))
        








        
        os.chdir(dataFilesPath)
