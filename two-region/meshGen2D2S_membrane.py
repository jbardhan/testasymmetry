
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

#def memMeshGen2D2S (d0,spacing,df,r,h):
def memMeshGen2D2S (distRange,r_ion,r_mem,h,particleCharge):

    #distRange = np.arange(d0,df+0.001,spacing)
    curDir = os.getcwd() #in the two-region folder
    meshFilesPath = os.path.join(curDir,'I-M-2D2S')
    try:
        os.makedirs(meshFilesPath)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    for root, dirs, files in os.walk(meshFilesPath):
        for f in files:
            os.unlink(os.path.join(root, f))
        for d in dirs:
            shutil.rmtree(os.path.join(root, d))

    os.chdir(meshFilesPath)
    #Membrane setup
    dielgeo = "diel%d.geo" % (r_mem)
    sterngeo = "stern%d.geo" % (r_mem)
    dielstl = "diel%d.stl" % (r_mem)
    sternstl = "stern%d.stl" % (r_mem)
    listDir = os.listdir(meshFilesPath)
    
    copy2(os.path.join(curDir,'template.geo'),meshFilesPath)
    # Generating *.geo file
    lc = 0.8 * np.sqrt(r_mem)
    lc_fine = 0.001 * r_mem
    lc_st = 0.9 * np.sqrt(r_mem+2)
    lc_fine_st = 0.001 * (r_mem+2)
    r_diel = r_mem;
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
    meshDiel, info = pymesh.collapse_short_edges(meshDiel, abs_threshold=0.2 , preserve_feature=False)
    elementsDiel = np.copy(meshDiel.elements)
    for i in range(elementsDiel.shape[0]):
        for j in range(elementsDiel.shape[1]):
            elementsDiel[i,j]+=1
    #Stern surface generation        
    meshStern = pymesh.load_mesh(sternstl)

    meshStern, info = pymesh.remove_duplicated_vertices(meshStern)
    meshStern, info = pymesh.remove_duplicated_faces(meshStern)
    meshStern, info = pymesh.collapse_short_edges(meshStern, abs_threshold=0.2 , preserve_feature=False)
    elementsStern = np.copy(meshStern.elements)
    #Because element numbers in PyMesh start from zero!
    for i in range(elementsStern.shape[0]):
        for j in range(elementsStern.shape[1]):
            elementsStern[i,j]+=1
        
    mem_Diel = "diel%d" % r_mem
    mem_Stern = "stern%d" % r_mem
    np.savetxt("%s.vert" % mem_Diel,meshDiel.vertices,fmt='%5.3f',delimiter=' ')
    np.savetxt("%s.face" % mem_Diel,elementsDiel,fmt='%5.0f',delimiter=' ')
    np.savetxt("%s.vert" % mem_Stern,meshStern.vertices,fmt='%5.3f',delimiter=' ')
    np.savetxt("%s.face" % mem_Stern,elementsStern,fmt='%5.0f',delimiter=' ')
    
    # Creating mem.srf and two empty files which run 
    f=open("mem.srf","w+")
    f.write("f\nf\n./%s\n./%s\n\n\n0\n" % (mem_Stern,mem_Diel))
    f.close()
    
    
    mem_Diel = "diel%d" % r_mem
    mem_Stern = "stern%d" % r_mem
    Path('%s/%s' % (meshFilesPath,mem_Diel)).touch()
    Path('%s/%s' % (meshFilesPath,mem_Stern)).touch()
    
    

    memMeshInfo = "%s%d\n%s%d\n%s%d\n%s%d\n%s%d%s" % ('# Elements in dielectric mesh = ',meshDiel.num_elements,
                                           '# Vertices in dielectric mesh = ',meshDiel.num_vertices,
                                           '# Elements in stern layer mesh = ',meshStern.num_elements,
                                           '# Vertices in stern layer mesh = ',meshStern.num_vertices,
                                           'Surface area = ',np.pi*(2 * r_mem ** 2 + 2 * r_mem * h),' AA')
    f = open('memMesh.info',"w+")
    f.write(memMeshInfo)
    f.close()
    os.remove(os.path.join(meshFilesPath,'template.geo'))

    fullMembraneChargeDist = 'po4_full.pqr'
    
    
    f = open(os.path.join(curDir,fullMembraneChargeDist))
    lines = f.readlines()
    f.close()
    list =[];
    for id in range(0,len(lines)):
        list.append(lines[id].split())
    df=pd.DataFrame(list,columns=["Type", "Atom_number", "Atom_type", "AA type", "Chain", "X", "Y", "Z", "q", "r"])
    df[["Atom_number", "Chain"]] = df[["Atom_number", "Chain"]].apply(pd.to_numeric)
    df[["X", "Y", "Z", "q", "r"]] = df[["X", "Y", "Z", "q", "r"]].apply(pd.to_numeric)
    df = df[((df.X)**2 + (df.Y)**2   < (r_mem-5)**2)].sort_values("Chain")
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
    np.savetxt('halfMempo4.pqr',df[0:int(len(df)/2)],fmt='%s')
    np.savetxt('fullMempo4.pqr', df.values, fmt='%s')

    count = 1
    
    for distance in distRange:

        molDiel = "mol_diel_%i.xyzr" % count
        molStern = "mol_stern_%i.xyzr" % count
        
        # Generating mesh input files
        
        # Cl mesh
        z_mol = np.float(h+2+np.float(distance)+2+r_ion)
        f = open(molDiel,"w+")
        mol_D_xyzr = "%10.5f%10.5f%10.5f%10.5f\n" % (0,0,z_mol,r_ion) # Na-Na or Na-Na(-)
        #Cl_D_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (-distance/2,0,0,2.312) # Cl-Cl or Na-Cl
        f.write(mol_D_xyzr)
        f.close()

        # Cl stern mesh
        f = open(molStern,"w+")   
        mol_S_xyzr = "%10.5f%10.5f%10.5f%10.5f\n" % (0,0,z_mol,r_ion+2) # Na-Na or Na-Na(-)
        #Cl_S_xyzr = "%8.5f%8.5f%8.5f%8.5f\n" % (-distance/2,0,0,4.312) # Cl-Cl or Na-Cl
        f.write(mol_S_xyzr)
        f.close()

        
        # generating vert and face files using msms
        # change the line below to the absolute path of msms in your system to generate meshes

        msmsPath = '/Users/Ali/Downloads/msms_MacOSX_2.6.1/msms.MacOSX.2.6.1'
        dielFaceVert = "mol_diel_%i" % count
        sternFaceVert = "mol_stern_%i" % count
        
        mol_D_msms = "%s%s%s%s%s" % (msmsPath,' -if ',molDiel,' -no_header -de 8.0 -prob 1.4 -of ',dielFaceVert)
        mol_S_msms = "%s%s%s%s%s" % (msmsPath,' -if ',molStern,' -no_header -de 2.0 -prob 1.4 -of ',sternFaceVert)
        os.system(mol_D_msms)
        os.system(mol_S_msms)
        
        
        molSrf = "molSrf_%i.xyzr" % count
        f=open(molSrf,"w+")
        f.write("f\nf\n./%s\n./%s\n\n\n0\n" % (sternFaceVert,dielFaceVert))
        f.close()

        
        Path('%s/%s' % (meshFilesPath,dielFaceVert)).touch()
        Path('%s/%s' % (meshFilesPath,sternFaceVert)).touch()
        
        # writing pqr files
            
        # Charge locations
        
        
        molPqr = "molPqr_%i.pqr" % count
        
        f = open(molPqr,"w+")
        z_mol = "%4.2f" % (h+2+distance+2+r_ion)
        f.write("%s %5d %4s %5s %5d %9.6f %9.6f %s %9.6f %9.6f" % ('ATOM',1,'NA','TMP',1,0,0,z_mol,particleCharge,1.0))
        #f.write("%s %5d %4s %5s %5d %s %9.6f %9.6f %9.6f %9.6f" % ('ATOM',1,'NA','TMP',1,x_Cl,0,0,-1,1.2979))
        f.close()

        # Removing the first 3 lines from msms output files
        count += 1
        '''listDir = os.listdir(meshFilesPath)
        threeFirstLines = []      
        for item in listDir:

            if ((item.endswith(".vert") and item.startswith("mol")) or (item.endswith(".face") and item.startswith("mol"))):
                with open(item) as f, open("temp.txt", "w") as out:
                    for x in range(3):
                        threeFirstLines.append(next(f))
                    for line in f:
                        out.write(line)
                os.remove(os.path.join(meshFilesPath,item))
                os.rename("temp.txt",item)'''

    listDir = os.listdir(meshFilesPath)
    for item in listDir:
        if item.endswith(".xyzr"):
            os.remove(os.path.join(meshFilesPath,item))
    os.chdir(curDir)
