import os
from pathlib import Path
import pymesh
import numpy as np
pwd = os.getcwd()

# R = Cylinder (membrane) radius
for R in range(1,11):
    
    # Each folder corresponds to one radius and has geometric information and
    # other information such as charge distribution for each charge location (*.pqr).
    # Geomtry file (*.geo) that is used to enerate each geometry is also included.
    # These folders are setup in a way that it needs minimal change in order to be read
    # by SLIC Matlab code
    
    Rdir = pwd + "/%d" % R
    os.chdir(Rdir)
    
    # Mesh files generated using Gmsh
    dielMeshFileName = 'test2_%d.stl' % R
    sternMeshFileName = 'test2_%d_stern.stl' % R
    
    #Dielectric surface generation
    meshDiel = pymesh.load_mesh(dielMeshFileName)
    meshDiel, info = pymesh.remove_duplicated_vertices(meshDiel)
    meshDiel, info = pymesh.remove_duplicated_faces(meshDiel)
    elementsDiel = np.copy(meshDiel.elements)
    for i in range(elementsDiel.shape[0]):
        for j in range(elementsDiel.shape[1]):
            elementsDiel[i,j]+=1

    #Stern surface generation        
    meshStern = pymesh.load_mesh(sternMeshFileName)
    meshStern, info = pymesh.remove_duplicated_vertices(meshStern)
    meshStern, info = pymesh.remove_duplicated_faces(meshStern)
    elementsStern = np.copy(meshStern.elements)

    #Because element numbers in PyMesh start from zero!
    for i in range(elementsStern.shape[0]):
        for j in range(elementsStern.shape[1]):
            elementsStern[i,j]+=1
        
    dielFileName = "cyl_%d" % R
    sternFileName = "cyl_%d_stern" % R

    np.savetxt("%s.vert" % dielFileName,meshDiel.vertices,fmt='%5.3f',delimiter=' ')
    np.savetxt("%s.face" % dielFileName,elementsDiel,fmt='%5.0f',delimiter=' ')
    np.savetxt("%s.vert" % sternFileName,meshStern.vertices,fmt='%5.3f',delimiter=' ')
    np.savetxt("%s.face" % sternFileName,elementsStern,fmt='%5.0f',delimiter=' ')
    
    # Creating test_2.srf and two empty files which run 
    f=open("test_2.srf","w+")
    f.write("f\nf\n./%s\n./%s\n\n\n0\n" % (dielFileName,sternFileName))
    f.close()
    
    
    test_2_Diel = "cyl_%d" % R
    test_2_Stern = "cyl_%d_stern" % R
    Path('%s/%s' % (Rdir,test_2_Diel)).touch()
    Path('%s/%s' % (Rdir,test_2_Stern)).touch()
    
    #Removing old pqr files (if any)
    for item in os.listdir(Rdir):
        if item.endswith(".pqr"):
            os.remove(os.path.join(Rdir, item))
    
    # writing pqr files in each folder
    lambda_0 = 0
    lambda_1 = 1
    lambdaSteps = 10
    for lambdaFEP in np.linspace(lambda_0,lambda_1,lambdaSteps+1):
        h_Cylinder = 50
        
        # Charge locations
        x_q = 0;
        y_q = 0;
        
        # Charge moves from the center of the membrane toward the top plane on z-axis
        z_q = h_Cylinder/2 * (1 + lambdaFEP)
        
        q = 1
        R_q = 2.27 # Ask Jay!
        pqrFileName = "test_lambda_%3.1f.pqr" % lambdaFEP
        f = open(pqrFileName,"w+")
        f.write("%s %5d %4s %5s %5d %15.6f %9.6f %9.6f %9.6f %9.6f" % ('ATOM',1,'NA','TMP',1,x_q,y_q,z_q,q,R_q))
        f.close()
    meshInfo = "%s%d\n%s%d\n%s%d\n%s%d\n%s%d%s" % ('# Elements in dielectric mesh = ',meshDiel.num_elements,
                                           '# Vertices in dielectric mesh = ',meshDiel.num_vertices,
                                           '# Elements in stern layer mesh = ',meshStern.num_elements,
                                           '# Vertices in stern layer mesh = ',meshStern.num_vertices,
                                           'Surface area = ',np.pi*(2 * R ** 2 + 2 * R * h_Cylinder),' AA')
    f = open('mesh.info',"w+")
    f.write(meshInfo)
    f.close()
    os.chdir(pwd)