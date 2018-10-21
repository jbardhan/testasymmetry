
import os
from pathlib import Path
import pymesh
import numpy as np
from shutil import copy2

def insert(originalfile, string, newname):
    with open(originalfile,'r') as f:
        with open('newfile.txt','w') as f2: 
            f2.write(string)
            f2.write(f.read())
    os.rename('newfile.txt',newname)

def membraneMeshGen (r_0, r_f, h):
	pwd = os.getcwd()
	
	# r = Cylinder (membrane) radius
	for r in range(r_0, r_f+1):

		curDir = pwd + "/%d" % r
		
		os.chdir(curDir)
		dielgeo = "diel%d.geo" % (r)
		sterngeo = "stern%d.geo" % (r)
		dielstl = "diel%d.stl" % (r)
		sternstl = "stern%d.stl" % (r)
		listDir = os.listdir(curDir)

		# Removing old files		
		for item in listDir:

			if (item.endswith(".stl") or item.endswith(".geo") or
				item.endswith(".vert") or item.endswith(".face") or
				item.endswith(".msh") or item.endswith(".pqr") or
				os.path.getsize(os.path.join(curDir, item)) == 0 or item.endswith(".info")):
				os.remove(os.path.join(curDir, item))

		copy2(os.path.join(pwd,'template.geo'),curDir)

		# Generating *.geo file
		
		# Defining mesh parameters
		# Parameter c is set by trial and error for now and chanegs with the radius
		c = 2
		lc = 0.5 * np.sqrt(c*r)
		lc_fine = 0.001 * r
		lc_st = 0.5 * np.sqrt(c*(r+2))
		lc_fine_st = 0.001 * (r+2)
		r_diel = r;
		

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
		dielMeshCommand = "%s%s%s%s" % ('gmsh ',dielgeo,' -2 -format stl ',dielstl)
		sternMeshCommand = "%s%s%s%s" % ('gmsh ',sterngeo,' -2 -format stl ',sternstl)

		os.system(dielMeshCommand)
		os.system(sternMeshCommand)


		# Each folder corresponds to one radius and has geometric information and
		# other information such as charge distribution for each charge location (*.pqr).
		# Geomtry file (*.geo) that is used to enerate each geometry is also included.
		# These folders are setup in a way that it needs minimal change in order to be read
		# by SLIC Matlab code
		
		
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
		
		# Creating test_2.srf and two empty files which run 
		f=open("test_2.srf","w+")
		f.write("f\nf\n./%s\n./%s\n\n\n0\n" % (dielFileName,sternFileName))
		f.close()
		
		
		test_2_Diel = "diel%d" % r
		test_2_Stern = "stern%d" % r
		Path('%s/%s' % (curDir,test_2_Diel)).touch()
		Path('%s/%s' % (curDir,test_2_Stern)).touch()

		
		# writing pqr files in each folder
		lambda_0 = 0
		lambda_1 = 1
		lambdaSteps = 10
		for lambdaFEP in np.linspace(lambda_0,lambda_1,lambdaSteps+1):
		    
		    # Charge locations
		    x_q = 0;
		    y_q = 0;
		    
		    # Charge moves from the center of the membrane toward the top plane on z-axis
		    z_q = h/2 * (1 + lambdaFEP)
		    
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
		                                       'Surface area = ',np.pi*(2 * r ** 2 + 2 * r * h),' AA')
		f = open('mesh.info',"w+")
		f.write(meshInfo)
		f.close()
		os.remove(os.path.join(curDir,'template.geo'))
		os.chdir(pwd)
