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

def membraneMeshGen (r_0,spacing, r_f, h,particleCharge):
	curDir = os.getcwd()
	radiiRange= np.arange(r_0,r_f+0.0001,spacing)
	meshDataPath = os.path.join(curDir,'mesh-membrane')
	try:
	    os.makedirs(meshDataPath)
	except OSError as e:
	    if e.errno != errno.EEXIST:
	        raise
	for root, dirs, files in os.walk(meshDataPath):
	    for f in files:
	        os.unlink(os.path.join(root,f))
	    for d in dirs:
	        shutil.rmtree(os.path.join(root,d))  	
	# r = Cylinder (membrane) radius
	for r_mem in radiiRange:

		dirPath = os.path.join(meshDataPath,"%3.1f" % r_mem)
		os.makedirs(dirPath)
		os.chdir(dirPath)

		dielgeo = "diel%d.geo" % (r_mem)
		sterngeo = "stern%d.geo" % (r_mem)
		dielstl = "diel%d.stl" % (r_mem)
		sternstl = "stern%d.stl" % (r_mem)
		listDir = os.listdir(dirPath)

		# Removing old files		
		for item in listDir:

			if (item.endswith(".stl") or item.endswith(".geo") or
				item.endswith(".vert") or item.endswith(".face") or
				item.endswith(".msh") or item.endswith(".pqr") or
				os.path.getsize(os.path.join(curDir, item)) == 0 or item.endswith(".info")):
				os.remove(os.path.join(dirPath, item))

		copy2(os.path.join(curDir,'template.geo'),dirPath)

		# Generating *.geo file
		lc = 0.65 * np.sqrt(r_mem)
		lc_fine = 0.00001 * r_mem
		lc_st = 0.8 * np.sqrt(r_mem+2)
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
		meshDiel, info = pymesh.collapse_short_edges(meshDiel, abs_threshold=0.4 , preserve_feature=False) ##### New
		meshDiel, info = pymesh.remove_duplicated_vertices(meshDiel)
		meshDiel, info = pymesh.remove_duplicated_faces(meshDiel)
		elementsDiel = np.copy(meshDiel.elements)
		for i in range(elementsDiel.shape[0]):
		    for j in range(elementsDiel.shape[1]):
		        elementsDiel[i,j]+=1

		#Stern surface generation        
		meshStern = pymesh.load_mesh(sternstl)
		meshStern, info = pymesh.collapse_short_edges(meshStern, abs_threshold=0.4 , preserve_feature=False) ##### New
		meshStern, info = pymesh.remove_duplicated_vertices(meshStern)
		meshStern, info = pymesh.remove_duplicated_faces(meshStern)
		
		elementsStern = np.copy(meshStern.elements)

		#Because element numbers in PyMesh start from zero!
		for i in range(elementsStern.shape[0]):
		    for j in range(elementsStern.shape[1]):
		        elementsStern[i,j]+=1
		    
		dielFileName = "diel%d" % r_mem
		sternFileName = "stern%d" % r_mem

		np.savetxt("%s.vert" % dielFileName,meshDiel.vertices,fmt='%5.3f',delimiter=' ')
		np.savetxt("%s.face" % dielFileName,elementsDiel,fmt='%5.0f',delimiter=' ')
		np.savetxt("%s.vert" % sternFileName,meshStern.vertices,fmt='%5.3f',delimiter=' ')
		np.savetxt("%s.face" % sternFileName,elementsStern,fmt='%5.0f',delimiter=' ')
		
		# Creating test_2.srf and two empty files which run 
		f=open("test_2.srf","w+")
		f.write("f\nf\n./%s\n./%s\n\n\n0\n" % (sternFileName,dielFileName))
		f.close()
		
		
		test_2_Diel = "diel%d" % r_mem
		test_2_Stern = "stern%d" % r_mem
		Path('%s/%s' % (dirPath,test_2_Diel)).touch()
		Path('%s/%s' % (dirPath,test_2_Stern)).touch()

		
		# writing membrane pqr files 

		#chargedMembrane = True

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

        
		lambda_z = [25.0,35.0,40.0,42.5,45.0,46.0,47.0,47.5,48.0,48.5]
		k=1

		for lambdaFEP in lambda_z:
		    
		    # Charge locations
		    x_q = 0;
		    y_q = 0;
		    
		    # Charge moves from the center of the membrane toward the top plane on z-axis
		    # z_q = h/2 * (1 + lambdaFEP)
		    
		    
		    R_q = 1 # Ask Jay!
		    fullpo4PqrFileName = "test_lambda_full_po4_%d.pqr" % k
		    halfpo4PqrFileName = "test_lambda_half_po4_%d.pqr" % k
		    pqrFileName = "test_lambda_%d.pqr" % k
		    f = open(pqrFileName,"w+")
		    pqr_ion = "%s %5d %4s %5s %5d %15.6f %9.6f %9.6f %9.6f %9.6f\n" % ('ATOM',1,'NA','TMP',1,x_q,y_q,lambdaFEP,particleCharge,R_q)
		    f.write(pqr_ion)
		    f.close()
		    k += 1

		    insert('fullMempo4.pqr',pqr_ion,fullpo4PqrFileName)
		    insert('halfMempo4.pqr',pqr_ion,halfpo4PqrFileName)

		meshInfo = "%s%d\n%s%d\n%s%d\n%s%d\n%s%d%s" % ('# Elements in dielectric mesh = ',meshDiel.num_elements,
		                                       '# Vertices in dielectric mesh = ',meshDiel.num_vertices,		                                       '# Elements in stern layer mesh = ',meshStern.num_elements,
		                                       '# Vertices in stern layer mesh = ',meshStern.num_vertices,
		                                       'Surface area = ',np.pi*(2 * r_mem ** 2 + 2 * r_mem * h),' AA')
		f = open('mesh.info',"w+")
		f.write(meshInfo)
		f.close()
		os.remove(os.path.join(dirPath,'template.geo'))
		os.chdir(dirPath)
