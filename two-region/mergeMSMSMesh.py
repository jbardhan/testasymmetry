
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import pymesh

with open('mol_stern.vert') as f_vert:
    lines_v = [line_v.rstrip('\n').split() for line_v in f_vert]
f_vert.close()
x_v = []
y_v =[]
z_v = []

for line_v in lines_v:
    x_v.append(float(line_v[0]))
    y_v.append(float(line_v[1]))
    z_v.append(float(line_v[2]))
    
vert = np.array([[x_v[i],y_v[i],z_v[i]] for i in range(len(lines_v))], dtype=np.float32)


with open('mol_stern.face') as f_face:
    lines_f = [line_f.rstrip('\n').split() for line_f in f_face]
f_face.close()
x_f = []
y_f =[]
z_f = []

for line_f in lines_f:
    x_f.append(int(line_f[0])-1)
    y_f.append(int(line_f[1])-1)
    z_f.append(int(line_f[2])-1)
    
face = np.array([[x_f[i],y_f[i],z_f[i]] for i in range(len(lines_f))], dtype=np.int)
mesh_ion = pymesh.form_mesh(vert, face)
mesh_mem = pymesh.load_mesh("stern15.stl")
output_mesh = pymesh.boolean(mesh_ion, mesh_mem,operation="union", engine = "cork")
pymesh.save_mesh("merged.stl", output_mesh);

with open('mol_diel.vert') as f_vert:
    lines_v = [line_v.rstrip('\n').split() for line_v in f_vert]
f_vert.close()
x_v = []
y_v =[]
z_v = []

for line_v in lines_v:
    x_v.append(float(line_v[0]))
    y_v.append(float(line_v[1]))
    z_v.append(float(line_v[2]))
    
vert1 = np.array([[x_v[i],y_v[i],z_v[i]] for i in range(len(lines_v))], dtype=np.float32)


with open('mol_diel.face') as f_face:
    lines_f = [line_f.rstrip('\n').split() for line_f in f_face]
f_face.close()
x_f = []
y_f =[]
z_f = []

for line_f in lines_f:
    x_f.append(int(line_f[0])-1)
    y_f.append(int(line_f[1])-1)
    z_f.append(int(line_f[2])-1)
face1 = np.array([[x_f[i],y_f[i],z_f[i]] for i in range(len(lines_f))], dtype=np.int)
diel_ion = pymesh.form_mesh(vert1, face1)
pymesh.save_mesh("ion_diel.stl", diel_ion);