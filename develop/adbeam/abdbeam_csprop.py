# Stupid fix for adding path
import os
base_path=r"C:\Users\oyvinpet\Anaconda3"
path=os.pathsep.join([os.path.join(base_path, i) for i in [r"", r"bin", r"Scripts", r"Library\mingw-w64\bin", r"Library\bin"]])
os.environ["PATH"]+=os.pathsep+path
# 
import abdbeam as ab
import numpy as np
#Stiffness
t=2.234e-02
E=2.100e+11
v=3.000e-01
t=float(t)
E=float(t)
v=float(v)
G=E/(2*(1+v))
sc=ab.Section()
# Create a materials dictionary:
mts=dict()
mts[1]=ab.Isotropic(t,E,v)
# Create a points dictionary based on Y and Z point coordinates:
pts=dict()
pts[1]=ab.Point(3.750e+00,0.000e+00)
pts[2]=ab.Point(9.750e+00,0.000e+00)
pts[3]=ab.Point(1.650e+01,1.250e+00)
pts[4]=ab.Point(1.475e+01,2.080e+00)
pts[5]=ab.Point(7.500e-01,2.500e+00)
pts[6]=ab.Point(0.000e+00,9.250e-01)
# Create a segments dictionary referencing point and material ids:
sgs=dict()
sgs[1]=ab.Segment(1,2,1)
sgs[2]=ab.Segment(2,3,1)
sgs[3]=ab.Segment(3,4,1)
sgs[4]=ab.Segment(4,5,1)
sgs[5]=ab.Segment(5,6,1)
sgs[6]=ab.Segment(6,1,1)
# Point the dictionaries to the section
sc.materials=mts
sc.points=pts
sc.segments=sgs
# Calculate and output section properties
sc.calculate_properties()
#sc.summary()
#ab.plot_section(sc, figsize=(6.4*0.8, 4.8*0.8))
#Create a single load case and calculate its internal loads
#sc.loads[1]=ab.Load(Vz_s=-100)
#sc.calculate_internal_loads()
#Plot internal loads
#ab.plot_section_loads(sc, 1, int_load_list=[''Nxy''],title_list=[''Abdbeam - Nxy (N/m)''],figsize=(6.4*0.8, 4.8*0.8))
EA=float(sc.p_c[0,0])
EIyy=float(sc.p_c[1,1])
EIzz=float(sc.p_c[2,2])
EIyz=float(sc.p_c[1,2])
GJ=float(sc.p_c[3,3])
A=EA/E
Iyy=EIyy/E
Izz=EIzz/E
Iyz=EIyz/E
J=GJ/G
ValueMatrix=np.array([ sc.yc , sc.zc , sc.ys , sc.zs , A , Iyy, Izz , Iyz, J , sc.principal_axis_angle])
dir_save='C:\\Cloud\\OD_OWP\\Work\\Abaqus\\Suspensionbridge\\crosssection_prop'
name_save='abdbeam_csprop_result'
np.savetxt((dir_save+ '\\' +name_save+ '.txt'), ValueMatrix , delimiter=',')
