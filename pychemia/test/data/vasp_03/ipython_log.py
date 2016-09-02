import pychemia
from mayavi.mlab import *

st = pychemia.code.vasp.read_poscar()
print(st)
vo = pychemia.code.vasp.VaspOutput()
vo.get_magnetization()

x = []
y = []
z = []
u = []
v = []
w = []

for i in range(st.natom):
    x.append(st.positions[i][0])
    y.append(st.positions[i][1])
    z.append(st.positions[i][2])
    u.append(vo.get_magnetization()['x']['total'][i])
    v.append(vo.get_magnetization()['y']['total'][i])
    w.append(vo.get_magnetization()['z']['total'][i])

quiver3d(x, y, z, u, v, w)
