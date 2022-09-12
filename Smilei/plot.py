execfile("/Users/sharong/Desktop/smilei/Smilei/scripts/Diagnostics.py")
import happi
S=happi.Open("/Users/sharong/Desktop/smilei/Smilei/Data3")
dt = S.namelist.Main.timestep
dx = S.namelist.Main.cell_length[0]
dy = S.namelist.Main.cell_length[1]
simulation_time = S.namelist.Main.simulation_time
FieldDiag = S.Field(0, "Rho_electron")

f = FieldDiag.getData()[0].T
x = FieldDiag.getAxis("x")/(16*dx)
y = FieldDiag.getAxis("y")/(8*dy)-17

print(x.max())
print(y.max())
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
fig=plt.figure(figsize=[11,5])
# plt.plot(y, f[:,18])
plt.imshow(f, cmap="jet", extent=[x[0],x[-1],y[0],y[-1]], aspect="auto")
plt.show()