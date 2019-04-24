import numpy as np
from mpl_toolkits import mplot3d
#plt.ioff()
plt.close("all")


arr = np.load('/Users/rsimons/Dropbox/rcs_foggie/anchor_files/all_anchors.npy')[()]

fig = plt.figure()
ax = plt.axes(projection='3d')

central_xyz_fit = np.load('/Users/rsimons/Dropbox/rcs_foggie/catalogs/center_natural.npy')[()]

xf = central_xyz_fit['x']
yf = central_xyz_fit['y']
zf = central_xyz_fit['z']



t = arange(0, shape(arr)[1]) - 1

DD = t + 44
central_x = xf[0] * DD**4. + xf[1] * DD**3. + xf[2] * DD**2. + xf[3] * DD + xf[4]
central_y = yf[0] * DD**4. + yf[1] * DD**3. + yf[2] * DD**2. + yf[3] * DD + yf[4]
central_z = zf[0] * DD**4. + zf[1] * DD**3. + zf[2] * DD**2. + zf[3] * DD + zf[4]


sat_n = 11

clrs = ['red', 'blue', 'green', 'orange', 'purple']



for s in arange(1,5):

    x = arr[s,t, sat_n, 0] - central_x
    y = arr[s,t, sat_n, 1] - central_y
    z = arr[s,t, sat_n, 2] - central_z
    ax.plot3D(x, y, z, clrs[s])

ax.set_xlim(-50, 50)
ax.set_ylim(-50, 50)
ax.set_zlim(-50, 50)

fig.savefig('/Users/rsimons/Dropbox/rcs_foggie/figures/sat_tracks/%i.png'%sat_n, dpi = 300)


















