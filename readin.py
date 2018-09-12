import astropy
from astropy.io import fits
import yt


if True:
	ds = yt.load('~/Dropbox/rcs_foggie/data/halo_008508/nref11n_nref10f_selfshield_z6/RD0018/RD0018')
	ds = yt.load('~/Dropbox/rcs_foggie/data/halo_008508/nref11n_selfshield_z15/RD0018/RD0018')
	ad = ds.all_data()


if True:
	x = histogram(ad['particle_position_x'].value, bins = linspace(0.45, 0.55, 10000))
	good_argmax_x = argmax(x[0])
	good_x = (x[1][good_argmax_x] + x[1][good_argmax_x+1])/2.

	y = histogram(ad['particle_position_y'].value, bins = linspace(0.45, 0.55, 20000))
	good_argmax_y = argmax(y[0])
	good_y = (y[1][good_argmax_y] + y[1][good_argmax_y+1])/2.

	z = histogram(ad['particle_position_z'].value, bins = linspace(0.45, 0.55, 10000))
	good_argmax_z = argmax(z[0])
	good_z = (z[1][good_argmax_z] + z[1][good_argmax_z+1])/2.


for i in arange(300):
	print i, 1.*cos(pi*(i)/100.),1*sin(pi*(i)/100.)

	L = [1*cos(pi*(i)/100.),0,1*sin(pi*(i)/100.)] # vector normal to cutting plane
	north_vector = [0,1,0]
	if i > 50:
		cn = 0.00002*(i - 50)
		x_w = max(0.002 - cn, 0.0005)
		y_w = max(0.002 - cn, 0.0005)
		z_w = max(0.002 - cn, 0.0005)
		W = [x_w, y_w, z_w]
	else:
		W = [0.002, 0.002, 0.002]

	c = [good_x, good_y, good_z]
	N = 512

	image = yt.off_axis_projection(ds, c, L, W, N, "density", north_vector =  north_vector)
	yt.write_image(np.log10(image), "./movie_1/%i_offaxis_projection.png" % i)	
