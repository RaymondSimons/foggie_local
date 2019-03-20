import astropy
from astropy.io import fits
import os
import yt




anchor_dir = '/Users/rsimons/Dropbox/rcs_foggie/anchor_files'


simnames = ['natural', 'nref11n_nref10f_selfshield_z6']

DD_arr = arange(44, 1100)
tracks_file = '/Users/rsimons/Dropbox/rcs_foggie/catalogs/satellite_tracks.npy'
dist_from_center_file = '/Users/rsimons/Dropbox/rcs_foggie/catalogs/satellite_dist_from_center.npy'
if False:
    tracks = nan * zeros((2, len(DD_arr), 12, 3))


    for sm, simname in enumerate(simnames):
        if simname == 'natural': nsats = 12
        elif simname == 'nref11n_nref10f_selfshield_z6': nsats = 11
        for d, DD in enumerate(DD_arr):
            print DD
            fits_file = anchor_dir + '/' + '%s_DD%.4i_anchorprops.fits'%(simname, DD)
            if os.path.isfile(fits_file):
                fits_data = fits.open(fits_file)
                for s, sat_n in enumerate(arange(nsats)):
                    if len(fits_data['SAT_%.2i'%sat_n].data['box_avg']) > 0:
                        tracks[sm, d, sat_n, 0:3] = fits_data['SAT_%.2i'%sat_n].data['box_avg'][0:3]                
            else:
                print 'no track file', DD
    np.save(tracks_file, tracks)
    
if True:
    tracks = np.load(tracks_file)
    mn_dists = nan * zeros((12, 12))
    mx_dists = nan * zeros((12, 12))

    for s1, sat_n1 in enumerate(arange(12)):
        for s2, sat_n2 in enumerate(arange(12)):
            dist = sqrt((tracks[0,:,sat_n1,0] - tracks[1,:,sat_n2,0])**2. +
                        (tracks[0,:,sat_n1,1] - tracks[1,:,sat_n2,1])**2. +
                        (tracks[0,:,sat_n1,2] - tracks[1,:,sat_n2,2])**2.)
            mn_dists[s1, s2] = nanmean(dist)
            mx_dists[s1, s2] = nanmax(dist)
        
        gd = where(mn_dists[s1,:] == nanmin(mn_dists[s1,:]))[0][0]
        print '%.1f  %.1f'%(mn_dists[s1,gd], mx_dists[s1, gd]), sat_n1, gd


if False:
    tracks = np.load(tracks_file)
    dist_from_center = nan * zeros((2, len(DD_arr), 12))
    for sm, simname in enumerate(simnames):
        if simname == 'natural': nsats = 12
        elif simname == 'nref11n_nref10f_selfshield_z6': nsats = 11
        central_xyz_fit = np.load('/Users/rsimons/Dropbox/rcs_foggie/catalogs/center_%s.npy'%simname)[()]

        xf = central_xyz_fit['x']
        yf = central_xyz_fit['y']
        zf = central_xyz_fit['z']

        for d, DD in enumerate(DD_arr):
            central_x = xf[0] * DD**4. + xf[1] * DD**3. + xf[2] * DD**2. + xf[3] * DD + xf[4]
            central_y = yf[0] * DD**4. + yf[1] * DD**3. + yf[2] * DD**2. + yf[3] * DD + yf[4]
            central_z = zf[0] * DD**4. + zf[1] * DD**3. + zf[2] * DD**2. + zf[3] * DD + zf[4]
            print DD
            fits_file = anchor_dir + '/' + '%s_DD%.4i_anchorprops.fits'%(simname, DD)

            for s, sat_n in enumerate(arange(nsats)):
                cenx, ceny, cenz = tracks[sm, d, sat_n, :] 
                cen_gal = yt.YTArray([cenx - central_x, ceny - central_y, cenz - central_z], 'kpc')
                dist_from_center[sm, d, sat_n] = sqrt(sum(cen_gal**2.))
    np.save(dist_from_center_file, dist_from_center)















