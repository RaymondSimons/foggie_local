#!/u/rcsimons/miniconda3/bin/python3.7
import os, sys, argparse
import numpy as np
from numpy import *
from PIL import Image
from joblib import Parallel, delayed
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def parse():
    '''
    Parse command line arguments
    ''' 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='''\
                                Generate the cameras to use in Sunrise and make projection plots
                                of the data for some of these cameras. Then export the data within
                                the fov to a FITS file in a format that Sunrise understands.
                                ''')
    parser.add_argument('-ax', '--ax', default=None, help='DD to use')
    parser.add_argument('-sat', '--sat', default=None, help='DD to use')
    parser.add_argument('-zoom', '--zoom', default=None, help='DD to use')

    args = vars(parser.parse_args())
    return args







DDs = arange(49, 1000)
simnames =  ['natural',
            'natural_v2',
            'natural_v3',
            'natural_v4',
            'nref11c_nref9f',
            'nref11n_nref10f']

'''
axs = ['x', 'y', 'z']
zooms = ['zoomin', 'zoomout', 'zoomoutfar']
sats = arange(6)



#DDs = arange(49, 50)
simnames =  ['natural',
            'natural_v2',
            'natural_v3',
            'natural_v4',
            'nref11c_nref9f',
            'nref11n_nref10f']
axs = ['z']
zooms = ['zoomin']
sats = arange(1)
'''

def combine_frames(sat, ax, zoom, DD, simnames):
    print (DD)
    imgs = [] 
    fname_out = '/nobackupp2/rcsimons/foggie_momentum/sat_figures/combined/%s/%s/%.4i_%.2i_%s_%s.png'%(ax, zoom, DD, sat, ax, zoom)

    for s, simname in enumerate(simnames):
        fname = '/nobackupp2/rcsimons/foggie_momentum/sat_figures/%s/%s/%s/%s_%.4i_%.2i_%s_%s.png'%(simname, ax, zoom, simname, DD, sat, ax, zoom)
        if os.path.isfile(im):
            im =  Image.open(fname)
            imgs.append(im)
        else:
            return
    imgs_12 = imgs[0:2]
    imgs_34 = imgs[2:4]
    imgs_56 = imgs[4:6]
    imgs_comb = []

    #for imgs in [imgs_12, imgs_34, imgs_56]:
    min_shape = sorted([(np.sum(i.size), i.size ) for i in imgs])[0][1]
    imgs_comb_temp = np.vstack((np.asarray( i.resize(min_shape) ) for i in imgs ))
    imgs_comb_temp = Image.fromarray( imgs_comb_temp)
    imgs_comb_temp.save(fname_out)


if __name__ == '__main__':
    args = parse()
    ax = args['ax'] 
    sat = int(args['sat'])
    zoom = args['zoom'] 
    Parallel(n_jobs = -1)(delayed(combine_frames)(sat, ax, zoom, DD, simnames) for d, DD in enumerate(DDs))
    png_names = '/nobackupp2/rcsimons/foggie_momentum/sat_figures/combined/%s/%s/'%(ax, zoom) + '%4d' + '_%.2i_%s_%s.png'%(sat, ax, zoom)
    os.system('ffmpeg -r 24 -f image2 -s 1920x1080 -start_number 49 -i %s -vframes 1000 -vcodec libx264 -crf 25  -pix_fmt yuv420p /nobackupp2/rcsimons/foggie_momentum/sat_figures/movies/%.2i_%s_%s.mp4'%(png_names, sat, ax, zoom))

'''
#Combine Frames
for sat in sats:
    for a, ax in enumerate(axs):
        for z, zoom in enumerate(zooms):
            Parallel(n_jobs = -1)(delayed(combine_frames)(sat, ax, zoom, DD, simnames) for d, DD in enumerate(DDs))

#Make Movie
for sat in sats:
    for a, ax in enumerate(axs):
        for z, zoom in enumerate(zooms):
            png_names = '/nobackupp2/rcsimons/foggie_momentum/sat_figures/combined/%s/%s/'%(ax, zoom) + '%4d' + '_%.2i_%s_%s.png'%(sat, ax, zoom)
            os.system('ffmpeg -r 24 -f image2 -s 1920x1080 -start_number 49 -i %s -vframes 1000 -vcodec libx264 -crf 25  -pix_fmt yuv420p /nobackupp2/rcsimons/foggie_momentum/sat_figures/movies/%.2i_%s_%s.mp4'%(png_names, sat, ax, zoom))
                
'''


















