from astroquery.skyview import SkyView
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy.wcs
import numpy as np
from astropy.io import ascii
import sys
import os
import pylab 
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from photutils.background import Background2D, MedianBackground
from astropy.stats import SigmaClip
import warnings
warnings.filterwarnings('ignore')


def make_finders(target_list, directory='finders', clobber=False):

    sv = SkyView()

    if not os.path.exists(directory):
        os.makedirs(directory)

    # Params for making finders
    params = {'radius': 14.0*u.arcmin,
              'pixels': '1000',
              'deedger': 'skyview.process.Deedger',
              'survey': ['2MASS-J']}

    target_table = ascii.read(target_list, names=('field','ra','dec','epoch'))

    for row in target_table:

        coord = SkyCoord(row['ra'], row['dec'], unit=(u.hour, u.deg))

        rah,decd = coord.to_string(style='hmsdms', sep=':', precision=3).split()
        name = row['field']

        outimagename = os.path.join(directory, str(name+'_finder.png'))
        if os.path.exists(outimagename) and not clobber:
            print(f'WARNING: {outimagename} exists and clobber=False.  Skipping...')
            continue

        images = sv.get_images(position=coord, **params)
        image = images[0][0]

        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
        bkg = Background2D(image.data, (9, 9), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

        image.data = image.data - bkg.background

        wcs = astropy.wcs.WCS(image.header)

        fig, ax = plt.subplots(1,1,figsize = (8,6), subplot_kw={'projection':wcs})
        plt.set_cmap('gray_r')

        smoothedimage = gaussian_filter(image.data, 1.3)
        ax.imshow(smoothedimage, origin='lower',vmin=np.percentile(image.data.flatten(), 10), \
                vmax=np.percentile(image.data.flatten(), 99.0))

        # Set limits to size of image
        ax.set_xlim([0,(image.data.shape[0])])
        ax.set_ylim([0,(image.data.shape[1])])

        # Plot compass
        plt.plot([(image.data.shape[0])-10,(image.data.shape[0]-200)],[10,10], 'r-', lw=2)
        plt.plot([(image.data.shape[0])-10,(image.data.shape[0])-10],[10,200], 'r-', lw=2)
        plt.annotate("N", xy=((image.data.shape[0])-20, 200), color='r',  xycoords='data',xytext=(-4,5), textcoords='offset points')
        plt.annotate("E", xy=((image.data.shape[0])-200, 20), color='r',  xycoords='data',xytext=(-12,-5), textcoords='offset points')

        # Set axis tics (not implemented correctly yet)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))

        plt.text(1.02, 0.95, name, transform=ax.transAxes, fontweight='bold')
        plt.text(1.02, 0.90,rah+"  "+decd, transform=ax.transAxes)
        plt.text(1.02, 0.85,'2MASS J-band', transform=ax.transAxes)

        plt.xlabel('RA (J2000)')
        plt.ylabel('Dec (J2000)')

        outimagename = os.path.join(directory, str(name+'_finder.png'))
        pylab.savefig(os.path.join(directory, str(name+'_finder.png')), bbox_inches = 'tight', dpi = 150)
        print("Saved to %s"%os.path.join(directory, str(name+'_finder.png')))

        plt.clf()

if __name__=="__main__":
    make_finders(sys.argv[1])




