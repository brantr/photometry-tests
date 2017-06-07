from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from sorting_routines import *
import numpy as np
import sep
import sys

#define the fits filename
fname = 'input.fits'
if(len(sys.argv)>1):
  fname = sys.argv[1]

#open the fits file
hdulist = fits.open(fname)

#get some info about the file
hdulist.info()

#get the header of the 0th element of hdulist
header = hdulist[0].header

#print the header
print header

#get the data from the 0th HDU,
#allowing for multiple layers in the data
if (header['NAXIS']==2):
  data = hdulist[0].data[:,:]
else:
  data = hdulist[0].data[0,:,:]

data = data.byteswap().newbyteorder()
image = data.copy()

#set pixels with nan to 0 
idx_nan  = np.where(np.isnan(image)==True)
image[idx_nan[0][:],idx_nan[1][:]] = 0

#print the shape of the image
print image.shape

#print the type of the data
print image.dtype.name

#copy the image
p_image = image.copy()

#find the minimum of the image
image_min = p_image.min()
#shift to zero 
p_image -= image_min
#scale the image to have
#a fixed dynamic range, measured
#from the max
image_max = np.mean(p_image) + 10*(np.std(p_image))
idx = np.where(p_image>image_max)
if(len(idx)>0):
  p_image[idx] = image_max
#dyn_range = 10.
#image_min_set = image_max/dyn_range
image_min = np.mean(p_image) - 1*(np.std(p_image))
if(image_min<1.0e-10):
  image_min = 1.0e-10
idx = np.where(p_image<image_min)
if(len(idx)>0):
  p_image[idx] = image_min
print image_min, np.median(p_image), image_max

#print image_min, image_max, image_min_set
#p_image = image.copy()

p_image = (np.log10(p_image)-np.log10(image_min))/(np.log10(image_max)-np.log10(image_min))
fig, ax = plt.subplots()
ax.imshow(p_image,origin="lower")
#fig.colorbar()

###########################
# Begin photometry
# see http://sep.readthedocs.io/en/v1.0.x/tutorial.html
###########################


#measure a spatially varying background
bkg = sep.Background(image)

#get a global mean and rms of the image bkg
print bkg.globalback
print bkg.globalrms

#get a background image
bkg_image = bkg.back()

#show the background
#plt.imshow(bkg_image, cmap='gray', origin='lower')
#plt.colorbar()
#plt.show()

#get a background rms image
bkg_rms_image = bkg.rms()

#show the background rms
#plt.imshow(bkg_rms_image, cmap='gray', origin='lower')
#plt.colorbar()
#plt.show()

#get a background subtracted image
data_sub = image - bkg


###########################
# Perform object detection
# see http://sep.readthedocs.io/en/v1.0.x/tutorial.html#Object-detection
###########################

#sigma_detection = 1.5 #1.5 sigma sources
sigma_detection = 3.0 #1.5 sigma sources

#identify objects based on sigma detection threshold
objects = sep.extract(data_sub, sigma_detection, err=bkg.globalrms)

#print the number of objects detected
print "Detected %d objects..." % len(objects)

#plot an ellipse for each object
for i in range(len(objects)):
    e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                width=6*objects['a'][i],
                height=6*objects['b'][i],
                angle=objects['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)

#show the figure
plt.show()


###########################
# Perform aperture photometry
# see http://sep.readthedocs.io/en/v1.0.x/tutorial.html#Aperture-photometry
###########################

r_aperture = 0.3 # arcseconds
pixel_scale = 0.0317
n_pixels = r_aperture / pixel_scale
#aperture correction for 0.3" radius = 0.204462939361

flux, fluxerr, flag = sep.sum_circle(data_sub, objects['x'], objects['y'], n_pixels, err=bkg.globalrms, gain=1.0)


#compare with star catalog
fname_cat = "data/star_simulations/stars_f090w.cat"
fp = open(fname_cat,"r")
fl = fp.readlines()
fp.close()
mags = np.zeros(len(fl))
for i in range(len(fl)):
  mags[i] = l[i].split()[5]

fi = sorted_index(-1.*flux)

for i in range(100):
  


#close the file
hdulist.close()