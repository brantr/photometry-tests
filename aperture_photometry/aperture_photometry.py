from astropy.io import fits
import matplotlib.pyplot as plt
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
plt.imshow(p_image,origin="lower")
plt.colorbar()
plt.show()

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
plt.imshow(bkg_image, cmap='gray', origin='lower')
plt.colorbar()
plt.show()

#get a background rms image
bkg_rms_image = bkg.rms()

#show the background rms
plt.imshow(bkg_rms_image, cmap='gray', origin='lower')
plt.colorbar()
plt.show()

#get a background subtracted image
data_sub = data - bkg

#close the file
hdulist.close()