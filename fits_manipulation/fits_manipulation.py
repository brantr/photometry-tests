from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
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

#get the data from the 0th HDU
data = hdulist[0].data

#print the shape of the image
print data.shape

#print the type of the data
print data.dtype.name

#plot the image
image = data.copy()
image_min = image.min()
image_max = image.max()
print image_min, image_max
image = np.log10(image)
plt.imshow(image,origin="lower")
plt.show()

#close the file
hdulist.close()