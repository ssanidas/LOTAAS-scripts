file_name = '/projects/0/lotaas/data/raw/L318124_red/stokes/SAP2/BEAM71/L318124_SAP2_BEAM71.fits'
subint_start = 1
subint_end = 4

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pyfits


fits = pyfits.open(file_name,memmap=True)
subint = fits['SUBINT'].data[subint_start:subint_end]['DATA']

#Downsample the data
N_channels = fits['SUBINT'].header['NCHAN']
N_spectra = (subint_end-subint_start)*fits['SUBINT'].header['NSBLK']
subint = subint.reshape(N_spectra,N_channels)
subint = subint[:,~np.all(subint == 0, axis=0)]  #Remove null lines from the spectrum

plt.figure(figsize=(20,10))
plt.xlabel('Time (ch)')
plt.ylabel('Frequency (ch)')
plt.imshow(subint.T,cmap=mpl.cm.hot_r,origin="lower",aspect='auto',interpolation='nearest')
plt.colorbar()

plt.savefig('DynSpect.png',format='png',bbox_inches='tight')



