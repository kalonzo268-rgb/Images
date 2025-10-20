from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

filename = "tess2025258040922-s0097-1-1-0294-s_ffic.fits"
hdul = fits.open(filename)

# See HDU structure
hdul.info()

# Get main image data and header
image_data = hdul[1].data
header = hdul[1].header

# Plot
plt.figure(figsize=(10, 10))
plt.imshow(image_data, cmap='gray', origin='lower',
           vmin=np.percentile(image_data, 5),
           vmax=np.percentile(image_data, 99))
plt.colorbar(label='Flux (e-/s)')
plt.title('TESS Full Frame Image')
plt.xlabel('Pixel X')
plt.ylabel('Pixel Y')
plt.show()

hdul.close()
