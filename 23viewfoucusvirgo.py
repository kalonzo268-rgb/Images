from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

filename = "tess2020084105920-s0023-1-4-0177-s_ffic.fits"
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

# Open the FITS file
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

filename = "tess2020084105920-s0023-1-4-0177-s_ffic.fits"  # Example Sector 23 filename
hdul = fits.open(filename)

image_data = hdul[1].data
header = hdul[1].header
wcs = WCS(header)

virgo_center = SkyCoord(ra=187.5*u.deg, dec=12*u.deg, frame='icrs')

px, py = wcs.world_to_pixel(virgo_center)
print(f"Virgo center pixel coordinates: x={px:.1f}, y={py:.1f}")

half_size = 300  # 50 pixels on each side
x0, x1 = int(px - half_size), int(px + half_size)
y0, y1 = int(py - half_size), int(py + half_size)

cutout = image_data[y0:y1, x0:x1]

plt.figure(figsize=(8, 8))
plt.imshow(cutout, cmap='gray', origin='lower',
           vmin=np.percentile(cutout, 5),
           vmax=np.percentile(cutout, 99))
plt.title('Region around Virgo Cluster (TESS Sector 23)')
plt.xlabel('Pixel X')
plt.ylabel('Pixel Y')
plt.colorbar(label='Flux (e-/s)')
plt.show()

print(image_data.shape)
