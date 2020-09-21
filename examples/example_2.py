import formido.fd as fd
from astropy.table import Table
from astropy.io import ascii as ascii
import numpy as np

import matplotlib.pyplot as plt

data = Table.read('data/boo1.fits')
mask = data['ZQUALITY'] > 2
data = data[mask]

iso = ascii.read('data/boo1_iso.txt')
q = fd.Formido()

# Add stars from our fits file
q.AddStars(data['RA'], data['DEC'], data['RMAG'], data['GMAG'])
# Store our magnitude adjustments based on the distance and E(B-V), as well as
# the associated band scalars for extinction.
adj = q.Adjustments(dist=66*1000, ebv=0.0151, rb_scalar=2.751, bb_scalar=3.793)
color = q.GenerateColor()
hf = q.GenerateHLProb(.175, {'ra': 210.0200348, 'dec': 14.5135002})
# Generate probabilities based on model isochrone which is adjusted
# by our extinction and distance modulus to match the star data.
idprob = q.GenerateIDProb(iso['gmag'] + adj['d_bmag'], iso['rmag']+ adj['d_rmag'], verbose=True)

cmap = plt.cm.get_cmap('plasma')

fig, ax = plt.subplots(2,2, figsize=(10,10))

fig.suptitle('Boötes 1 Membership Probability')

ax[0,0].set_title('Boötes 1')
ax[0,0].scatter(data['RA'], data['DEC'])
ax[0,0].set_xlabel('RA')
ax[0,0].set_ylabel('Dec')

ax[0,1].set_title('Half-light Radius Probability')
hl_plot = ax[0,1].scatter(data['RA'], data['DEC'], c=hf, cmap=cmap, s=hf*50)
ax[0,1].set_xlabel('RA')
ax[0,1].set_ylabel('Dec')
cb = plt.colorbar(hl_plot, ax=ax[0,1], fraction=0.05, orientation='vertical')
cb.set_label('Probability')

stars = q.get()

ax[1,0].set_title('Color Magnitude Diagram')
ax[1,0].scatter(color, stars['rb'])
ax[1,0].set_ylim(23,16)
ax[1,0].scatter(iso['gmag'] - iso['rmag'] + adj['d_color'], iso['rmag'] + adj['d_rmag'], label='Isochrone', c='r', s=5)
ax[1,0].legend()

ax[1,1].set_title('Membership Probability based on Isochrone Distance')
id_plot = ax[1,1].scatter(color, stars['rb'], c=idprob, cmap=cmap)
ax[1,1].set_ylim(23,16)
ax[1,1].scatter(iso['gmag'] - iso['rmag'] + adj['d_color'], iso['rmag'] + adj['d_rmag'], label='Isochrone', c='r', s=5)
ax[1,1].legend()
cb = plt.colorbar(id_plot, ax=ax[1,1], fraction=0.05, orientation='vertical')
cb.set_label('Probability')

plt.savefig('plots/example_2.pdf')
