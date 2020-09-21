import formido.fd as fd
from astropy.table import Table
import matplotlib.pyplot as plt

data = Table.read('data/hercules.fits')
mask = data['ZQUALITY'] > 2
data = data[mask]

q = fd.Formido()

# Add the stars from our fits file
q.AddStars(data['RA'], data['DEC'], data['RMAG'], data['GMAG'])
# Print magnitude adjustments based on extinction and dist. modulus
print(q.Adjustments(dist=132*1000, ebv=0.0549, rb_scalar=2.751, bb_scalar=3.793))
# Calculate the half-light radius probability. 
hf = q.GenerateHLProb(5.83/60, center={'ra': 247.7722015, 'dec': 12.7851944})

cmap = plt.cm.get_cmap('plasma')

fig, ax = plt.subplots(2, figsize=(10,10))
fig.suptitle('Hercules with Half-light Radius Distance Probability')
ax[0].scatter(data['RA'], data['DEC'])
ax[0].set_xlabel('RA')
ax[0].set_ylabel('Dec')
p = ax[1].scatter(data['RA'], data['DEC'], c=hf, s=hf*50, cmap=cmap)
ax[1].set_xlabel('RA')
ax[1].set_ylabel('Dec')

cb = plt.colorbar(p, ax=ax[1], fraction=0.05, orientation='horizontal')
cb.set_label('Probability')

plt.savefig('plots/example_1.pdf')