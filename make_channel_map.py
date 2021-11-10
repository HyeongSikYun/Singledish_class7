import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata
from astropy.io import fits
from modules_v06 import *
from astropy.wcs import WCS

#################################################
dirct = data_directory
min_v = minimum_velocity; max_v = maximum_velocity; bin_v = velocity_bin_size # [km/s]
cent_v = system_velocity # [km/s]
int_max = max_val_of_img; int_min = 0-(int_max/3.0)
fits_file = 'cube_file_name.fits'
###################################################
hdulist = read_fits(dirct,fits_file)

def read_coord_info(hdr):
  ra_npx = hdr['NAXIS1']
  de_npx = hdr['NAXIS2']
  ra_ref = hdr['CRVAL1']
  de_ref = hdr['CRVAL2']
  ra_rpx = hdr['CRPIX1']
  de_rpx = hdr['CRPIX2']
  ra_del = hdr['CDELT1']
  de_del = hdr['CDELT2']
  ve_npx = hdr['NAXIS3']
  ve_ref = hdr['CRVAL3']
  ve_rpx = hdr['CRPIX3']
  ve_del = hdr['CDELT3']
  coord_info = [ra_npx,ra_ref,ra_rpx,ra_del,
                de_npx,de_ref,de_rpx,de_del,
                ve_npx,ve_ref,ve_rpx,ve_del]
  return coord_info

coords = read_coord_info(hdulist[0].header)
min_vel = (coords[9]-(coords[10]-1)*coords[11])*1.0e-3
max_vel = (coords[9]+(coords[8]-coords[10])*coords[11])*1.0e-3
vspace = np.linspace(min_vel,max_vel,coords[8])

new_hdr = copy.deepcopy(hdulist[0].header)
new_hdr['NAXIS'] = 2
del new_hdr['NAXIS3']
del new_hdr['NAXIS4']
del new_hdr['CTYPE3']
del new_hdr['CRVAL3']
del new_hdr['CDELT3']
del new_hdr['CRPIX3']
del new_hdr['CROTA3']
del new_hdr['CTYPE4']
del new_hdr['CRVAL4']
del new_hdr['CDELT4']
del new_hdr['CRPIX4']
del new_hdr['CROTA4']
wcs = WCS(new_hdr)

data = hdulist[0].data[0,:,:,:]
print 'spectral cube data is readed'
min_v = np.float(min_v); max_v = np.float(max_v)
bin_v = np.float(bin_v); cent_v = np.float(cent_v)

fig_num1 = np.int((max_v - cent_v + bin_v/2.0)/bin_v)
fig_num2 = np.int((cent_v - min_v - bin_v/2.0)/bin_v)
fig_num = fig_num1 + fig_num2 + 1
layout_x = 2
layout_y = (fig_num)/layout_x
if fig_num%layout_x != 0:
  layout_y = layout_y +1
fig = plt.figure(101, figsize=(2*layout_x,2*layout_y), dpi=150)
fig_layout = gridspec.GridSpec(layout_y,layout_x)
left_end_vr = cent_v - bin_v/2.0 - (fig_num2*bin_v)
channel_images = []
dv = np.abs(coords[11]*1.0e-3)
for fig_idx in range(fig_num):
  v_range = [left_end_vr+fig_idx*bin_v,left_end_vr+(fig_idx+1)*bin_v]
  c_vr = np.mean(v_range)
  print 'start to generate the channel image for '+ str(v_range)
  vr_idx1 = vspace > v_range[0]; vr_idx2 = vspace < v_range[1]
  vr_idx = np.logical_and(vr_idx1, vr_idx2)
  
  channel_cube = data[vr_idx,:,:]
  channel_img = np.sum(channel_cube, axis=0)*dv
  layout_x_id = fig_idx % layout_x; layout_y_id = fig_idx / layout_x
  channel_map = fig.add_subplot(fig_layout[layout_y_id,layout_x_id],projection=wcs)
  plt.imshow(channel_img, origin='lower', vmin=int_min,
             vmax=int_max, cmap='jet', aspect=1.0)
  channel_images = np.append(channel_images,channel_img)
  c_vr_str = "{0:3.1f} km/s".format(c_vr)
  ax = plt.gca()
  plt.text(0.1,0.9,c_vr_str,fontsize=10,transform=ax.transAxes,
           ha = 'left', va = 'top')
  
fig.savefig('channel_map.eps')
print 'maximum integrated intensity in channel map is '+ str(np.max(channel_images))
print 'channel map is stored'
print 'finished'









