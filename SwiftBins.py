from astropy.io import ascii
import numpy as np

batdata = ascii.read("BAT_105m_catalog_07jul2019.txt",delimiter='|',names=["ID", "BAT_NAME",     "RA",    "DEC",    "SNR",  "COUNTERPART_NAME",  "OTHER_NAME", "CTPT_RA","CTPT_DEC", "FLUX", "FLUX_LO", "FLUX_HI", "CONTA", "GAMM", "GAMM_LO", "GAMM_HI",   "CHI_SQ_R", "REDSHIFT",  "LUM", "ASSOC_STREN",  "CL2", "TYPE"])

brightsort = np.argsort(batdata['FLUX'])

dirs = open("DirectoryList.txt",'r').readlines()

ra_list = []
dec_list = []



for this_dir in dirs:
    this_dir= this_dir.rstrip()
    this_ra = int(this_dir[:3])
    this_dec = 90 - int(this_dir[3:])

    ra_list.append(this_ra)
    dec_list.append(this_dec)

ra_arr = np.array(ra_list)
dec_arr = np.array(dec_list)


for src_ra, src_dec, src_flux, src_name in zip(batdata['RA'][brightsort], batdata['DEC'][brightsort], batdata['FLUX'][brightsort], batdata['COUNTERPART_NAME'][brightsort]):

    ra_mask1 = ra_arr - src_ra < 1.5
    ra_mask2 = ra_arr - src_ra > -1.5
    ra_mask = np.logical_and(ra_mask1,ra_mask2)
    
    
    dec_mask1 = np.abs(dec_arr - src_dec) < 1.5
    dec_mask2 = -1.*np.abs(dec_arr - src_dec) > -1.5

    dec_mask = np.logical_and(dec_mask1,dec_mask2)
    
    pos_mask = np.logical_and(ra_mask,dec_mask)

    ra_filtered = ra_arr[pos_mask]
    dec_filtered = dec_arr[pos_mask]
    #print(ra_filtered,dec_filtered)
    
    
    if ra_filtered.shape[0] == 1:
        print("%s %03i%03i  %.3f  %.3f   %.3f" %(src_name, ra_filtered[0],90-dec_filtered[0], src_ra, src_dec, src_flux))

