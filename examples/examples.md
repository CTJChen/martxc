## Example for using the martxc codes
Prerequisite - please install martxc following the README.md file's instructions. 
Please also download the sample event file and attitude file.

Two different attitude files are included, first is a single-pointing, 100ks exposure pointing at the center of North Ecliptic Pole. The second attitude file is a slice of the erosita all-sky survey based on the [sample attitude file](https://www.sternwarte.uni-erlangen.de/research/sixte/data/eRASS_Pc87M55_3dobi_att_remeis.fits.bz2) provided by the eROSITA team.

Several event files, simulated using the [SIXTE](https://www.sternwarte.uni-erlangen.de/research/sixte/) simulation software are included. Most of the file names should be self-explanatory. For any questions, please raise a github issue.

### Example 1 - make an exposure map using only the attitude file. 
This example uses martexpmap.py to generate an exposure map. This would be useful for finding out the depth of ART-XC all-sky survey coverage of a certain sky region. Since the PSF half-power-diameter of ART-XC is ~45 arcsec, the exposuer map is calculated using a square pixel with 45.0 arcsec on each side.

```bash
python martexpmap.py -att ./examples/att_survey_nep_4deg_by_4deg.fits -out expmap_survey.fits -fov 18.0 -vig artxc_vignetting.fits -rasize 45.0 -decsize 45.0 -time 100000
```
This would make an exposure map showing the survey coverage of the first 100 ks.
Note that the field-of-view radius of ART-XC is 18 arcmin, which has to be manually set for now. 
Realistic detector mask for each module would be included in the future releases. 

### Example 2 - make an image and an exposure map from a simulated event list.
Here we will use the event list simulated based on observing a single luminous source (5-30 keV flux ~ 1e-9 erg/s/cm^2).
After getting the image, we can then compute the exposure map for each image pixel. This could be used to obtain the count rates for each pixel.

```bash
python martevt2img.py -evt ./examples/evt_single_source_5ks.fits -out img_single_source_5ks.fits -rasize 45.0 -decsize 45.0 -box 268.0 272.0 63.5 68.5
python martexpmap.py -att ./examples/att_60ks_4_pointings.fits -out exp_4_pointings_5ks.fits -fov 18.0 -img img_single_source_5ks.fits -time 5000.0 -vig artxc_vignetting.fits
```

For simplicity, we did not use SIXTE's background modules when simulating our event lists. Alternatively, we simulate expected background count rate per pixel given our pixel size, using the NuSTAR background spectra, ART-XC response matricies, and ART-XC PSF. To the first order, ART-XC's background count-rate per pixel is 1e-3 count/s. We can use this information to add Poissonian noise to the image:
```bash
python martpoisson.py -img img_single_source_5ks.fits -out img_noise_single_source_5ks.fits -combine True -nrate 1e-3 -exp exp_4_pointings_5ks.fits 
```
We can also define a DS9 region on the source and obtain the auxillary response matrix file (ARF) based on the ART-XC on-axis ARF:
```bash
python martmkarf.py -arf art-xc_v0.0.arf -vig artxc_vignetting.fits -psf artxc_psf_eef.fits -out offaxis_arf.arf -img img_single_source_5ks.fits -region ds9.reg
