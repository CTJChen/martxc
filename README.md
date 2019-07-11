# martxc [![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
martxc is a collection of data analysis and reduction tools for the X-ray telescope
[ART-XC](https://wwwastro.msfc.nasa.gov/artxc/). 

## Installation and updates
In a terminal, under the path in which you want to have the software installed, 
run the following commands:
```bash
git clone https://github.com/CTJChen/martxc
cd martxc
python setup.py install
```

## Instructions and list of scripts

The current version works with "level 2" event lists generated with [SIXTE](https://www.sternwarte.uni-erlangen.de/research/sixte/). 

martevt2img.py - converts an event list fits file to an image with designated RA and DEC pixel size and optional RA/DEC range parameters. 

martexpmap.py - computes an exposure map based on an attitude file with designated RA and DEC pixel size and optional RA/DEC range parameters. Also works with products of martevt2img.py to make expsure maps with shared fits header keywords. 

martmkarf.py - generate aperture and off-axis corrected auxillary response file.

martpoisson.py - a simplified script for adding background noise to a fits image. Requires prior knowledge on the background count rate per pixel (not a part of martxc).

The python scripts should be executable, you can add the martxc path to your system PATH environment if necessary.
```bash
export PATH=/PATH_TO_MARTXC/:$PATH
```
After doing so, the following command will display the help doc of martevt2img.py:
```bash
martevt2img.py -h
```

See the [example document](https://github.com/CTJChen/martxc/blob/master/examples/examples.md) for further instructions.

For any questions/issues/bugs, please report on this github repository or send an email to ctchen dot git @ gmail.

## martxc uses the following python packages:

* python 2.7+ or 3.2+ (with argparse included in the standard library).
* numpy, scipy, astropy

## License
The is a free, open-source software under the AGPLv3.0 license. martxc is developed and maintained by [Chien-Ting Chen](https://github.com/CTJChen) at the NASA [MSFC X-ray group](https://wwwastro.msfc.nasa.gov/).

