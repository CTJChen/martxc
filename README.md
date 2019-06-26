# martxc [![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
martxc is a collection of data analysis and reduction tools for the X-ray telescope
ART-XC. 

## List of tools

martevt2img.py - convert X-ray event list to an image.

martexpmap.py - compute exposure map based on an attitude file.

martpoisson.py - simulate poisson noise based on an exposure map or a given exposure value.


## Installation
In a terminal, under the path in which you want to have the software installed, 
run the following commands:
```bash
git clone https://github.com/CTJChen/martxc
cd martxc
python setup.py install
```
## martxc uses the following python packages:

* numpy, scipy, astropy

## License
The repository is developed by Chien-Ting Chen (C) 2019- at the NASA [MSFC X-ray group](https://wwwastro.msfc.nasa.gov/).

