try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

long_description = ""
with open('README.md') as f:
    long_description = f.read()


setup(name='martxc_local',
      version='0.1',
      author='Chien-Ting Chen',
      author_email='ct.chen@nasa.gov',
      license='MIT',
      description='A collection of ART-XC data analysis and reduction scripts',
      packages=['martxclib'],
      scripts=['martexpmap.py', 'martevt2img.py', 'martspecextract.py','martspecextract_raw.py',\
      'martpoisson.py', 'martmkarf.py', 'martfilter.py', 'martpha2pi.py', 'marteef.py', 
      'martmerge.py', 'martfindtile.py'],
      long_description=long_description,
      install_requires=[
          "scipy",
          "numpy",
          "astropy>=3.1",
      ],
      )
