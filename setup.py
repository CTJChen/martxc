try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

long_description = ""
with open('README.rst') as f:
    long_description = f.read()


setup(name='martxc',
      version='0.1',
      author='Chien-Ting Chen',
      author_email='ct.chen@nasa.gov',
      license='MIT',
      description='A collection of ART-XC data analysis and reduction scripts',
      packages=['martxc'],
      scripts=['martxcexp.py'],
      description='ART-XC exposuremap script',
      long_description=long_description,
      install_requires=[
          "scipy",
          "numpy",
          "astropy",
          "tqdm",
          "argparse",
          "joblib",
          "healpy",
      ],
      )
