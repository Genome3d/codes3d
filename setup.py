#!/usr/bin/env python

from distutils.core import setup

setup(name='CoDeS3D',
	  version='0.1',
	  description='Contextualizing Developmental SNPs in 3 Dimensions',
	  author='Cameron Ekblad',
	  author_email='cekb635@aucklanduni.ac.nz',
	  url='https://github.com/alcamerone/codes3d',
	  packages=['codes3d'],
	  install_requires=[
	  	'configobj',
	  	'pandas',
	  	'pybedtools',
	  	'requests'
	  ],
	  dependency_links=['https://github.com/wikipathways/wikipathways-api-client-py',
	  				    'https://github.com/matplotlib/matplotlib']

	 )