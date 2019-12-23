from setuptools import setup, find_packages
from os.path import join, dirname
import grav
setup(
    name='grav',
    version='0.0.1',
    packages=find_packages(),
    long_description=open(join(dirname(__file__), 'README.md')).read(),
    install_requires=['numpy', 'scipy', 'matplotlib', 'astropy', 'lmfit'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: MIT License',
        'Topic :: Education',
        'Programming Language :: Python :: 3',
    ],
    entry_points={
		'console_scripts':['grav_fly = grav.grav_fly:main']
	},
    test_suite='tests',
    include_package_data=True,
)
