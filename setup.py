from setuptools import setup, find_packages
from os.path import join, dirname
import grav_fly
setup(
    name='grav_fly',
    version='0.0.1',
    packages=find_packages(),
    long_description=open(join(dirname(__file__), 'README.txt')).read(),
    install_requires=['numpy', 'scipy', 'matplotlib', 'astropy', 'lmfit'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: MIT License',
        'Topic :: Education',
        'Programming Language :: Python :: 3',
    ],
    include_package_data=True,
)
