#!/usr/bin/env python

from setuptools import find_packages, setup

from tmber import __version__

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='tmber',
    version=__version__,
    author='dceoy',
    author_email='dnarsil+github@gmail.com',
    description='Tumor Mutational Burden Analyzer',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/dceoy/tmber.git',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['docopt', 'pandas', 'psutil', 'pyyaml'],
    entry_points={
        'console_scripts': ['tmber=tmber.cli:main']
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development'
    ],
    python_requires='>=3.6',
)
