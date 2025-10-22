from setuptools import setup, find_packages
import glob
import os

data_files = [f for f in glob.glob('data/*') if os.path.isfile(f)]

setup(
    name='swineotype',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[],
    entry_points={
        'console_scripts': [
            'swineotype = swineotype.main:main',
        ],
    },
    data_files=[('share/swineotype/data', data_files)]
)
