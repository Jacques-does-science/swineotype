from setuptools import setup, find_packages
import glob

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
    data_files=[('share/swineotype/data', glob.glob('data/*'))]
)
