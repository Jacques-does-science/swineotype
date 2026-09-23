import os
import re

from setuptools import setup, find_packages


def read_version() -> str:
    """Single source of truth: swineotype/__init__.py."""
    src = open(os.path.join('swineotype', '__init__.py')).read()
    return re.search(r'^__version__\s*=\s*"([^"]+)"', src, re.M).group(1)


setup(
    name='swineotype',
    version=read_version(),
    packages=find_packages(),
    # The reference data lives inside the package and is found with
    # importlib.resources, so an install needs no data_files/share/ layout.
    package_data={'swineotype': ['data/*']},
    # These were previously declared nowhere, so `pip install swineotype`
    # produced an environment that imported click and yaml and failed. They
    # only appeared to be satisfied because scripts/install_swineotype.sh
    # conda-installs them separately.
    python_requires='>=3.10',
    install_requires=[
        'click>=8.0',
        'PyYAML>=5.4',
    ],
    extras_require={
        'test': ['pytest>=7.0'],
    },
    entry_points={
        'console_scripts': [
            'swineotype = swineotype.main:main',
        ],
    },
)
