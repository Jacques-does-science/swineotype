import os
import re

from setuptools import setup, find_packages

# Ship exactly the reference data the tool reads at runtime -- including the
# resolver manifest, which run records cite for the reference-data version.
# A bare data/* glob also swept up stray indexes and scratch files.
DATA_FILES = [
    'data/suis_wzxwzy_whitelist.fasta',
    'data/suis_resolver_refs.fasta',
    'data/suis_resolver_refs.manifest.json',
]
data_files = [f for f in DATA_FILES if os.path.isfile(f)]
missing = [f for f in DATA_FILES if not os.path.isfile(f)]
if missing:
    raise SystemExit(f'missing reference data: {missing}')


def read_version() -> str:
    """Single source of truth: swineotype/__init__.py."""
    src = open(os.path.join('swineotype', '__init__.py')).read()
    return re.search(r'^__version__\s*=\s*"([^"]+)"', src, re.M).group(1)


setup(
    name='swineotype',
    version=read_version(),
    packages=find_packages(),
    # These were previously declared nowhere, so `pip install swineotype`
    # produced an environment that imported click/pandas/yaml and failed. They
    # only appeared to be satisfied because scripts/install_swineotype.sh
    # conda-installs them separately.
    python_requires='>=3.10',
    install_requires=[
        'click>=8.0',
        'pandas>=1.3',
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
    data_files=[('share/swineotype/data', data_files)]
)
