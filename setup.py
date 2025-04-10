import re
from setuptools import setup, find_packages

with open('README.md') as fh:
    long_description = fh.read()

with open('sarand/__init__.py') as fh:
    info = fh.read()
    version = re.search('^__version__\s*=\s*"(.*)"', info, re.M).group(1)

setup(
    name='sarand',
    packages=find_packages(),
    version=version,
    license="GPLv3",
    description="Tool to extract the neighborhood of the target Antimicrobial Resistance (AMR) genes from the assembly graph.",
    author='Somayeh Kafaie',
    author_email='so.kafaie@gmail.com',
    url="https://github.com/beiko-lab/sarand",
    keywords=["Metagenomic Assembly graph", "Antimicrobial resistance", "Context extraction"],
    python_requires='>=3.7',
    long_description=long_description,
    long_description_content_type="text/markdown",
    # what about the ones in the requirement file?????? should I include them here?????
    install_requires=['dna_features_viewer', 'numpy', 'pillow', 'matplotlib', 'gfapy', 'pandas', 'biopython'],
    # not sure how to set this correctly???????????????????
    # package_data={
    #     'sarand': ['data/CARD_AMR_seq.fasta'],
    # },
    # data_files=[('data', ['data/CARD_AMR_seq.fasta']), ('test', ['test/metagenome_data/Ecoli_NC_010488.fna',\
    #                 'test/metagenome_data/klebsiella_NC_009650.fna', 'test/metagenome_data/staphylo_NC_002758.fna'])],
    # How to include test files???? similar to the above or in manifest.in
    include_package_data=True,
    package_data={'': ['data/*.fasta']},
    entry_points={
        'console_scripts': [
            'sarand = sarand.__main__:main'
        ],
    },
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    zip_safe=True,
)
