import setuptools
import re
from distutils.core import setup

with open('README.md') as fh:
    long_description = fh.read()

with open('sarand/__init__.py') as fh:
    info = fh.read()
    version = re.search('^__version__\s*=\s*"(.*)"',
                        info, re.M).group(1)

setup(
    name='sarand',
    packages=['sarand'],
    version=version,
    license="GPLv3",
    description="Tool to extract the neighborhood of the target Antimicrobial Resistance (AMR) genes from the assembly graph.",
    author='Somayeh Kafaie',
    author_email='so.kafaie@gmail.com',
    url="https://github.com/somayeh-aut/AMR_context",
    #download_url=f"https://github.com/pha4ge/archive/v{version}.tar.gz",
    keywords=["Metagenomic Assembly graph", "Antimicrobial resistance", "Context extraction"],
    python_requires='>=3.6.10',
    long_description=long_description,
    long_description_content_type="text/markdown",
    #what about the ones in the requirement file?????? should I include them here?????
    install_requires=['dna_features_viewer', 'numpy',
                      'matplotlib','gfapy', 'pandas'],
    #not sure how to set this correctly???????????????????
    # package_data={
    #     'sarand': ['data/CARD_AMR_seq.fasta'],
    # },
    # data_files=[('data', ['data/CARD_AMR_seq.fasta']), ('test', ['test/metagenome_data/Ecoli_NC_010488.fna',\
    #                 'test/metagenome_data/klebsiella_NC_009650.fna', 'test/metagenome_data/staphylo_NC_002758.fna'])],
    #How to include test files???? similar to the above or in manifest.in
    include_package_data=True,
    package_data={'': ['data/*.fasta']},
    entry_points={
        'console_scripts': [
            'sarand = sarand.sarand_main:main'
            ],
        },
    classifiers=[
        #"Development Status :: 5 - Production/Stable",
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    zip_safe=True,
    )
