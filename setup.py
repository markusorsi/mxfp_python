from setuptools import setup, find_packages

VERSION = '1.1.2' 
DESCRIPTION = 'MXFP'
LONG_DESCRIPTION = 'Open-source version of MXFP based on the RDKit.'

# Setting up
setup(
        name="mxfp", 
        version=VERSION,
        author="Markus Orsi",
        author_email="<markus.orsi@unibe.ch>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], 
        keywords=['cheminformatics', 'mxfp', 'fingerprint'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)