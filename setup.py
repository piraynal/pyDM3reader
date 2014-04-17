from setuptools import setup, find_packages
setup(
    name = "DM3lib",
    version = "1.0",
    packages = ['dm3_lib'],

    install_requires = ['pillow>=2.3.1', 'numpy', 'scipy'],

    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst'],
        '': ['demo.py', 'utilities.py'],
    },

    # metadata for upload to PyPI
    author = "Pierre-Ivan Raynal",
    author_email = "raynal@univ-tours.fr",
    description = "Python module for parsing GATAN DM3 (DigitalMicrograph) files",
    license = "GPLv3",
    keywords = "DigitalMicrograph Gatan Microscopy",
    url = "http://microscopies.med.univ-tours.fr/",   # project home page, if any

    # could also include long_description, download_url, classifiers, etc.
)
