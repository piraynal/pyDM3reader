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
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Visualization'],
    # List of classifiers:
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    # could also include long_description, download_url, classifiers, etc.
)
