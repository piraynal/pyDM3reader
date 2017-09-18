from setuptools import setup, find_packages
setup(
    name = "dm3_lib",
    version = "1.5dev",
    packages = ['dm3_lib'],

    install_requires = ['pillow>=2.3.1', 'numpy'],

    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst'],
        'dm3_lib': ['demo/demo.py', 'demo/utilities.py'],
    },

    # metadata for upload to PyPI
    author = "Pierre-Ivan Raynal",
    author_email = "raynal@univ-tours.fr",
    description = "Python module for parsing GATAN DM3|DM4 (DigitalMicrograph) files",
    license = "MIT",
    keywords = "GATAN DigitalMicrograph DM3 DM4 Transmission Electron Microscopy",
    url = "http://microscopies.med.univ-tours.fr/pydm3reader/",   # project home page, if any
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License (MIT)'
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Visualization'],
    # List of classifiers:
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    # could also include long_description, download_url, classifiers, etc.
)
