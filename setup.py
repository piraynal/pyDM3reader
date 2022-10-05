from setuptools import setup, find_packages

setup(
    name = "dm3_lib",
    version = "2.0beta",
    packages = ['dm3_lib'],

    install_requires = ['pillow>=2.3.1', 'numpy'],

    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst'],
        'dm3_lib': ['demo/demo.py', 'demo/utilities.py'],
    },

    author = "Pierre-Ivan Raynal",
    author_email = "raynal@univ-tours.fr",
    description = "Python module for parsing GATAN DM3|DM4 (DigitalMicrograph) files",
    license = "MIT",
    keywords = "GATAN DigitalMicrograph DM3 DM4 Transmission Electron Microscopy",
    url = "https://microscopies.med.univ-tours.fr/pydm3reader/",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License (MIT)'
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Visualization'],
)
