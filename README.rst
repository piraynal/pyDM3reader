=================
Python DM3 Reader
=================


Author:
`Pierre-Ivan Raynal <mailto:raynal@univ-tours.fr>`_
(`https://microscopies.med.univ-tours.fr/
<https://microscopies.med.univ-tours.fr/>`_)

Contributions by
`Philippe Mallet-Ladeira <mailto:philippe.mallet@cemes.fr>`_, 
`Jordan Chess <mailto:jchess@uoregon.edu>`_
and `Chris Harris <mailto:chris.harris@kitware.com>`_

This Python module is an adaptation and extension of the `DM3_Reader`
`ImageJ plug-in <http://rsb.info.nih.gov/ij/plugins/DM3_Reader.html>`_ by
`Greg Jefferis <mailto:jefferis@gmail.com>`_.

It allows to extract image data, metadata (specimen, operator, HV, MAG, etc.)
and thumbnail from `GATAN DigitalMicrograph 3` and `4` files. In particular,
it can dump all metadata (“Tags”), pass thumbnail and image data as PIL Images
or Numpy arrays, as well as save the included thumbnail in a PNG file.
(It was initially meant to be called by a script indexing electron microscope
images in a database.)


Dependencies
============

 - Python 3
 - Numpy
 - Pillow (fork of Python Imaging Library)

Optional Dependencies
---------------------

 - Matplotlib
 - SciPy 
 - IPython IDE


Usage
=====

Use in your own script would typically require lines such as::

    import dm3_lib as dm3
    import matplotlib.pyplot as plt
    # parse DM3 file
    dm3f = dm3.DM3("sample.dm3")
    # print some useful image information
    print dm3f.info
    print "pixel size = %s %s"%dm3f.pxsize
    # display image
    plt.ion()
    plt.matshow(dm3f.imagedata, vmin=dm3f.cuts[0], vmax=dm3f.cuts[1])
    plt.colorbar(shrink=.8)

A more detailed example is located in the ``site-packages/dm3_lib/demo`` directory
under the name ``demo.py``.

Known Issues
============

Not all data types are implemented yet.
