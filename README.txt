#######################
## Python DM3 Reader ##
#######################


Author: 
Pierre-Ivan Raynal <raynal@univ-tours.fr>
(http://microscopies.med.univ-tours.fr/)


This Python module is an transposition and adaptation of the DM3_Reader ImageJ plug-in (http://rsb.info.nih.gov/ij/plugins/DM3_Reader.html) by Greg Jefferis <jefferis@stanford.edu>.

It allows to extract thumbnail, image data and metadata (specimen, operator, HV, MAG, etc.) from GATAN DigitalMicrograph 3 files (it was initially meant to be called by a script indexing electron microscope images in a database). In particular, it can dump all metadata (“Tags”), pass thumbnail and image data as PIL Images or Numpy arrays, as well as save the thumbnail view in a PNG file.


Dependencies
------------

    Python Imaging Library
    NumPy
    SciPy

Optional:

	IPython IDE
	Matplotlib


Usage
-----

Use in your own script would typically require lines such as:

	import DM3lib as dm3
	import matplotlib.pyplot as plt
	# parse DM3 file
	dm3f = dm3.DM3("sample.dm3")
	# print some useful image information
	print dm3f.info
	print "pixel size = %s %s"%dm3f.pxsize
	# display image
	plt.matshow(dm3f.imagedata, vmin=dm3f.cuts[0], vmax=dm3f.cuts[1])
	plt.colorbar()


Known Issues
------------

This DM3 parser has only been tested on single-image files (i.e. not stacks) and not all data types are implemented.
