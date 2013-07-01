#!/usr/bin/python

################################################################################
## Python script for parsing GATAN DM3 (DigitalMicrograph) files
## and extracting various metadata
## --
## warning: *tested on single-image files only*
## --
## based on the DM3_Reader plug-in (v 1.3.4) for ImageJ by Greg Jefferis <jefferis@stanford.edu>
## http://rsb.info.nih.gov/ij/plugins/DM3_Reader.html
## --
## Python adaptation: Pierre-Ivan Raynal <raynal@univ-tours.fr>
## http://microscopies.med.univ-tours.fr/
################################################################################

import sys, os, time
import struct
from PIL import Image
import numpy
import scipy.misc

__all__ = ["DM3","version"]

version='1.0.dev'

debugLevel = 0   # 0=none, 1-3=basic, 4-5=simple, 6-10 verbose


### utility fuctions ###
# Image to Array
def im2ar( im ):
	if im.mode in ('L','I','F'):
	    # Warning: only works with PIL.Image.Image whose mode is 'L', 'I' or 'F'
    	#          => error if mode == 'I;16' for instance
	    a = scipy.misc.fromimage( im )
    	return a
#	else:
#		return False
		
## Array to image file
def ar2imfile(filename, a):
    scipy.misc.imsave(filename, a)


### binary data reading functions ###

def readLong(f):
	'''Read 4 bytes as integer in file f'''
	read_bytes = f.read(4)
	return struct.unpack('>l', read_bytes)[0]

def readShort(f):
	'''Read 2 bytes as integer in file f'''
	read_bytes = f.read(2)
	return struct.unpack('>h', read_bytes)[0]

def readByte(f):
	'''Read 1 byte as integer in file f'''
	read_bytes = f.read(1)
	return struct.unpack('>b', read_bytes)[0]

def readBool(f):
	'''Read 1 byte as boolean in file f'''
	read_val = readByte(f)
	return (read_val!=0)

def readChar(f):
	'''Read 1 byte as char in file f'''
	read_bytes = f.read(1)
	return struct.unpack('c', read_bytes)[0]

def readString(f, len=1):
	'''Read len bytes as a string in file f'''
	read_bytes = f.read(len)
	str_fmt = '>'+str(len)+'s'
	return struct.unpack( str_fmt, read_bytes )[0]

def readLEShort(f):
	'''Read 2 bytes as *little endian* integer in file f'''
	read_bytes = f.read(2)
	return struct.unpack('<h', read_bytes)[0]

def readLELong(f):
	'''Read 4 bytes as *little endian* integer in file f'''
	read_bytes = f.read(4)
	return struct.unpack('<l', read_bytes)[0]

def readLEUShort(f):
	'''Read 2 bytes as *little endian* unsigned integer in file f'''
	read_bytes = f.read(2)
	return struct.unpack('<H', read_bytes)[0]

def readLEULong(f):
	'''Read 4 bytes as *little endian* unsigned integer in file f'''
	read_bytes = f.read(4)
	return struct.unpack('<L', read_bytes)[0]

def readLEFloat(f):
	'''Read 4 bytes as *little endian* float in file f'''
	read_bytes = f.read(4)
	return struct.unpack('<f', read_bytes)[0]

def readLEDouble(f):
	'''Read 8 bytes as *little endian* double in file f'''
	read_bytes = f.read(8)
	return struct.unpack('<d', read_bytes)[0]


## constants for encoded data types ##
SHORT = 2
LONG = 3
USHORT = 4
ULONG = 5
FLOAT = 6
DOUBLE = 7
BOOLEAN = 8
CHAR = 9
OCTET = 10
STRUCT = 15
STRING = 18
ARRAY = 20

# - association data type <--> reading function
readFunc = {
	SHORT: readLEShort,
	LONG: readLELong,
	USHORT: readLEUShort,
	ULONG: readLEULong,
	FLOAT: readLEFloat,
	DOUBLE: readLEDouble,
	BOOLEAN: readBool,
	CHAR: readChar,
	OCTET: readChar,    # difference with char???
}

## list of image DataTypes ##
dataTypes = {
	0:	'NULL_DATA',
	1:	'SIGNED_INT16_DATA',
	2:	'REAL4_DATA',
	3:	'COMPLEX8_DATA',
	4:	'OBSELETE_DATA',
	5:	'PACKED_DATA',
	6:	'UNSIGNED_INT8_DATA',
	7:	'SIGNED_INT32_DATA',
	8:	'RGB_DATA',
	9:	'SIGNED_INT8_DATA',
	10: 'UNSIGNED_INT16_DATA',
	11: 'UNSIGNED_INT32_DATA',
	12: 'REAL8_DATA',
	13: 'COMPLEX16_DATA',
	14: 'BINARY_DATA',
	15: 'RGB_UINT8_0_DATA',
	16: 'RGB_UINT8_1_DATA',
	17: 'RGB_UINT16_DATA',
	18: 'RGB_FLOAT32_DATA',
	19: 'RGB_FLOAT64_DATA',
	20: 'RGBA_UINT8_0_DATA',
	21: 'RGBA_UINT8_1_DATA',
	22: 'RGBA_UINT8_2_DATA',
	23: 'RGBA_UINT8_3_DATA',
	24: 'RGBA_UINT16_DATA',
	25: 'RGBA_FLOAT32_DATA',
	26: 'RGBA_FLOAT64_DATA',
	27: 'POINT2_SINT16_0_DATA',
	28: 'POINT2_SINT16_1_DATA',
	29: 'POINT2_SINT32_0_DATA',
	30: 'POINT2_FLOAT32_0_DATA',
	31: 'RECT_SINT16_1_DATA',
	32: 'RECT_SINT32_1_DATA',
	33: 'RECT_FLOAT32_1_DATA',
	34: 'RECT_FLOAT32_0_DATA',
	35: 'SIGNED_INT64_DATA',
	36: 'UNSIGNED_INT64_DATA',
	37: 'LAST_DATA',
	}

## other constants ##
IMGLIST = "root.ImageList."
OBJLIST = "root.DocumentObjectList."
MAXDEPTH = 64

## END constants ##


class DM3(object):
	## utility functions
	def __makeGroupString(self):
		tString = self.__curGroupAtLevelX[0]
		for i in range( 1, self.__curGroupLevel+1 ):
			tString += '.' + self.__curGroupAtLevelX[i]		
		return tString

	def __makeGroupNameString(self):
		tString = self.__curGroupNameAtLevelX[0]
		for i in range( 1, self.__curGroupLevel+1 ):
			tString += '.' + str( self.__curGroupNameAtLevelX[i] )
		return tString

	def __readTagGroup(self):	
		# go down a level
		self.__curGroupLevel += 1
		# increment group counter
		self.__curGroupAtLevelX[self.__curGroupLevel] += 1
		# set number of current tag to -1 --- readTagEntry() pre-increments => first gets 0
		self.__curTagAtLevelX[self.__curGroupLevel] = -1
		if ( debugLevel > 5):
			print "rTG: Current Group Level:", self.__curGroupLevel
		# is the group sorted?
		sorted = readByte(self.__f)
		isSorted = (sorted == 1)
		# is the group open?
		opened = readByte(self.__f)
		isOpen = (opened == 1)
		# number of Tags
		nTags = readLong(self.__f)
		if ( debugLevel > 5):
			print "rTG: Iterating over the", nTags, "tag entries in this group"
		# read Tags
		for i in range( nTags ):
			self.__readTagEntry()
		# go back up one level as reading group is finished
		self.__curGroupLevel += -1
		return 1

	def	__readTagEntry(self):
		# is data or a new group?
		data = readByte(self.__f)
		isData = (data == 21)
		self.__curTagAtLevelX[self.__curGroupLevel] += 1
		# get tag label if exists
		lenTagLabel = readShort(self.__f)
		if ( lenTagLabel != 0 ):
			tagLabel = readString(self.__f, lenTagLabel)
		else:
			tagLabel = str( self.__curTagAtLevelX[self.__curGroupLevel] )
		if ( debugLevel > 5):
			print str(self.__curGroupLevel)+"|"+__makeGroupString()+": Tag label = "+tagLabel
		elif ( debugLevel > 0 ):
			print str(self.__curGroupLevel)+": Tag label = "+tagLabel
		if isData:
			# give it a name
			self.__curTagName = self.__makeGroupNameString()+"."+tagLabel
			# read it
			self.__readTagType()
		else:
			# it is a tag group
			self.__curGroupNameAtLevelX[self.__curGroupLevel+1] = tagLabel
			self.__readTagGroup()  # increments curGroupLevel
		return 1

	def __readTagType(self):
		delim = readString(self.__f, 4)
		if ( delim != "%%%%" ):
			raise Exception, hex( self.__f.tell() )+": Tag Type delimiter not %%%%"
		nInTag = readLong(self.__f)
		self.__readAnyData()
		return 1

	def __encodedTypeSize(self, eT):
		# returns the size in bytes of the data type	
		if eT == 0:
			width = 0
		elif eT in (BOOLEAN, CHAR, OCTET):
			width = 1
		elif eT in (SHORT, USHORT):
			width = 2
		elif eT in (LONG, ULONG, FLOAT):
			width = 4
		elif eT == DOUBLE:
			width = 8
		else:
			# returns -1 for unrecognised types
			width=-1	
		return width

	def __readAnyData(self):
		## higher level function dispatching to handling data types to other functions	
		# - get Type category (short, long, array...)
		encodedType = readLong(self.__f)
		# - calc size of encodedType
		etSize = self.__encodedTypeSize(encodedType)
		if ( debugLevel > 5):
			print "rAnD, " + hex( f.tell() ) + ": Tag Type = " + str(encodedType) +  ", Tag Size = " + str(etSize)
		if ( etSize > 0 ):
			self.__storeTag( self.__curTagName, self.__readNativeData(encodedType, etSize) )
		elif ( encodedType == STRING ):
			stringSize = readLong(self.__f)
			self.__readStringData(stringSize)
		elif ( encodedType == STRUCT ):
			# does not store tags yet
			structTypes = self.__readStructTypes()
			self.__readStructData(structTypes)
		elif ( encodedType == ARRAY ):
			# does not store tags yet
			# indicates size of skipped data blocks
			arrayTypes = self.__readArrayTypes()
			self.__readArrayData(arrayTypes)
		else:
			raise Exception, "rAnD, " + hex(self.__f.tell()) + ": Can't understand encoded type"		
		return 1
	
	def __readNativeData(self, encodedType, etSize):
		# reads ordinary data types
		if encodedType in readFunc.keys():
			val = readFunc[encodedType](self.__f)
		else:
			raise Exception, "rND, " + hex(self.__f.tell()) + ": Unknown data type " + str(encodedType)		
		if ( debugLevel > 3 ):
			print "rND, " + hex(self.__f.tell()) + ": " + str(val)
		elif ( debugLevel > 0 ):
			print val
		return val

	def __readStringData(self, stringSize):
		# reads string data
		if ( stringSize <= 0 ):
			rString = ""
		else:
			if ( debugLevel > 3 ):
				print "rSD @ " + str(f.tell()) + "/" + hex(f.tell()) +" :",
			## !!! *Unicode* string (UTF-16)... convert to Python unicode str
			rString = readString(self.__f, stringSize)
			rString = unicode(rString, "utf_16_le")
			if ( debugLevel > 3 ):
				print rString + "   <"  + repr( rString ) + ">"
		if ( debugLevel > 0 ):
			print "StringVal:", rString
		self.__storeTag( self.__curTagName, rString )
		return rString
	
	def __readArrayTypes(self):
		# determines the data types in an array data type
		arrayType = readLong(self.__f)
		itemTypes=[]
		if ( arrayType == STRUCT ):
			itemTypes = self.__readStructTypes()
		elif ( arrayType == ARRAY ):
			itemTypes = self.__readArrayTypes()
		else:
			itemTypes.append( arrayType )
		return itemTypes

	def __readArrayData(self, arrayTypes):
		# reads array data
		
		arraySize = readLong(self.__f)
		
		if ( debugLevel > 3 ):
			print "rArD, " + hex( f.tell() ) + ": Reading array of size = " + str(arraySize)
		
		itemSize = 0
		encodedType = 0
		
		for i in range( len(arrayTypes) ):
			encodedType = int( arrayTypes[i] )
			etSize = self.__encodedTypeSize(encodedType)
			itemSize += etSize
			if ( debugLevel > 5 ):
				print "rArD: Tag Type = " + str(encodedType) + ", Tag Size = " + str(etSize)
			##! readNativeData( encodedType, etSize ) !##
		
		if ( debugLevel > 5 ):
			print "rArD: Array Item Size = " + str(itemSize)
	
		bufSize = arraySize * itemSize
		
		if ( (not self.__curTagName.endswith("ImageData.Data"))
				and  ( len(arrayTypes) == 1 )
				and  ( encodedType == USHORT )
				and  ( arraySize < 256 ) ):
			# treat as string
			val = self.__readStringData( bufSize )
		else:
			# treat as binary data
			# - store data size and offset as tags 
			self.__storeTag( self.__curTagName + ".Size", bufSize )
			self.__storeTag( self.__curTagName + ".Offset", self.__f.tell() )
			# - skip data w/o reading
			self.__f.seek( self.__f.tell() + bufSize )
		
		return 1

	def __readStructTypes(self):
		# analyses data types in a struct
		
		if ( debugLevel > 3 ):
			print "Reading Struct Types at Pos = " + hex(self.__f.tell())
	
		structNameLength = readLong(self.__f)
		nFields = readLong(self.__f)
	
		if ( debugLevel > 5 ):
			print "nFields = ", nFields
	
		if ( nFields > 100 ):
			raise Exception, hex(self.__f.tell())+": Too many fields"
			
		fieldTypes = []	
		nameLength = 0
		for i in range( nFields ):
			nameLength = readLong(self.__f)
			if ( debugLevel > 9 ):
				print i + "th namelength = " + nameLength
			fieldType = readLong(self.__f)
			fieldTypes.append( fieldType )
	
		return fieldTypes

	def __readStructData(self, structTypes):
		# reads struct data based on type info in structType
		for i in range( len(structTypes) ):
			encodedType = structTypes[i]
			etSize = self.__encodedTypeSize(encodedType)
	
			if ( debugLevel > 5 ):
				print "Tag Type = " + str(encodedType) + ", Tag Size = " + str(etSize)
	
			# get data
			self.__readNativeData(encodedType, etSize)
		
		return 1
	
	def __storeTag(self, tagName, tagValue):
		# NB: all tag values (and names) stored as unicode objects;
		#     => can then be easily converted to any encoding
		# - /!\ tag names may not be ascii char only (e.g. '\xb5', i.e. MICRO SIGN)
		tagName = unicode(tagName, 'latin-1')
		# - convert tag value to unicode if not already unicode object (as for string data)
		tagValue = unicode(tagValue)
		# store Tags as list and dict
		self.__storedTags.append( tagName + " = " + tagValue )
		self.__tagDict[tagName] = tagValue
		
	### END utility functions ###
	
	def __init__(self, filename, dump=False, dump_dir='/tmp', debug=0):
		'''DM3 object: parses DM3 file and extracts Tags; dumps Tags in a txt file if dump==True.'''
		
		## initialize variables ##
		self.debug = debug
		self.__filename = filename
		self.__chosenImage = 1
		# - track currently read group
		self.__curGroupLevel = -1
		self.__curGroupAtLevelX = [ 0 for x in range(MAXDEPTH) ]
		self.__curGroupNameAtLevelX = [ '' for x in range(MAXDEPTH) ]
		# - track current tag
		self.__curTagAtLevelX = [ '' for x in range(MAXDEPTH) ]
		self.__curTagName = ''
		# - open file for reading
		self.__f = open( self.__filename, 'rb' )
		# - create Tags repositories
		self.__storedTags = []
		self.__tagDict = {}
		
		if self.debug>0:
			t1 = time.time()
		isDM3 = True
		## read header (first 3 4-byte int)
		# get version
		fileVersion = readLong(self.__f)
		if ( fileVersion != 3 ):
			isDM3 = False
		# get indicated file size
		fileSize = readLong(self.__f)
		# get byte-ordering
		lE = readLong(self.__f)
		littleEndian = (lE == 1)
		if not littleEndian:
			isDM3 = False
		# check file header, raise Exception if not DM3
		if not isDM3:
			raise Exception, "%s does not appear to be a DM3 file."%os.path.split(self.__filename)[1]
		elif self.debug > 0:
			print "%s appears to be a DM3 file"%(self.__filename)
			
		if ( debugLevel > 5 or self.debug > 1):
			print "Header info.:"
			print "- file version:", fileVersion
			print "- lE:", lE
			print "- file size:", fileSize, "bytes"

		# set name of root group (contains all data)...
		self.__curGroupNameAtLevelX[0] = "root"
		# ... then read it
		self.__readTagGroup()
		if self.debug > 0:
			print "-- %s Tags read --"%len(self.__storedTags)
		
		if self.debug>0:
			t2 = time.time()
			print "| parse DM3 file: %.3g s"%(t2-t1)
				
		# dump Tags in txt file if requested
		if dump:
			dump_file = os.path.join(dump_dir, os.path.split(self.__filename)[1]+".tagdump.txt")
			try:
				dumpf = open( dump_file, 'w' )
			except:
				print "Warning: cannot generate dump file."
			else:
				for tag in self.__storedTags:
					dumpf.write( tag.encode('latin-1') + "\n" )
				dumpf.close

	def getFilename(self):
		return self.__filename
	filename = property(getFilename)

	def getTags(self):
		return self.__tagDict
	tags = property(getTags)

	def getInfo(self, info_charset='latin1'):
		'''Extracts useful experiment info from DM3 file and 
		exports thumbnail to a PNG file if 'make_tn' set to 'True'.'''
		
		# define useful information
		info_keys = {
			'descrip': 'root.ImageList.1.Description',
			'acq_date': 'root.ImageList.1.ImageTags.DataBar.Acquisition Date',
			'acq_time': 'root.ImageList.1.ImageTags.DataBar.Acquisition Time',
			'name': 'root.ImageList.1.ImageTags.Microscope Info.Name',
			'micro': 'root.ImageList.1.ImageTags.Microscope Info.Microscope',
			'hv': 'root.ImageList.1.ImageTags.Microscope Info.Voltage',
			'mag': 'root.ImageList.1.ImageTags.Microscope Info.Indicated Magnification',
			'mode': 'root.ImageList.1.ImageTags.Microscope Info.Operation Mode',
			'operator': 'root.ImageList.1.ImageTags.Microscope Info.Operator',
			'specimen': 'root.ImageList.1.ImageTags.Microscope Info.Specimen',
		#	'image_notes': 'root.DocumentObjectList.10.Text' # = Image Notes 		
			}

		# get experiment information
		infoDict = {}
		for key, tag_name in info_keys.items():
			if self.tags.has_key(tag_name):
				# tags supplied as Python unicode str; convert to chosen charset (typ. latin-1)
				infoDict[key] = self.tags[tag_name].encode(info_charset)

		# return experiment information
		return infoDict

	info = property(getInfo)
	
	def getThumbnail(self, asDict=False):
		'''Returns thumbnail as Image or dict.'''
		# get thumbnail
		tn_size = int( self.tags['root.ImageList.0.ImageData.Data.Size'] )
		tn_offset = int( self.tags['root.ImageList.0.ImageData.Data.Offset'] )
		tn_width = int( self.tags['root.ImageList.0.ImageData.Dimensions.0'] )
		tn_height = int( self.tags['root.ImageList.0.ImageData.Dimensions.1'] )
			
		if self.debug > 0:
			print "Notice: tn data in %s starts at %s"%(os.path.split(self.__filename)[1], hex(tn_offset))
			print "Notice: tn size: %sx%s px"%(tn_width,tn_height)
			
		sizeError = False
		if (tn_width*tn_height*4) != tn_size:
			raise Exception, "Cannot extract thumbnail from %s"%os.path.split(self.__filename)[1]
		else:
			self.__f.seek( tn_offset )			
			rawdata = self.__f.read(tn_size)	
			# - read as 16-bit LE unsigned integer
			tn = Image.fromstring( 'F', (tn_width,tn_height), rawdata, 'raw', 'F;32' )
			# - rescale and convert px data
			tn = tn.point(lambda x: x * (1./65536) + 0)
			tn = tn.convert('L')
		
		if asDict:
			# - fill tnDict
			tnDict = {}
			tnDict['size'] = tn.size
			tnDict['mode'] = tn.mode
			tnDict['rawdata'] = tn.tostring()
			return tnDict
		else:
			return tn

	thumbnail = property(getThumbnail)

	def getThumbnailData(self):
		'''Returns thumbnail data as numpy.array'''
		return im2ar(self.thumbnail)

	thumbnaildata = property(getThumbnailData)
				
	def makePNGThumbnail(self, tn_file=''):
		'''Save thumbnail as PNG file.'''
		# - cleanup name
		if tn_file == '':
			tn_path = os.path.join('./', os.path.split(self.filename)[1]+'.tn.png')
		else:
			if os.path.splitext(tn_file)[1] != '.png':
				tn_path = os.path.splitext(tn_file)[0] + '.png'
			else:
				tn_path = tn_file
		# - save tn file
		try:
			self.thumbnail.save(tn_path, 'PNG')
			if self.debug > 0:
				print "Thumbnail saved as '%s'."%tn_path
		except:
			print "Warning: could not save thumbnail."

	def getImage(self):
		'''Extracts image data as Image'''
		
		# PIL "raw" decoder modes for the various image dataTypes
		dataTypesDec = {
			1: 'F;16S',    #16-bit LE signed integer
			2: 'F;32F',    #32-bit LE floating point
			6: 'F;8',      #8-bit unsigned integer
			7: 'F;32S',    #32-bit LE signed integer
			9: 'F;8S',     #8-bit signed integer
			10: 'F;16',    #16-bit LE unsigned integer
			11: 'F;32',    #32-bit LE unsigned integer
			#14: 'F;8',     #binary
			}
		
		# get relevant Tags			
		data_offset = int( self.tags['root.ImageList.1.ImageData.Data.Offset'] )
		data_size = int( self.tags['root.ImageList.1.ImageData.Data.Size'] )
		data_type = int( self.tags['root.ImageList.1.ImageData.DataType'] )
		im_width = int( self.tags['root.ImageList.1.ImageData.Dimensions.0'] )
		im_height = int( self.tags['root.ImageList.1.ImageData.Dimensions.1'] )

		if self.debug>0:
			print "Notice: image data in %s starts at %s"%(os.path.split(self.__filename)[1], hex(data_offset))
			print "Notice: image size: %sx%s px"%(im_width,im_height)
						
		# check if image DataType is implemented, then read
		if data_type in dataTypesDec.keys():
			decoder = dataTypesDec[data_type]
			if self.debug>0:
				print "Notice: image data read as %s"%decoder
				t1 = time.time()
			self.__f.seek( data_offset )
			rawdata = self.__f.read(data_size)
			im = Image.fromstring( 'F', (im_width,im_height), rawdata, 'raw', decoder )
			if self.debug>0:
				t2 = time.time()
				print "| read image data: %.3g s"%(t2-t1)
		else:	
			raise Exception, "Cannot extract image data from %s: unimplemented DataType (%s:%s)."%(os.path.split(self.__filename)[1],data_type,dataTypes[data_type])
			
		return im
		
	image = property(getImage)

	def getImageData(self):
		'''Extracts image data as numpy.array'''
		return im2ar(self.image)

	imagedata = property(getImageData)
	
	def getDisplayCuts(self):
		'''Returns display level limits.'''
		display_min = int( float( self.tags['root.DocumentObjectList.0.ImageDisplayInfo.LowLimit'] ) )
		display_max = int( float( self.tags['root.DocumentObjectList.0.ImageDisplayInfo.HighLimit'] ) )
		cuts = (display_min, display_max)
		return cuts
	
	cuts = property(getDisplayCuts)

	def getPixelSize(self):
		'''Returns pixel size and unit.'''
		pixel_size = float( self.tags['root.ImageList.1.ImageData.Calibrations.Dimension.0.Scale'] )
		unit = self.tags['root.ImageList.1.ImageData.Calibrations.Dimension.0.Units']
		if unit == u'\xb5m':
			unit = 'micron'
		else:
			unit = unit.encode('ascii')
		if self.debug>0:
			print "pixel size = %s %s"%(pixel_size,unit)		
		return (pixel_size,unit)

	pxsize = property(getPixelSize)
			
if __name__ == '__main__':
	print "DM3lib v.%s"%version

