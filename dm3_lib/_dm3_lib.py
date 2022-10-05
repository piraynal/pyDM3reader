#!/usr/bin/env python3
"""Python module for parsing GATAN DM3 and DM4 files"""

################################################################################
## Python script for parsing GATAN DM3/DM4 (DigitalMicrograph) files
## --
## based on the DM3_Reader plug-in for ImageJ
## by Greg Jefferis <jefferis@gmail.com>
## https://imagej.nih.gov/ij/plugins/DM3_Reader.html
## --
## Python version: Pierre-Ivan Raynal <raynal@univ-tours.fr>
## https://microscopies.med.univ-tours.fr/
################################################################################

import os.path
import struct
import numpy
from PIL import Image

__all__ = ["DM3", "VERSION", "SUPPORTED_DATA_TYPES"]

VERSION = '2.0dev'

debugLevel = 0   # 0=none, 1-3=basic, 4-5=simple, 6-10 verbose

## 'str' function renamed for readability
unicode_str = str

### utility fuctions ###

### binary data reading functions ###

def readByte(f):
    """Read 1 byte as integer in file f"""
    read_bytes = f.read(1)
    return struct.unpack('>b', read_bytes)[0]

def readShort(f):
    """Read 2 bytes as BE integer in file f"""
    read_bytes = f.read(2)
    return struct.unpack('>h', read_bytes)[0]

def readLong(f):
    """Read 4 bytes as BE integer in file f"""
    read_bytes = f.read(4)
    return struct.unpack('>l', read_bytes)[0]

def readLongLong(f):
    """Read 8 bytes as BE integer in file f"""
    read_bytes = f.read(8)
    return struct.unpack('>q', read_bytes)[0]

def readBool(f):
    """Read 1 byte as boolean in file f"""
    read_val = readByte(f)
    return (read_val!=0)

def readChar(f):
    """Read 1 byte as char in file f"""
    read_bytes = f.read(1)
    return struct.unpack('c', read_bytes)[0]

def readString(f, len_=1):
    """Read len_ bytes as a string in file f"""
    read_bytes = f.read(len_)
    str_fmt = '>'+str(len_)+'s'
    return struct.unpack( str_fmt, read_bytes )[0]

def readLEShort(f):
    """Read 2 bytes as *little endian* integer in file f"""
    read_bytes = f.read(2)
    return struct.unpack('<h', read_bytes)[0]

def readLELong(f):
    """Read 4 bytes as *little endian* integer in file f"""
    read_bytes = f.read(4)
    return struct.unpack('<l', read_bytes)[0]

def readLELongLong(f):
    """Read 8 bytes as *little endian* integer in file f"""
    read_bytes = f.read(8)
    return struct.unpack('<q', read_bytes)[0]

def readLEUShort(f):
    """Read 2 bytes as *little endian* unsigned integer in file f"""
    read_bytes = f.read(2)
    return struct.unpack('<H', read_bytes)[0]

def readLEULong(f):
    """Read 4 bytes as *little endian* unsigned integer in file f"""
    read_bytes = f.read(4)
    return struct.unpack('<L', read_bytes)[0]

def readLEULongLong(f):
    """Read 8 bytes as *little endian* unsigned integer in file f"""
    read_bytes = f.read(8)
    return struct.unpack('<Q', read_bytes)[0]

def readLEFloat(f):
    """Read 4 bytes as *little endian* float in file f"""
    read_bytes = f.read(4)
    return struct.unpack('<f', read_bytes)[0]

def readLEDouble(f):
    """Read 8 bytes as *little endian* double in file f"""
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
LONGLONG = 11    # DM4 only
BELONGLONG = 12    # DM4 only
STRUCT = 15
STRING = 18
ARRAY = 20

# - association data type <--> reading function
readFunc = {
    SHORT: readLEShort,
    LONG: readLELong,
    LONGLONG: readLELongLong,    # DM4
    BELONGLONG: readLongLong,    # DM4
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
    0:  'NULL_DATA',
    1:  'SIGNED_INT16_DATA',
    2:  'REAL4_DATA',
    3:  'COMPLEX8_DATA',
    4:  'OBSELETE_DATA',
    5:  'PACKED_DATA',
    6:  'UNSIGNED_INT8_DATA',
    7:  'SIGNED_INT32_DATA',
    8:  'RGB_DATA',
    9:  'SIGNED_INT8_DATA',
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

## supported Data Types
dT_supported = [1, 2, 6, 7, 9, 10, 11, 14]
SUPPORTED_DATA_TYPES = {i: dataTypes[i] for i in dT_supported}

## other constants ##
IMGLIST = "root.ImageList."
OBJLIST = "root.DocumentObjectList."
MAXDEPTH = 64

DEFAULTCHARSET = 'utf-8'
## END constants ##


class DM3(object):
    """DM3 object. """

    ## utility functions
    def _makeGroupString(self):
        tString = str(self._curGroupAtLevelX[0])
        for i in range( 1, self._curGroupLevel+1 ):
            tString += '.{}'.format(self._curGroupAtLevelX[i])
        return tString

    def _makeGroupNameString(self):
        tString = self._curGroupNameAtLevelX[0]
        for i in range( 1, self._curGroupLevel+1 ):
            tString += '.' + str( self._curGroupNameAtLevelX[i] )
        return tString

    def _readIntValue(self):
        if (self._fileVersion == 4):
            Val = readLongLong(self._f)
        else:
            Val = readLong(self._f)
        return Val

    def _readTagGroup(self):
        # go down a level
        self._curGroupLevel += 1
        # increment group counter
        self._curGroupAtLevelX[self._curGroupLevel] += 1
        # set number of current tag to -1
        # --- readTagEntry() pre-increments => first gets 0
        self._curTagAtLevelX[self._curGroupLevel] = -1
        if ( debugLevel > 5):
            print("rTG: Current Group Level:", self._curGroupLevel)
        # is the group sorted?
        sorted_ = readByte(self._f)
        isSorted = (sorted_ == 1)
        # is the group open?
        opened = readByte(self._f)
        isOpen = (opened == 1)
        # number of Tags
        nTags = self._readIntValue()
        if ( debugLevel > 5):
            print("rTG: Iterating over the", nTags, "tag entries in this group")
        # read Tags
        for i in range( nTags ):
            self._readTagEntry()
        # go back up one level as reading group is finished
        self._curGroupLevel += -1
        return 1

    def _readTagEntry(self):
        # is data or a new group?
        data = readByte(self._f)
        isData = (data == 21)
        self._curTagAtLevelX[self._curGroupLevel] += 1
        # get tag label if exists
        lenTagLabel = readShort(self._f)
        if ( lenTagLabel != 0 ):
            tagLabel = readString(self._f, lenTagLabel).decode('latin-1')
        else:
            tagLabel = str( self._curTagAtLevelX[self._curGroupLevel] )
        if ( debugLevel > 5):
            print("{}|{}:".format(self._curGroupLevel, self._makeGroupString()),
                  end=' ')
            print("Tag label = "+tagLabel)
        elif ( debugLevel > 1 ):
            print(str(self._curGroupLevel)+": Tag label = "+tagLabel)
        # if DM4 file, get tag data size
        if (self._fileVersion == 4):
            lenTagData = readLongLong(self._f)
            if ( debugLevel > 1 ):
                print(str(self._curGroupLevel)+": Tag data size = "+str(lenTagData)+" bytes")
        if isData:
            # give it a name
            self._curTagName = self._makeGroupNameString()+"."+tagLabel
            # read it
            self._readTagType()
        else:
            # it is a tag group
            self._curGroupNameAtLevelX[self._curGroupLevel+1] = tagLabel
            self._readTagGroup()  # increments curGroupLevel
        return 1

    def _readTagType(self):
        delim = readString(self._f, 4).decode('latin-1')
        if ( delim != '%%%%' ):
            raise Exception(hex( self._f.tell() )
                            + ": Tag Type delimiter not %%%%")
        nInTag = self._readIntValue()
        self._readAnyData()
        return 1

    def _encodedTypeSize(self, eT):
        # returns the size in bytes of the data type
        if eT == 0:
            width = 0
        elif eT in (BOOLEAN, CHAR, OCTET):
            width = 1
        elif eT in (SHORT, USHORT):
            width = 2
        elif eT in (LONG, ULONG, FLOAT):
            width = 4
        elif eT in (DOUBLE, LONGLONG, BELONGLONG):
            width = 8
        else:
            # returns -1 for unrecognised types
            width = -1
        return width

    def _readAnyData(self):
        ## higher level function dispatching to handling data types
        ## to other functions
        # - get Type category (short, long, array...)
        encodedType = self._readIntValue()
        # - calc size of encodedType
        etSize = self._encodedTypeSize(encodedType)
        if ( debugLevel > 5):
            print("rAnD, " + hex( self._f.tell() ) + ":", end=' ')
            print("Tag Type = " + str(encodedType) + ",", end=' ')
            print("Tag Size = " + str(etSize))
        if ( etSize > 0 ):
            self._storeTag( self._curTagName,
                            self._readNativeData(encodedType, etSize) )
        elif ( encodedType == STRING ):
            stringSize = self._readIntValue()
            self._readStringData(stringSize)
        elif ( encodedType == STRUCT ):
            # does not store tags yet
            structTypes = self._readStructTypes()
            self._readStructData(structTypes)
        elif ( encodedType == ARRAY ):
            # does not store tags yet
            # indicates size of skipped data blocks
            arrayTypes = self._readArrayTypes()
            self._readArrayData(arrayTypes)
        else:
            raise Exception("rAnD, " + hex(self._f.tell())
                            + ": Can't understand encoded type")
        return 1

    def _readNativeData(self, encodedType, etSize):
        # reads ordinary data types
        if encodedType in readFunc:
            val = readFunc[encodedType](self._f)
        else:
            raise Exception("rND, " + hex(self._f.tell())
                            + ": Unknown data type " + str(encodedType))
        if ( debugLevel > 3 ):
            print("rND, " + hex(self._f.tell()) + ": " + str(val))
        elif ( debugLevel > 1 ):
            print(val)
        return val

    def _readStringData(self, stringSize):
        # reads string data
        if ( stringSize <= 0 ):
            rString = ""
        else:
            if ( debugLevel > 3 ):
                print("rSD @ " + str(self._f.tell()) + "/" + hex(self._f.tell()) +" :", end=' ')
            rString = readString(self._f, stringSize)
            # /!\ UTF-16 unicode string => convert to Python unicode str
            rString = rString.decode('utf-16-le')
            if ( debugLevel > 3 ):
                print(rString + "   <"  + repr( rString ) + ">")
        if ( debugLevel > 1 ):
            print("StringVal:", rString)
        self._storeTag( self._curTagName, rString )
        return rString

    def _readArrayTypes(self):
        # determines the data types in an array data type
        arrayType = self._readIntValue()
        itemTypes = []
        if ( arrayType == STRUCT ):
            itemTypes = self._readStructTypes()
        elif ( arrayType == ARRAY ):
            itemTypes = self._readArrayTypes()
        else:
            itemTypes.append( arrayType )
        return itemTypes

    def _readArrayData(self, arrayTypes):
        # reads array data
        arraySize = self._readIntValue()
        
        if ( debugLevel > 3 ):
            print("rArD, " + hex( self._f.tell() ) + ":", end=' ')
            print("Reading array of size = " + str(arraySize))

        itemSize = 0
        encodedType = 0

        for i in range( len(arrayTypes) ):
            encodedType = int( arrayTypes[i] )
            etSize = self._encodedTypeSize(encodedType)
            itemSize += etSize
            if ( debugLevel > 5 ):
                print("rArD: Tag Type = " + str(encodedType) + ",", end=' ')
                print("Tag Size = " + str(etSize))
            ##! readNativeData( encodedType, etSize ) !##

        if ( debugLevel > 5 ):
            print("rArD: Array Item Size = " + str(itemSize))

        bufSize = arraySize * itemSize

        if ( (not self._curTagName.endswith("ImageData.Data"))
                and  ( len(arrayTypes) == 1 )
                and  ( encodedType == USHORT )
                and  ( arraySize < 256 ) ):
            # treat as string
            val = self._readStringData( bufSize )
        else:
            # treat as binary data
            # - store data size and offset as tags
            self._storeTag( self._curTagName + ".Size", bufSize )
            self._storeTag( self._curTagName + ".Offset", self._f.tell() )
            # - skip data w/o reading
            self._f.seek( self._f.tell() + bufSize )

        return 1

    def _readStructTypes(self):
        # analyses data types in a struct

        if ( debugLevel > 3 ):
            print("Reading Struct Types at Pos = " + hex(self._f.tell()))
        
        structNameLength = self._readIntValue()
        nFields = self._readIntValue()

        if ( debugLevel > 5 ):
            print("nFields = ", nFields)

        if ( nFields > 100 ):
            raise Exception(hex(self._f.tell())+": Too many fields")

        fieldTypes = []
        nameLength = 0
        for i in range( nFields ):
            nameLength = self._readIntValue()
            if ( debugLevel > 9 ):
                print("{}th nameLength = {}".format(i, nameLength))
            fieldType = self._readIntValue()
            fieldTypes.append( fieldType )

        return fieldTypes

    def _readStructData(self, structTypes):
        # reads struct data based on type info in structType
        for i in range( len(structTypes) ):
            encodedType = structTypes[i]
            etSize = self._encodedTypeSize(encodedType)

            if ( debugLevel > 5 ):
                print("Tag Type = " + str(encodedType) + ",", end=' ')
                print("Tag Size = " + str(etSize))

            # get data
            self._readNativeData(encodedType, etSize)

        return 1

    def _storeTag(self, tagName, tagValue):
        # store Tags as list and dict
        # NB: all tag values (and names) stored as unicode objects;
        #     => can then be easily converted to any encoding
        if ( debugLevel == 1 ):
            print(" - storing Tag:")
            print("  -- name:  ", tagName)
            print("  -- value: ", tagValue, type(tagValue))
        # - convert tag value to unicode if not already unicode object
        self._storedTags.append( tagName + " = " + unicode_str(tagValue) )
        self._tagDict[tagName] = unicode_str(tagValue)

    ### END utility functions ###

    def __init__(self, file_, debug=0):
        """DM3 object: parses DM3|DM4 file (either a file path or a file-like object)."""

        ## initialize variables ##
        self._debug = debug
        self._outputcharset = DEFAULTCHARSET
        self._chosenImage = 1
        # - track currently read group
        self._curGroupLevel = -1
        self._curGroupAtLevelX = [ 0 for x in range(MAXDEPTH) ]
        self._curGroupNameAtLevelX = [ '' for x in range(MAXDEPTH) ]
        # - track current tag
        self._curTagAtLevelX = [ '' for x in range(MAXDEPTH) ]
        self._curTagName = ''
        # - create Tags repositories
        self._storedTags = []
        self._tagDict = {}

        ## get binary data
        if isinstance(file_, str):
        # - open file at file path for reading
            self._filename = file_
            self._f = open( self._filename, 'rb' )
        else:
        # - assume we have been passed a file-like object
            if file_.mode == 'rb':
                self._f = file_
                self._filename = file_.name
            else:
                raise Exception("File should be *binary* data.")
            
        ## parse header
        isDM3,isDM4 = (False, False)
        # get version
        fileVersion = readLong(self._f)
        if (fileVersion == 3):
            isDM3 = True
        elif (fileVersion == 4): 
            isDM4 = True
        # get size of root tag directory, check consistency
        fileSize = os.path.getsize(self._filename)
        sizeOK = True
        if isDM3:
            rootLen = readLong(self._f)
            if (rootLen != fileSize - 16):
                sizeOK = False
        elif isDM4:
            rootLen = readLongLong(self._f)
            if (rootLen != fileSize - 24):
                sizeOK = False
        # get byte-ordering
        lE = readLong(self._f)
        littleEndian = (lE == 1)
        if not littleEndian:
            isDM3,isDM4 = (False, False)
            
        # raise Exception if not DM3 or DM4
        if not (isDM3 or isDM4):
            raise Exception("'%s' does not appear to be a DM3/DM4 file."
                            % os.path.split(self._filename)[1])
        elif self._debug > 0:
            print("'%s' appears to be a DM%s file" % (self._filename, fileVersion))

        if ( debugLevel > 5 or self._debug > 1):
            print("Header info. found:")
            print("- file version:", fileVersion)
            print("- byte order:", lE)
            print("- root tag dir. size:", rootLen, "bytes")
            print("- file size:", fileSize, "bytes")
            if not sizeOK:
                msg = "Warning: file size and root tag dir. size inconsistent"
                print("+ %s"%msg)
        
        self._fileVersion = fileVersion
        
        # set name of root group (contains all data)...
        self._curGroupNameAtLevelX[0] = "root"
        # ... then read it
        self._readTagGroup()
        if self._debug > 0:
            print("-- %s Tags read --" % len(self._storedTags))

        # fetch image characteristics
        tag_root = 'root.ImageList.1'
        self._data_type = int( self.tags["%s.ImageData.DataType" % tag_root] )
        self._im_width = int( self.tags["%s.ImageData.Dimensions.0" % tag_root] )
        self._im_height = int( self.tags["%s.ImageData.Dimensions.1" % tag_root] )
        try:
            self._im_depth = int( self.tags['root.ImageList.1.ImageData.Dimensions.2'] )
        except KeyError:
            self._im_depth = 1        

        if self._debug > 0:
            print("Notice: image size: %sx%s px" % (self._im_width, self._im_height))
            if self._im_depth>1:
                print("Notice: %s image stack" % (self._im_depth))

    @property
    def file_version(self):
        """Returns file format version (i.e., 3 or 4)."""
        return self._fileVersion
    
    @property
    def data_type(self):
        """Returns image DataType."""
        return self._data_type

    @property
    def data_type_str(self):
        """Returns image DataType string."""
        return dataTypes[self._data_type]

    @property
    def width(self):
        """Returns image width (px)."""
        return self._im_width

    @property
    def height(self):
        """Returns image height (px)."""
        return self._im_height

    @property
    def depth(self):
        """Returns image depth (i.e. number of images in stack)."""
        return self._im_depth

    @property
    def size(self):
        """Returns image size (width,height[,depth])."""
        if self._im_depth > 1:
            return (self._im_width, self._im_height, self._im_depth)
        else:
            return (self._im_width, self._im_height)

    @property
    def outputcharset(self):
        """Returns Tag dump/output charset."""
        return self._outputcharset

    @outputcharset.setter
    def outputcharset(self, value):
        """Set Tag dump/output charset."""
        self._outputcharset = value

    @property
    def filename(self):
        """Returns full file path."""
        return self._filename

    @property
    def tags(self):
        """Returns all image Tags."""
        return self._tagDict

    def dumpTags(self, dump_dir='/tmp'):
        """Dumps image Tags in a txt file."""
        dump_file = os.path.join(dump_dir,
                                 os.path.split(self._filename)[1]
                                 + ".tagdump.txt")
        try:
            dumpf = open( dump_file, 'w' )
        except:
            print("Warning: cannot generate dump file.")
        else:
            for tag in self._storedTags:
                dumpf.write( "{}\n".format(tag.encode(self._outputcharset)))
            dumpf.close

    @property
    def info(self):
        """Extracts useful experiment info from DM3 file."""
        # define useful information
        tag_root = 'root.ImageList.1.ImageTags'
        info_ = {
            'gms_v': "GMS Version.Created",
            'gms_v_': "GMS Version.Saved",
            'device': "Acquisition.Device.Name",
            'acq_date': "DataBar.Acquisition Date",
            'acq_time': "DataBar.Acquisition Time",
            'binning': "DataBar.Binning",
            'hv': "Microscope Info.Voltage",
            'hv_f': "Microscope Info.Formatted Voltage",
            'mag': "Microscope Info.Indicated Magnification",
            'mag_f': "Microscope Info.Formatted Indicated Mag",
            'mode': "Microscope Info.Operation Mode",
            'micro': "Session Info.Microscope",
            'operator': "Session Info.Operator",
            'specimen': "Session Info.Specimen",         
            'name_old': "Microscope Info.Name",
            'micro_old': "Microscope Info.Microscope",
            'operator_old': "Microscope Info.Operator",
            'specimen_old': "Microscope Info.Specimen",
        #    'image_notes': "root.DocumentObjectList.10.Text' # = Image Notes
            }
        # get experiment information
        infoDict = {}
        for key in info_.keys():
            tag_name = "%s.%s" % (tag_root, info_[key])
            if tag_name in self.tags:
                # tags supplied as Python unicode str; convert to chosen charset
                # (typically latin-1 or utf-8)
                infoDict[key] = self.tags[tag_name].encode(self._outputcharset)
        # return experiment information
        return infoDict

    @property
    def imagedata(self):
        """Extracts image data as numpy.array"""

        # numpy dtype strings associated to the various image dataTypes
        dT_str = {
            1: '<i2',     #16-bit LE signed integer
            2: '<f4',     #32-bit LE floating point
            6: 'u1',      #8-bit unsigned integer
            7: '<i4',     #32-bit LE signed integer
            9: 'i1',      #8-bit signed integer
            10: '<u2',    #16-bit LE unsigned integer
            11: '<u4',    #32-bit LE unsigned integer
            14: 'u1',     #binary
            }

        # get relevant Tags
        tag_root = 'root.ImageList.1'
        data_offset = int( self.tags["%s.ImageData.Data.Offset" % tag_root] )
        data_size = int( self.tags["%s.ImageData.Data.Size" % tag_root] )
        data_type = self._data_type
        im_width = self._im_width
        im_height = self._im_height
        im_depth = self._im_depth

        if self._debug > 0:
            print("Notice: image data in %s starts at %s" % (
                os.path.split(self._filename)[1], hex(data_offset)
                ))

        # check if image DataType is implemented, then read
        if data_type in dT_str:
            np_dt = numpy.dtype( dT_str[data_type] )
            if self._debug > 0:
                print("Notice: image data type: %s ('%s'), read as %s" % (
                    data_type, dataTypes[data_type], np_dt
                    ))
            self._f.seek( data_offset )
            # - fetch image data
            rawdata = self._f.read(data_size)
            # - convert raw to numpy array w/ correct dtype
            ima = numpy.fromstring(rawdata, dtype=np_dt)
            # - reshape to matrix or stack
            if im_depth > 1:
                ima = ima.reshape(im_depth, im_height, im_width)
            else:
                ima = ima.reshape(im_height, im_width)
        else:
            raise Exception(
                "Cannot extract image data from %s: unimplemented DataType (%s:%s)." %
                (os.path.split(self._filename)[1], data_type, dataTypes[data_type])
                )

        # if image dataType is BINARY, binarize image
        # (i.e., px_value>0 is True)
        if data_type == 14:
            ima[ima>0] = 1

        return ima


    @property
    def Image(self):
        """Returns image data as PIL Image"""

        # define PIL Image mode for the various (supported) image dataTypes,
        # among:
        # - '1': 1-bit pixels, black and white, stored as 8-bit pixels
        # - 'L': 8-bit pixels, gray levels
        # - 'I': 32-bit integer pixels
        # - 'F': 32-bit floating point pixels
        dT_modes = {
            1: 'I',     # 16-bit LE signed integer
            2: 'F',     # 32-bit LE floating point
            6: 'L',     # 8-bit unsigned integer
            7: 'I',     # 32-bit LE signed integer
            9: 'I',     # 8-bit signed integer
            10: 'I',    # 16-bit LE unsigned integer
            11: 'I',    # 32-bit LE unsigned integer
            14: 'L',    # "binary"
            }
        
        # define loaded array dtype if has to be fixed to match Image mode
        dT_newdtypes = {
            1:  'int32',      # 16-bit LE integer to 32-bit int
            2:  'float32',    # 32-bit LE float to 32-bit float
            9:  'int32',      # 8-bit signed integer to 32-bit int
            10: 'int32',      # 16-bit LE u. integer to 32-bit int
            }   

        # get relevant Tags
        data_type = self._data_type
        im_width = self._im_width
        im_height = self._im_height    
        im_depth = self._im_depth

        # fetch image data array
        ima = self.imagedata
        # assign Image mode
        mode_ = dT_modes[data_type]

        # reshape array if image stack
        if im_depth > 1:
            ima = ima.reshape(im_height*im_depth, im_width)

        # load image data array into Image object (recast array if necessary)
        if data_type in dT_newdtypes:
            im = Image.fromarray(ima.astype(dT_newdtypes[data_type]),mode_)
        else:
            im = Image.fromarray(ima,mode_)

        return im


    @property
    def contrastlimits(self):
        """Returns display range (cuts)."""
        tag_root = 'root.DocumentObjectList.0'
        low = int(float(self.tags["%s.ImageDisplayInfo.LowLimit" % tag_root]))
        high = int(float(self.tags["%s.ImageDisplayInfo.HighLimit" % tag_root]))
        cuts = (low, high)
        return cuts

    @property
    def cuts(self):
        """Returns display range (cuts)."""
        return self.contrastlimits

    @property
    def pxsize(self):
        """Returns pixel size and unit."""
        tag_root = 'root.ImageList.1'
        pixel_size = float(
            self.tags["%s.ImageData.Calibrations.Dimension.0.Scale" % tag_root])
        unit = self.tags["%s.ImageData.Calibrations.Dimension.0.Units" %
                         tag_root]
        if unit == '\xb5m':
            unit = 'micron'
        else:
            unit = unit.encode('ascii')
        if self._debug > 0:
            print("pixel size = %s %s" % (pixel_size, unit))
        return (pixel_size, unit)


    @property
    def tnImage(self):
        """Returns thumbnail as PIL Image."""
        # get thumbnail
        tag_root = 'root.ImageList.0'
        tn_size = int( self.tags["%s.ImageData.Data.Size" % tag_root] )
        tn_offset = int( self.tags["%s.ImageData.Data.Offset" % tag_root] )
        tn_width = int( self.tags["%s.ImageData.Dimensions.0" % tag_root] )
        tn_height = int( self.tags["%s.ImageData.Dimensions.1" % tag_root] )

        if self._debug > 0:
            print("Notice: tn data in %s starts at %s" % (
                os.path.split(self._filename)[1], hex(tn_offset)
                ))
            print("Notice: tn size: %sx%s px" % (tn_width, tn_height))

        if (tn_width*tn_height*4) != tn_size:
            raise Exception("Cannot extract thumbnail from %s"
                            % os.path.split(self._filename)[1])
        else:
            self._f.seek( tn_offset )
            rawdata = self._f.read(tn_size)
            # - read as 32-bit LE unsigned integer
            tn = Image.frombytes( 'F', (tn_width, tn_height), rawdata,
                                   'raw', 'F;32' )
            # - rescale and convert px data
            tn = tn.point(lambda x: x * (1./65536) + 0)
            tn = tn.convert('L')
        # - return image
        return tn

    @property
    def thumbnaildata(self):
        """Fetch thumbnail image data as numpy.array"""
 
        # get useful thumbnail Tags
        tag_root = 'root.ImageList.0'
        tn_size = int( self.tags["%s.ImageData.Data.Size" % tag_root] )
        tn_offset = int( self.tags["%s.ImageData.Data.Offset" % tag_root] )
        tn_width = int( self.tags["%s.ImageData.Dimensions.0" % tag_root] )
        tn_height = int( self.tags["%s.ImageData.Dimensions.1" % tag_root] )

        if self._debug > 0:
            print("Notice: tn data in %s starts at %s" % (
                os.path.split(self._filename)[1], hex(tn_offset)
                ))
            print("Notice: tn size: %sx%s px" % (tn_width, tn_height))

        # get thumbnail data
        if (tn_width*tn_height*4) == tn_size:
            self._f.seek(tn_offset)
            rawtndata = self._f.read(tn_size)
            print('## rawdata:', len(rawtndata))
           # - read as 32-bit LE unsigned integer
            np_dt_tn = numpy.dtype('<u4')
            tndata = numpy.fromstring(rawtndata, dtype=np_dt_tn)
            print('## tndata:', len(tndata))
            tndata = tndata.reshape(tn_height, tn_width)
            # - rescale and convert to integer
            tndata = tndata/65536. + 0.
            tndata = tndata.astype(int)
            # - return thumbnail data
            return tndata
        else:
            raise Exception("Cannot extract thumbnail from %s"
                            % os.path.split(self._filename)[1])

    def makePNGThumbnail(self, tn_file=''):
        """Save thumbnail as PNG file."""
        # - cleanup name
        if tn_file == '':
            tn_path = os.path.join('./',
                                   os.path.split(self.filename)[1]+'.tn.png')
        else:
            if os.path.splitext(tn_file)[1] != '.png':
                tn_path = os.path.splitext(tn_file)[0] + '.png'
            else:
                tn_path = tn_file
        # - save tn file
        try:
            self.thumbnail.save(tn_path, 'PNG')
            if self._debug > 0:
                print("Thumbnail saved as '%s'." % tn_path)
        except:
            print("Warning: could not save thumbnail.")


## MAIN ##
if __name__ == '__main__':
    print("dm3_lib %s" % VERSION)

