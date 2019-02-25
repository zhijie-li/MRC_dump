#!/usr/bin/env python


import numpy as np
import struct
import os
import argparse

NUMPY_MODE = {
        0: np.dtype(np.int8),
        1: np.dtype(np.int16),
        2: np.dtype(np.float32),
        6: np.dtype(np.uint16),
        np.dtype(np.int8):    0,
        np.dtype(np.int16):   1,
        np.dtype(np.float32): 2,
        np.dtype(np.uint16):  6
       }

DSIZE_TABLE={
                0 : 1,
                1 : 2,
                2 : 4,  #normally should be this
                3  : 4,
                4 : 8,
                6 : 2
  }
MODE_TABLE_HUMAN={  0 : "8-bit signed integer (range -128 to 127)",
                1 : "16-bit Int (2 Bytes)",
                2 : "32-bit Real (4 Bytes)",  #normally should be this
                3  : "Complex 16-bit (4 Bytes)",
                4 : "Complex 32-bit  (8 Bytes)",
                5 : "mode 5: unkown",
                6 : "16-bit unsigned Int (2 Bytes)"
  }
###RELION recently introduced mode 101:4-bit int for movies



def parse_MRC2014_header(header_data,filesize):
  '''
  Interpreting MRC2014 file (.mrc) header. This one is specific for MRC2014.
  http://www.ccpem.ac.uk/mrc_format/mrc2014.php
  '''
  nx,ny,nz, \
  mapmode,    \
  nxstart, nystart, nzstart,    \
  mx, my, mz,   \
  cella, cellb, cellc, cellalpha,cellbeta,cellgamma,    \
  mapc,mapr,maps,    \
  dmin,dmax,dmean,    \
  ispg, nsymbt \
  = struct.unpack("<iiiiiiiiiiffffffiiifffii", header_data[:96])

  _lskflg         = struct.unpack("<i",   header_data[96:100] )#ccp4 Xtal map only. lskflg:Flag for skew transformation, =0 none, =1 if foll
  _skwmat_string  =                       header_data[100:136] #ccp4 Xtal map only
  _skwtrn_string  =                       header_data[136:148] #ccp4 Xtal map only
  _futureuse_str  =\
  extra           =                       header_data[96:96+100]
  exttyp          =                       header_data[104:104+4]
  nversion        =                       header_data[108:108+4]
  orix,oriy,oriz, = struct.unpack("<fff", header_data[196:196+12]) # this is defined only in MRC 2014, not in ccp4 maps
  map_str         =                       header_data[208:208+4]   #'MAP '
  machst_string   =                       header_data[212:212+4]
  machst,         = struct.unpack("<i",   header_data[212:212+4] ) #"in practice it is safe to use 0x44 0x44 0x00 0x00 for little endian machines",  68*256+68=17476, or 'DD' as string. Sometimes (serialEM) it can be 0x44 0x41:'DA',16708 
  rms,            = struct.unpack( "f" ,  header_data[216:216+4])
  nlabl,          = struct.unpack( "i" ,  header_data[220:220+4])  #need to use i, 'l' will cause problem on some 64bit system(linux 64)
  label           =                       header_data[224:224+800]



  OK = 1
  #some checks
  #1 "MAP "
  if   map_str[:3].decode() != "MAP": 
    #it really should be 'MAP ' but some programs such as motioncor2 uses \x00 instead of \x20 as the 4th character..
    print("Warning: the \"MAP(4D 41 50)\" keyword is not found at Bytes 208-210: [" + _map_str + "] \n")
    print (map_str)
    OK =0
  if (exttyp.decode()=='CCP4'):
    print("Warning: the EXTTYPE indicats that the extended header is for a CCP4 map ")
    print (exttype)
    OK =0

  #2 size
  #size of map is determined by column * row * sections * bytes/point.
  #bytes/point:
  #13-16  MODE  0 8-bit signed integer (range -128 to 127)
  #1 16-bit signed integer
  #2 32-bit signed real
  #3 transform : complex 16-bit integers
  #4 transform : complex 32-bit reals
  #6 16-bit unsigned integer  2

  _bytes_per_point = DSIZE_TABLE[mapmode]
  _expected_data_size = nx * ny * nz * _bytes_per_point
  _expected_file_size = _expected_data_size + 1024 + nsymbt #_nsymbt size of extended header (which follows main header) in bytes  7

  #s= "\nfile size: " + str(filesize)  + " data block size: " + str(_expected_data_size) + "\nheader + SYMM: "+ str(1024+_nsymbt) + "\n"  + "headersize: 1024  SYMM blocks (80 Bytes each): " + str(_nsymbt) +" Bytes"
  #print(s)

  if (nsymbt % 80) != 0:
    print ( "Warning! The extended header (NSYMBT) <{} bytes> is not multiple of 80 Bytes.\n".format(nsymbt) )
    #OK = 0 #not everyone follows 80-byte rule. e.g. EPU.
  if _expected_file_size != filesize:
    print ( "Error! The calculated file size is different from the actual size \n" )
    OK = 0
    #exit()
  if _expected_file_size == filesize and (nsymbt % 80) == 0:
    print ("File size check OK")


  #3 Spacegroup
  #if _ispg <= 0 or _ispg >=400:
  #  print ( "This file does not contain a crystallographic space group\n")
  
  #4 endianess
  #Note 11 Bytes 213 and 214 contain 4 `nibbles' (half-bytes) indicating the
  #representation of float, complex, integer and character datatypes. Bytes 215 and
  #216 are unused. The CCP4 library contains a general representation of datatypes,
  #but in practice it is safe to use 0x44 0x44 0x00 0x00 for little endian
  #machines, and 0x11 0x11 0x00 0x00 for big endian machines. The CCP4 library uses
  #this information to automatically byte-swap data if appropriate, when
  #tranferring data files between machines.

  endianness='LE' #little endian by default
  if(bytearray(machst_string)[0] & ord('\x44')):
   #== 68 and (machst_string[1] & ord('\x44') == 68  or machst_string[1] & ord('\x41') == 65)): 
#  if(ord(machst_string[0]) & ord('\x44') == 68 and (ord(machst_string[1]) & ord('\x44') == 68  or ord(machst_string[1]) & ord('\x41') == 65)): 
    endianness='LE'
    #print("MACHST indicates LE")
  if(bytearray(machst_string)[0] & ord('\x11') ==  ord('\x11')): 
#  if(ord(machst_string[0]) & ord('\x11') ==  ord('\x11')): 
    endianness='BE'
  if (machst_string[:2] !='DD' and machst_string[:2] !='DA' ):
    print("Warning: the MACHST string check failed. This file may not be in little endian! (should be 0x44 0x44 0x00 0x00 \"DD\") or 0x44 0x41 0x00 0x00 \"DA\")")
    endianness='BE'
    OK =0




  header_inf={
    'MRC_NX'                : nx             ,
    'MRC_NY'                : ny             ,
    'MRC_NZ'                : nz             ,
    'c'                 : nx             , #for backward compatibility with ccp4, everything has two keys. MRC2014 items always uppercase and start with MRC_
    'r'                 : ny             ,
    's'                 : nz             ,
    'MRC_MAPMODE'           : mapmode       ,
    'mapmode'           : mapmode       ,
    'MRC_NXSTART'           : nxstart       ,
    'MRC_NYSTART'           : nystart       ,
    'MRC_NZSTART'           : nzstart       ,
    'nxstart'           : nxstart       ,
    'nystart'           : nystart       ,
    'nzstart'           : nzstart       ,
    'MRC_MX'                : mx            ,
    'MRC_MY'                : my            ,
    'MRC_MZ'                : mz            ,
    'nx'                : mx            ,
    'ny'                : my            ,
    'nz'                : mz            ,
    'MRC_CELL_A'            : cella             ,
    'MRC_CELL_B'            : cellb             ,
    'MRC_CELL_C'            : cellc             ,
    'MRC_CELL_ALPHA'        : cellalpha         ,
    'MRC_CELL_BETA'         : cellbeta          ,
    'MRC_CELL_GAMMA'        : cellgamma         ,
    'x'                 : cella        ,
    'y'                 : cellb        ,
    'z'                 : cellc        ,
    'alpha'             : cellalpha   ,
    'beta'              : cellbeta    ,
    'gamma'             : cellgamma   ,
    'MRC_MAPC'              : mapc          ,
    'MRC_MAPR'              : mapr          ,
    'MRC_MAPS'              : maps          ,
    'mapc'              : mapc          ,
    'mapr'              : mapr          ,
    'maps'              : maps          ,
    'MRC_DMIN'              : dmin          ,
    'MRC_DMAX'              : dmax          ,
    'MRC_DMEAN'             : dmean         ,
    'amin'              : dmin          ,
    'amax'              : dmax          ,
    'amean'             : dmean         ,
    'MRC_ISPG'              : ispg          ,
    'MRC_NSYMBT'            : nsymbt        ,
    'ispg'              : ispg          ,
    'nsymbt'            : nsymbt        ,

    'MRC_EXTRA'             : extra         ,
    'lskflg'            : _lskflg        ,
    'skwmat_string'     : _skwmat_string ,
    'skwtrn_string'     : _skwtrn_string ,
    'futureuse_str'     : _futureuse_str ,
    'MRC_EXTTYP'            : exttyp        ,
    'exttyp'            : exttyp        ,

    'MRC_NVERSION'          : nversion      ,
    'nversion'          : nversion      ,
    'MRC_ORIX'              : orix          ,
    'MRC_ORIY'              : oriy          ,
    'MRC_ORIZ'              : oriz          ,
    'orix'              : orix          ,
    'oriy'              : oriy          ,
    'oriz'              : oriz          ,
    'MRC_MAP_STR'           : map_str       ,
    'MRC_MACHST_STRING'     : machst_string ,
    'MRC_MACHST_INT'        : machst       ,
    'map_str'           : map_str       ,
    'machst_string'     : machst_string ,
    'machst_int'        : machst       ,
    'MRC_rms'               : rms          ,
    'arms'               : rms          ,
    'MRC_NLABL'             : nlabl         ,
    'MRC_LABEL'             : label         ,
    'nlabl'             : nlabl         ,
    'label'             : label         ,
    'OK'                : OK             ,
    'endianness'        : endianness ,
    'ori_1024'          :header_data[:1024]
  }

  import re
  p = re.compile(b'\x00')
  labelshort=p.sub('', label)

  outlist= \
  'nx, ny, nz:              [{}\t{}\t{}]\n'.format(nx,ny,nz) +\
  'start x, y, z:           [{}\t{}\t{}]\n'.format(nxstart,nystart,nzstart)  +\
  'map c, r, s:             [{}\t{}\t{}]\n'.format(mapc,mapr,maps) +\
  'mx, my, mz:              [{}\t{}\t{}]\n'.format(mx,my,mz) +\
  'origins x, y, z:         [{}\t{}\t{}]\n'.format(orix,oriy,oriz) +\
  'a, b, c:                 [{}\t{}\t{}]\n'.format(cella,cellb,cellc) +\
  'dmin, dmax:              [{}\t{}]   \n'.format(dmin,dmax,dmean) + \
  'dmean, rms:              [{:.2f}\t{}]   \n'.format(dmean,rms) + \
  'ispg, nsymbt:            [{}\t{}]\n'.format(ispg,nsymbt) + \
  'map_str:                 [{}]\n'.format(map_str) +\
  'machst_string:           [{:}]\n'.format(machst_string) +\
  'nlabl                    [{}]\n'.format(nlabl) +\
  'label                    [{}]\n'.format(labelshort) +\
  'Mapmode:                 [{} => {}]\n'.format(str(mapmode),MODE_TABLE_HUMAN[mapmode]) +\
  'Endianness:              [{}]\n'.format(endianness)
  
  print(outlist)
  if(OK==1):
    print ("file check OK")
  else:
    print ("file check finished with warning")

  return      header_inf
  
  
def read_MRC2014(filename):
    file_size=os.path.getsize(filename)
    header_inf={}
    exthdr=''
    data_1d=np.ndarray(0)
    print("\nReading MRC file <{}>\n".format (filename))
    try:
        f=open(filename,'rb')
        f.seek(0)
        header = f.read(1024)
        header_inf=parse_MRC2014_header(header,file_size)                   
        
        if header_inf['nsymbt']>0:
          exthdr = f.read(header_inf['nsymbt'])                       
        x,y,z=header_inf['MRC_NX'],header_inf['MRC_NY'],header_inf['MRC_NZ']            
        data_1d=np.fromfile(f,dtype=NUMPY_MODE[header_inf['MRC_MAPMODE']],count=x*y*z)               
    finally:                 
      f.close()
    return header_inf,exthdr,data_1d

def proc(filename,args):
    file_size=os.path.getsize(filename)
    
    header_inf,exthdr,data_1d=read_MRC2014(filename)

  
def main():
    print('MRC lib and tools       <zhijie.li@utoronto.ca>')
    parser = argparse.ArgumentParser(description='MRC file reading and manipulation')
    parser.add_argument("mrc_files",metavar='.mrc files',nargs='+',type=str, help="The .mrc files to process.")
    args=parser.parse_args()                           
    
    for filename in args.mrc_files:
        proc(filename,args)

if __name__ == '__main__':
  main()
  

