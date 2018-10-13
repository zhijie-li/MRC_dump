#!/usr/bin/env python
from __future__ import print_function

import struct
import os
import sys
import yaml
import csv
import numpy as np
from EM_tiff import save_tiff8, save_tiff16_no_rescaling
import mrc
#import scipy.misc
#import skimage
#from skimage import io
#zlib is used by the write PNG function



def parse_ccp4_map_header(header_data,filesize):
  '''
  Interpreting MRC/CCP4 map file (.map .mrc) header, 1024 Bytes. This will be followed by Symm records, 80 bytes each, and then the data block to the EOF.
  Based on:
  http://www.ccpem.ac.uk/mrc_format/mrc2014.php
  '''
  #still need to write some check: spacegroup vs nxyz, mapxrs... endianness (machst)
  #but EM maps are quite simple

  #print               (header_data[:100])
  
  _c, _r, _s, \
  _mapmode,    \
  _ncstart, _nrstart, _nsstart,    \
  _nx, _ny, _nz,   \
  _x, _y, _z,_alpha,_beta,_gamma,    \
  _mapc,_mapr,_maps,    \
  _amin,_amax,_amean,    \
  _ispg, _nsymbt,_lskflg \
  = struct.unpack("<iiiiiiiiiiffffffiiifffiii", header_data[:100])
  #= struct.unpack("llllllllllfffffflllffflll", header_data[:100]) #need to use i, 'l' will cause problem on some 64bit system(linux 64)
  #print (_c, _r, _s,  _mapmode,  _ncstart, _nrstart, _nsstart,  _nx, _ny, _nz,  _x, _y, _z,_alpha,_beta,_gamma,  _mapc,_mapr,_maps,  _amin,_amax,_amean,  _ispg, _nsymbt,_lskflg)


  _skwmat_string =                     header_data[100:136]
  _skwtrn_string =                     header_data[136:148]
  _futureuse_str =                     header_data[96:96+100]
  _exttyp  			 =										 header_data[104:104+4]
  _nversion			 =										 header_data[108:108+4]
  _orix,_oriy,_oriz,= struct.unpack("<fff", header_data[196:196+12]) # this is defined only in MRC 2014, not in ccp4 maps
  _map_str       =                     header_data[208:208+4]   #'MAP '
  _machst_string =                     header_data[212:212+4]
  _arms,         =struct.unpack( "f" ,header_data[216:216+4])
  _nlabl,        =struct.unpack( "i" ,header_data[220:220+4])  #need to use i, 'l' will cause problem on some 64bit system(linux 64)
  _label         =                     header_data[224:224+800]



  OK = 1
  #some checks
  #1 "MAP "
  if   _map_str[:3].decode() != "MAP": #the .decode() part: because _map_str is byte literal... b'MAP' 
    #it really should be 'MAP ' but some programs such as motioncor2 uses \x00 instead of \x20 as the 4th character..
    print("Warning: the \"MAP(4D 41 50)\" keyword is not found at Bytes 208-210: [" + _map_str + "] \n")
    
    print (_map_str)
    OK =0
#    sys.exit()
  #2 size
  #size of map is determined by column * row * sections * bytes/point.
  #bytes/point:
  #13-16  MODE  0 8-bit signed integer (range -128 to 127)
  #1 16-bit signed integer
  #2 32-bit signed real
  #3 transform : complex 16-bit integers
  #4 transform : complex 32-bit reals
  #6 16-bit unsigned integer  2
  _mode_table={  0 : 1,
                1 : 2,
                2 : 4,  #normally should be this
                3  : 4,
                4 : 8,
                5 : 0,
                6 : 2
  }
  _mode_table_human={  0 : "envelope, signed 8-bit Int -128..127 (1 Byte)",
                1 : "16-bit Int (2 Bytes)",
                2 : "32-bit Real (4 Bytes)",  #normally should be this
                3  : "Complex 16-bit Int x 2 (4 Bytes)",
                4 : "Complex 32-bit Real x 2 (8 Bytes)",
                5 : "mode 5: unkown",
                6 : "16-bit unsigned Int (2 Bytes)"
  }
  
  _bytes_per_point = _mode_table[_mapmode]
  _expected_data_size = _c * _r * _s * _bytes_per_point
  _expected_file_size = _expected_data_size + 1024 + _nsymbt #_nsymbt size of extended header (which follows main header) in bytes  7

  s= "\nfile size: " + str(filesize)  + " data block size: " + str(_expected_data_size) + "\nheader + SYMM: "+ str(1024+_nsymbt) + "\n" \
    + "headersize: 1024  SYMM blocks (80 Bytes each): " + str(_nsymbt) +" Bytes"
  print(s)

  if (_nsymbt % 80) != 0:
    print ( "Warning! The total size of symbol blocks is not multiples of 80 Bytes \n" )
  if _expected_file_size != filesize:
    print ( "Warning! The calculated file size is different from the actual size \n" )
  if _expected_file_size == filesize and (_nsymbt % 80) == 0:
    print ("File size check OK")
  else:
    OK = 0
  #3 Spacegroup
  if _ispg <= 0 or _ispg >=400:
    print ( "This file does not contain a crystallographic space group\n")
  #4 endianess
  #Note 11 Bytes 213 and 214 contain 4 `nibbles' (half-bytes) indicating the
  #representation of float, complex, integer and character datatypes. Bytes 215 and
  #216 are unused. The CCP4 library contains a general representation of datatypes,
  #but in practice it is safe to use 0x44 0x44 0x00 0x00 for little endian
  #machines, and 0x11 0x11 0x00 0x00 for big endian machines. The CCP4 library uses
  #this information to automatically byte-swap data if appropriate, when
  #tranferring data files between machines.
  endianness='LE' #little endian by default
  if ord(_machst_string[0]) & ord('\x40') == 1: #will only check for float
    endianness='LE'
  if ord(_machst_string[0]) & ord('\x10') == 1:
    endianness='BE'
    
      


  header_inf={
    'c'                 : _c             ,
    'r'                 : _r             ,
    's'                 : _s             ,
    'mapmode'           : _mapmode       ,
    'ncstart'           : _ncstart       ,
    'nrstart'           : _nrstart       ,
    'nsstart'           : _nsstart       ,
    'nx'                : _nx            ,
    'ny'                : _ny            ,
    'nz'                : _nz            ,
    'x'                 : _x             ,
    'y'                 : _y             ,
    'z'                 : _z             ,
    'alpha'             : _alpha         ,
    'beta'              : _beta          ,
    'gamma'             : _gamma         ,
    'mapc'              : _mapc          ,
    'mapr'              : _mapr          ,
    'maps'              : _maps          ,
    'amin'              : _amin          ,
    'amax'              : _amax          ,
    'amean'             : _amean         ,
    'ispg'              : _ispg          ,
    'nsymbt'            : _nsymbt        ,
    'lskflg'            : _lskflg        ,
    'skwmat_string'     : _skwmat_string ,
    'skwtrn_string'     : _skwtrn_string ,
    'futureuse_str'     : _futureuse_str ,
    'exttyp'						: _exttyp        ,
    'nversion'					: _nversion			 ,
    'orix'							: _orix          ,
    'oriy'							: _oriy          ,
    'oriz'							: _oriz          ,
    'map_str'           : _map_str       ,
    'machst_string'     : _machst_string ,
    'arms'              : _arms          ,
    'nlabl'             : _nlabl         ,
    'label'             : _label         ,
    'OK'                : OK             ,
    'endianness'        : endianness
  }


  
  outlist= \
  'c, r, s:                 [{}\t{}\t{}]\n'.format(_c,_r,_s) +\
  'nstart c, r, s:          [{}\t{}\t{}]\n'.format(_ncstart,_nrstart,_nsstart)  +\
  'map c, r, s:             [{}\t{}\t{}]\n'.format(_mapc,_mapr,_maps) + '\n'+\
  'mx, my, mz:              [{}\t{}\t{}]\n'.format(_nx,_ny,_nz) +\
  'origins x, y, z:         [{}\t{}\t{}]\n'.format(_orix,_oriy,_oriz) +'\n' +\
  'a, b, c:                 [{}\t{}\t{}]\n'.format(_x,_y,_z) +\
  'alpha, beta, gamma:      [{}\t{}\t{}]\n'.format(_alpha,_beta,_gamma)+ \
  'amin, amax:              [{}\t{}]   \n'.format(_amin,_amax,_amean) + \
  'amean, arms:             [{}\t{}]   \n'.format(_amean,_arms) + \
  'ispg, nsymbt,lskflg:     [{}\t{}\t{}]\n'.format(_ispg,_nsymbt,_lskflg) + \
  'skwmat_string:           [{}]\n'.format(_skwmat_string) + \
  'skwtrn_string    :       [{}]\n'.format(_skwtrn_string) + \
  'futureuse_str:           [{}]\n'.format(_futureuse_str) + \
  'map_str:                 [' + _map_str + "]\n" +\
  'machst_string:           [' + str(bytearray(_machst_string)).encode('hex') + "]\n" +\
  'nlabl                    ['+ str(_nlabl)        + "]\n"    +\
  'label                    ['+ _label  + "]\n" +\
  'Mapmode:                 [' +str(_mapmode) +"] "+ _mode_table_human[_mapmode] + "\n" +\
  'Endianness:              [' + endianness + "]\n"

  
  print (outlist)
  

  return      header_inf
#endo of parse_ccp4_map_header


        
#################################3
def unpack_data(header_inf,f):

    le_be={
    'LE' : '<',
    'BE' : '>'
    }    

    data_type_table={
                0 : 'b'  ,   #  "envelope, signed 8-bit Int -128..127 (1 Byte)",
                1 : 'h'  ,   #  "16-bit Int (2 Bytes)",
                2 : 'f'  ,   #  "32-bit Real (4 Bytes)",  #normally should be this
                3 : 'hh' ,    #  "Complex 16-bit Int x 2 (4 Bytes)",
                4 : 'ff' ,    #  "Complex 32-bit Real x 2 (8 Bytes)",
                5 : ''   ,   #  "mode 5: unkown",
                6 : 'H'     #  "16-bit unsigned Int (2 Bytes)"
    }
    
    f_str=data_type_table[header_inf['mapmode']]

    _mode_table={  0 : 1,
                1 : 2,
                2 : 4,  #normally should be this
                3  : 4,
                4 : 8,
                5 : 0,
                6 : 2
    }

    data_len=_mode_table[header_inf['mapmode']]
    
    format_str= le_be[header_inf['endianness']]+ f_str * header_inf['c'] #format string for a column
    
    data= [[[0.0 for x in xrange(header_inf['c'])] for y in xrange(header_inf['r']) ]for z in xrange(header_inf['s']) ] 
    data1d = [0.0 for x in xrange(header_inf['c']*header_inf['r']*header_inf['s'])]
    count=0

    for s in xrange(header_inf['s']):
        for r in xrange(header_inf['r']):
                data[s][r]=struct.unpack(format_str,f.read(data_len * header_inf['c']))
                start=s*header_inf['r']*header_inf['c']+r*header_inf['c']
                data1d[  start: start+header_inf['c'] ] = data[s][r]
                #print (start )
                count+= header_inf['c'] * data_len

    print (str(count) + " Bytes of map data read." )
    #print( data1d)
    return data, data1d
#########################################3
def write_map_float4LE(data): #takes 1D float array

    fmt='<f'
    data_str=''
    for f in data:
        data_str+=struct.pack(fmt,f)
    return data_str

####end of reading data block############    

def mapping_xyz(header_inf):
    lut={1:'x', 2:'y',3:'z'} #lookup table for maping mapc,mapr,maps to c r s
    mapping={
    'c':lut[header_inf['mapc']],
    'r':lut[header_inf['mapr']],
    's':lut[header_inf['maps']]
    }
    inv_mapping={v :k for k, v in mapping.items()}

    print (mapping)
    #print (inv_mapping)
    return inv_mapping

#####################
def save_yaml(data):
    with open('dump.yml', 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)
#####################
def save_csv(data,header_inf):
    with open('dump.csv', 'wb') as outfile:
        spamwriter = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for s in xrange(header_inf['s']):
            for r in xrange(header_inf['r']):
                spamwriter.writerow(data[s][r])
def bin(array,x,y): #bin 1d array

    binned=0

    return binned

############################main##################
def main():
    filename="testmap.map"
    #pngname=filename+".slice1.png"
    if len(sys.argv)>1 and sys.argv[1]:
      filename=sys.argv[1]
      pngname=filename+".slice1.png"
      if len(sys.argv)>2 and sys.argv[2]:
          pngname=sys.argv[2]
    elif len(sys.argv)>3: 
      print ("This program only accepts 1-2 parameters") #the 1st is the input, 2nd is the output png filename if given
      sys.exit(0)
    
    
    file_size=os.path.getsize(filename)
    header_inf,exthdr,data_1d=mrc.read_MRC2014(filename)
    x,y,z=header_inf['c'],header_inf['r'],header_inf['s']
    data_3d=data_1d.reshape(z,y,x)
    if(len(exthdr)>0):
      print("Extended header found: {} Bytes, using {} Bytes".format(len(exthdr),header_inf['nsymbt']))
    for s in xrange(0,header_inf['MRC_NZ']):
      tif8_name=filename+'{:03d}'.format(s)+'.uint8.tif'
      save_tiff8(data_3d[s],exthdr,tif8_name,header_inf['amax'],header_inf['amin'],y,x)
      tif16_name=filename+'{:03d}'.format(s)+'.tif'
      save_tiff16_no_rescaling(data_3d[s],exthdr,tif16_name,header_inf['amax'],header_inf['amin'],y,x)
    
 #   f=open(filename,'rb')
 #   try:
 #       f.seek(0)
 #       header = f.read(1024)
 #       if header:
 #           header_inf=parse_ccp4_map_header(header,file_size)
 #           if header_inf['OK'] == 1: #header check OK
 #               #save_yaml(header_inf) #optional
 #               if header_inf['nsymbt']>0:
 #                   print ("SYMM blocks: (" + str(header_inf['nsymbt']/80) + ")")
 #                   sym_blc = f.read(header_inf['nsymbt'])
 #                   for i in xrange(0,header_inf['nsymbt']/80):
 #                     print ("[" + sym_blc[i*80:i*80+80]+ "]" )
 #               #data=f.read(file_size - 1024 - header_inf['nsymbt'])
 #               #data_1d=unpack_data(header_inf,f) 
 #               x,y,z=header_inf['c'],header_inf['r'],header_inf['s']
 #               data_1d=np.fromfile(f,dtype='float32',count=x*y*z)
 #               mapping_xyz(header_inf)
 #               
 #               #save_yaml(data) #optional
 #               #save_csv(data_3dmatrix,header_inf)#optional
##        sum=0;
##        for d in data_1d:
##            sum+=d;
##            if header_inf['amin']>d:
##                header_inf['amin']=d
##            if header_inf['amax']<d:
##                header_inf['amax']=d
 #               if header_inf['arms']==0:
 #                   header_inf['amean']=new_amean=np.mean(data_1d)
 #                   header_inf['amin']=new_amin=np.min(data_1d)
 #                   header_inf['amax']=new_amax=np.max(data_1d)
 #                   
 #                   header_inf['arms']=new_arms= np.sqrt(np.mean(np.square(data_1d-new_amean)))
 #                   
 #                   print ("updated min, max, mean, rms: {} {} {} {}".format(new_amin,new_amax,new_amean, new_arms))
 #               x,y,z=header_inf['c'],header_inf['r'],header_inf['s']
 #               #data_slice=np.reshape(data,(z,y,x))
 #               
 #               for s in xrange(0,header_inf['s']):
 #                   
 #                   start=s*header_inf['c']*header_inf['r']
 #                   end=start+ header_inf['c']*header_inf['r']
 #                   print ("Saving slice <{:03d}>".format(s),     "datablock:          [{}\t{}]\n".format(start,end))
 #
                    
    #save_PNG(filename+'{:03d}'.format(s)+'.png',data_1d[start:end],header_inf,8,0,rmscut=6,slice_num=s)
    #finally:
    #    f.close()

    #alldata=np.array(data_3dmatrix)
    #the following is for a 2-fold average of map
    #from operator import add
    #def cond_add(x,y):
    #    if x==0.0 or y == 0.0:
    #        return 0.0
    #    else:
    #        return x+y
    #            
    #data2=map(cond_add,data_1d, reversed(data_1d))
    #
    #outdata= header + write_map_float4LE(data2)
    #fo=open("ave.mrcs",'wb')
    #fo.write(outdata)
    #fo.close
    
    #some stats
    #import matplotlib.pyplot as plt
    #hist, bin_edges = np.histogram(alldata,  bins=np.arange(-0.5,1.5,0.1))
    #print (hist)
    #plt.hist(alldata, bins='auto')  # arguments are passed to np.histogram
    #plt.title("Histogram with 'auto' bins")
    #plt.show()
    
    
#####################3
main()