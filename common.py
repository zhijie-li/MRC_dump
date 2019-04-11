from __future__ import print_function
from __future__ import division
import sys,os
import numpy as np
print('###############################################################################################')
print( "Python:",os.path.realpath(sys.executable), sys.version.split('\n')[0])
print( "numpy:", np.version.version)
PY3 = sys.version_info > (3,)
if(PY3):
  print("\n   !!!   Warning: Python3 detected     !!!")
  
  print("   !!! This program is written python2 !!!\n")

print("\n")


def save_safe_yaml(data,filename):
    import yaml

    yamltext=yaml.safe_dump(data , default_flow_style=False)
    with open(filename, 'w') as outfile:
        yaml.safe_dump(data, outfile, default_flow_style=False)
    return yamltext
def save_yaml(data,filename):
    import yaml

    with open(filename, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)
def print_safe_yaml(data):
  import yaml

  print (yaml.safe_dump(data , default_flow_style=False))
def print_yaml(data):
  import yaml

  print (yaml.dump(data , default_flow_style=False))


def print_hex(astr):
  hexs=" ".join("{:02x}".format(ord(c)) for c in astr)
  
  print ("<{}>".format( hexs), end='\t')


def ndary_tobytes(ndary):
  '''
  equivalent to ndarry.tobytes
  this is for older version of numpy
  '''
  fmtstr={
    np.dtype(np.int8):     'b'       ,   #signed 1-byte integer
    np.dtype(np.int16):    'h'       ,   #Signed 2-byte integer
    np.dtype(np.int32):    'i'       ,   #Signed 4-byte integer
    np.dtype(np.float32):  'f'       ,   #4-byte float
    np.dtype(np.uint16):   'd'       ,   #8-byte float
  }


  #fmt=fmtstr[ndary.dtype]
  #onedary=ndary.reshape(ndary.size)
  #a=onedary.tolist()
  #bstring= struct.pack(fmt * len(a),*a)
  #===========method2========
  #based on https://stackoverflow.com/questions/43925624/fastest-method-to-dump-numpy-array-into-string
  tmp = io.BytesIO()
  np.save(tmp, ndary)
  #ndary.tofile(tmp)

  bstring=tmp.getvalue()
  #the npy header length is saved in bytes 9 and 10 in the header, total headerlength = 6(0x93NUMPY)+2(vv)+2(LL)+header_length
  start=10+ord(bstring[9])*256+ord(bstring[8]) #normally 80
  #print( start)
  return bstring[start:]





def find_dir(filename):
  if(os.path.isfile(filename)):
      file_dir = os.path.dirname(os.path.abspath(filename))
      mrcname=os.path.basename(filename)
      (root,ext)=os.path.splitext(mrcname)
  return (root,ext,mrcname,file_dir)
