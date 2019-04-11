from __future__ import print_function
from __future__ import division

import numpy as np

def doc():
	'''
<libtiff_ctypes.py>

class TIFF(ctypes.c_void_p):
    """ Holds a pointer to TIFF object.

    To open a tiff file for reading, use

      tiff = TIFF.open (filename, more='r')
      
    To read an image from a tiff file, use

      image = tiff.read_image()

    where image will be a numpy array.

    To read all images from a tiff file, use

      for image in tiff.iter_images():
          # do stuff with image

    To creat a tiff file containing numpy array as image, use

      tiff = TIFF.open(filename, mode='w')
      tiff.write_image(array)
      tiff.close()

    To copy and change tags from a tiff file:

      tiff_in =  TIFF.open(filename_in)
      tiff_in.copy (filename_out, compression=, bitspersample=,
      sampleformat=,...)
    """

	'''
	pass
	
def read_tiff(filename):
  from PIL import Image
  im = Image.open(filename)
  import numpy
  data = numpy.array(im)
  
  return data

  
def bindata(dataori,binning): #dataori needs to be 2d
  binned1=np.zeros(dataori.shape,dtype=np.float32)
  bin_rtn=binned1
  y,x=dataori.shape
  
  if(binning > 1):
      
      for i in range (0,binning):
          binned1 += np.roll(dataori,-i,axis=1)
       #   print(binned1.mean(),i,binned1.shape)      

      binned2=np.copy(binned1[: ,  0:x//binning*binning:binning])
      binned=np.copy(binned2)
      #print(binned2.mean(),i,binned2.shape)      
      for i in range (1,binning):
          binned += np.roll(binned2,-i,axis=0)
          #print(binned.mean(),i,binned.shape)      
     
  
      bin_rtn=np.copy(binned[0:y//binning*binning:binning,:])
      
      bin_rtn/=(binning**2)
  #print (bin_rtn.shape)      
  return bin_rtn


def save_tiff8(buff,desc_data,tif_name,amax,amin,ysize,xsize,udate=False,sigma=None):
  '''
  save numpy ndarray as 8bit image
  will rescale from min to max
  '''
  try:
    from libtiff import TIFF
  except ImportError as e:
    print ("Need <libtiff> for saving TIFF files:\n  pip install libtiff  ")
    return

  data = np.asarray(buff,dtype=np.float32) 
  scale=1.0
  d0=np.zeros(ysize*xsize)
  if (udate): #update data
    amax=data.max()
    amin=data.min()
  if (sigma is not None): #scaling using a sigma cutoff
    amean=data.mean()
    arms=data.std()
    
    amax_t=amean+sigma*arms
    amin_t=amean-sigma*arms
    if(amax_t<amax):
      amax=amax_t
    if(amin_t>amin):
      amin=amin_t
      
  if (amax-amin)>0 :
    scale=255/(amax-amin)
    d0=(data-amin)*scale
  if (amax-amin)<0 :
    scale=255/(amax-amin)
    d0=(data-amax)*scale

  d1=d0.reshape(ysize,xsize)
  d2=d1.astype(np.uint8,casting='unsafe')
  #print("{} {} {}".format(d2[0][0], d2.shape,d2.dtype))
  
  #tiff = libtiff.TIFFimage(d2,description=desc_data)
  #tiff.write_file(tif_name,compression='lzw')
  #del tiff
  tiff = TIFF.open(tif_name, mode='w')
  tiff.SetField('ImageDescription', desc_data)
  tiff.write_image(d2,compression='lzw')
  tiff.close()

def save_tiff16(buff,desc_data,tif_name,amax,amin,ysize,xsize):
  try:
    from libtiff import TIFF
  except ImportError as e:
    print ("Need <libtiff> for saving TIFF files:\n  pip install libtiff  ")
    return

  data = np.asarray(buff,dtype=np.float32) 
  d0=data*65535/(amax-amin)
  d1=d0.reshape(ysize,xsize)
  d2=d1.astype(np.uint16)
  #print("{} {} {}".format(d2[0][0], d2.shape,d2.dtype))
  tiff = TIFF.open(tif_name, mode='w')
  tiff.SetField('ImageDescription', desc_data)
  tiff.write_image(d2,compression='lzw')
  tiff.close()

def save_tiff16_no_rescaling(buff,desc_data,tif_name,amax,amin,ysize,xsize): 
  '''save buffer as-is (dtype), will reshape to xsize x ysize'''
  try:
    from libtiff import TIFF
  except ImportError as e:
    print ("Need <libtiff> for saving TIFF files:\n  pip install libtiff  ")
    return

  #data = np.asarray(buff,dtype=np.int16) 
  d2=buff.reshape(ysize,xsize)
  tiff = TIFF.open(tif_name, mode='w')
  tiff.SetField('ImageDescription', desc_data)
  tiff.write_image(d2,compression='lzw')
  tiff.close()

def save_tiff16_positive(buff,desc_data,tif_name,amax,amin,ysize,xsize): #shift values (x-amin) to change all data points to positive
  try:
    from libtiff import TIFF
  except ImportError as e:
    print ("Need <libtiff> for saving TIFF files:\n  pip install libtiff  ")
    return

  data = buff-amin
  
  data1 = np.asarray(data,dtype=np.uint16) 
  #print('{} {} {} {}'.format(data1.min(),data1.mean(),data1.max(),data1.std()))
  d2=data1.reshape(ysize,xsize)
  tiff = TIFF.open(tif_name, mode='w')
  tiff.SetField('ImageDescription', desc_data)
  tiff.write_image(d2,compression='lzw')
  tiff.close()

#####################
def save_PNG(outfile,buf,header_inf,bit,mode,rmscut=0):
    amean,arms,amin,amax,x,y,z=header_inf['amean'],header_inf['arms'],header_inf['amin'],header_inf['amax'],header_inf['c'],header_inf['r'],header_inf['s']
    
    if rmscut >0:
        cutoff=rmscut*arms #normally 6rms will cover >99.99% pixels
        print("Generating PNG using RMS cutoff = {} ".format(rmscut))
        zeropoint=amin if (amean-cutoff)<amin else (amean-cutoff) #use amin if 6 rms is smaller than amin
        maxpoint=amax if (amean+cutoff)>amax else (amean+cutoff) #same, to maximize constrast
    if rmscut <=0:
        print("Generating PNG with {} = 0 {} = 255".format(amin,amax))
        zeropoint=amin
        maxpoint=amax
    if zeropoint!=maxpoint:
        scale=255/(maxpoint-zeropoint)
    else: 
        scale=0
    


    buf1=''
    for p in buf:
        if p<=zeropoint:
           buf1+=(chr(0)) 
        elif p>=maxpoint:
           buf1+=(chr(255)) 
        else:
            buf1+=(chr (int (scale*(p-zeropoint))))
    
    imdata = png_data(buf1, x, y,bit,mode)
    with open(outfile, 'wb') as fd:        
        fd.write(imdata)        
        fd.close
#################
def png_data(buf, width, height,bit,mode):
    """ buf: must be bytes or a bytearray in Python3.x,
        a regular string in Python2.x.
    """
    import zlib, struct
    bpp=0 #byte per pixel
    if  mode == 0:
        bpp = 1 * bit / 8  #GrayScale
    if  mode == 6:
        bpp = 4 * bit / 8                                                                        #RGBA
    if  mode == 2:
        bpp = 3 * bit / 8                                                                        #RGB

    # reverse the vertical line order and add null bytes at the start
    width_byte = width *bpp 
    
    s=(height) * width_byte
    raw_data = b''.join(
        chr(0) + buf[span : span + width_byte]  for span in range(0,s-1, width_byte)
    )

    def png_pack(png_tag, data):
        chunk_head = png_tag + data
        return (struct.pack("!I", len(data)) +
                chunk_head +
                struct.pack("!I", 0xFFFFFFFF & zlib.crc32(chunk_head)))
    #tEXt_data=b''.join([chr(10),chr(10),'##################not-so-secret message########################',chr(10),chr(10),chr(0),msg,chr(10),chr(10),"##end of message##",chr(10)])
    
    
    return b''.join([
        chr(137),b'PNG',chr(13),chr(10),chr(26),chr(10),
        png_pack(b'IHDR', struct.pack("!2I5B", width, height, bit, mode, 0, 0, 0)), # bit 0 palette  bit 1 color bit 2 alpha: 00000111
        #png_pack(b'tEXt', tEXt_data),
        png_pack(b'IDAT', zlib.compress(raw_data, 9)),
        png_pack(b'IEND', b'')])

