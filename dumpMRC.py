#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

#import struct
import os
#import sys
#import yaml
#import csv
import numpy as np
#import scipy.stats
import EMAN2star #needs to copy EMAN2star.py from EMAN2 folder to program folder

import argparse

from EM_tiff import save_tiff8, save_tiff16_no_rescaling, bindata
import mrc

#import scipy.misc
#import skimage
#from skimage import io
#zlib is used by the write PNG function


############################main##################



def modify_image(data3d):
    pass

def bin_inf(data_2d,binning,sigma):
  binned=bindata(data_2d,binning)
  if(binning <=1):
    binned =data_2d
  y,x=binned.shape
 # print(y,x)
  newmax=binmax=  binned.max()
  newmin=binmin=  binned.min()
  binrms=  binned.std()
  binmean=  binned.mean()

  if(binmean+sigma*binrms<binmax):
    newmax=binmean+sigma*binrms
  if(binmean-sigma*binrms>binmin):
    newmin=binmean-sigma*binrms
#  print(newmax,newmin,binrms,binmean,y,x)
  return binned,newmax,newmin,binrms,binmean,y,x
    
def bin_inf_3d(data_3d,binning,sigma,header_inf):
  z,y,x=data_3d.shape
  if(binning <=1):
    binned =data_3d
    newmax=header_inf['amax']
    newmin=header_inf['amin']
    binmean=header_inf['amean']
    binrms=header_inf['arms']
    return binned,newmax,newmin,binrms,binmean,z,y,x
    
  else:
    binned=np.zeros((z,int(y/binning),int(x/binning)),dtype=data_3d.dtype)
    for i in range(z):
        
        binned[i]=bindata(data_3d[i],binning)
    bz,by,bx=binned.shape

    newmax=binmax=  binned.max()
    newmin=binmin=  binned.min()
    binrms=  binned.std()
    binmean=  binned.mean()
    
    if(binmean+sigma*binrms<binmax):
      newmax=binmean+sigma*binrms
    if(binmean-sigma*binrms>binmin):
      newmin=binmean-sigma*binrms
#    print(newmax,newmin,binrms,binmean,y,x)
    return binned,newmax,newmin,binrms,binmean,bz,by,bx

def display_image(data_2d,binning=2,sigma=3): #needs 2d data
  from PIL import Image

  (binned,newmax,newmin,binrms,binmean,y,x)=bin_inf(data_2d,binning,sigma)
  
  img = Image.fromarray(binned)
  #img.save('my.png')
  img.show()

def matplot_image(data_2d,binning=2,sigma=3): #needs 2d data

  (binned,newmax,newmin,binrms,binmean,y,x)=bin_inf(data_2d,binning,sigma)

  from matplotlib import pyplot as plt
  plt.imshow(binned, interpolation='nearest')
  plt.show()



def proc(filename,args,ctf,binning=4,sigma=3,dump=True):
    if args.bin is not None:
      binning=args.bin
      dump=True
    if args.sigma is not None:
      sigma=args.sigma
      dump=True
    if args.show:
      dump=False
    if args.info:
      dump=False
      
##READ AND GET SCALING DONE    
    file_size=os.path.getsize(filename)
    header_inf,exthdr,data_1d=mrc.read_MRC2014(filename)
    x,y,z=header_inf['c'],header_inf['r'],header_inf['s']
    data_3d=data_1d.reshape(z,y,x)
    
      
    if header_inf['arms']<=0 or args.update_stats_frame0: #unreliable header inf, update
        header_inf['amin']=data_3d[0].min()
        header_inf['amax']=data_3d[0].max()
        header_inf['amean']=data_3d[0].mean()
        header_inf['arms']=data_3d[0].std()
        print ("Updated based on frame 0: min {} max {} mean {} rms {}".format(header_inf['amin'],header_inf['amax'],header_inf['amean'],header_inf['arms']))
    if args.update_stats: #force update stats with all data, this is for maps
        header_inf['amin']=data_1d.min()
        header_inf['amax']=data_1d.max()
        header_inf['amean']=data_1d.mean()
        header_inf['arms']=data_1d.std()
        print ("Updated based on all data: min {} max {} mean {} rms {}".format(header_inf['amin'],header_inf['amax'],header_inf['amean'],header_inf['arms']))
      
      #print(data_3d.dtype)
    newmax=header_inf['amax']
    newmin=header_inf['amin']
    if(header_inf['amean']+6*header_inf['arms']<header_inf['amax']):
        newmax=header_inf['amean']+6*header_inf['arms']
    if(header_inf['amean']-6*header_inf['arms']>header_inf['amin']):
        newmin=header_inf['amean']-6*header_inf['arms']
        
    if(len(exthdr)>0):
        print("Extended header found: {} Bytes, using {} Bytes".format(len(exthdr),header_inf['nsymbt']))

    if args.csv:
      import csv
      import numpy as np
      with open(filename+'.csv', 'w') as csv_f:
        data_2d=data_1d.reshape(z*y,x)
        np.savetxt(csv_f, data_2d, '%s', ',')



###SHOW frame1

    if args.show:
      display_image(data_3d[0],binning=binning,sigma=sigma)
      #matplot_image(data_3d[0],binning=binning,sigma=sigma)
      #exit()
      
###dump frames
    if dump==True:
      print("Dumping {} {} bytes {} x {} x {}".format(filename,file_size,x,y,z))
      (binned,newmax,newmin,binrms,binmean,bz,by,bx)=bin_inf_3d(data_3d,binning,sigma,header_inf)

      for s in range(0,header_inf['MRC_NZ']):
        tif8_name=filename+'.'
        if(header_inf['nz']>1):
          tif8_name+='{:03d}_'.format(s)
        tif8_name+='m{:1.1f}_'.format(header_inf['amean'])
        if(filename in ctf):
          tif8_name+=ctf[filename]+'_'
        tif8_name+='uint8.tif'
        
        save_tiff8(binned[s],exthdr,tif8_name,newmax,newmin,by,bx)
#        tif16_name=filename+'{:03d}'.format(s)+'.tif'

        #Saving in tif is not really useful for float data. But if the data are in int8 or 16 or 32, compression can be good.
        #save_tiff16_no_rescaling(data_3d[s],exthdr,tif16_name,header_inf['amax'],header_inf['amin'],y,x)


    
def readctfstar():
  if(os.path.isfile('micrographs_all_gctf.star')):
    ctfstar = EMAN2star.StarFile('micrographs_all_gctf.star')
    ctfstar.readfile()
    ctf={}
    for idx, val in enumerate(ctfstar['rlnMicrographName']):
        ctf[val]= '{:2.1f}um'.format(   ctfstar['rlnDefocusU'][idx]/10000)
        print (val,ctf[val])
    return ctf
  else: return {}
def main():
    print(os.path.realpath(__file__))
    print('dumpMRC v2018-10-29       <zhijie.li@utoronto.ca>')
    print('###############################################################################################')

    
    parser = argparse.ArgumentParser(description='Dump MRC file to tiff slices')
    parser.add_argument("mrc_files",metavar='.mrc files',nargs='+',type=str, help="The .mrc files to process.")
    parser.add_argument("--sigma", help="rescale tiff8 image by sigma",                        type=float, metavar="N")
    parser.add_argument("-show",action='store_true',default=False, help="Display the first frame only")
    parser.add_argument("-info",action='store_true',default=False, help="Display the header info only")
    parser.add_argument("-csv",action='store_true',default=False, help="save the map's data to a csv file")
    parser.add_argument("-update_stats",action='store_true',default=False, help="update stats")
    parser.add_argument("-update_stats_frame0",action='store_true',default=False, help="update stats based on first 2D slice")
    parser.add_argument("--gain", help="Gain reference (FEI raw format)", type=str, metavar="N")
    parser.add_argument("--bin", help="Binning factor for uint8.tif, needs to be positive integer", type=int,default=2, metavar="N")
    

    args=parser.parse_args()                           
    ctf=readctfstar()

    for filename in args.mrc_files:
      proc(filename,args,ctf)

      
main()