#!/usr/bin/env python
from __future__ import print_function

#import struct
import os
#import sys
#import yaml
#import csv
#import numpy as np
#import scipy.stats

import argparse

from EM_tiff import save_tiff8, save_tiff16_no_rescaling
import mrc

#import scipy.misc
#import skimage
#from skimage import io
#zlib is used by the write PNG function


############################main##################
def proc(filename,args):
    
    
    file_size=os.path.getsize(filename)
    header_inf,exthdr,data_1d=mrc.read_MRC2014(filename)
    x,y,z=header_inf['c'],header_inf['r'],header_inf['s']
    data_3d=data_1d.reshape(z,y,x)
    print("Dumping {} {} bytes {} x {} x {}".format(filename,file_size,x,y,z))

    if header_inf['arms']<=0: #unreliable header inf, update
      header_inf['amin']=data_3d[0].min()
      header_inf['amax']=data_3d[0].max()
      header_inf['amean']=data_3d[0].mean()
      header_inf['arms']=data_3d[0].std()
      print ("Updated based on frame 0: min {} max {} mean {} rms {}".format(header_inf['amin'],header_inf['amax'],header_inf['amean'],header_inf['arms']))

    #print(data_3d.dtype)
    newmax=header_inf['amax']
    newmin=header_inf['amin']
    if(header_inf['amean']+6*header_inf['arms']<header_inf['amax']):
      newmax=header_inf['amean']+6*header_inf['arms']
    if(header_inf['amean']-6*header_inf['arms']>header_inf['amin']):
      newmin=header_inf['amean']-6*header_inf['arms']
      
    if(len(exthdr)>0):
      print("Extended header found: {} Bytes, using {} Bytes".format(len(exthdr),header_inf['nsymbt']))

    for s in xrange(0,header_inf['MRC_NZ']):
      tif8_name=filename+'{:03d}'.format(s)+'.uint8.tif'
      save_tiff8(data_3d[s],exthdr,tif8_name,newmax,newmin,y,x)
      tif16_name=filename+'{:03d}'.format(s)+'.tif'
      #Saving in tif is not really useful for float data. But if the data are in int8 or 16 or 32, compression can be good.
      #save_tiff16_no_rescaling(data_3d[s],exthdr,tif16_name,header_inf['amax'],header_inf['amin'],y,x)


def main():
    print('dumpMRC v2018-10-29       <zhijie.li@utoronto.ca>')
    parser = argparse.ArgumentParser(description='Dump MRC file to tiff slices')
    parser.add_argument("mrc_files",metavar='.mrc files',nargs='+',type=str, help="The .mrc files to process.")
    args=parser.parse_args()


    for filename in args.mrc_files:
      proc(filename,args)

      
main()