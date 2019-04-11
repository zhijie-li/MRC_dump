#!/usr/bin/env python2.7
from __future__ import print_function

import argparse

import re
import os
import sys
import yaml
import csv
import numpy as np
import EMAN2star 
import mrc
from EM_tiff import save_tiff8, bindata

def save_safe_yaml(data,filename):
    yamltext=yaml.safe_dump(data , default_flow_style=False)
    with open(filename, 'w') as outfile:
        yaml.safe_dump(data, outfile, default_flow_style=False)
    return yamltext
def save_yaml(data,filename):
    with open(filename, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)
        

def draw(data_3d,u,a,step=1,size=100,binning=1):
  
  r=size/2/binning;
  y,x=data_3d.shape
  #print(data_3d.shape,x,y)
  u[0]=int(u[0]/binning)
  u[1]=int(u[1]/binning)

  x1,x2=0,x-1
  y1,y2=0,y-1
  if(u[0]-r>0):
    x1=u[0]-r
  if(u[0]+r<x):
    x2=u[0]+r
  if(u[1]-r>0):
    y1=u[1]-r
  if(u[1]+r<y):
    y2=u[1]+r
  for i in range(x1,x2):
    if(i%step==0 and abs(u[0]-i)>r/2):
     # data_3d[u[1]][i]=a
      data_3d[y1][i]=a
      #data_3d[y1+1][i]=a
      data_3d[y2][i]=a
      #data_3d[y2-1][i]=a
  for i in range(y1,y2):
    if(i%step==0 and abs(u[1]-i)>r/2):
      #data_3d[i][u[0]]=a
      data_3d[i][x2]=a
      #data_3d[i][x1+1]=a
      #data_3d[i][x2-1]=a
      data_3d[i][x1]=a
    
        

        

def plot_mrc(filename,particle_used,args, binning=4,sigma=3,marker_size=120):
    if args.bin is not None:
      binning=args.bin
    if args.sigma is not None:
      sigma=args.sigma
    if args.marker is not None:
      marker_size=args.marker
  
    file_size=os.path.getsize(filename)
    header_inf,exthdr,data_1d=mrc.read_MRC2014(filename)
    x,y,z=header_inf['c'],header_inf['r'],header_inf['s']
    data_3d=data_1d.reshape(z,y,x)
          

    
    basename=re.search('[^/]+(?=.mrc$)',filename).group(0)
    mrcname=re.search('[^/]+.mrc$',filename).group(0)
    numb=len(particle_used[filename])
    tif8_name=('{:04d}'.format(numb))+basename+'.'

#    print (filename,basename,mrcname)

    batch='links_{:03d}-{:03d}/'.format(int(numb/20)*20,int(numb/20)*20+20)
    if(not os.path.isdir(batch)):
        os.mkdir(batch)
    os.symlink(filename,batch+mrcname) #make symbolic link to the original mrc


#    if(header_inf['nz']>1):
#      tif8_name+='{:03d}_'.format(s)
    #tif8_name+='m{:1.1f}_'.format(header_inf['amean'])
    tif8_name_all=tif8_name+'_a.'
    tif8_name_used=tif8_name+'_u.'
    
    tif8_name+='tif'
    tif8_name_all+='tif'
    tif8_name_used+='tif'
#binning
    binned=data_3d[0]  
    if(binning!=1):
      binned=bindata(data_3d[0],binning)

 #   if header_inf['arms']<=0: #unreliable header inf, update
 #     header_inf['amin']=data_3d[0].min()
 #     header_inf['amax']=data_3d[0].max()
 #     header_inf['amean']=data_3d[0].mean()
 #     header_inf['arms']=data_3d[0].std()
 #     print ("Updated based on frame 0: min {} max {} mean {} rms {}".format(header_inf['amin'],header_inf['amax'],header_inf['amean'],header_inf['arms']))
    newmax=binmax=  binned.max()
    newmin=binmin=  binned.min()
    binrms=  binned.std()
    binmean=  binned.mean()
    if(binmean+sigma*binrms<binmax):
      newmax=binmean+sigma*binrms
    if(binmean-sigma*binrms>binmin):
      newmin=binmean-sigma*binrms

    ybin,xbin=binned.shape
    save_tiff8(binned,exthdr,batch+tif8_name,newmax,newmin,ybin,xbin)
    
    for k in particle_used[filename]:
        
        draw(binned,k,newmax,step=1,size=marker_size,binning=binning)

    save_tiff8(binned,exthdr,batch+tif8_name_all,newmax,newmin,ybin,xbin)




def proc(filename,args):
  import sqlite3
  conn = sqlite3.connect(filename)
  c = conn.cursor()
  
  t = ('IMAGE_ASSET_ID')

  #c.execute('SELECT * FROM {} WHERE {}=?'.format('IMAGE_ASSETS', t),('',))
  c.execute('SELECT * FROM {}'.format('IMAGE_ASSETS'))
  rows = c.fetchall()
 
  for row in rows:
        print(row)
   
  conn.close()
  

def main():
  print('read_cisTEM_db v2019-01-12       <zhijie.li@utoronto.ca>')
  parser = argparse.ArgumentParser(description='process cisTEM sqlite3 db file to get particle positions')
  parser.add_argument("db_files",metavar='.db files',nargs='+',type=str, help="The .db file to process.")
  #parser.add_argument("--marker", help="Marker size, in pixels",                      type=int, metavar="N")
  #parser.add_argument("--bin", help="Binning factor, needs to be positive integer",                      type=int, metavar="N")
  #parser.add_argument("--sigma", help="Rescaling, how many sigmas",                      type=float, metavar="N")
  #parser.add_argument("-csv",action='store_true',default=False, help="export csv file")
  #parser.add_argument("-yaml",action='store_true',default=False, help="export yaml file ")
  #parser.add_argument("-txt",action='store_true',default=False, help="export txt file")
  #parser.add_argument("-tif",action='store_true',default=False, help="export txt file")
  
  args=parser.parse_args()                           

  #if ((args.marker is not None ) or (  args.sigma is not None ) or ( args.bin is not None )):
  #  args.tif=True
  for filename in args.db_files:
    proc(filename,args)

#########################################################  
main()                                                   


