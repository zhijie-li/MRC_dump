#!/usr/bin/env python2.7

from __future__ import print_function
import common

import argparse
import subprocess
import re
import os
import sys
import yaml
import csv
import numpy as np
import EMAN2star 
import mrc
from EM_tiff import save_tiff8, bindata
        

def shapeary(x,y,a,u,line=False):
  ary0=[[-x/2,b] for b in range(int(-y/2),int(y/2))]
  ary0.extend(  [    [c,y/2] for c in range(int(-y/2),int(y/2))])
  ary0.extend(  [    [x/2,d] for d in range(int(-y/2),int(y/2))])
  if(line):
    ary0.extend(  [    [0,e] for e in range(0,int(y/2))])
    
  #print(ary0)
  rot=np.complex(np.cos(a),np.sin(a))
  ary=[]
  for p in ary0:
    rotp=np.complex(p[0],p[1])*rot
    ary.append([int(np.real(rotp)+u[0]),int(np.imag(rotp)+u[1])])
  return ary
  
def draw2(data_3d,u,a,step=1,sizex=100,sizey=100,binning=1,line=False):
  
  
  rx=sizex/2/binning;
  ry=sizey/2/binning;
  y,x=data_3d.shape
  #print(data_3d.shape,x,y)
  px=int(u[0]/binning)
  py=int(u[1]/binning)
  angle=u[2]
  
  shape_ary=shapeary(rx,ry,angle,[px,py],line=line)
  for p in shape_ary:
    if(0<=p[1]<y and 0<=p[0]<x):
      data_3d[p[1]][p[0]]=a
  
  
  
  
  
def draw(data_3d,u,a,step=1,sizex=100,sizey=100,binning=1,angle=0):
  
  rx=sizex/2/binning;
  ry=sizey/2/binning;
  r=int(rx+ry)/2
  y,x=data_3d.shape
  #print(data_3d.shape,x,y)
  px=int(u[0]/binning)
  py=int(u[1]/binning)

  x1,x2=0,x-1
  y1,y2=0,y-1
  if(px-r>0):
    x1=px-r
  if(px+r<x):
    x2=px+r
  if(py-r>0):
    y1=py-r
  if(py+r<y):
    y2=py+r
  for i in range(x1,x2):
    if(i%step==0 and abs(px-i)>r/2):
      data_3d[y1][i]=a
      data_3d[y2][i]=a
  for i in range(y1,y2):
    if(i%step==0 and abs(py-i)>r/2):
      data_3d[i][x2]=a
      data_3d[i][x1]=a
    
        


def write_relion_manual_pick_star(fn,XYary):
    #print(XYary)
    header='''
#cshack.py
data_
loop_ 
_rlnCoordinateX #1 
_rlnCoordinateY #2 
_rlnClassNumber #3 
_rlnAnglePsi #4 
_rlnAutopickFigureOfMerit #5 
'''
    body=''
    for xy in XYary:
      body+='{:12.6f}   {:12.6f}         -999   -999.00000   -999.00000\n'.format(xy[0],xy[1])
    file = open(fn,'w') 
    file.write(header) 
    file.write(body) 
    file.close



def find_dir(filename):
  work_dir = os.getcwd()
  if(not( os.path.isfile(filename) or os.path.isdir(filename))):
    print ("not found",filename)
  else: 
    filename=os.path.abspath(filename)
#    print (filename)

  file_dir=filename
  mrcname=''
  job_dir    =None
  project_dir=None
  (root,ext)=('','')

  if(os.path.isfile(filename)):
    file_dir = os.path.dirname(os.path.abspath(filename))
    mrcname=os.path.basename(filename)
    (root,ext)=os.path.splitext(mrcname)
  
  job_dir    =re.search('^.+\/P\d+\/J\d+(?=\/*)',file_dir) #20190107_1230.mrc
  if(job_dir != None and os.path.isdir(job_dir.group(0))):
    job_dir=job_dir.group(0)
    project_dir=re.search('^.+\/P\d+(?=\/)',job_dir) #20190107_1230.mrc
    if(project_dir != None and os.path.isdir(project_dir.group(0))):
      project_dir = project_dir.group(0)
  titan_datetime=''
  FoilHole=''
  exposure_pos=''
  shortroot=''
  if(ext=='.mrc' or ext=='.mrcs' or ext=='.tif' or ext=='.MRC' or ext=='.MRCS' or ext=='.TIF'):
    titan_datetime=''
    FoilHole=''
    exposure_pos=''
    shortroot=root

    titan_datetime=re.search('201\d\d\d\d\d_\d\d\d\d(?=.Fractions)',root)
    if(titan_datetime is not None): 
      titan_datetime=titan_datetime.group(0) #20190107_1230.mrc
      FoilHole==''
      exposure_pos=''
      #print (root)
      F=re.search('(?<=FoilHole_)\d+(?=_)',root)
      if(F is not None):
        FoilHole=F.group(0)
      else:
        print ("!!!!",root)
      e=re.search('(?<=Data_)\d+(?=_)',root)
      if(e is not None):
        exposure_pos=e.group(0)
        
      else:
        print ("!!!!>>",root)
      shortroot=titan_datetime+'_'+FoilHole[4:]+exposure_pos[4:]
  return (mrcname,titan_datetime,root,ext,work_dir,job_dir,project_dir,file_dir,FoilHole,exposure_pos,shortroot)


def gen_manualpick_star(filename,particle_used):
        

      (mrcname,titan_datetime,file_root,file_ext,work_dir,job_dir,proj_dir,filedir,FoilHole,exposure_pos,shortroot)=find_dir(filename)
      fn= shortroot +'_manualpick.star'
      (batch,numb)=make_batchdir(filename,particle_used) #make subdir

      write_relion_manual_pick_star(batch+fn,particle_used)

def make_batchdir(filename,particle_used,b=10):
    numb=len(particle_used)
    batch='links_{:03d}-{:03d}/'.format(int(numb/b)*b,int(numb/b)*b+b)

    if(not os.path.isdir(batch)):
        os.mkdir(batch)
    
    return (batch,numb)
  
def plot_mrc(filename,particle_used,reduced_particles,args, binning=4,sigma=3,marker_size=120):
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
    

    
    (mrcname,titan_datetime,file_root,file_ext,work_dir,job_dir,proj_dir,filedir,FoilHole,exposure_pos,shortroot)=find_dir(filename)

    (batch,numb)=('',0)
    if(len(reduced_particles)==0):
      (batch,numb)=make_batchdir(filename,particle_used)
    else:
      (batch,numb)=make_batchdir(filename,reduced_particles)
    

    
    shortname=shortroot+'.mrc'
    if(args.makelink_to_micrograph):
      os.symlink(filename,batch+shortname) #make symbolic link to the original mrc
    
    tif8_name=('{:04d}_'.format(numb))+shortroot+'.'
    tif8_name_all=tif8_name+'_a.'
    #tif8_name_used=tif8_name+'_u.'
   
    tif8_name+='tif'
    tif8_name_all+='tif'
    #tif8_name_used+='tif'
#binning
    binned=data_3d[0]  
    if(binning!=1):
      binned=bindata(data_3d[0],binning) #this will change the original image

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
    
    for k in particle_used:
        
        draw2(binned,k,newmax,step=1,sizex=marker_size*2,sizey=marker_size*2,binning=binning)
    if(len(reduced_particles)>0):
      for k in reduced_particles:
        draw2(binned,k,newmax,step=1,sizex=marker_size*2,sizey=marker_size*2,binning=binning,line=True)

    save_tiff8(binned,exthdr,batch+tif8_name_all,newmax,newmin,ybin,xbin)


  
def print_content(a,r=1):
  dt=a.dtype  
  
  for j in range(len(a[0].dtype)):  
    print ('{:16s} {:30s}'.format(dt[j],dt.names[j]),sep='\t',end='\t')
    
    for i in range (r):  
      print(a[i][dt.names[j]],sep='\t',end='\t')
    
    print()    

def print_content_file(a,fn,r=1):
  dt=a.dtype  
  log = open(fn, "w")

  for j in range(len(a[0].dtype)):  
    print ('{:16s} {:30s}'.format(dt[j],dt.names[j]),sep='\t',end='\t',file = log)
    
    for i in range (r):  
      print(a[i][dt.names[j]],sep='\t',end='\t',file = log)
    
    print(file = log)    
  log.close()
  
def guess_job_type(filename='job.log'):
  jobtype=None
  if(os.path.isfile(filename)):
    with open(filename,'r') as log:
      for line in log:
        end=line.find('cryosparc2_compute.jobs.jobregister')
        if end>0:
          jobtype=re.search('[^\s]+run[^\s]*',line[:end]).group(0)
          break

      log.close()
  if(jobtype==None): #test if it is Gctf result
    t=0
    if(os.path.isfile(filename)):
      with open(filename,'r') as log:
        for line in log:
          if(line.find('bin/Gctf')>0):
            t+=1
          if(line.find('Kai Zhang@MRC Laboratory of Molecular Biology')>0):
            t+=1
          if(line.find('kzhang@mrc-lmb.cam.ac.uk')>0):
            t+=1
          if(t==3):
            jobtype='Gctf'
            break
        log.close()
  return jobtype

def analyze_job_structure(filename):
  
  (mrcname,titan_datetime,file_root,file_ext,work_dir,job_dir,proj_dir,filedir,FoilHole,exposure_pos,shortroot)=find_dir(filename)
  job_type=None
  J=None
  P=None
  if(job_dir != None):
    job_type=guess_job_type(job_dir+'/job.log')
    print('Job type: <{}>'.format(job_type))
    J=os.path.basename(job_dir)
  if(proj_dir != None):
    P=os.path.basename(proj_dir)  
  return (mrcname,titan_datetime,file_root,file_ext,work_dir,job_dir,proj_dir,filedir,job_type, P,J,shortroot)


def get_dir_size(directory):
  p = subprocess.Popen(['du', '-hs',directory], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output, err = p.communicate(b"")
  rc = p.returncode
  return output

def clean_job(dirname):
  
  (mrcname,titan_datetime,file_root,file_ext,work_dir,job_dir,proj_dir,filedir,job_type, P,J,shortroot)=analyze_job_structure(dirname)
  if(job_dir==None):
    return
  print("\n",get_dir_size(job_dir).split()[0],job_dir,'\tType: <{}>'.format(job_type) )
  files = [f for f in os.listdir( job_dir ) if os.path.isfile(job_dir+'/'+ f) ]

  if(job_type =='class2D.run' or job_type =='refine.run' or job_type=='local_refine.run'):
    prefix='(?<=cryosparc_{:s}_{:s}_)'.format(P,J)
    counts=set()
    for f in files:
      m=re.search(prefix+'\d+(?=_)',f)
      if m is not None:
        counts.add(m.group(0))
    maxi=None
    if len(counts)>0:
      maxi=max(counts)
      #print (counts,maxi)      
      for f in files:
        m=re.search(prefix+'\d+(?=_)',f)
        if m is not None:
          if (m.group(0) != maxi):
            print(job_dir+'/'+f)
            os.remove(job_dir+'/'+f)
            pass
          else: 
            #print(f)
            pass
  if(job_type =='Gctf'):
    img_dir= job_dir+'/motioncorrected/'
    if(os.path.isdir( job_dir+'/motioncorrected/')):
      img_dir= job_dir+'/motioncorrected/'
    if(os.path.isdir( job_dir+'/imported/')):
      img_dir= job_dir+'/imported/'  
    files = [f for f in os.listdir( img_dir) if os.path.isfile(img_dir+'/'+ f) ]
    for f in files:
      m=re.search('.ctf$',f)
      if m is not None:
        os.remove(img_dir+'/'+f)
        print(img_dir+'/'+f)
     
        
##################################################    

def proc(filename,args,particle_used):
      
  a=np.load(filename)
  
  
  (mrcname,titan_datetime,file_root,file_ext,work_dir,job_dir,proj_dir,filedir,FoilHole,exposure_pos,shortroot)=find_dir(filename)
  print("\ncurrent dir: " + work_dir)
  
  prefix=os.path.basename(proj_dir)+'_'+os.path.basename(job_dir)+'_'+file_root
  
  print("Project dir: [" + os.path.basename(proj_dir)+"]\tJob: [" + os.path.basename(job_dir)+ "]\tBasename: [" + mrcname + "]",file_root,file_ext)

  
  print("\nInput file: " + os.path.abspath(filename),'{} x {}'.format(len(a[0].dtype),len(a)),'Stride: {}\n'.format(a.strides),sep='\t')

  
  np.set_printoptions(threshold=100000000)


  print_content(a)
  print_content_file(a,file_root+'_cs.txt')

  print(a.shape)
  #exit()
  
  particle_file_set=set();
                              

  
  
  if (args.delCTF): #a very specific usage. Sometimes csparc keeps deleted ctf job in .cs file, causing new jobs to fail
    toextract=[]
    dt=a.dtype  
    
    for j in range(len(a[0].dtype)):  
      name=dt.names[j]
      
      if(name[:3]=='ctf'):
        print (j,name,"deleted")
      else:
        toextract.append(name)
    print(toextract)
    b=a[toextract]
    print('==============modified======================')
    print(a.shape,b.shape)
    print_content(b)
    np.save(prefix+'modified.cs',b)
    
    
      
  
  if (args.locate_images):
    dt=a.dtype  
    
    if('location/micrograph_path' not in dt.names):
      print("Input .cs file does not contain ['location/micrograph_path'] ")
    else:  
      for idx in range(len(a)):

          
            
        #proj_dir=os.path.realpath('../')
        
        symbolic_path=proj_dir + '/' + a[idx]['location/micrograph_path']
        realp=None
        existing_micrograph=None
        if(os.path.isfile (symbolic_path)):
          existing_micrograph=realp=os.path.realpath(symbolic_path)
          
                  
        #FoilHole_1110033_Data_1109469_1109470_20190102_1555_Fractions_rigid_aligned.mrc
        mrc_root=''
        mrc_root_re    =re.search('FoilHole_\w+_Fractions',a[idx]['location/micrograph_path'])
        if (mrc_root_re is not None):
          mrc_root = mrc_root_re.group(0)
          #print(mrc_root)
        
        
        
        fract_x=a[idx]['location/center_x_frac']
        fract_y=a[idx]['location/center_y_frac']
        [dimension_y,dimension_x]=a[idx]['location/micrograph_shape']
        
        coordx=int(dimension_x*fract_x)
        coordy=int(dimension_y*fract_y)
        
        angle=0
        ncc_score=0
        pick_stats_power =0
        if('pick_stats/angle_rad' in dt.names):
          angle=a[idx]['pick_stats/angle_rad']
          ncc_score =a[idx]['pick_stats/ncc_score']
          pick_stats_power =a[idx]['pick_stats/power']
          pick_stats_template_idx=a[idx]['pick_stats/template_idx']
          
        #if('alignments2D/alpha' in dt.names):
        #  angle=a[idx]['alignments2D/alpha']      
      
      
        if( mrc_root!='' and not mrc_root in particle_used):
            particle_used[mrc_root]={}
            
            particle_used[mrc_root]['particles']=[]
            particle_used[mrc_root]['particles_cs_line']=[]
            
            particle_used[mrc_root]['reduced_particles']=[]
            particle_used[mrc_root]['reduced_particles_cs_line']=[]
            particle_used[mrc_root]['micrograph']=None


        if( mrc_root!=''):
          particle_used[mrc_root]['particles'].append([coordx,coordy,angle,ncc_score,pick_stats_power,pick_stats_template_idx,realp,job_dir])
          particle_used[mrc_root]['particles_cs_line'].append(a[idx])

#          if(args.cs1_as_ref):
#              if(args.cs_files[0]==filename):
#                ref_micrograph_list[mrc_root]=[a[idx]['location/micrograph_path'],a[idx]['location/micrograph_uid']]
#                #print (a[idx]['location/micrograph_path'])
#                pass

#              else:
#                if(mrc_root in ref_micrograph_list.keys()):
#                  print (ref_micrograph_list[mrc_root])
#                  a[idx]['location/micrograph_path']=ref_micrograph_list[mrc_root][0]
#                  a[idx]['location/micrograph_uid']=ref_micrograph_list[mrc_root][1]
#                  print (a[idx]['location/micrograph_path'])
#

        else:
          print ("!!!!!!",a[idx]['location/micrograph_path'])
        if((existing_micrograph is not None) and (particle_used[mrc_root]['micrograph'] is None)):
          particle_used[mrc_root]['micrograph']=existing_micrograph


             
    


  if args.txt:
    with open(prefix+'.txt','w') as txt:
      print(a,'\n',a.dtype,file=txt)
      pass
  if args.csv:
    with open(prefix+'.csv', 'w') as csv_f:
        csv_f.write(','.join(a.dtype.names) + '\n')
        np.savetxt(csv_f, a, '%s', ',')


##################################################  
  
def main():
  print('cryoSPARC Tools v2019-02-13       <zhijie.li@utoronto.ca>')
  parser = argparse.ArgumentParser(description='Analyzing cryoSPARC results. Need passthrough.cs')
  parser.add_argument("cs_files",metavar='.cs files',nargs='+',type=str, help="The .cs files to process.")
  parser.add_argument("--marker", help="Marker size, in pixels",
                      type=int, metavar="N")
  parser.add_argument("--bin", help="Binning factor, needs to be positive integer",
                      type=int, metavar="N")
  parser.add_argument("--sigma", help="Rescaling, how many sigmas",
                      type=float, metavar="N")
  parser.add_argument("-delCTF",action='store_true',default=False, help="remove CTF columns and save 'modified.cs'")

  parser.add_argument("-csv",action='store_true',default=False, help="export csv file")
  parser.add_argument("-yaml",action='store_true',default=False, help="export yaml file ")
  parser.add_argument("-txt",action='store_true',default=False, help="export txt file")
  parser.add_argument("-tif",action='store_true',default=False, help="export tif file")
  parser.add_argument("-locate_images",action='store_true',default=False, help="Try to locate images")
  parser.add_argument("-cleanjobs",action='store_true',default=False, help="delete intermediate files in 2D classification runs")
  
  parser.add_argument("-reduce_particles",action='store_true',default=False, help="analyze multiple .cs files and summarize particle locations")
  parser.add_argument("-makelink_to_micrograph",action='store_true',default=False, help="make symbolic links to the micrographs in the plot directory")
#  parser.add_argument("-cs1_as_ref",action='store_true',default=False, help="make new cs file, using the first supplied reference passthrough.cs as reference")
  parser.add_argument("--hack_pp",type=str, metavar="xxx.cs", help="consolidate particles coordinates from multiple cs files, then save them in the supplied particle picking cs file")

  args=parser.parse_args()
  ####plot particles on images###
  
  if(args.hack_pp is not None): #depends on these data
    args.reduce_particles=True
    args.locate_images=True


  if ((args.marker is not None ) or (  args.sigma is not None ) or ( args.bin is not None )):
    args.tif=True
  
  
  particle_used={}
  #ref_micrograph_list={}

  for filename in args.cs_files:
    if args.cleanjobs:
      clean_job(filename)
      
    else:
      
      proc(filename,args,particle_used)  #,  ref_micrograph_list)

      
      

  if args.yaml:
    common.save_safe_yaml(particle_used,'particles.yaml')

      
  if( args.reduce_particles):
    acount=0
    for mrc_root in particle_used.keys():
      reject=np.zeros(len(particle_used[mrc_root]['particles']))
      for i,u in enumerate(particle_used[mrc_root]['particles']):
        if(reject[i]!=1):
          for j in range(i+1,len(particle_used[mrc_root]['particles'])):
            if(reject[j]!=1):
              a=particle_used[mrc_root]['particles'][i]
              b=particle_used[mrc_root]['particles'][j]
              if((a[0]-b[0])**2+(a[1]-b[1])**2 < 5000): #70 pixel cutoff
                reject[j]=1
      for i,u in enumerate(particle_used[mrc_root]['particles']):
        if(reject[i]!=1):
          acount+=1
          particle_used[mrc_root]['reduced_particles'].append(particle_used[mrc_root]['particles'][i])
          
          particle_used[mrc_root]['reduced_particles_cs_line'].append(particle_used[mrc_root]['particles_cs_line'][i])
          #print( len(particle_used[mrc_root]['reduced_particles_cs_line']))
    print("Aggregated particles: {}".format(acount))  
    
  if(args.hack_pp is not None):
    a=np.load(args.hack_pp)
    microlist={}
    ml=[]
    for l in a:
      if(l['location/micrograph_uid'] not in microlist.keys()):
        ml.append(l['location/micrograph_uid'])
        microlist[l['location/micrograph_uid']]={}
        
        #microlist[l['location/micrograph_uid']]['standard']=l
        #microlist[l['location/micrograph_uid']]['path']=l['location/micrograph_path']
        
        microlist[l['location/micrograph_uid']]['mrc_root']=''
        r=re.search('FoilHole_\w+_Fractions',l['location/micrograph_path'])
        if(r is not None):
          microlist[l['location/micrograph_uid']]['mrc_root']=r.group(0)
        else:
          print(l['location/micrograph_path'])
        #  print (microlist[l['location/micrograph_uid']]['mrc_root'])
        microlist[l['location/micrograph_uid']]['particle_record']=[]
      
      microlist[l['location/micrograph_uid']]['particle_record'].append(l)
  
    
    all=0
    for k in ml:
      all+=len(microlist[k]['particle_record'])
      #print(k,len(microlist[k]['particle_record']))
    print ("{} micrographs in picking job, {} particle records.".format(len(microlist.keys()),all))
    print ("{} micrographs in aggregated particle cs. ".format(len(particle_used.keys())))
    b=a[0]
    
    #print (b.shape,b)
    

    
    
    c=0
    for micrograph_uid in ml:
      idx=-1

      mrc_root=microlist[micrograph_uid]['mrc_root']
      #print (mrc_root)
      if(mrc_root in particle_used.keys()):
        
        for p in particle_used[mrc_root]['reduced_particles_cs_line']:
          
          idx+=1
          blank=microlist[micrograph_uid]['particle_record'][idx] #assume initial picking, each image should have thaousands of particles
          #modify only the particle picking data:
          blank['location/center_x_frac'] =p['location/center_x_frac']
          blank['location/center_y_frac'] =p['location/center_y_frac']
          blank['pick_stats/ncc_score'	 ]=p['pick_stats/ncc_score'	 ] 
          blank['pick_stats/power'	     ]=p['pick_stats/power'	     ] 
          blank['pick_stats/template_idx']=p['pick_stats/template_idx']
          blank['pick_stats/angle_rad'   ]=p['pick_stats/angle_rad'   ]
          if(c==0):
            b= blank
          else:
            b= np.append(b,blank)
            
          
          #print (c,idx,b.shape,b)
          c+=1
      if(idx==-1)        :
          print("No particles in file {}".format(mrc_root))
      #if(c>=3):        break        
            
    print_content(b,r=1)
    print (b.shape)
    np.save('reduced.cs',b)


 

  if args.tif:
    for mrc_root in particle_used.keys():
      
      if(particle_used[mrc_root]['micrograph'] is not None):
        micrograph=particle_used[mrc_root]['micrograph']
        plot_mrc(micrograph,particle_used[mrc_root]['particles'],particle_used[mrc_root]['reduced_particles'],args)

      if( args.reduce_particles):    
        gen_manualpick_star(micrograph,particle_used[mrc_root]['reduced_particles'])
      else:
        gen_manualpick_star(micrograph,particle_used[mrc_root]['particles'])

      
  #for micrograph in particle_used.keys():
  #    print(len(particle_used[micrograph]),micrograph)
  
#########################################################  
main()                                                   


