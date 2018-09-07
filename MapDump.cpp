#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <time.h>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <dirent.h>

#define PNG_DEBUG 3
#include <png.h>



int findparam(int argc,char *argv[], char *infile[], char *outdir, char *prefix,  int * statframe, int * bin_factor, int * vbin,float * sigma); //find the input file list


void statarray(float * ary, int size,float *min,float *max,float *rms,float *mean);


int readmap_header(const char *filename0, char *header, int *nsymbt, int *c, int *r, int *s , float *amin, float *amax, float *amean, float *arms, int * mode, int * datasize);

float ** ReadMap(const char *filename0, int nsymbt,int c, int r,int s ,int mode); 

int read_slice(const char *filename0, float *array,int nsymbt,int c, int r,int s ,int mode) ;
//int readmap(const char *filename0, map_data *array,int nsymbt,int c, int r, int s, int datasize)

//int read_whole_array(const char *filename0, float *array,int nsymbt,int c, int r,int s, int datasize) ;//read the whole array (1.6G for K2 in real4)
//template <class map_data>
//int read_slice      (const char *filename0, map_data *array,int nsymbt,int c, int r,int slice, int datasize) ;//slice start from 0

int dumparray(const char *filename0, char *png, int arraysize);

int bin(float * array, float * array_bin,int x, int y, int bin_factor);

//int addtoheader(char * out_header, float max,float min,float mean,float rms, int num_map);
//int addtoheader(char * out_header, float max,float min,float mean,float rms, char * label);

//void write_png_file(char* file_name,int width,int height,png_byte color_type, png_byte bit_depth, png_bytep * row_pointers);

int rescale(float *, float,float,int);
int writeImage(const char* filename, int width, int height, float *buffer);

struct headerinf{
    int c;
    int r;
    int s;
    int nsymbt;
    };



int main(int argc,char *argv[])
{ 
  float t1,t2;
  t1=clock();
//  printf ("\nMapDump: dump CCP4/MRC map/movie/image files to PNG.\n Input must have .mrcs extension\n\n");

//parse input
  char **infile = new char*[2048]; 
  char *outdir=new char[1024];
  char *prefix=new char[1024];
  for (int i = 0; i < 2048; ++i) {    infile[i] = new char[1024];  }
	int statframe=0;
	int bin_factor=0;
	int vbin=0;
	float sigma=1.5; //6-sigma cutoff 
  int num_infile= findparam(argc,argv, infile, outdir, prefix, &statframe, &bin_factor, &vbin, &sigma);
  std::cout<<"Num of files to dump: "<<num_infile<<'\n';
	for(int f=0;f<num_infile;f++){
		printf("%d %s\n",f+1, infile[f]);
	}
  std::cout<<"Num of slices used for statistics: "<<statframe<<'\n'<<'\n';
  std::cout<<"Binning: "<<bin_factor<<'\n'<<'\n';

//
	for(int f=0;f<num_infile;f++){

  char *header=new char[1024];
  char *symmop=new char[10240];//more than enough
  
  
    int   c=0,r=0,s=0;
    int   nsymbt=0;

    char * fn=new char[1024];
    strcpy(fn, infile[f]);
		std::cout<<"\n\n====Processing file ["<<fn<<"]===========\n";

    float amin=0,amax=0,amean=0,arms=0;
    int datasize=4; //default real4
    int mode=2; //default
    readmap_header(const_cast<char*>(fn), header, &nsymbt, &c,&r,&s,&amin,&amax,&amean,&arms, &mode, &datasize);//get size of the first file

    int arraysize=c*r*s;
    int slice_size=c*r;


    float **  array= ReadMap(const_cast<char*>(fn), nsymbt,c,r,s, mode); //read one slice



		if (arms <=0 && statframe==0){
			printf("\nThe RMS found in header is zero or negative. updatinging with frame 0\n");
			statframe=1;
			}

		
    if(statframe>0){
              float amin_n=0;
              float amax_n=0;
              float arms_n=0;
              float amean_n=0;
		
		int statcont=0;
    for(int z=0;z<s && z<statframe;z++){ //statframe is upper limit of the frames to quant

              printf ("\n\nUpdating stats using slice %d....\n\n",z);
              float amin_s;
              float amax_s;
              float arms_s;
              float amean_s;
              statarray(array[z],c*r,&amin_s,&amax_s,&arms_s,&amean_s);
              amin_n+=amin_s;
              amax_n+=amax_s;
              arms_n+=arms_s;
              amean_n+=amean_s;
              statcont++;

            }
              amin_n/=statcont;
              amax_n/=statcont;
              arms_n/=statcont;
              amean_n/=statcont;
			arms=arms_n;
			amin=amin_n;
			amax=amax_n;
			amean=amean_n;			
			printf("\nUpdated min %8.3f max %8.3f mean %8.3f RMS %8.3f\n\n", amin, amax,  amean, arms);
	}
    for(int z=0;z<s-vbin+1;z++){     
        char  png_name[1024];

        float highcut=amean+sigma*arms; //default sigma=6
        if(highcut>amax){highcut=amax;}
            
        float scale=128/(sigma*arms)/bin_factor/bin_factor/vbin;
        float shift=-amean*bin_factor*bin_factor*vbin;
        
        sprintf(png_name,"%s.slice_%04d.png",fn,z);

        
        
        
            
         char title[1024]; 
         
         printf ("saving image to file <%s> scaling factor %f\n", png_name, scale);
         
         float * outframe=array[z];
         if(vbin>1 ){
         	float * Vbin_fp=new float  [slice_size];
         	outframe= Vbin_fp;
         	
         	for(int i=0;i<slice_size;i++){
         		*(Vbin_fp+i)=0;
         		for(int j=0;j<vbin;j++){
         			*(Vbin_fp+i)+=*(array[z+j]+i);
         		}
         	}
         		
         	z+=vbin-1;
         	}
         	

        if(bin_factor<=0){
					rescale(array[z],scale,shift,slice_size);
					writeImage(png_name, c,r, outframe);
					}
        	else{
		          float * array_bin=new float [(c/bin_factor)*(r/bin_factor)];      

	        		bin(outframe,array_bin,c,r,bin_factor); 
	        		//std::cout<<"ok\n";
							rescale(array_bin,scale,shift,(c/bin_factor)*(r/bin_factor));
							//std::cout<<(c/bin_factor)*(r/bin_factor)<<"\n";
							writeImage(png_name, c/bin_factor,r/bin_factor, array_bin);
							delete array_bin;
							delete outframe;
        		}
        
        }


	//    dumparray("dump.ary", png,  arraysize);
    
			for(int z=0;z<s;z++){delete array[z];}

	  delete array;
	  delete symmop;
	  delete header;

	}//end of file cycle
	
  t2 = clock() - t1;
  printf ("%f seconds.\n",((float)t2)/CLOCKS_PER_SEC);
	return 0;
}

int readmap_header(const char *filename0, char *header, int *nsymbt, int *c, int *r, int *s , float *amin, float *amax, float *amean, float *arms, int *mode, int *dtsz)
{
  int fsize=0;
  std::ifstream in (filename0,std::ios::in|std::ios::binary|std::ios::ate);
	if(!in)
	{
	  printf ("file does not exist\n");
  }
  if (in.is_open()){
    fsize = in.tellg();
    in.seekg (0, std::ios::beg);
    in.read(( char*)c,4); 
    in.read(( char*)r,4);
    in.read(( char*)s,4);
    in.read(( char*)mode,4);

		if(*mode == 0)    {*dtsz=1;}
		if(*mode == 1)    {*dtsz=2;}
		if(*mode == 2)    {*dtsz=4;}
		if(*mode == 3)    {*dtsz=4;}
		if(*mode == 4)    {*dtsz=8;}
		if(*mode == 6)    {*dtsz=2;}
			
    printf("x y z: %d %d %d\n",*c,*r,*s);
    printf("mode: %d,  %d Bytes per pixel\n",*mode,*dtsz);

    int arraysize=(*c)*(*r)*(*s);

    in.seekg (76, std::ios::beg);
    in.read(( char*)amin,4); 
    in.read(( char*)amax,4);
    in.read(( char*)amean,4);
    in.seekg (216, std::ios::beg);
    in.read(( char*)arms,4); 
		
    printf("min %8.3f  max %8.3f  mean %8.3f  RMS %8.3f\n",*amin,*amax,*amean,*arms);

    
    in.seekg (92, std::ios::beg);
    

    in.read((char*)nsymbt,4);
    printf("SYMM bytes: %d\n",*nsymbt);


    printf( "Array size: %d \nSYMM block bytes: nsymbt_str %d \n\n",arraysize,*nsymbt);
		int arraysize_dt=arraysize * (*dtsz);
    if( arraysize_dt +1024+*nsymbt != fsize){      printf ("error! file size %d does not equal to data array size (%d) + 1024 + %d \n!!", fsize, arraysize_dt,*nsymbt);      }

  }
  in.close();
  return fsize;
}
/*

int readmap(const char *filename0, map_data *array,int nsymbt,int c, int r, int s, int datasize)
{
  int fsize=0;
  int arraysize;
  char *header= new char[1024];
  char *symmop= new char[nsymbt];
  std::ifstream in (filename0,std::ios::in|std::ios::binary|std::ios::ate);
	if(!in)	{	  printf ( "file does not exist\n");    return 0;  }

  if (in.is_open()){

    fsize = in.tellg();
    arraysize=c*r*s;
    if(arraysize*4+1024+nsymbt != fsize){      printf ("error! file size %d does not equal to arraysize (%d) + 1024 + %d \n!!", fsize, arraysize,nsymbt);      }

    in.seekg (0, std::ios::beg);
    in.read(( char *)header,1024);
    in.read(( char *)symmop,nsymbt);
    in.read(( char *)array,(arraysize*datasize));
  }

  in.close();
  delete header;
  delete symmop;
  return fsize;
}


int read_whole_array(const char *filename0, float *array,int nsymbt,int c, int r,int s, int datasize) //slice start from 0
{
  int fsize=0;

  std::ifstream in (filename0,std::ios::in|std::ios::binary|std::ios::ate);
	if(!in)	{	  printf ( "file does not exist\n");    return 0;  }

  if (in.is_open()){
    fsize = in.tellg();
    int start=1024+nsymbt;
    int end=1024+nsymbt+s*c*r*datasize;
    if(end>fsize){      printf ("error! file size %d <  \n!!", fsize);      }
        printf ("reading the data array of (%d)bytes x %d x %d, %d frames at bytes %12d - %12d ... ",datasize, c,r,s,start,end);

    in.seekg (start, std::ios::beg);
    in.read(( char *)array,(c*r*s*datasize));
  }

  in.close();
  return c*r*s*datasize;
}
*/
float ** ReadMap(const char *fn, int nsymbt,int c, int r,int s ,int mode) //slice start from 0
{
    float **array = new float * [s];
  	for(int z=0;z<s;z++){ 
	  	 	array[z]    =new float [c*r];
	      read_slice(const_cast<char*>(fn), array[z],nsymbt,c,r,z, mode); //read one slice
	  	}
	return array;
}


int read_slice(const char *filename0, float *array,int nsymbt,int c, int r,int s ,int mode) //slice start from 0
{
  long int fsize=0;
	int datasize=4;
		if(mode==0){datasize=1;}
		if(mode==2){datasize=4;}

  std::ifstream in (filename0,std::ios::in|std::ios::binary|std::ios::ate);
	if(!in)	{	  printf ( "file does not exist\n");    return 0;  }

  if (in.is_open()){
    fsize = in.tellg();
    long int start=1024+nsymbt+(s)*c*r*datasize;
    long int end=1024+nsymbt+(s+1)*c*r*datasize;

    if(end>fsize){      printf ("error! file size %ld <  \n!!", fsize);      }
        //printf ("reading slice at bytes %12ld - %12ld ... \n",start,end);

    in.seekg (start, std::ios::beg);
    if(mode==2){in.read(( char *)array,(c*r*datasize));}
    if(mode==0){
    	char * buffer = new char [c*r*datasize];
    	in.read(( char *)buffer,(c*r*datasize));
    	for (int i=0;i<(c*r*datasize);i++){
    		*(array+i)=(float)*(buffer+i);
    		
    		}
    		delete buffer;
    	}
    
    
  }
  in.close();
  return c*r*datasize;
}


int dumparray(const char *filename0, char *png, int arraysize)
{	
  std::ofstream out;
  out.open(filename0,std::ios::binary);
	if(!out)
	{
	  printf ("file does not exist\n");

  }
  if (out.is_open()){
    out.write(( char *)png,(arraysize));
  }

  out.close();
  return 0;
}


void statarray(float * ary, int size,float *min,float *max,float *rms,float *mean)
{	
    *min=*max=*(ary);
    *mean=*rms=0;

    for (int i=0;i<size; i++){
        if(*min>*(ary+i)){*min=*(ary+i);}
        if(*max<*(ary+i)){*max=*(ary+i);}
        *mean+=*(ary+i);
    }
    *mean=*mean/size;
    
    for (int i =0;i<size; i++){

       *rms+=(*(ary+i)-*mean)*(*(ary+i)-*mean)/size;
    }
    *rms=sqrt(*rms); // a lot less divisions but may overflow
    
    printf("Updated  min %8.3f max %8.3f mean %8.3f RMS %8.3f  \n\n",*min,*max,*mean,*rms);
    
}




//void write_png_file(char* file_name,int width,int height,png_byte color_type, png_byte bit_depth, png_bytep * row_pointers)
//{
//
//				int number_of_passes;
//				png_structp png_ptr;
//				
//				png_infop info_ptr;
//
//
//
//        /* create file */
//        FILE *fp = fopen(file_name, "wb");
//        if (!fp)
//                abort_("[write_png_file] File %s could not be opened for writing", file_name);
//
//
//        /* initialize stuff */
//        png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
//
//        if (!png_ptr)
//                abort_("[write_png_file] png_create_write_struct failed");
//
//        info_ptr = png_create_info_struct(png_ptr);
//        if (!info_ptr)
//                abort_("[write_png_file] png_create_info_struct failed");
//
//        if (setjmp(png_jmpbuf(png_ptr)))
//                abort_("[write_png_file] Error during init_io");
//
//        png_init_io(png_ptr, fp);
//
//
//        /* write header */
//        if (setjmp(png_jmpbuf(png_ptr)))
//                abort_("[write_png_file] Error during writing header");
//
//        png_set_IHDR(png_ptr, info_ptr, width, height,
//                     bit_depth, color_type, PNG_INTERLACE_NONE,
//                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
//
//        png_write_info(png_ptr, info_ptr);
//
//
//        /* write bytes */
//        if (setjmp(png_jmpbuf(png_ptr)))
//                abort_("[write_png_file] Error during writing bytes");
//
//        png_write_image(png_ptr, row_pointers);
//
//
//        /* end write */
//        if (setjmp(png_jmpbuf(png_ptr)))
//                abort_("[write_png_file] Error during end of write");
//
//        png_write_end(png_ptr, NULL);
//
//        /* cleanup heap allocation */
//        for (y=0; y<height; y++)
//                free(row_pointers[y]);
//        free(row_pointers);
//
//        fclose(fp);
//}
//


int bin(float * array, float * array_bin,int x, int y, int bin_factor)
{

	int idx=0;
    for (int i=0;i<y-y%bin_factor;i=i+bin_factor){
        for (int j=0;j<x-x%bin_factor;j=j+bin_factor){ //each element in new array
        	*(array_bin+idx)=0;
            	for (int yy=i;yy<i+bin_factor;yy++)
            		{
            			for (int xx=j;xx<j+bin_factor;xx++)
            			{
            				*(array_bin+idx)	+= *(array+xx+x*(yy));
            			}
								}
	       	idx++;
        }
    }
    //std::cout<<idx<<"\n";
} 
int rescale(float * ary, float scale,float shift, int size){
    //printf(" rescaling with scale = %f ...\n",scale);
    for(int i=0;i<size;i++)
    {
        
        *(ary+i)=(*(ary+i)+shift )* scale+127;
      	if(*(ary+i)<0){*(ary+i)=0 ;}
      	if(*(ary+i)>255){*(ary+i)=255 ;}
        
        }
    }

inline void setRGB(png_byte *ptr, float val)
{
    
    ptr[0]= round(val);
	
}

int writeImage(const char* filename, int width, int height, float *buffer)
{
	int code = 0;
	FILE *fp = NULL;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep row = NULL;
	
	// Open file for writing (binary mode)
	fp = fopen(filename, "wb");
	if (fp == NULL) {
		fprintf(stderr, "Could not open file %s for writing\n", filename);
		code = 1;
		goto finalise;
	}

	// Initialize write structure
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL) {
		fprintf(stderr, "Could not allocate write struct\n");
		code = 1;
		goto finalise;
	}

	// Initialize info structure
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL) {
		fprintf(stderr, "Could not allocate info struct\n");
		code = 1;
		goto finalise;
	}

	// Setup Exception handling
	if (setjmp(png_jmpbuf(png_ptr))) {
		fprintf(stderr, "Error during png creation\n");
		code = 1;
		goto finalise;
	}

	png_init_io(png_ptr, fp);

	// Write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, width, height,
			8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
/*
	// Set title
	if (title != NULL) {
		png_text title_text;
		title_text.compression = PNG_TEXT_COMPRESSION_NONE;
		
		strcpy(title_text.key ,"Title");
		title_text.text = title;
		png_set_text(png_ptr, info_ptr, &title_text, 1);
	}
*/
	png_write_info(png_ptr, info_ptr);

	// Allocate memory for one row (3 bytes per pixel - RGB)
	row = new png_byte[width];

	// Write image data
	int x, y;
	for (y=0 ; y<height ; y++) {
		for (x=0 ; x<width ; x++) {
			setRGB(&(row[x]), buffer[y*width + x]);
		}
		png_write_row(png_ptr, row);
	}

	// End write
	png_write_end(png_ptr, NULL);

	finalise:
	if (fp != NULL) fclose(fp);
	if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
	if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
	if (row != NULL) delete row;

	return code;
}






//modified from map_average_gcc.cpp



int findparam(int argc,char *argv[], char *infile[], char *outdir, char *prefix, int * statframe, int * bin_factor, int * vbin, float * sigma) //find the input file list
{
  char *dir=NULL;
  char *filelist=NULL;
  char *ext=new char[1024]; strcpy(ext,".mrc");
  int file_count=0;
  int dir_sig=0;  
  int fl_sig=0;
  
  
 	for (int i=1;i<argc;i++)
	{
		//printf ("%d %d %s\n",argc,i,argv[i]);
	  std::string current=argv[i];
	  
 	  if(strcmp(argv[i],"-d")==0  && i<argc-1) { dir_sig=1;      dir=argv[i+1]; i+=1;continue;}
 	  if(strcmp(argv[i],"-l")==0  && i<argc-1) { fl_sig=1;  filelist=argv[i+1]; i+=1;continue;}
 	  if(strcmp(argv[i],"-o")==0  && i<argc-1) {      strcpy(outdir,argv[i+1]); i+=1;continue;}
 	  if(strcmp(argv[i],"-prefix")==0  && i<argc-1) { strcpy(prefix,argv[i+1]); i+=1;continue;}
 	  if(strcmp(argv[i],"-statframe")==0  && i<argc-1) { 
 	  	char *tmp= new char[100]; 
 	  	strcpy(tmp,argv[i+1]);
 	  	*statframe=std::__cxx11::stoi(tmp);
// 	  	printf ("stat : %s", tmp); 
			i+=1;continue;
		}
 	  if(strcmp(argv[i],"-bin")==0  && i<argc-1) { 
 	  	char *tmp= new char[100]; 
 	  	strcpy(tmp,argv[i+1]);
 	  	*bin_factor=std::__cxx11::stoi(tmp);
			i+=1;continue;
		}
 	  if(strcmp(argv[i],"-vbin")==0  && i<argc-1) { 
 	  	char *tmp= new char[100]; 
 	  	strcpy(tmp,argv[i+1]);
 	  	*vbin=std::__cxx11::stoi(tmp);
			i+=1;continue;
		}
 	  if(strcmp(argv[i],"-sigma")==0  && i<argc-1) { 
 	  	char *tmp= new char[100]; 
 	  	strcpy(tmp,argv[i+1]);
 	  	*sigma=std::__cxx11::stof(tmp);
			i+=1;continue;
		}

 	  if ( current.length()>4 && (current.length()-current.rfind(ext) ==4 || current.length()-current.rfind(ext) ==5 )) 
 	  { 
       strcpy ( infile[file_count] ,argv[i] );
       file_count++;
//       if(file_count>=2048){printf ("cannot handle more than 2048 files\n");return file_count;}
    }
  }

  if(fl_sig==1)
    {
      std::string  line;
      std::ifstream in (filelist);
      if (in.is_open())
      {
        while ( getline (in,line) )
        {
 	        if ( line.length()>4 && line.length()-line.rfind(ext) ==4) 
 	        { 
            strcpy ( infile[file_count] , (char*)line.c_str() );
            file_count++;
        if(file_count>=2048){printf ("cannot handle more than 2048 files\n");return file_count;}
          }
        }
      in.close();
      }
    }


  if(dir_sig==1)
    {
      DIR *dir1;
      struct dirent *ent;
      if ((dir1 = opendir (dir)) != NULL) {
        while ((ent = readdir (dir1)) != NULL) {
          std::string line=ent->d_name;
 	        if ( line.length()>4 && line.length()-line.rfind(ext) ==4) 
 	        { 
            strcpy ( infile[file_count] , (char*)line.c_str() );
            file_count++;
            if(file_count>=2048){printf ("cannot handle more than 2048 files\n");return file_count;}
          }
        }
      closedir (dir1);
      } 
    }

  return file_count;
}
