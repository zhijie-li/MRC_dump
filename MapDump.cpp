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



void statarray(float * ary, int size,float *min,float *max,float *rms,float *mean);

int readmap_header(const char *filename0, char *header, int *nsymbt, int *c, int *r, int *s , float *amin, float *amax, float *amean, float *arms);


int readmap(const char *filename0, float *array, int nsymbt, int c, int r, int s);

int read_slice(const char *filename0, float *array,int nsymbt,int c, int r,int slice) ;//slice start from 0
int dumparray(const char *filename0, char *png, int arraysize);

int bin(float * array, float * array_bin,int x, int y);

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
  printf ("\nMapDump: dump CCP4/MRC map/movie/image files to PNG.\n\n");

  char *header=new char[1024];
  char *symmop=new char[10240];//more than enough
  
    int   c=0,r=0,s=0;
    int   nsymbt=0;

    char * fn=new char[1024];
    strcpy(fn,argv[1]);


    float amin=0,amax=0,amean=0,arms=0;
    readmap_header(const_cast<char*>(fn), header, &nsymbt, &c,&r,&s,&amin,&amax,&amean,&arms);//get size of the first file


    int arraysize=c*r*s;
    int slice_size=c*r;
    float *array    =new float[slice_size];

    float *array_bin =new float[slice_size/4];
    float *array_bin2 =new float[slice_size/16];
    float *array_bin3 =new float[slice_size/64];
    
    





    for(int z=0;z<s;z++){

        read_slice(const_cast<char*>(fn), array,nsymbt,c,r,z); //read one slice

        if(z==0 ){ //need to updata rms based on slice 1
                printf ("\n\nUpdating stats using slice %d....\n\n",z);
              statarray(array,c*r,&amin,&amax,&arms,&amean);
            }


        
        char  png_name[1024];

        float lowcut=amean-6*arms; //default
        if(lowcut<0){lowcut=0;}
        if(lowcut<amin){lowcut=amin;}
        float highcut=amean+6*arms; //default
        if(highcut>amax){highcut=amax;}
            
        float scale=255/(highcut-lowcut);
        
        sprintf(png_name,"%s.slice_%04d.png",fn,z);

        
        
        rescale(array,scale,lowcut,c*r);
        
            
        //  array_bin=array;      
        //bin(array,array_bin,c,r); 
         //bin(array_bin,array_bin2,c/2,r/2); 
         //bin(array_bin2,array_bin3,c/4,r/4);
         char title[1024]; 
         
         printf ("saving image to file <%s>\n", png_name);
        writeImage(png_name, c,r, array);
        

        }

//    dumparray("dump.ary", png,  arraysize);
    

  delete array;


  t2 = clock() - t1;
  printf ("%f seconds.\n",((float)t2)/CLOCKS_PER_SEC);
	return 0;
}

int readmap_header(const char *filename0, char *header, int *nsymbt, int *c, int *r, int *s , float *amin, float *amax, float *amean, float *arms)
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
    printf("x y z: %d %d %d\n",*c,*r,*s);

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


    printf( "Array size: %d * 4 (float32) \nSYMM block bytes: nsymbt_str %d \n\n",arraysize,*nsymbt);

    if(arraysize*4+1024+*nsymbt != fsize){      printf ("error! file size %d does not equal to arraysize (%d) + 1024 + %d \n!!", fsize, arraysize,*nsymbt);      }

  }
  in.close();
  return fsize;
}

int readmap(const char *filename0, float *array,int nsymbt,int c, int r, int s)
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
    in.read(( char *)array,(arraysize*4));
  }

  in.close();
  delete header;
  delete symmop;
  return fsize;
}

int read_slice(const char *filename0, float *array,int nsymbt,int c, int r,int slice) //slice start from 0
{
  int fsize=0;

  std::ifstream in (filename0,std::ios::in|std::ios::binary|std::ios::ate);
	if(!in)	{	  printf ( "file does not exist\n");    return 0;  }

  if (in.is_open()){
    fsize = in.tellg();
    int start=1024+nsymbt+(slice)*c*r*4;
    int end=1024+nsymbt+(slice+1)*c*r*4;
    if(end>fsize){      printf ("error! file size %d <  \n!!", fsize);      }
        printf ("reading slice at bytes %12d - %12d ... ",start,end);

    in.seekg (start, std::ios::beg);
    in.read(( char *)array,(c*r*4));
  }

  in.close();
  return c*r*4;
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
        *mean+=*ary;
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


int bin(float * array, float * array_bin,int x, int y)
{
    for (int i=0;i<y/2;i++){
        for (int j=0;j<x/2;j++){
            int idx=i*x/2+j;
            *(array_bin+idx)=(*(array+i*2*x+j*2)+*(array+(i*2+1)*x+j*2+1)+*(array+(i*2+1)*x+j*2)+*(array+(i*2)*x+j*2+1))*4;
            
        
        }
    }
} 
int rescale(float * ary, float scale,float min, int size){
    printf(" rescaling with min = %f scale = %f ...\n",min,scale);
    for(int i=0;i<size;i++)
    {
        
        *(ary+i)=(*(ary+i)-min )* scale;
        
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
