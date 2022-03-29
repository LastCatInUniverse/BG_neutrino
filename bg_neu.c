#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define NAMESIZE 512
#define STRINGSIZE 512

// global variables
char inname[NAMESIZE],
     outname[NAMESIZE];
FILE *infile=NULL,
     *outfile=NULL;
char key_xyz=0,
     key_tbl=0,
     key_png=0,
     key_gal=0;
int n_ra, n_dec, n_rmag;

double flux_sun=60, // billions, actually flux_sun=6e+10,
       rmag_sun=-26.7,
       ra_box=0,
       dec_box=0,
       box=0;
int    sigma=0;

struct csvdata{
  char name[NAMESIZE];
  int n;
  double *value;
};
struct csvdata *data;
int n_par=0;


// short information
void print_help(){
  printf("~~~~~~~~~~ %s ~~~~~~~~~~\n","programname");
  printf("Usage:\n");
  printf("$> bg_neu -xyz <input file> -o <output file>\n");
  printf("  -xyz   convert decart coordinates to spherical\n");
  printf("  -tbl   also print table of objects in tex and pdf formats\n");
  printf("  -gal   use galactic coordinates\n");
  printf("  -png   use PNG format for output image\n");
  printf("  nra=   number of column with RA\n");
  printf("  ndec=  number of column with Dec\n");
  printf("  nrmag= number of column with Rmag\n");
  printf("  ra=    center of selected area in degree, default whole map\n");
  printf("  dec=   center of selected area in degree, default whole map\n");
  printf("  box=   side of selected area in degree, default whole map\n");
  printf("  sigma= use smoothing for output image\n");
  exit(1);
}

// read external arguments
void read_args(char **arg){
  while(*++arg!=NULL){
    // read keys
    if(!strcmp(*arg,"--help"))
      print_help();
    if(!strcmp(*arg,"-xyz")){
      key_xyz=1;
      continue;
    }
    if(!strcmp(*arg,"-tbl")){
      key_tbl=1;
      continue;
    }
    if(!strcmp(*arg,"-gal")){
      key_gal=1;
      continue;
    }
    if(!strcmp(*arg,"-png")){
      key_png=1;
      continue;
    }
    // read name for output file
    if(!strcmp(*arg,"-o")){
      arg++;
      strncpy(outname,*arg,NAMESIZE);
      continue;
    }
    // read parameters
    if(!strncmp(*arg,"nra=",4)){
      n_ra=atoi(*arg+4)-1;
      continue;
    }
    if(!strncmp(*arg,"ndec=",5)){
      n_dec=atoi(*arg+5)-1;
      continue;
    }
    if(!strncmp(*arg,"nrmag=",6)){
      n_rmag=atoi(*arg+6)-1;
      continue;
    }
    if(!strncmp(*arg,"ra=",3)){
      ra_box=atof(*arg+3);
      continue;
    }
    if(!strncmp(*arg,"dec=",4)){
      dec_box=atof(*arg+4);
      continue;
    }
    if(!strncmp(*arg,"box=",4)){
      box=atof(*arg+4);
      continue;
    }
    if(!strncmp(*arg,"sigma=",6)){
      sigma=atoi(*arg+6);
      continue;
    }
    // read name for inpiut file
    if(*inname==0){
      strncpy(inname,*arg,NAMESIZE);
      continue;
    }
    fprintf(stderr,"Unknown parameter '%s'\n",*arg);
    exit(1);
  }
  if(*inname==0){
    fprintf(stderr,"Required name of input file\n");
    exit(1);
  }
  if(*outname==0){
    strcpy(outname,inname);
    sprintf(strrchr(outname,'.'),"_result.csv");
  }
}

void open_files(){
  if((infile=fopen(inname,"r"))==NULL){
    fprintf(stderr,"Can not open file '%s'\n",inname);
    exit(1);
  }
  if((outfile=fopen(outname,"w"))==NULL){
    fprintf(stderr,"Can not write file '%s'\n",outname);
    exit(1);
  }
}

double read_coord(char *str){
  double res;
  char sep=0;
  int sign=1;
  while(*str==' ') str++;  // shift to first non-space letter
  if(strchr(str,' ')!=NULL) sep=' ';
  if(strchr(str,':')!=NULL) sep=':';
//   double min,sec;
//   sprintf(str,"%f%f%f\n",&res,&min,&sec);
//   if(res<0) sign=-1;
//   else sign=1;
//   return res+sign*min/60+sign*sec/3600;
  res=atof(str);                // deg
  if(sep==0) return res; // if format with one value
  if(res<0) sign=-1;
  str=strchr(str,sep)+1;
  res=res+sign*atof(str)/60;    // min
  str=strchr(str,sep)+1;
  res=res+sign*atof(str)/3600;  // sec  
 return res;
}

void read_data(){
  char buf[STRINGSIZE],*p;
  int n=0;
  // calc number of parameters
  fgets(buf,STRINGSIZE,infile);
  p=buf;
  while(strchr(p,',')!=NULL){
    p=strchr(p,',')+1;
    n_par++;
  }
  // calc number of objects in data
  while(fgets(buf,STRINGSIZE,infile)){
    n++;
  }
  rewind(infile);
  // allocate memory
  data=calloc(sizeof(*data),(n_par+1));
  for(int i=0;i<n_par;i++)
    data[i].value=calloc(sizeof(*data[i].value),(n+1));
  // read header
  p=fgets(buf,STRINGSIZE,infile);
  int i=0;
  while(strchr(p,',')!=NULL){
    data[i].n=n;
    strcpy(data[i].name,p);
    *strchr(data[i].name,',')=0;
    //detect RA, Dec & Rmag
    if(strstr(data[i].name,"RA")!=NULL && n_ra==0) n_ra=i;
    if(strstr(data[i].name,"Dec")!=NULL && n_dec==0) n_dec=i;
    if(strstr(data[i].name,"Rmag")!=NULL && n_rmag==0) n_rmag=i;
    p=strchr(p,',')+1;
    i++;
  }
  // read main data
  for(int i=0;(p=fgets(buf,STRINGSIZE,infile))!=NULL;i++)
    for(int j=0;j<n_par;j++){
      if(j==n_ra || j==n_dec)
        data[j].value[i]=read_coord(p);
      else
        data[j].value[i]=atof(p);
      p=strchr(p,',')+1;
    }
}

void plot_sky(){
  FILE *pyfile=NULL;
  char pyname[]="plot.py",
       cmd[512],
       outepsname[512],
       figsize[512];
  // figure size and ratio
  if(box==0)
    strcpy(figsize,"14,7");
  else
    strcpy(figsize,"7,7");
  // name for eps file
  strcpy(outepsname,outname);
  if(key_png==0)
    sprintf(strrchr(outepsname,'.'),".eps");
  else
    sprintf(strrchr(outepsname,'.'),".png");
  // create python script file
  if((pyfile=fopen(pyname,"w"))==NULL){
    fprintf(stderr,"Can't create file '%s'. Exit.\n",pyname);
    exit(1);
  }
  // print script to file
  fprintf(pyfile,"import matplotlib.pyplot as plt\n");
//  fprintf(pyfile,"import matplotlib.cbook as cbook\n");
  fprintf(pyfile,"import numpy as np\n");
  fprintf(pyfile,"import matplotlib.cm as cm\n");
  fprintf(pyfile,"import pandas as pd\n");
  fprintf(pyfile,"import math\n");
  fprintf(pyfile,"from scipy.ndimage.filters import gaussian_filter\n");
  fprintf(pyfile,"\n");
  fprintf(pyfile,"def deg2rad(x):\n");
  fprintf(pyfile,"  return ((x+180)%360-180)/180*math.pi\n");
  fprintf(pyfile,"def myplot(x, y, s, bins=1000):\n");
  fprintf(pyfile,"  heatmap, xedges, yedges = np.histogram2d(x,y,bins=bins)\n");
  fprintf(pyfile,"  heatmap = gaussian_filter(heatmap,sigma=s)\n");
  fprintf(pyfile,"  extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]\n");
  fprintf(pyfile,"  return heatmap.T, extent\n");
  fprintf(pyfile,"data=pd.read_csv('%s')\n",outname);
  fprintf(pyfile,"fig = plt.figure(figsize=(%s))\n",figsize);
  if(box==0)
    fprintf(pyfile,"fig.add_subplot(projection=\"aitoff\")\n");
  fprintf(pyfile,"plt.title(\"%s\"  )\n",inname);
  fprintf(pyfile,"plt.grid(True)\n");
  if(box!=0){
    fprintf(pyfile,"fig,ax=plt.subplots(figsize=(%s))\n",figsize);
    fprintf(pyfile,"ax.set_aspect(1)\n");
    if(key_gal==0){
      fprintf(pyfile,"ax.set_xlabel('RA,degrees')\n");
      fprintf(pyfile,"ax.set_ylabel('Dec,degrees')\n");
    }
    else{
      fprintf(pyfile,"ax.set_xlabel('Galactic longitude,degrees')\n");
      fprintf(pyfile,"ax.set_ylabel('Galactic latitude, degrees')\n");
    }
  }
  // read csv file
  if(box==0){
    fprintf(pyfile,"ra=deg2rad(data.RA)\n");
    fprintf(pyfile,"dec=deg2rad(data.Dec)\n");
  }
  else{
    fprintf(pyfile,"ra=data.RA\n");
    fprintf(pyfile,"dec=data.Dec\n");
  }
  fprintf(pyfile,"\n");
  if(sigma==0){
    fprintf(pyfile,"image=plt.scatter(ra,dec,data.Flux/5,color=\"k\")\n");
  }
  else{
//    fprintf(pyfile,"fig=plt.figure(figsize=(%s))\n",figsize);
//    fprintf(pyfile,"ax=fig.subplots(1,1,1)\n");
//    fprintf(pyfile,"fig = plt.subplots(2, 2)\n");
//    fprintf(pyfile,"ax=fig.flatten()\n");
    fprintf(pyfile,"img, extent = myplot(ra,dec,%d)\n",sigma);
    fprintf(pyfile,"plt.imshow(img,aspect=abs((extent[1]-extent[0])/(extent[3]-extent[2])),extent=extent,origin='lower',cmap=cm.jet)\n");
//    fprintf(pyfile,"plt.set_title(\"Title\")\n");
  }
//  fprintf(pyfile,"msft.plot(\"RA\",\"Dec\",subplots=True)\n");
  fprintf(pyfile,"\n");
  fprintf(pyfile,"plt.savefig(\"%s\")\n",outepsname);
  sprintf(cmd,"python %s\n",pyname);
  fclose(pyfile);
  system(cmd);
}

void result_xyz(){
  double teta, fi, flux, x,y,z;
  int i,n_x,n_y,n_z;
  // find x,y,z in data
  for(i=0;i<n_par;i++){
    if(!strcmp(data[i].name,"\"x\"")) n_x=i;
    if(!strcmp(data[i].name,"\"y\"")) n_y=i;
    if(!strcmp(data[i].name,"\"z\"")) n_z=i;
  }
  // print header
  fprintf(outfile,"\"id\",\"RA\",\"Dec\",\"Flux\"\n");
  // find min and max
  double min_x,max_x,min_y,max_y,min_z,max_z;
  min_x=max_x=data[n_x].value[0];
  min_y=max_y=data[n_y].value[0];
  min_z=max_z=data[n_z].value[0];
  for(i=0;i<data[0].n;i++){
    if(data[n_x].value[i]<min_x) min_x=data[n_x].value[i];
    if(data[n_x].value[i]>max_x) max_x=data[n_x].value[i];
    if(data[n_y].value[i]<min_y) min_y=data[n_y].value[i];
    if(data[n_y].value[i]>max_y) max_y=data[n_y].value[i];
    if(data[n_z].value[i]<min_z) min_z=data[n_z].value[i];
    if(data[n_z].value[i]>max_z) max_z=data[n_z].value[i];
  }
  // results
  for(i=0;i<data[0].n;i++){
    x=data[n_x].value[i]-(max_x-min_x)/2;
    y=data[n_y].value[i]-(max_y-min_y)/2;
    z=data[n_z].value[i]-(max_z-min_z)/2;
    teta=atan(z/sqrt(x*x+y*y));;
    if(x>0 && y>=0) fi=atan(y/x);
    if(x<0) fi=M_PI+atan(y/x);
    if(x>0 && y<0) fi=2*M_PI+atan(y/x);
    if(x==0 && y>0) fi=M_PI/2;
    if(x==0 && y<0) fi=3*M_PI/2;
    teta=teta*180/M_PI;
    fi=fi*180/M_PI;
    flux=data[n_rmag].value[i];
    fprintf(outfile,"%d,%f,%f,%f\n",i+1,fi,teta,flux);
  }
}

double calc_flux(double rmag){
//    return rmag;
  double flux;
  flux = flux_sun*pow(10,rmag)/pow(10,rmag_sun);
  return flux;
}

void result_sky(){
  if(n_ra==0){
    fprintf(stderr,"Can't find RA in data\n");
    exit(1);}
  if(n_dec==0){
    fprintf(stderr,"Can't find Dec in data\n");
    exit(1);}
  if(n_rmag==0){
    fprintf(stderr,"Can't find Rmag in data\n");
    exit(1);}
  if(ra_box<0-box || ra_box>360+box ||dec_box<-90-box ||dec_box>90+box){
    fprintf(stderr,"Incorrect box parameters\n");
    exit(1);}
  // print header
  fprintf(outfile,"\"id\",\"RA\",\"Dec\",\"Flux\"\n");
  int i;
  double ra, dec,flux, glat,glon,glon2;
  double ra_gal,dec_gal,glon_gal;
//          ra_gal=read_coord("12:49:00")*360/24;
//          dec_gal=read_coord("+27.4");
         ra_gal=192.85833*M_PI/180;
         dec_gal=27.12833*M_PI/180;
         glon_gal=32.93192+90;
  // print result data
  for(i=0;i<data[0].n;i++){
    ra=data[n_ra].value[i]*360/24; // hours to degree
    dec=data[n_dec].value[i];
    if(key_gal){
       ra=ra*M_PI/180;
       dec=dec*M_PI/180;
       glat=sin(dec)*sin(dec_gal)+cos(dec)*cos(dec_gal)*cos(ra-ra_gal);
         glat=asin(glat);
       glon=cos(dec)*sin(ra-ra_gal)/cos(glat);
         glon=asin(glon);
       glon2=(cos(dec_gal)*sin(dec)-sin(dec_gal)*cos(dec)*cos(ra-ra_gal))/cos(glat);
       if(glon2<=0)
         glon=M_PI-glon;
       glon=glon_gal-glon*180/M_PI;
       if(glon<=0)
         glon=glon+360;
       glat=glat*180/M_PI;
       ra=glon;
       dec=glat;
    }
//    flux=data[n_rmag].value[i];
    flux=calc_flux(data[n_rmag].value[i]);
    if(box==0)
      fprintf(outfile,"%d,%f,%f,%f\n",i+1,ra,dec,flux);
    else
      if(ra>ra_box-box/2 && ra<ra_box+box/2 && dec>dec_box-box/2 && dec<dec_box+box/2)
        fprintf(outfile,"%d,%f,%f,%f\n",i+1,ra,dec,flux);
  }
  fclose(outfile);
}

int main(int argc, char **argv){
  read_args(argv);
  open_files();
  read_data();
  if(key_xyz) result_xyz();
    else result_sky();
  plot_sky();
  // print data information
  fprintf(stdout,"# data=%s result=%s\n# columns ra=%d dec=%d r_mag=%d\n",inname,outname,n_ra,n_dec,n_rmag);
  if(box!=0)
    fprintf(stdout,"# box side=%g center ra=%g dec=%g sigma=%d\n",box,ra_box,dec_box,sigma);
}
