

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



void freadl(FILE *stream, char **string)
{    
    unsigned long counter = 0;
    char *line = NULL;
    int next = fgetc(stream);

    do {
        next = fgetc(stream);
        if (next == EOF) {
            free(line);
            break;
        }
        ++counter;
        line = (char*)realloc(line, counter + 1);
        if (line == NULL) {
            puts("line == NULL");
            exit(EXIT_FAILURE);
        }
        line[counter - 1] = (char)next;
    } while (next != '\n');
    line[counter - 1] = '\0';
    *string = (char *)malloc(strlen(line) + 1);
    if (*string == NULL) {
        puts("*string == NULL");
    } else {
        strcpy(*string, line);
    }
    free(line);

}

int writebinaryarray(char *filename, char *headerinfo, real *data)
{
	int status=0;
	
	return status;

}

int readbinaryarray(char *filename, char *headerinfo, real *data)
{
	int status=0;
	int i;
	
	return status;
}

int writevacconfig(char *name,int n,params p, meta md, real *w, real *wd, state st)
{
  int status=0;
  int i1,j1,k1,ifield;
  int ni,nj,nk;
  
  char tcfg[300];  
  char configfile[300];
  
  char buffer[800];
  double dbuffer[12];
  int ibuffer[5];
  char ext[3];

  ni=p.n[0];
  nj=p.n[1];
    #ifdef USE_SAC_3D
  nk=p.n[2];
    #endif

 
   sprintf(configfile,"%s",name);
   
   #ifdef USE_MPI
   
   char *pch1,*pch2;
   pch1 = strtok (configfile,".");
   sprintf(tcfg,"%s",pch1);
   pch2 = strtok (NULL,".");

   printf("here1 %s %d %d \n",tcfg,ni,nj);
   //sprintf(ext,"%s",pch2);

      //set the input filename corresponding to proc id
 
     #ifdef USE_SAC3D
	      if(p.ipe>99)
		sprintf(configfile,"%s%d_np0%d0%d0%d_%d.out",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.pnpe[2],p.ipe);
	      else if(p.ipe>9)
		sprintf(configfile,"%s%d_np0%d0%d0%d_0%d.out",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.pnpe[2],p.ipe);
	      else
		sprintf(configfile,"%s%d_np0%d0%d0%d_00%d.out",tcfg,p.pnpe[0],p.pnpe[1],p.pnpe[2],p.ipe);  	     
     #else
	      if(p.ipe>99)
		sprintf(configfile,"%s%d_np0%d0%d_%d.out",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.ipe);
	      else if(p.ipe>9)
		sprintf(configfile,"%s%d_np0%d0%d_0%d.out",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.ipe);
	      else
		sprintf(configfile,"%s%d_np0%d0%d_00%d.out",tcfg,st.it,p.pnpe[0],p.pnpe[1],p.ipe);  	     	     
     #endif

 
   #else


         //save file containing current data
      sprintf(configfile,"%s_%d.out",name,st.it);
  #endif
   

     // sprintf(configfile,"%s",name);
      printf("write vac check dims %d %d %d %lf\n",ni,nj,st.it,st.t);
//printf("here3 %s \n",configfile); 
      FILE *fdt=fopen(configfile,"w");
      //FILE *fdt=fopen("out/test.out","w");


      fwrite(md.name,sizeof(char)*79,1,fdt);

      // ibuffer[0]=24;
       //fwrite(ibuffer,sizeof(int),1,fdt);
      //fwrite(md.name,sizeof(char)*79,1,fdt);
      //*line2:
      //*   it          - timestep (integer)
      //*   t           - time     (real)
      //*   ndim        - dimensionality, negative sign for gen. coord (integer)
      //*   neqpar      - number of equation parameters (integer)
      //*   nw          - number of flow variables (integer)
      //sprintf(buffer,"%ld %lg %ld %ld %ld\n",st.it,st.t,3,4,8);
      //it,time,ndim,neqpar,nw
      printf("st.it=%d\n",st.it);
      ibuffer[0]=st.it;
      dbuffer[0]=st.t;
      //ibuffer[0]=1;
      //dbuffer[0]=10.98; st.it
      fwrite(ibuffer,sizeof(int),1,fdt);
      fwrite(dbuffer,sizeof(double),1,fdt);

    #ifdef USE_SAC
      ibuffer[0]=2;
      ibuffer[1]=6;
      ibuffer[2]=10;
      //ibuffer[0]=st.it;
      //ibuffer[1]=(float)st.t;


      //ibuffer[2]=2;  //report from distribution routine to reset the data
      //ibuffer[3]=6;
      //ibuffer[4]=10;
   #endif
    #ifdef USE_SAC_3D
      ibuffer[0]=3;
      ibuffer[1]=7;
      ibuffer[2]=13;
      //ibuffer[2]=3;
    #endif
      fwrite(ibuffer,sizeof(int)*3,1,fdt);
      //fwrite(ibuffer,sizeof(int)*5,1,fdt);
      //line3:
      //*   nx()        - the grid dimensions      (ndim integers)
      //sprintf(buffer,"%ld %ld\n",ni,nj);
      ibuffer[0]=ni;
      ibuffer[1]=nj;

      //ibuffer[0]=128;
      //ibuffer[1]=256;
    #ifdef USE_SAC
      fwrite(ibuffer,sizeof(int)*2,1,fdt);
    #endif
    #ifdef USE_SAC_3D
      ibuffer[2]=nk;
      fwrite(ibuffer,sizeof(int)*3,1,fdt);
    #endif
      //*line4:
      //*   eqpar()     - equation parameters from filenameini (neqpar reals)
      //sprintf(buffer,"%lg %lg %lg %lg %lg %lg\n",p.gamma,p.eta,p.g[0],p.g[1],0,0);
      dbuffer[0]=p.gamma;
      dbuffer[1]=p.eta;
      dbuffer[2]=p.g[0];
      dbuffer[3]=p.g[1];

    #ifdef USE_SAC
      dbuffer[4]=0;
      dbuffer[5]=0;
     fwrite(dbuffer,sizeof(double)*6,1,fdt);
    #endif
    #ifdef USE_SAC_3D
      dbuffer[4]=p.g[2];
      dbuffer[5]=0;
      dbuffer[6]=0;
     fwrite(dbuffer,sizeof(double)*7,1,fdt);
    #endif
 

      //*line5:
      //*   varnames    - names of the coordinates, variables, equation parameters
      //*                 eg. 'x y rho mx my e bx by  gamma eta' (character*79)

    #ifdef USE_SAC_3D
      sprintf(buffer,"x y z rho mx my mz e bx by bz gamma eta g1 g2 g3");
    #else
      sprintf(buffer,"x y rho mx my e bx by  gamma eta g1 g2");
    #endif

      fwrite(buffer,sizeof(char)*79,1,fdt);

    #ifdef USE_SAC_3D
      for(ifield=0;ifield<16;ifield++)   
   #else
       for(ifield=0;ifield<12;ifield++)   
   #endif 

    #ifdef USE_SAC_3D
   for( k1=0;k1<nk;k1++)
    #endif
for( j1=0;j1<nj;j1++)
  
	{
//energyb,rhob,b1b,b2b         
       for( i1=0;i1<ni;i1++)     
      {
          //if(ifield==2) printf("density %lG ",w[(j1*ni+i1)]);
               if(ifield==0)
               dbuffer[0]=wd[(j1*ni+i1)+(ni*nj*(pos1))];
               else if(ifield==1)
               dbuffer[0]=wd[(j1*ni+i1)+(ni*nj*(pos2))];
    #ifdef USE_SAC_3D
               if(ifield==0)
               dbuffer[0]=wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos1))];
               else if(ifield==1)
               dbuffer[0]=wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos2))];
               else if(ifield==2)
               dbuffer[0]=wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos3))];
               
    #else
    
    #endif
               else
    #ifdef USE_SAC_3D
                dbuffer[0]=w[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(ifield-3))];
    #else
                dbuffer[0]=w[(j1*ni+i1)+(ni*nj*(ifield-2))];              
    #endif

                fwrite(dbuffer,sizeof(double),1,fdt);
                //printf("%g ",dbuffer[0]);		

        }     
      }
      //printf("density nowt \n \n ");
      buffer[0]='\n';
      fwrite(buffer,sizeof(char),1,fdt);
     // printf("\n");
      fclose(fdt);

  return status;
}




int readbinvacconfig(char *name,params p, meta md, real *w,real *wd, state st)
{

  int status=0;
  int i1,j1,k1,ifield;
  int ni,nj,nk;
  
  char tcfg[300];  
  char configfile[300];
  
  char buffer[800];
  double dbuffer[12];
  int ibuffer[5];
  char ext[3];
  long lsize;
  size_t result;

  //char *bigbuf;
  ni=p.n[0];
  nj=p.n[1];
    #ifdef USE_SAC_3D
  nk=p.n[2];
    #endif






      //save file containing current data
      //sprintf(configfile,"out/v%s.out",name);
      sprintf(configfile,"%s",name);
      printf("check dims %d %d \n",ni,nj);
      FILE *fdt=fopen(configfile,"r");


      fread(buffer,sizeof(char),80,fdt);
      for(i1=0;i1<81;i1++)
         putchar(buffer[i1]);
      printf("starting %s\n ",buffer);
       strcpy(md.name,buffer);

      //*line2:
      //*   it          - timestep (integer)
      //*   t           - time     (real)
      //*   ndim        - dimensionality, negative sign for gen. coord (integer)
      //*   neqpar      - number of equation parameters (integer)
      //*   nw          - number of flow variables (integer)
      //sprintf(buffer,"%ld %lg %ld %ld %ld\n",st.it,st.t,3,4,8);
      //it,time,ndim,neqpar,nw

      fread(ibuffer,sizeof(int),1,fdt);
      fread(dbuffer,sizeof(double),1,fdt);
      st.it=ibuffer[0]=st.it;
      st.t=dbuffer[0]=st.t;
      printf("st.it=%d st.t=%f\n",st.it,st.t);
      fread(ibuffer,sizeof(int),3,fdt);

      //fread(ibuffer,sizeof(int)*3,1,fdt);

      //line3:
      //*   nx()        - the grid dimensions      (ndim integers)
      //sprintf(buffer,"%ld %ld\n",ni,nj);
      //fread(ibuffer,sizeof(int)*2,1,fdt);
    #ifdef USE_SAC
      fread(ibuffer,sizeof(int)*2,1,fdt);
    #endif
    #ifdef USE_SAC_3D
      
      fread(ibuffer,sizeof(int)*3,1,fdt);
      //nk=ibuffer[2];
    #endif
      //ni=ibuffer[0];
      //nj=ibuffer[1];







      //*line4:
      //*   eqpar()     - equation parameters from filenameini (neqpar reals)
      //sprintf(buffer,"%lg %lg %lg %lg %lg %lg\n",p.gamma,p.eta,p.g[0],p.g[1],0,0);
      //fread(dbuffer,sizeof(double)*6,1,fdt);



    #ifdef USE_SAC
      fread(dbuffer,sizeof(double)*6,1,fdt);
      //dbuffer[4]=0;
      //dbuffer[5]=0;
     
    #endif
    #ifdef USE_SAC_3D
      fread(dbuffer,sizeof(double)*7,1,fdt);
      p.g[2]=dbuffer[4];
      //dbuffer[5]=0;
      //dbuffer[6]=0;     
    #endif

      p.gamma=dbuffer[0];
      p.eta=dbuffer[1];
      p.g[0]=dbuffer[2];
      p.g[1]=dbuffer[3];
      printf("%f %f %f %f\n",dbuffer[0],dbuffer[1],dbuffer[2],dbuffer[3]);



      //*line5:
      //*   varnames    - names of the coordinates, variables, equation parameters
      //*                 eg. 'x y rho mx my e bx by  gamma eta' (character*79)
      sprintf(buffer,"x y rho mx my mz e bx by bz gamma eta g1 g2 g3\n");
      fread(buffer,sizeof(char)*79,1,fdt);



      // printf("read buffer%s\n",buffer);



    #ifdef USE_SAC_3D
      for(ifield=0;ifield<16;ifield++)   
   #else
       for(ifield=0;ifield<12;ifield++)   
   #endif 

    #ifdef USE_SAC_3D
   for( k1=0;k1<nk;k1++)
    #endif


         for( i1=0;i1<ni;i1++)   
	{
//energyb,rhob,b1b,b2b         
     
for( j1=0;j1<nj;j1++)  
      {



               fread(dbuffer,sizeof(double),1,fdt);
               if(ifield==0)
               	wd[(j1*ni+i1)+(ni*nj*(pos1))]=dbuffer[0];
               else if(ifield==1)
               	wd[(j1*ni+i1)+(ni*nj*(pos2))]=dbuffer[0];
    #ifdef USE_SAC_3D
               if(ifield==0)
                 wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos1))]=dbuffer[0];
               else if(ifield==1)
                 wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos2))]=dbuffer[0];
               else if(ifield==2)
                 wd[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(pos3))]=dbuffer[0];
               
    #else
    
    #endif
               else
    #ifdef USE_SAC_3D
                w[(k1*ni*nj+j1*ni+i1)+(ni*nj*nk*(ifield-3))]=dbuffer[0];
    #else
                w[(j1*ni+i1)+(ni*nj*(ifield-2))]=dbuffer[0]; 
              //if(ifield==2)
              //  printf("%g ",w[(j1*ni+i1)+(ni*nj*(ifield-2))]);             
    #endif

               






               // fread(dbuffer,12*sizeof(double),1,fdt);		
                //i1*p.dx[0]=dbuffer[0];
                //j1*p.dx[1]=dbuffer[1];
               //if(j1==2 || j1==3)
                //  printf("%d %d %d %d %lg\n", ni,nj,i1,j1,dbuffer[2]);
 


        } 
        // if(ifield==2)
       // printf("\n");    
      }

      // printf("read bin vac read fields\n");
      fclose(fdt);
      //free(bigbuf);
  return status;


}

