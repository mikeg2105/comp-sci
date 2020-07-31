/*The patient data inflamation study

See
http://swcarpentry.github.io/python-novice-inflammation/

and
http://swcarpentry.github.io/matlab-novice-inflammation/

Illustrates
  1. performance comparison for different ODE solvers,
  2. difference between euler, predictor corrector methods and fourth oder runge kutta methods for solution of ODE's
  3. use of OpenGL and Glut.

The vanderpol oscillator is an oscillator with a non linear damping term.
One of its applications has been the modelling of action potentials for neorons. 
http://en.wikipedia.org/wiki/Van_der_Pol_oscillator

To compile this program use

  g++ -o vanderpol vdpodesolvers.c

Or use 

  make all

OpenGLUT is the OpenGL utility toolkit for developing grapjics and visualisation based applications.

To enable the GLUT uncomment the correct lines from the file make_include_windows
http://openglut.sourceforge.net/

*/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#define NUM_POINTS    50000
#include <stdlib.h>



#ifdef USE_GLUT
	#include<GL/gl.h>
	#include<GL/glu.h>
	#include<GL/freeglut.h>


	#include <GL/freeglut.h>


	#define drawOneLine(x1,y1,x2,y2)  glBegin(GL_LINES);  \
   		glVertex2f ((x1),(y1)); glVertex2f ((x2),(y2)); glEnd();
#endif

float fh=0.001;
float fR=2.0;
float ps0=0.5;
float v0=0.3;
float fom=1.0;
float fgam=0.1;
int st=0;  //st is solver type rk4=0,  adams=2, euler=3, eulerpredictorcorrector=4
float vcurrent=v0;
float pscurrent=ps0;
float vnext=v0;
float psnext=ps0;
float t=0;
float **ps;
float *v;
float *ta;
//float ps[ NUM_POINTS ][ 2 ];
//float v[ NUM_POINTS ];
//float ta[ NUM_POINTS ];

/*
 * A handle on the main window.
 */
int main_window;

char *readLine(FILE *file, char *line) {

    if (file == NULL) {
        printf("Error: file pointer is null.");
        exit(1);
    }

    int maximumLineLength = 128;
    char *lineBuffer = (char *)malloc(sizeof(char) * maximumLineLength);

    if (lineBuffer == NULL) {
        printf("Error allocating memory for line buffer.");
        exit(1);
    }

    char ch = getc(file);
    int count = 0;

    while ((ch != '\n') && (ch != EOF)) {
        if (count == maximumLineLength) {
            maximumLineLength += 128;
            lineBuffer = (char *)realloc(lineBuffer, maximumLineLength);
            if (lineBuffer == NULL) {
                printf("Error reallocating space for line buffer.");
                exit(1);
            }
        }
        lineBuffer[count] = ch;
        count++;

        ch = getc(file);
    }

    lineBuffer[count] = '\0';
    /*char line[count + 1];*/
    line=(char *)calloc(count+1,sizeof(char));
    strncpy(line, lineBuffer, (count + 1));
    free(lineBuffer);
    char *constLine = line;
    return constLine;
}

/*----------------------*/ 
const char* getfield(char* line, int num)
{
    const char* tok;
    for (tok = strtok(line, ";");
            tok && *tok;
            tok = strtok(NULL, ";\n"))
    {
        if (!--num)
            return tok;
    }
    return NULL;
}






float timer()
{

   /*REAL secs;
   clock_t Time;
   Time = clock();

   secs = (REAL) Time / (REAL) CLOCKS_PER_SEC;
   return secs;*/
   double retval;
	static long zsec=0;
	static long zusec=0;
	double esec;
	
	struct timeval tp;
	struct timezone tzp;
	
	gettimeofday(&tp, &tzp);
	
	if(zsec==0) zsec=tp.tv_sec;
	if(zusec==0) zusec=tp.tv_usec;
	
	retval=(tp.tv_sec - zsec)+(tp.tv_usec-zusec)*0.000001;
	return retval;

}

float second(){ return timer();}

float f(float t,float ps,float v,float R, float gam0, float om) {
	
  float fqr=(ps/R);
  float f1= v*gam0*(1-(fqr*fqr))-om*om*ps;
	
  return f1;
}

//function [ps,v,ta,time]=rk4(ps,v,ta,ha,R,gam0,om)
float rk4(float **ps, float *v, float *ta, float h, float R, float gam0, float om){
  float tinitial=timer();
  //float fqr;
  //float rk1,rk2,rk3,rk4;
  
  for( int j=0;j<NUM_POINTS-2;j++)
  {
    t=j*h;
    ps[j][0]=t;
    //ta(i,j)=t;
    float rk1=f(t,ps[j][1],v[j],R,gam0,om);
    float rk2=f(t+(h/2),ps[j][1],v[j]+(rk1*h/2),R,gam0,om);
    float rk3=f(t+(h/2),ps[j][1],v[j]+(rk1*h/2),R,gam0,om);
    float rk4=f(t+(h/2),ps[j][1],v[j]+(rk1*h/2),R,gam0,om);
    //fqr=(ps[j][1]/R);
    //float f1= v*gam0*(1-(fqr*fqr))-om*om*ps[j][1];
  
    //rk1=v[j]*gam0*(1-(fqr*fqr))-om*om*ps[j][1];
    //rk2=(v[j]+(rk1*h/2))*gam0*(1-(fqr*fqr))-om*om*ps[j][1];
    //rk3=(v[j]+(rk1*h/2))*gam0*(1-(fqr*fqr))-om*om*ps[j][1];
    //rk4=(v[j]+(rk1*h/2))*gam0*(1-(fqr*fqr))-om*om*ps[j][1];
    
    v[j+1]=v[j]+(h/6)*(rk1+2*rk2+2*rk3+rk4);

    rk1=v[j];
    rk2=v[j]+(rk1*h/2);
    rk3=v[j]+(rk2*h/2);
    rk4=v[j]+(h*rk3);
    ta[j+1]=t;
    ps[j+1][1]=ps[j][1]+(h/6)*(rk1+2*rk2+2*rk3+rk4);
    //plot2d(ps(t+1),v(t+1),-1,rect=[-4,-6,4,6]);
    //xset('wshow');
    //xset('wwpc');
    
  }//end;
  float  tfinal=timer();
     //plot(ta(i,:),ps(i,:))
     return(tfinal-tinitial);
   
}//endfunction

float adams(float **ps, float *v, float *ta, float h, float R, float gam0, float om){
  float tinitial=timer();
  
  
  for( int j=0;j<4;j++)
  {
    t=j*h;
    ps[j][0]=t;
    //ta(i,j)=t;
    float rk1=f(t,ps[j][1],v[j],R,gam0,om);
    float rk2=f(t+(h/2),ps[j][1],v[j]+(rk1*h/2),R,gam0,om);
    float rk3=f(t+(h/2),ps[j][1],v[j]+(rk1*h/2),R,gam0,om);
    float rk4=f(t+(h/2),ps[j][1],v[j]+(rk1*h/2),R,gam0,om);
    v[j+1]=v[j]+(h/6)*(rk1+2*rk2+2*rk3+rk4);

    rk1=v[j];
    rk2=v[j]+(rk1*h/2);
    rk3=v[j]+(rk2*h/2);
    rk4=v[j]+(h*rk3);
    ta[j+1]=t;
    ps[j+1][1]=ps[j][1]+(h/6)*(rk1+2*rk2+2*rk3+rk4);
    //plot2d(ps(t+1),v(t+1),-1,rect=[-4,-6,4,6]);
    //xset('wshow');
    //xset('wwpc');
    
  }//end;
  
  float vam0=f(t-(4*h),ps[0][1],v[0],R,gam0,om);
  float vam1=f(t-(3*h),ps[1][1],v[1],R,gam0,om);
  float vam2=f(t-(2*h),ps[2][1],v[2],R,gam0,om);
  float vam3, fvt,pse,fpse,vt;
  
  
  for( int j=4;j<NUM_POINTS-2;j++)
  {
    t=j*h;
	
	  //run adams moulton predictor-corrector steps
        ps[j][0]=t;

	    ta[j]=t;
	    vam3=f(t-h,ps[j-1][1],v[j-1],R,gam0,om);
	    vt=v[j]-(h/24)*(9*vam0-37*vam1+59*vam2-55*vam3);
	    fvt=f(t,ps[j][1],vt,R,gam0,om);
	    //corrector step
	    v[j+1]=v[j]+(h/24)*(vam1-5*vam2+19*vam3+9*fvt);
	
	    fvt=f(t,ps[j][1],v[j],R,gam0,om);
	
	    vam0=vam1;
	    vam1=vam2;
	    vam2=vam3;
	
	
	
	    pse=ps[j][1]-(h/24)*(9*v[j-4]-37*v[j-3]+59*v[j-2]-55*v[j-1]);
	    fpse=(pse-ps[j][1])/h; //numerical derivative using the estimated value
	    //fpse=0;
	    ps[j+1][1]=ps[j][1]+(h/24)*(v[j-3]-5*v[j-2]+19*v[j-1]+9*fpse);    
  }//end;
  float  tfinal=timer();
     //plot(ta(i,:),ps(i,:))
     return(tfinal-tinitial);
   
}//endfunction

float euler(float **ps, float *v, float *ta, float h, float R, float gam0, float om){
  float tinitial=timer();
  for( int j=0;j<NUM_POINTS-2;j++)
  {
    t=j*h;
    ta[j]=t;
    ps[j][0]=t;
    
    float fvt=f(t,ps[j][1],v[j],R,gam0,om);
    //euler step
    v[j+1]=v[j]+h*fvt;
    //Euler step
    ps[j+1][1]=ps[j][1]+h*v[j+1];
    ta[j+1]=t;
  }//end;
  float  tfinal=timer();
  //plot(ta(i,:),ps(i,:))
  return(tfinal-tinitial);
   
}//endfunction

float eulerpc(float **ps, float *v, float *ta, float h, float R, float gam0, float om){
  float tinitial=timer();
  for( int j=0;j<NUM_POINTS-2;j++)
  {
    t=j*h;
    ta[j]=t;
    ps[j][0]=t;
    float fvt=f(t,ps[j][1],v[j],R,gam0,om);
    //euler pc step
    float v1=v[j]+h*fvt;
    float ps1=ps[j][1]+v1*h;
    v[j+1]=v[j]+h*(fvt+f(t,ps1,v1,R,gam0,om))/2;
 
    //Euler pc step
    float p1=ps[j][1]+h*v[j];
    ps[j+1][1]=ps[j][1]+h*(p1+v[j+1])/2;
    ps[j+1][1]=ps[j][1]+h*v[j+1];

    ta[j+1]=t;   
  }//end;
  float  tfinal=timer();
     //plot(ta(i,:),ps(i,:))
     return(tfinal-tinitial);
   
}//endfunction


void initmodel()
{
	int i=0;
	float time;
	FILE* stream = fopen("input", "r");
	
	
	/*another line reader*/
	       FILE * fp;
       char * line = NULL;
       size_t len = 0;
       ssize_t read;

       fp = fopen("/etc/motd", "r");
       if (fp == NULL)
           exit(EXIT_FAILURE);

       while ((line = readLine(fp, line)) != '\0') {
           printf("Retrieved line of length %zu :\n", read);
           printf("%s", line);
       }

       fclose(fp);
       if (line)
           free((void *)line);

	/*read data from csv file*/
    //char line[1024];
    while (fgets(line, 1024, stream))
    {
        char* tmp = strdup(line);
        printf("Field 3 would be %s\n", getfield(tmp, 3));
        // NOTE strtok clobbers tmp
        free(tmp);
    }
	
	
	
	ps=(float **)calloc(NUM_POINTS,sizeof(float *));
	for(i=0; i<NUM_POINTS; i++)
	  ps[i]=(float *)calloc(2,sizeof(float));
    v=(float *)calloc(NUM_POINTS,sizeof(float));
    ta=(float *)calloc(NUM_POINTS,sizeof(float));
//float ps[ NUM_POINTS ][ 2 ];
//float v[ NUM_POINTS ];
//float ta[ NUM_POINTS ];
	for(i=0; i<NUM_POINTS; i++)
	{
		//ps[i][0]=20+i/2;//-(NUM_POINTS/2);
		//ps[i][1]=80+50*sin((double)(4*3.1412*i/NUM_POINTS));
		ps[i][0]=0;
		ps[i][1]=0;
		ta[i]=0;
	}
 		ps[0][1]=ps0; 
 		v[0]=v0;


    


		
}

void runmodel()
{
  float time;
  int i;
  FILE *mfptr; /* mfptr = file pointer for results file*/
  switch(st)
   {
   	case 0:
   		time=rk4( ps, v, ta, fh, fR, fgam, fom);
   	break;
   	
   	case 2:
   	   	time=adams(ps, v, ta, fh, fR, fgam, fom);   	
   	break;
   	
   	case 3:
  		time=euler(ps, v, ta, fh, fR, fgam, fom);    	
   	break;
   	
   	case 4:
   		time=eulerpc(ps, v, ta, fh, fR, fgam, fom);
   	
   	break;
   	
   	default:
   		time=rk4(ps, v, ta, fh, fR, fgam, fom);   	
   	break;
   }
   printf("run time is %f seconds\n\n",20*time);


  //write results to output file
        if((mfptr=fopen("models.dat", "w"))==NULL)
		printf("The file could not be opened!\n");
	else
	{

			fprintf(mfptr, "Time \t Amplitude\n");



                 for(i=0; i<NUM_POINTS;i++)
   				fprintf(mfptr,"%g %g\n",ps[i][0],ps[i][1]);

		fclose(mfptr);
	}


 //  for(i=0; i<NUM_POINTS;i++)
 //  				printf("%f %f\n",ps[i][0],ps[i][1]);	
}

void freemodel()
{
    int i=0;
	for(i=0; i<NUM_POINTS; i++)
	  free(ps[i]);
	free(ps);
    free(v);
    free(ta);
	
}






#ifdef USE_GLUT

void init(void) 
{
   glClearColor (0.0, 0.0, 0.0, 0.0);
   glShadeModel (GL_FLAT);
}

void display(void)
{
   int i;
   glClear (GL_COLOR_BUFFER_BIT);
   glColor3d( 1.0, 0.0, 0.0 );  /* Red */

 
   for(i=0; i<NUM_POINTS-1; i++)
   {
	  drawOneLine (4*ps[ i ][0], 5*(ps[ i ][1]+12), 4*ps[ i+1 ][0], 5*(ps[ i+1 ][1]+12));
    }    

   glFlush ();   //glEnd( ); 
 
}



void reshape (int w, int h)
{
   glViewport (0, 0, (GLsizei) w, (GLsizei) h);
   glMatrixMode (GL_PROJECTION);
   glLoadIdentity ();
   
    gluOrtho2D (0.0, (GLdouble) w, 0.0, (GLdouble) h);
}

/* ARGSUSED1 */
void keyboard(unsigned char key, int x, int y)
{
   switch (key) {
      case 27:     /* The escape key*/
         exit(0);
         break;
   }
}

#endif



int main(int argc, char** argv)
{


#ifdef USE_GLUT
   glutInit(&argc, argv);
   glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
   glutInitWindowSize (400, 150); 
   glutInitWindowPosition (100, 100);
   
   main_window = glutCreateWindow( "Vanderpol Oscillator" );
   init ();
#endif


   initmodel();
   runmodel();
 
#ifdef USE_GLUT
   glutDisplayFunc(display); 
   glutReshapeFunc(reshape);
   glutKeyboardFunc(keyboard);
   glutMainLoop();
#endif
   freemodel();

   return 0;
}
