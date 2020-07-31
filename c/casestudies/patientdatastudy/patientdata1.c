/*program21.c*/
/*Example of Reading data from a sequential file */
/*fileread.c*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define NROW 60
#define NCOL 40



char *readLine(FILE *file, char *line) {

    if (file == NULL) {
        printf("Error: file pointer is null.");
        exit(1);
    }

    int maximumLineLength = 1280;
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
    line=(char *)calloc(count+1,sizeof(char));
    strncpy(line, lineBuffer, (count + 1));
	
    free(lineBuffer);
    char *constLine = line;
    return constLine;
}



main()
{

	
 	FILE * fp;
    char * line = NULL;
	int rcount,ccount,col,row;
	char *tokenptr;
    char *tokenseparators = ",";
   
    char patientfilename[] ="data/inflammation-01.csv";
	printf("Simulation log reader.\n");
	
	/*Each file has 60 rows and 40 columns*/
	float patientdata[NROW][NCOL];
	
	
	/*initialize the patient data*/
	for(row=0; row<NROW; row++)
		for(col=0; col<NCOL; col++)
			 patientdata[row][col]=0.0;
	
	       fp = fopen(patientfilename, "r");
       if (fp == NULL)
           exit(EXIT_FAILURE);

	   rcount=0;
	   /*first count number of numbers in the line*/
       while (rcount<NROW) {
		   
		   line = readLine(fp, line);
		   
           /*printf("Retrieved line of length %s :\n", line);
           printf("%s\n", line);*/
       

	   
	         /*use tokenizer to readline*/			tokenptr = strtok(line, tokenseparators);
		    ccount=0;
			while(tokenptr != NULL)
			{
						
						if(ccount<NCOL  && rcount<NROW)
						    /*printf("%f\n",atof(tokenptr));*//*This line for checking during debug*/
							patientdata[rcount][ccount]=atof(tokenptr);					     
					    else
							printf("Overflow %d %d\n",rcount,ccount); /*Warning in case of file overflow error*/
						tokenptr = strtok(NULL, tokenseparators);
						ccount++;
			}  
			rcount++;
	   
	   }
	   
	   	for(row=0; row<NROW; row++)
		{
			for(col=0; col<NCOL; col++)
				printf("%f ",patientdata[row][col]);
			printf("\n\n");
		}

	   
       fclose(fp);
       if (line)
           free((void *)line);
	
	
       
    

	return 0;
}



