#include <stdio.h>

/*Matrix transpose example*/


#define MAX_ROW      50
#define MAX_COL      50


int main(void) {

int i, j, m, n, temp;
int mat[MAX_ROW][MAX_COL];


    // variable dim  is set to smaller value of defined
    // maximal number of rows and columns


    int dim = (MAX_ROW < MAX_COL)? MAX_ROW : MAX_COL;



    // storing matrix size

    do {

        printf("Input number of rows < %d   :", dim);
        scanf("%d", &m);
        printf("Input number of columns < %d:", dim);
        scanf("%d", &n);

    }   while (m < 1 || m > dim || n < 1 || n > dim);



    // storing matrix elements

    printf("\nInput of matrix elements :\n");
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {

         printf("Input element [%d][%d] : ", i, j);
           scanf("%d", &mat[i][j]);

        }
    }  



    // printing matrix before transposing

    printf("\n\nMatrix before transposing:\n");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {

         printf("%3d", mat[i][j]);

      }
        printf("\n");
    }




   // transposing

    for ( i=0; i<m; ++i ) {
   // looping must start from i+1 !!
        for ( j=i+1; j<n; ++j ) {

               temp = mat[i][j];
                mat[i][j] = mat[j][i];
                mat[j][i] = temp;
        }
    }



    // print after transposing
    // number of rows becomes number of columns ...

    printf("\nMatrix after transposing:\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {

              printf("%3d", mat[i][j]);
        }
        printf("\n");
    }
    
    //system("PAUSE");
} // main


