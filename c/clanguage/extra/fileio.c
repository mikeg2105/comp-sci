#include<stdio.h>
/*Differerent examples of writing data to a file 
  als illustrates how to use structures

  Compile using
   gcc -o fileio filio.c

*/



	//C++ sequential file access example
	// Create a structure to hold the data we want
	struct Item
	{
	char szName[20];
	char szType[20];
	char szPrice[5];
	} ;




int main()
{
	FILE* examplefile;
	examplefile=fopen("Example.txt","w");
	fprintf(examplefile,"Example of some text to write to file");
	fclose(examplefile);

	char szStringFromFile[100];
	int iCounter = 0;
	examplefile=fopen("Example.txt","r");
	while(!feof(examplefile) && iCounter < 99)
	{
		szStringFromFile[iCounter++]=(char)fgetc(examplefile);
	}

	szStringFromFile[iCounter] = '\0';
	fclose(examplefile);
	puts(szStringFromFile);


	// Create some variables & data
	struct Item one={"Apple","Fruit","0.20"};
	struct Item two={"Bannana","Fruit","0.15"};
	struct Item three={"Peas","Veg","0.10"};
	struct Item four={"Digestive","Biscuit","0.30"};
	struct Item five={"Cornflakes","Breakfast","0.05"};

	// Empty variables that we will use to read data back in
	struct Item ToRead1={};
	struct Item ToRead4={};
	struct Item ToRead3={};

	// Open a file
	FILE* file;
	file = fopen("Items.txt","w+");

	// Write the data to the file
	fwrite(&one,sizeof(struct Item),1,file);
	fwrite(&two,sizeof(struct Item),1,file);
	fwrite(&three,sizeof(struct Item),1,file);
	fwrite(&four,sizeof(struct Item),1,file);
	fwrite(&five,sizeof(struct Item),1,file);

	// Move to the start of the file and read the first item
	fseek(file,0,SEEK_SET);
	fread(&ToRead1,sizeof(struct Item),1,file);

	// Read the 4th item
	fseek(file,(sizeof(struct Item)*3),SEEK_SET);
	fread(&ToRead4,sizeof(struct Item),1,file);

	// Read the 3rd item
	fseek(file,(sizeof(struct Item)*2),SEEK_SET);
	fread(&ToRead3,sizeof(struct Item),1,file);

	// Output to the screen
	printf("Item 1 is: %s\n",one.szName);
	printf("Item 4 is: %s\n",four.szName);
	printf("Item 3 is: %s\n",three.szName);

	// Close the file
	fclose(file);


	return 0;
}
