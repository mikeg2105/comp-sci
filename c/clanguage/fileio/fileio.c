#include<stdio.h>


//C++ sequential file access example
// Create a structure to hold the data we want
struct Item
{
char szName[20];
char szType[20];
char szPrice[5];
};



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
Item one={"Apple","Fruit","0.20"};
Item two={"Bannana","Fruit","0.15"};
Item three={"Peas","Veg","0.10"};
Item four={"Digestive","Biscuit","0.30"};
Item five={"Cornflakes","Breakfast","0.05"};

// Empty variables that we will use to read data back in
Item ToRead1={};
Item ToRead4={};
Item ToRead3={};

// Open a file
FILE* file;
file = fopen("Items.txt","w+");

// Write the data to the file
fwrite(&one,sizeof(Item),1,file);
fwrite(&two,sizeof(Item),1,file);
fwrite(&three,sizeof(Item),1,file);
fwrite(&four,sizeof(Item),1,file);
fwrite(&five,sizeof(Item),1,file);

// Move to the start of the file and read the first item
fseek(file,0,SEEK_SET);
fread(&ToRead1,sizeof(Item),1,file);

// Read the 4th item
fseek(file,(sizeof(Item)*3),SEEK_SET);
fread(&ToRead4,sizeof(Item),1,file);

// Read the 3rd item
fseek(file,(sizeof(Item)*2),SEEK_SET);
fread(&ToRead3,sizeof(Item),1,file);

// Output to the screen
printf("Item 1 is: %s\n",one.szName);
printf("Item 4 is: %s\n",four.szName);
printf("Item 3 is: %s\n",three.szName);

// Close the file
fclose(file);


return 0;
}
