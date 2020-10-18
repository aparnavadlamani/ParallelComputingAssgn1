#include<iostream>
#include<fstream>
#include<limits>

using namespace std;

int main(int argc, char** argv)
{ 
    int v = atoi(argv[1]);
    int** graph = new int* [v];
    for(int i=0;i<v;i++)
    {
        graph[i] = new int[v];
    }
    for(int i=0;i<v;i++)
    {
        for(int j=i+1;j<v;j++)
        {
            if(rand()%2==1)
            {
                graph[i][j] = 1;
                graph[j][i] = 1;
            }
        }
    }
    fstream myfile;
    myfile.open("matrix.txt",fstream::out);
    //myfile<<v<<endl;
    for(int i=0;i<v;i++)
    {
        for(int j=0;j<v;j++)
        {
            myfile<<graph[i][j]<<" ";

        }
        myfile<<endl;
    }
    myfile.close();
    return 0;
}