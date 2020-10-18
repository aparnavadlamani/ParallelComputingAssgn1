#include<iostream>
#include<mpi.h>

using namespace std;

int main(int argc, char** argv)
{
    int rank, np;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int n = atoi(argv[1]);
    int s = atoi(argv[2]);
    s = s-1;
    int d = atoi(argv[3]);
    d = d-1;
    const char* filename = argv[4];

    int dims[2];
    dims[0] = 0;
    dims[1] = 0;
    int periods[2];
    int coordinates[2];
    int left, right, up, down;
    int *a, *b, *c;
    periods[0] = 1;
    periods[1] = 1;
    MPI_Dims_create(np, 2, dims);

    int len = 1;
    //cout<<dims[0]<<" "<<dims[1]<<endl;
    if(dims[0]!=dims[1])
    {
        if(rank==0)
        {
            cout<<"Number of processors must be a square"<<endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    MPI_Comm Grid_Comm;
    MPI_Status status;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &Grid_Comm);

    int new_n = n/dims[0];
    int row_index = s%new_n;
    int col_index = d%new_n;

    a = new int[new_n*new_n];
    b = new int[new_n*new_n];
    c = new int[new_n*new_n];

    FILE* ftr;
    int num;
    if((ftr = fopen(filename,"r"))==NULL)
    {
        cout<<"Error opening file"<<endl;
        return 0;
    }
    int row = (rank/dims[0])*new_n;
    int col = (rank%dims[0])*new_n;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            fscanf(ftr, "%d", &num);
            if(i>=row&&i<row+new_n&&j>=col&&j<col+new_n)
            {
                a[(i-row)*new_n+(j-col)] = num;
                b[(i-row)*new_n+(j-col)] = num;
                c[(i-row)*new_n+(j-col)] = 0;
            }
        }
    }
    fclose(ftr);
    // for(int i=0;i<new_n;i++)
    // {
    //      for(int j=0;j<new_n;j++)
    //      {
    //          cout<<a[i*new_n+j]<<" ";
    //      }
    //      cout<<endl;
    // }
    int duration = 0;
    int start = MPI_Wtime();
    int pr = (s/new_n)*dims[0] + (d/new_n); 
    int path_found = 0;
    if(rank==pr)
    {
        if(a[(s%new_n)*new_n+d%new_n]==1)
        {
            path_found = 1;
            len=1;
        }
    }
    MPI_Bcast(&path_found, 1, MPI_INT, pr, Grid_Comm);
    MPI_Cart_coords(Grid_Comm, rank, 2, coordinates);

    while(path_found==0&&len<=n)
    {
        len+=1;
        MPI_Cart_coords(Grid_Comm, rank, 2, coordinates);
        MPI_Cart_shift(Grid_Comm, 1, coordinates[0], &left, &right);
        MPI_Cart_shift(Grid_Comm, 0, coordinates[1], &up, &down);
        MPI_Sendrecv_replace(a, new_n*new_n, MPI_INT, left, 0, right, 0, Grid_Comm, &status);
        MPI_Sendrecv_replace(b, new_n*new_n, MPI_INT, up, 0, down, 0, Grid_Comm, &status);
        for(int i=0;i<new_n;i++)
        {
            for(int j=0;j<new_n;j++)
            {
                for(int k=0;k<new_n;k++)
                {
                    c[i*new_n+k] += a[i*new_n+j]*b[j*new_n+k];
                }
            }
        }
        for(int x=0;x<dims[0];x++)
        {    
            MPI_Cart_shift(Grid_Comm, 1, 1, &left, &right);
            MPI_Cart_shift(Grid_Comm, 0, 1, &up, &down);
            MPI_Sendrecv_replace(a, new_n*new_n, MPI_INT, left, 0, right, 0, Grid_Comm, &status);
            MPI_Sendrecv_replace(b, new_n*new_n, MPI_INT, up, 0, down, 0, Grid_Comm, &status);

            for(int i=0;i<new_n;i++)
            {
                for(int j=0;j<new_n;j++)
                {
                    for(int k=0;k<new_n;k++)
                    {
                        c[i*new_n+k] += a[i*new_n+j]*b[j*new_n+k];
                    }
                }
            }
        }
        MPI_Barrier(Grid_Comm);
        if(rank==pr)
        {
            if(c[(s%new_n)*new_n+d%new_n]!=0)
            {
                path_found = 1;
            }
            // for(int z1=0;z1<new_n;z1++)
            // {
            //     for(int z2=0;z2<new_n;z2++)
            //     {
            //         cout<<c[z1*new_n+z2]<<" ";
            //     }
            //     cout<<endl;
            // }
        }
        MPI_Bcast(&path_found, 1, MPI_INT, pr, Grid_Comm);
        for(int i=0;i<new_n;i++)
        {
            for(int j=0;j<new_n;j++)
            {
                a[i*new_n+j] = c[i*new_n+j];
                c[i*new_n+j] = 0;
            }
        }        
    }
    int end = MPI_Wtime();
    duration = end - start;
    if(rank==0)
    {
        cout<<"Time: "<<duration<<endl;
        cout<<"Path Length: "<<len<<endl;
    }
    free(a);
    free(b);
    free(c);
    MPI_Finalize();
    return 0;
}