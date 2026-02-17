#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define SERVER_NODE 0
#define MAX_ITER 50

// not very optimized :X
int fib(int i){
    if(i<1){
        return i;
    }else{
        return fib(i-2) + fib(i-1);
    }
}

// rapid calculation on fibionacci
int main(int argc, char ** argv) {
    int rang, world_size;

    int tab_all[4][MAX_ITER];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rang);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);

    for(int i=0; i<MAX_ITER; ++i){
        tab_all[&rang%4][i] = fib(i);
    }

    MPI_Finalize();

    printf("Final result: \n");
    for(int i=0; i<MAX_ITER; ++i){
        printf("%d\t%d\t%d\t%d\n",tab_all[0][i],tab_all[1][i],tab_all[2][i],tab_all[3][i]);
    }

    return 0;
}