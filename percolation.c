#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<stdbool.h>
#include<omp.h>
#include "/usr/include/mpi/mpi.h"
#include<string.h>
#include<limits.h>

// define constants
#define L 46    
#define N (L*L)
#define originProcess 0
#define amountOfProcesses 2
#define amountOfThreads 2
#define split (N/amountOfProcesses)

/* structure of coordinates */
typedef struct coord{
    int row;
    int col;
} COORD;

/* structure of cluster, with the id, and coordinates within it */
typedef struct cluster{
    int id;
    COORD pos[N];
    int length;
} CLUSTER;

/* function declarations */
void printGrid(int[L][L]);
void printClusters(CLUSTER[], int);
void printCoords(COORD[], int);
void printArr(int[], int);
    
void percolate(double, int);
void dfs(int, int, int, int, int[L][L], int[L][L], CLUSTER[], int);
CLUSTER* getCluster(int, CLUSTER[], int);
void mergeClusters(CLUSTER*, CLUSTER*, CLUSTER[], int*, int[L][L]);

void doesPercolate(CLUSTER[], int, bool*, bool*);
int maxCluster(CLUSTER[], int);

bool contains(int[N], int, int);
bool containsAll(int[], int);
int mod (int, int);
double randomProb();


/* program begins here */
int main() {
    double prob;
    MPI_Init(NULL, NULL);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // get the probability of filling a site with '1', and call percolate() for all processes.
    if(world_rank == originProcess){
        printf("Enter a number for 0 to 1 from the seeding probability:\n");
        scanf("%lf", &prob);
        MPI_Send(&prob, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        percolate(prob, world_rank);
    }
    else{
        MPI_Recv(&prob, 1, MPI_DOUBLE, originProcess, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        percolate(prob, world_rank);
    }
    MPI_Finalize();
    return 0;
}


/* this performs all functionality of this application, it finds the amount of clusters in the lattice
*  the size of greatest cluster, and whether row and/or column percolation occurs */
void percolate(double prob, int world_rank) {
    bool row_perco = false;
    bool col_perco = false; 

    /* Depth First Search function used to explore the lattice and build the clusters */
    void dfs(int row, int col, int minEdge, int maxEdge, int lattice[L][L], int seen[L][L], CLUSTER cList[N], int idx){
        if(row < minEdge || row >= maxEdge || seen[row][col] > 0 || lattice[row][col] == 0){
            return;
        }
        // use local pointer to address of cluster in array to improve readability
        CLUSTER* currCluster = &cList[idx];

        // assign cluster ID to site
        seen[row][col] = currCluster->id;

        // adding site to cluster
        int nxt = currCluster->length;
        currCluster->pos[nxt] = (COORD){
            .row = row,
            .col = col
        };
        currCluster->length++;

        // continue to explore the lattice recursively and build the cluster
        dfs(mod(row+1, L), col, minEdge, maxEdge, lattice, seen, cList, idx);
        dfs(mod(row-1, L), col, minEdge, maxEdge, lattice, seen, cList, idx);
        dfs(row, mod(col+1, L), minEdge, maxEdge, lattice, seen, cList, idx);
        dfs(row, mod(col-1, L), minEdge, maxEdge, lattice, seen, cList, idx);
    }

    // set the number of threads per process
    omp_set_num_threads(amountOfThreads);

    // for the origin process...
    if(world_rank == originProcess){
        // set up the lattice
        int (*lattice)[L] = malloc(sizeof(int[L][L]));
        double seeding_prob = prob;
        srand(time(NULL));
  
        //populate lattice randomly
        for(int i=0; i < L; ++i){
            for(int j=0; j < L; ++j){
                if (randomProb() < seeding_prob){
                     lattice[i][j] = 1;
                }
                else{
                     lattice[i][j] = 0;
                }       
            }
        }

        // send info to other processes
        for(int i=1; i<amountOfProcesses; i++){
            // send the lattice
            printf("\n%dx%d lattice.\n ",L, L);
            MPI_Send(lattice, N, MPI_INT, i, 0, MPI_COMM_WORLD);
            // calculate the split of the lattice and send it
            int begin = 0;
            int end = N;
            int vals[2] = {begin, end};
            MPI_Send(vals, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        //free(lattice);
    }
    // for all other processes...
    else{
        int (*lattice)[L] = malloc(sizeof(int[L][L]));
        // receive the lattice from the original process
        MPI_Recv(lattice, N, MPI_INT, originProcess, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // receive split info
        int vals[2];
        MPI_Recv(vals, 2, MPI_INT, originProcess, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int begin = vals[0];
        int end = vals[1];

        // initialise values
        CLUSTER (*clusterList) = calloc(N, sizeof(int[N]));
        int (*seen)[L] = calloc(N, sizeof(int[L][L]));
        int amountOfClusters = 1;
        int listTop = 0;
        int nthreads = 2;

        /* split that process into threads, and perfrom this block in parallel
        *  and share the variables in shared() across all threads */
        #pragma omp parallel shared(seen, amountOfClusters, clusterList, listTop)
        {
            // define the start and end of the rows each thread covers in the lattice
            int thr = omp_get_thread_num();
            int beginRow = thr*(end/nthreads)/L;
            int endRow = (thr+1)*(end/nthreads)/L;

            /* statically split the for loop for the two threads
            *  each thread does N/nthreads iterations */
            #pragma omp for schedule(static, N/nthreads)
            for(int i=begin; i<end; i++) {
                int row = i/L;
                int col = i%L;
                //if a site is occupied, perform DFS
                if(lattice[row][col] == 1 && seen[row][col] == 0){
                    // add cluster to clusterList
                    #pragma omp critical
                    {
                        CLUSTER fst = { .id = amountOfClusters++, .length = 0};
                        clusterList[listTop] = fst;
                    }
                    dfs(row, col, beginRow, endRow, lattice, seen, clusterList, listTop++);
                }
            }

            int idList[N];
            int top_idList = 0;

            /* Perform this code block such that only one thread performs it at a tim
            *  Stitch the segments of the lattice that each thread worked on together
            *  Go along the edges of the segments, and compare with the edges of the adjacent segments
            *  Merge clusters together accordingly */
            #pragma omp critical
            {
                // for each column in a row
                for(int i=0; i<L; i++){
                    // if the current site and the site of the previous row are not empty, and are not part of the same cluster, then merge their clusters together
                    if(lattice[beginRow][i] != 0 && lattice[mod(beginRow-1, L)][i] != 0 && (seen[beginRow][i] != seen[mod(beginRow-1, L)][i])
                        && (!contains(idList, top_idList, seen[beginRow][i]) || !contains(idList, top_idList, seen[mod(beginRow-1, L)][i]))) {
                        //get the clusters that have these sites that meet on the edge of the segments
                        CLUSTER* c1 = getCluster(seen[mod(beginRow-1, L)][i], clusterList, listTop);
                        CLUSTER* c2 = getCluster(seen[beginRow][i], clusterList, listTop);
                        // the id of c1 will be the id of the merged cluster, so add it to the 'seen' array of cluster ids
                        idList[top_idList++] = c1->id;
                        // merge the clusters together in the list
                        mergeClusters(c1, c2, clusterList, &listTop, seen);
                    }
                }
            }

        }
        printGrid(lattice);
        //printGrid(seen);
        free(lattice);
        free(seen);
        
        doesPercolate(clusterList, listTop, &row_perco, &col_perco);
        int maxSize = maxCluster(clusterList, listTop);
        free(clusterList);

        if(maxSize < 0){
            fprintf(stderr, "ERROR: No max cluster found");
        }

        printf("\nNumber of clusters: %d\n", listTop);
        printf("Number of sites in biggest cluster: %d\n", maxSize);
        printf("Row percolation: %s\n", row_perco ? "true" : "false");
        printf("Col percolation: %s\n", col_perco ? "true" : "false");
    }
}


/* assign boolean value to row_perco and col_perco depending if the lattice is row percolating or column percolating */
void doesPercolate(CLUSTER list[], int top, bool* row_perco, bool* col_perco){
    for(int i=0; i<top; i++){
        int rowp[L] = {0};
        int colp[L] = {0};
        for(int j=0; j<list[i].length; j++){
            if(rowp[list[i].pos[j].row] == 0){
                rowp[list[i].pos[j].row] = 1;
            }
            if(colp[list[i].pos[j].col] == 0){
                colp[list[i].pos[j].col] = 1;
            }
        }
        (*row_perco) = containsAll(rowp,1);
        (*col_perco) = containsAll(colp,1);
    }
}

/* check if every element of an integer array is equal to the integer num */
bool containsAll(int arr[L], int num){
    for(int i=0; i<L; i++){
        if(arr[i] != num){
            return false;
        }
    }
    return true;
}

/* return the size of the cluster with the most sites */
int maxCluster(CLUSTER list[N], int top){
    int currMax = -1;
    for(int i=0; i<top; i++){
        if(list[i].length > currMax){
            currMax = list[i].length;
        }
    }
    return currMax;
}

/* return whether an integer array contains an integer or not 
*  list[N]: integer array
*  top: the pointer for the array
*  val: the integer to check for
*/
bool contains(int list[N], int top, int val){
    for(int i=0; i<top; i++){
        if(list[i] == val){
            return true;
        }
    }
    return false;
}

/* get the cluster with the id */
CLUSTER* getCluster(int id, CLUSTER list[], int size){
    for(int i=0; i<size; i++){
        if(list[i].id == id){
            return &list[i];
        }
    }
    fprintf(stderr, "ERROR: No cluster exists.\n");
    exit(EXIT_FAILURE);
}

/* merge two clusters together */
void mergeClusters(CLUSTER* c1, CLUSTER* c2, CLUSTER list[], int* size, int seen[L][L]){
    //assign all values of cluster c2 to cluster c1, and upate seen grid
    for(int i=0; i<c2->length; i++){
        seen[c2->pos[i].row][c2->pos[i].col] = c1->id;
        c1->pos[c1->length++] = c2->pos[i];
    }

    //mark the index of the cluster to remove
    int to_remove = -1;
    for(int i=0; i<(*size); i++){
        if(list[i].id == c2->id){
            to_remove = i;
        }
    }

    //remove the cluster by decreasing the size of the array
    //and assigning its address to the last element.
    if(to_remove >= 0 && to_remove < (*size)){
        list[to_remove] = list[(*size)-1];
        (*size)--;
    }
    else{
        fprintf(stderr, "ERROR: index out of range.\n");
        exit(EXIT_FAILURE);
    }
}

/* find the modulus and handle negatives */
int mod (int a, int b)
{
   if(b < 0)
     return -mod(-a, -b);   
   int ret = a % b;
   if(ret < 0)
     ret+=b;
   return ret;
}

/* return a random number between 0 and 1 */
double randomProb(){
    return (double)rand()/(double)RAND_MAX;
}

/* prints the grid in ASCII to the terminal*/
void printGrid(int lattice[L][L]){
    for(int i = 0; i < L; ++i){
        printf("\n");
        for(int j = 0; j < L; ++j){
            printf("%d ", lattice[i][j]);
        }
    }
    printf("\n");
}

/* prints each cluster and the coordinates in each cluster */
void printClusters(CLUSTER list[N], int listTop){
    for(int i=0; i<listTop; i++){
        printf("cluster_id: %d\n", list[i].id);
        printCoords(list[i].pos, list[i].length);
        printf("\n");
    }
}

/* prints each coordinate in a array of coordinates */
void printCoords(COORD pos[N], int length){
    for(int i=0; i<length; i++){
        printf("(%d,%d) | ", pos[i].row, pos[i].col);
    }
    printf("\n");
}

/* print the values of an array */
void printArr(int arr[], int top){
    printf("[");
    for(int i=0; i<top; i++){
        printf("%d", arr[i]);
        if(i+1 != top){
            printf(", ");
        }
    }
    printf("]: %d \n", top);
}


/* use these commands to compile and execute the application */
// mpicc -o percolation percolation.c -fopenmp
// mpirun -np 2 ./percolation