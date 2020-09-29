#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<stdbool.h>
#include<omp.h>
#include<mpi.h>
#include<string.h>
#include<limits.h>

#define L 512
#define N (L*L)
#define M ((2*L)*L)
#define EMPTY (-N-1)
#define numthreads 4
#define originProcess 0
#define finalProcess 65
#define amountProcesses 64
#define split (N/amountProcesses)


//structure for the cluster information
struct cluster{
    int id;
    int cols[L];
    int rows[L];
    int coltop;
    int rowtop;
    int size;
};

//structure for the edge information
struct edge{
    int pos;
    int id;
    int size;
};

void percolate(char, double, int);
double rando();
void makeGrid(double, int[L][L]);
void makeBond(double, char[2*L][L]);
void printGrid(int[L][L]);
void printBond(char[2*L][L]);
bool isvalueinarray(int, int[split], int);
bool isvaluein2darraystruct(int, struct cluster *, int);
int returnvaluein3darraystruct(int, struct cluster *, int);

void main()
{
    char option;
    double prob;
   
    MPI_Init(NULL, NULL);
    //Find out rank, size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_rank <= amountProcesses || world_rank == finalProcess){
        if (world_rank == originProcess){
            printf("Enter a number for 0 to 1 from the seeding probability:\n");
            scanf("%lf", &prob);
            printf("Type 's' for site percolation:\n");
            scanf("%c", &option);
            scanf("%c", &option);
             for(int i=1;i<=amountProcesses;i++){
                MPI_Send(&prob, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                MPI_Send(&option, 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
            }
            MPI_Send(&prob, 1, MPI_DOUBLE, finalProcess, 0, MPI_COMM_WORLD);
            MPI_Send(&option, 1, MPI_CHAR, finalProcess, 0, MPI_COMM_WORLD);
            percolate(option, prob, world_rank);
        }
        else{
            MPI_Recv(&prob, 1, MPI_DOUBLE, originProcess, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&option, 1, MPI_CHAR, originProcess, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);            
            percolate(option, prob, world_rank);
        }
    }
    MPI_Finalize();
   }


void percolate(char option, double prob, int world_rank) {
    
    omp_set_num_threads(numthreads);

    if (option == 's'){
        
        
        int nthreads;
      
        //for original process, populate the lattice and send it to other processes
        if (world_rank == originProcess){
            int (*occ)[L] = malloc(sizeof(int[L][L]));

            double seeding_prob = prob;
            srand(time(NULL));
            
            //populate lattice
            for(int i=0;i<L;i++){
                for(int j=0;j<L;j++){
                    if (rando() < seeding_prob){
                         occ[i][j] = 1;
                    }
                    else{
                         occ[i][j] = 0;
                    }       
                }
            }
            
            //send to every other process
            for(int i=1;i<=amountProcesses;i++){
                MPI_Send(occ, N, MPI_INT, i, 0, MPI_COMM_WORLD);
                int begin = split*(i-1);
                int end = split*i;
                int vals[2] = {begin, end};
                MPI_Send(vals, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            free(occ);
        }

        else if (world_rank !=originProcess && world_rank <=amountProcesses && world_rank!=finalProcess) {
            //allocate arrays on the heap
            int (*nn)[4] = malloc(sizeof(int[N][4]));
            int (*occ)[L] = malloc(sizeof(int[L][L]));

            //receive the lattice from original process 
            MPI_Recv(occ, N, MPI_INT, originProcess, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            /* NEAREST NEIGHBOURS */
            for (int i=0; i<N; i++) {
                nn[i][0] = (i+1)%N;  //EAST
                nn[i][1] = ((i+N)-1)%N; //WEST
                nn[i][2] = (i+L)%N;  //SOUTH
                nn[i][3] = ((i+N)-L)%N; //NORTH
                if (i%L==0) nn[i][1] = i+L-1;
                if ((i+1)%L==0) nn[i][0] = i-L+1;   
            }


            /* Receive where to search in the lattice from the original process */
            int vals[2];
            MPI_Recv(vals, 2, MPI_INT, originProcess, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int begin = vals[0];
            int end = vals[1];
            int lim[numthreads], amount_upper_edge[numthreads], amount_lower_edge[numthreads], global_max=-1, perc_type;
            struct cluster (*global_cluster_struct)[split] = malloc(sizeof(struct cluster [numthreads][split])); 
            struct edge (*global_lower_edge)[L] = malloc(sizeof(struct edge [numthreads][L]));
            struct edge (*global_upper_edge)[L] = malloc(sizeof(struct edge [numthreads][L]));
        
            //printf("Enter '0' or 1' or '2' for row or coloumn or both percolation\n");
            //scanf("%d", &perc_type);

            //if (perc_type == 1) row_flag = true;
            //else if (perc_type == 2) col_flag = true;
            #pragma omp parallel
                {
                    //initialise data structures
                    int (*keeptrack) = malloc(sizeof(int[split])); 

                    int ktop = 0, ctop=0, max=0, lower_etop=0, upper_etop=0, id = omp_get_thread_num();
                  
                    struct cluster cluster_struct[L];
                    struct edge lower_edge_struct[L];
                    struct edge upper_edge_struct[L];
                   
                    nthreads = omp_get_num_threads();
                    
                    #pragma omp for schedule(static, split/omp_get_num_threads())
                    for(int i=begin; i<end; i++) {
                        //printf("world_rank: %d | i: %d===================================\n", world_rank, i);
                        if (isvalueinarray(i, keeptrack, ktop))
                            continue;

                        bool row_percolate = false, col_percolate = false;
                        int current_node=0, neighbour_node=0, vtop=0, rowtop=0, coltop=0, qtop=0, tempi=i/L, tempj=i%L;

                        int row_arr[L];
                        int col_arr[L];
                        int queue[split];
                        int visited[split];

                        queue[qtop++] = -1; queue[qtop] = i;
                        visited[vtop++] = i;
                        keeptrack[ktop++] = i;
                        row_arr[rowtop++] = tempi;
                        col_arr[coltop++] = tempj;

                        //if a site is occupied
                        if (occ[tempi][tempj] == 1){
                           
                            //while the stack is not empty
                            while(qtop!=0){
                                //printf("qtop: %d\n",qtop);
                                current_node = queue[qtop]; 
                                //set current node, and the neighbour nodes
                                int upper_row_boundary = ((id+1)*(split/nthreads) + split*(world_rank-1));
                                int lower_row_boundary = ((id)*(split/nthreads) + split*(world_rank-1));

                                int dir;
                                for (dir = 0; dir < 4; dir++){
                                    if (occ[nn[current_node][dir]/L][nn[current_node][dir]%L] == 1 && !isvalueinarray(nn[current_node][dir], visited, vtop)
                                        && (nn[current_node][dir]/L < upper_row_boundary/L) && (nn[current_node][dir]/L >= lower_row_boundary/L)){
                                            //if next row hasn't been visited before then add to percolation tracking
                                            if (!isvalueinarray(nn[current_node][dir]/L, row_arr, rowtop)) 
                                                row_arr[rowtop++] = nn[current_node][dir]/L;
                                            if (!isvalueinarray(nn[current_node][dir]%L, col_arr, coltop))
                                                col_arr[coltop++] = nn[current_node][dir]%L;
                                            //add the next node to array of visited nodes
                                            visited[vtop++] = nn[current_node][dir];
                                            //keep track of all visited nodes 
                                            keeptrack[ktop++] = nn[current_node][dir];            
                                            //add next node to top of the queue
                                            queue[++qtop] = nn[current_node][dir];
                                            break;
                                    }
                                    //if there is no next_node to go to, decrement the queue
                                    else if (dir == 3)
                                        queue[qtop--] = 0;
                                }                    
                            }

                            if (rowtop == L/nthreads) row_percolate = true;
                            if (coltop == L/nthreads) col_percolate = true;
                            if (vtop>max) max = vtop;
                            //iterates through the visited sites
                            for(int v=0; v<vtop; v++){
                                if (!isvaluein2darraystruct(i, cluster_struct, ctop)){
                                	//ctop is the current cluster
    								cluster_struct[ctop].id = i;
                                    cluster_struct[ctop].size = vtop;
    								//copies the columns and rows of the cluster to the struct array
    								for(int k=0;k<coltop;k++) cluster_struct[ctop].cols[k] = col_arr[k];
    								for(int k=0;k<rowtop;k++) cluster_struct[ctop].rows[k] = row_arr[k];
    								//keeps track of how many rows and columns the cluster covers
    								cluster_struct[ctop].rowtop = rowtop;
    								cluster_struct[ctop++].coltop = coltop;
                                }
                                //boundaries of grid
                                int lrb = ((id+1)*(split/nthreads) + split*(world_rank-1))/L;
                                int urb = ((id)*(split/nthreads) + split*(world_rank-1))/L;
                
                                if (visited[v]/L == (lrb - 1)){
                                    lower_edge_struct[lower_etop].pos = visited[v];
                                    lower_edge_struct[lower_etop].size = vtop;
                                    lower_edge_struct[lower_etop++].id = i;
                                }
                                else if (visited[v]/L == (urb)){
                                    upper_edge_struct[upper_etop].pos = visited[v];
                                    upper_edge_struct[upper_etop].size = vtop;
                                    upper_edge_struct[upper_etop++].id = i;
                                }
                            }
                        }
                    } //end of for loop
            
                    
                    #pragma omp critical
                    {
                        //converts data structures in threads to global data structures
                        for(int i=0; i<lower_etop; i++) 
                            global_lower_edge[id][i] = lower_edge_struct[i];

                        for(int i=0; i<upper_etop; i++) 
                            global_upper_edge[id][i] = upper_edge_struct[i];

                        for(int i=0; i<ctop; i++) 
                            global_cluster_struct[id][i] = cluster_struct[i];

                        if (max > global_max)
                            global_max = max;

                        //keeps track of amount of edges and clusters for global access
                        lim[id] = ctop;
                        amount_lower_edge[id] = lower_etop;
                        amount_upper_edge[id] = upper_etop;
                        
                    }
                    free(keeptrack); 
                }

                free(occ);
                free(nn);
               
    			struct cluster (*new_cluster_struct) = malloc(sizeof(struct cluster [split]));
                
                // iterate through each thread and assign cluster information to global structure
    			int cc=0;
                for(int i=0;i<numthreads;i++){
                    for (int j=0; j<lim[i]; j++){
    					new_cluster_struct[cc++] = global_cluster_struct[i][j];
    				}
                }
                free(global_cluster_struct);
                

                //send the amount of edges for each thread to the final process
                MPI_Send(amount_lower_edge, numthreads, MPI_INT, finalProcess, 0, MPI_COMM_WORLD);
                MPI_Send(amount_upper_edge, numthreads, MPI_INT, finalProcess, 0, MPI_COMM_WORLD);

                //for each thread send the information about the edges, e.g. position, id etc
                for(int n=0; n<numthreads; n++){
                    struct edge* ltemp = global_lower_edge[n];
                    struct edge* utemp = global_upper_edge[n];

                    MPI_Send(ltemp, (3*amount_lower_edge[n]), MPI_INT, finalProcess, 0, MPI_COMM_WORLD);    
                    MPI_Send(utemp, (3*amount_upper_edge[n]), MPI_INT, finalProcess, 0, MPI_COMM_WORLD);
                }

                free(global_lower_edge);
                free(global_upper_edge);


                //send the cluster for this process to the final process e.g. id, size, rows/columns
                MPI_Send(new_cluster_struct, (cc*(4+2*L)), MPI_INT, finalProcess, 0, MPI_COMM_WORLD);
                free(new_cluster_struct);
            
        }
        /* FOR THE FINAL PROCESS COMPILE ALL THE DATA FROM THE OTHER PROCESSES AND JOIN CLUSTERS */
        else if (world_rank == finalProcess){

            int visited_cluster[split];
            int border_cluster[split];
            int vc=0, bc=0;
            struct cluster (*final_cluster) = malloc(sizeof(struct cluster [N]));

            struct edge (*global_upper)[L] = malloc(sizeof(struct edge [numthreads*amountProcesses][L]));
            struct edge (*global_lower)[L] = malloc(sizeof(struct edge [numthreads*amountProcesses][L]));

            int lamount[(numthreads*amountProcesses)];
            int uamount[(numthreads*amountProcesses)];
            int lc = 0, uc = 0, fc = 0, la=0, ua=0;


            //JOIN ALL THE INFORMATION FROM OTHER PROCESSES
            for(int i=1;i<=amountProcesses;i++){
                int lower_edge[numthreads];
                int upper_edge[numthreads];
                struct cluster (*recv) = malloc(sizeof(struct cluster [split]));
            
                MPI_Status status;

                int num;
            
                MPI_Recv(lower_edge, numthreads, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(upper_edge, numthreads, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for(int n=0; n<numthreads; n++){
                    struct edge ltest[lower_edge[n]];
                    struct edge utest[upper_edge[n]];
                    MPI_Recv(ltest, (3*lower_edge[n]), MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(utest, (3*upper_edge[n]), MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    
                    for(int j=0; j<lower_edge[n];j++){
                        global_lower[lc][j] = ltest[j];
                    }
                    for(int j=0; j<upper_edge[n];j++){
                        global_upper[uc][j] = utest[j];
                    }
                    lc++;
                    uc++;
                    lamount[la++] = lower_edge[n];
                    uamount[ua++] = upper_edge[n];

                }
                
                MPI_Recv(recv, split*(4+2*L), MPI_INT, i, 0, MPI_COMM_WORLD, &status);
               
                MPI_Get_count(&status, MPI_INT, &num);

                for(int j=0;j<(num/(4+2*L));j++){
        
                    final_cluster[fc++] = recv[j];
                }
                free(recv);              
            }


            //update amount of threads to be the amount of threads multiplied by the amount of processes
            nthreads = uc;
            for (int id = 0; id < nthreads; id++){
             
                //iterates through bottom edges and finds sites that connect in the next thread
                for (int j=0; j<lamount[id]; j++){

                    //sets the next position in the thread under the current thread
                    int next_pos;
                    if ((global_lower[id][j].pos + L) > (split*world_rank)){
                        next_pos = (global_lower[id][j].pos + L)%(split*world_rank) + (split*(world_rank-1));
                    }
                    else{
                        next_pos = (global_lower[id][j].pos + L)%(split*world_rank);
                    }

                    for(int k=0; k<uamount[(id+1)%nthreads]; k++){
                       

                       //CHECKS IF THERE IS A SITES DIRECTLY UNDERNEATH CURRENT SITE IN NEXT THREAD
                        if (global_upper[(id+1)%nthreads][k].pos == next_pos){ 
                          
                            int next_cols[split], nc=0, next_rows[split], nr=0;

                            //ADDING THE CLUSTER DIRECTLY BELOW
                            if(!isvalueinarray(global_upper[(id+1)%nthreads][k].id, visited_cluster, vc)){
                            
                                visited_cluster[vc++] = global_upper[(id+1)%nthreads][k].id;

                                if (!isvalueinarray(global_lower[id][j].id, visited_cluster, vc))
                                    visited_cluster[vc++] = global_lower[id][j].id;

                                int r = returnvaluein3darraystruct(global_lower[id][j].id, final_cluster, fc);
                                int r2 = returnvaluein3darraystruct(global_upper[(id+1)%nthreads][k].id, final_cluster, fc);
                                
                                if (final_cluster[r2].id==final_cluster[r].id) 
                                    continue;

                                final_cluster[r].size += final_cluster[r2].size; //add the size of the clusters together



								for(int i=0;i<fc;i++){
									if (final_cluster[i].id == final_cluster[r2].id){
										nc = final_cluster[i].coltop;
										nr = final_cluster[i].rowtop;
										for(int j=0; j<nc;j++) 
                                            next_cols[j] = final_cluster[i].cols[j];
										for(int k=0; k<nr;k++)
											next_rows[k] = final_cluster[i].rows[k];
									}
								}
								for(int i=0;i<fc;i++){
									if (final_cluster[i].id == final_cluster[r].id){
										for(int j=0; j<nc; j++){
											if(!isvalueinarray(next_cols[j], final_cluster[i].cols, final_cluster[i].coltop))
												final_cluster[i].cols[final_cluster[i].coltop++] = next_cols[j];  //one not in list
										}
										for(int q=0; q<nr; q++){
											if(!isvalueinarray(next_rows[q], final_cluster[i].rows, final_cluster[i].rowtop))
												final_cluster[i].rows[final_cluster[i].rowtop++] = next_rows[q];
										}
									}
								}


                                for(int i=0; i<lamount[(id+1)%nthreads]; i++){
                                    if (global_lower[(id+1)%nthreads][i].id == global_upper[(id+1)%nthreads][k].id)
                                        global_lower[(id+1)%nthreads][i].id = final_cluster[r].id;
                                }
                                for(int h=0; h<uamount[(id+1)%nthreads]; h++){
                                    if (global_upper[(id+1)%nthreads][h].id == global_upper[(id+1)%nthreads][k].id && global_upper[(id+1)%nthreads][h].pos != global_upper[(id+1)%nthreads][k].pos)
                                        global_upper[(id+1)%nthreads][h].id = final_cluster[r].id;
                                }
                                global_upper[(id+1)%nthreads][k].id = final_cluster[r].id;
                            }

                            //ADDING A CLUSTER THAT DOES NOT PERCOLATE THAT IS ABOVE
                            else if (isvalueinarray(global_upper[(id+1)%nthreads][k].id, visited_cluster, vc)
                                && !isvalueinarray(global_lower[id][j].id, visited_cluster, vc) && (id != numthreads-1)){

                                visited_cluster[vc++] = global_lower[id][j].id;
                                

                                int r = returnvaluein3darraystruct(global_lower[id][j].id, final_cluster, fc);
                                int r2 = returnvaluein3darraystruct(global_upper[(id+1)%nthreads][k].id, final_cluster, fc);

                                if (final_cluster[r2].id==final_cluster[r].id)
                                    continue;
                                
                                final_cluster[r2].size += final_cluster[r].size;
								
								for(int i=0;i<fc;i++){
									if (final_cluster[i].id == final_cluster[r].id){
										nc = final_cluster[i].coltop;
										nr = final_cluster[i].rowtop;
										for(int j=0; j<nc;j++)
											next_cols[j] = final_cluster[i].cols[j];
										for(int j=0; j<nr;j++)
											next_rows[j] = final_cluster[i].rows[j];
									}
								}
								for(int i=0;i<fc;i++){
									if (final_cluster[i].id == final_cluster[r2].id){
										for(int j=0; j<nc; j++){
											if(!isvalueinarray(next_cols[j], final_cluster[i].cols, final_cluster[i].coltop))
												final_cluster[i].cols[final_cluster[i].coltop++] = next_cols[j];	
										}
										for(int j=0; j<nr; j++){
											if(!isvalueinarray(next_rows[j], final_cluster[i].rows, final_cluster[i].rowtop))
												final_cluster[i].rows[final_cluster[i].rowtop++] = next_rows[j];
										}
									}
								}
        
                                 for(int i=0; i<lamount[id]; i++){
                                    if (global_lower[id][i].id == global_upper[(id+1)%nthreads][k].id);
                                        global_lower[id][i].id = final_cluster[r2].id;
                                }
                            }

                            //ADDING THE VERY BORDERS AT THE TOP AND BOTTOM TOGETHER PROPERLY
                            else if (isvalueinarray(global_upper[(id+1)%nthreads][k].id, visited_cluster, vc) && id==(nthreads-1) && !isvalueinarray(global_lower[id][j].id, border_cluster, bc)){

								border_cluster[bc++] = global_lower[id][j].id;

                                int r = returnvaluein3darraystruct(global_lower[id][j].id, final_cluster, fc);
                                int r2 = returnvaluein3darraystruct(global_upper[(id+1)%nthreads][k].id, final_cluster, fc);

                                if (final_cluster[r2].id==final_cluster[r].id)
                                    continue;

                                final_cluster[r2].size += final_cluster[r].size;

                                for (int i=0; i<lamount[id]; i++){
                                    if (global_lower[id][i].id == global_upper[(id+1)%nthreads][k].id)
                                        global_lower[id][i].id = final_cluster[r2].id; 
                                }
                            }
                        }
                    }
                }
            }
            free(global_lower);
            free(global_upper);
            

            /* WORK OUT THE CLUSTER SIZE, BY TAKING THE MAX FROM THE STRUCTURE,
                AND WORK OUT THE PERCOLATIONS BASED ON WHETHER THE AMOUNT OF ROWS/COLUMNS
                IS EQUAL TO THE TOTAL AMOUNT OF ROWS/COLUMNS */

            bool row_perco = false, col_perco = false;
            int row_perco_id=-1, col_perco_id=-1, fmax=-1;
            
            for (int i=0; i<fc; i++){
				if(final_cluster[i].coltop == L){
					col_perco = true;
					col_perco_id = final_cluster[i].id;
				}
				if(final_cluster[i].rowtop == L){
					row_perco = true;
					row_perco_id = final_cluster[i].id;
				}
                if (final_cluster[i].size>fmax){
                    fmax = final_cluster[i].size;
                }
            }
            free(final_cluster);

            printf("The biggest cluster is %d sites\n", fmax);
            printf("Row Percolation: %s\n", (row_perco ? "true" : "false"));
            printf("Column Percolation: %s\n", (col_perco ? "true" : "false"));
        }
        
    }

        //int (*aa)[6] = malloc(sizeof(int[M][6]));
    //char (*bond)[L] = malloc(sizeof(char[2*L][L]));

    //makeBond(prob, bond);

    /*NEAREST NEIGHBOURS*/
    /*For each site record the nearest neighbours in an array*/
    
    /*
    for (int i=0; i<M; i++){
        if((i/L)%2==1){            
            //VERTICAL BOND 
            aa[i][0] = (i+2*L)%M;  //SOUTH SOUTH  
            aa[i][1] = (i+L)%M;  //SOUTH EAST
            aa[i][2] = (i-1+L)%M;  //SOUTH WEST              
            aa[i][3] = ((i-2*L)+M)%M;  //NORTH NORTH         
            aa[i][4] = ((i-L)+M)%M;  //NORTH EAST                
            aa[i][5] = (i-1-L);  //NORTH WEST-  
            
            //vertical ends
            if(i%L == 0){
                aa[i][2] = ((i-1)+2*L)%M;
                aa[i][5] = (i-1);
            }
        }

        else if ((i/L)%2 == 0){
            //HORIZONTAL BOND
            aa[i][0] = (i+1);    //EAST EAST   
            aa[i][1] = (((i+1-L)+M)%M);  //EAST NORTH  MAKE IT WRAP AROUND
            aa[i][2] = (i+1+L);  //EAST SOUTH
            aa[i][3] = (i-1);    //WEST WEST      
            aa[i][4] = (((i-L)+M)%M);    //WEST NORTH     
            aa[i][5] = ((i+L)%M);    //WEST SOUTH 
            
            //horizontal ends
            if((i+1)%L == 0){
                aa[i][0] = (i+1-L);
                aa[i][1] = (((i+1)-2*L)+M)%M;
                aa[i][2] = i+1;
            }
            if(i%(2*L) == 0){
                aa[i][3] = (i-1+L);
            }
        }
    }
    */
    /*
    if (world_rank == finalProcess){
    if (option == 'b'){
    
        //int tempi, tempj;
        int rpcount = 0;
        int cpcount = 0;
        int ktop =0 ;
        int keeptrack[M];
        ktop++;

        bool row_flag = false;
        bool col_flag = false;
        bool both_flag = true;
        
        int perc_type;
        printf("Enter '0' or '1' or '2' for row or coloumn or both percolation\n");
        scanf("%d", &perc_type);
        
        if (perc_type == 0) row_flag = true;
        else if (perc_type == 1) col_flag = true;
        else if (perc_type == 2) both_flag = true;
        
        int max = 0;
        
        
        int k;
        for (k = 0; k < M; k++){
            if (isvalueinarray(k, keeptrack, ktop)){
                continue;
            }


            int current_node = 0;
            int vtop2=0, qtop2=0, rowtop=0, coltop=0;
            int queue2[M], visited2[M];

            int tempi = k/L;
            int tempj = k%L;

            int row_arr[L], col_arr[L];

            queue2[qtop2] = -1;
            qtop2++;
            queue2[qtop2] = k;
            
            keeptrack[ktop] = k;
            ktop++;
            
            visited2[vtop2] = k;
            vtop2++;
            
            if (tempi%2 == 1){
                row_arr[rowtop] = tempi;
                rowtop++;
            }
            
            if (tempi%2 == 0){
                col_arr[coltop] = tempj;
                coltop++;
            }

            if (bond[tempi][tempj] == '-') {
                
                while (qtop2 != 0){
                    current_node = queue2[qtop2];

                    int dir;
                    for (dir = 0; dir < 6; dir++){
                        
                        if (bond[aa[current_node][dir]/L][aa[current_node][dir]%L] == '-' && !isvalueinarray(aa[current_node][dir], visited2, vtop2)) {
                            
                            //if next bond is verticle
                            if((aa[current_node][dir]/L)%2==1){
                                if (!isvalueinarray(aa[current_node][dir]/L, row_arr, rowtop)) {
                                row_arr[rowtop] = aa[current_node][dir]/L;
                                rowtop++;
                                }
                            }
                            
                            //if next bond is horizontal
                            if((aa[current_node][dir]/L)%2==0){
                                if (!isvalueinarray(aa[current_node][dir]%L, col_arr, coltop)) {
                                col_arr[coltop] = aa[current_node][dir]%L;
                                coltop++;
                                }
                            }
                        
                            //add the next node to array of visited nodes
                            visited2[vtop2] = aa[current_node][dir];
                            vtop2++;
                            
                            //keep track of all visited nodes 
                            keeptrack[ktop] = aa[current_node][dir];
                            ktop++;
                            
                            //add next node to top of the queue
                            qtop2++;
                            queue2[qtop2] = aa[current_node][dir];
                            break;
                        }
                        //if there is no next_node to go to, decrement the queue
                        else if (dir == 5) {
                            queue2[qtop2] = 0;
                            qtop2--;
                        }
                    }
                }
                if (rowtop == L) rpcount++;
                if(coltop == L) cpcount++;
            }

            if (vtop2 > max) max = vtop2;
        }

        if (row_flag == true){
             if (rpcount > 0) printf("ROW PERCOLATES");
             else printf("NO ROW PERCOLATION");
        }
        else if (col_flag == true){
            if (cpcount > 0) printf("COLUMN PERCOLATES");
            else printf("NO COLUMN PERCOLATION");
        }
        else if (both_flag == true){
            if (cpcount > 0 && rpcount) printf("ROW AND COLUMN PERCOLATES");
            else printf("NO ROW AND COLUMN PERCOLATION");
        }

        free(bond);

        printf("\nThe biggest cluster is %d bonds.", max);

    }
    }
    */
}


/*Checks if a value is an in array, give those and the size of the array*/
bool isvalueinarray(int val, int arr[split], int size){
    int i;
    for (i=0; i < size; i++) {
        if (arr[i] == val)
            return true;
    }    
    return false;
}

bool isvaluein2darraystruct(int val, struct cluster cluster_struct[split], int size){
    int i;
    for (i=0; i < size; i++) {
        if (cluster_struct[i].id == val)
            return true;  
    }
    return false;
}

int returnvaluein3darraystruct(int val, struct cluster cluster_struct[split], int size){

    int t;
    for (t=0; t<size; t++){
        if (cluster_struct[t].id == val)
            return t;
    }
    return -1;
}


/*Returns a random number between 0 and 1*/
double rando(){
    return (double)rand()/(double)RAND_MAX;
}

/*Occupies the grid using the given seeding probability*/
void makeBond(double seeding_prob, char bond[2*L][L]){
    int i, j;
    for(i=0;i<2*L;i++){
        for(j=0;j<L;j++){
            if (rando() < seeding_prob){
                bond[i][j] = '-';
            }
            else{
                bond[i][j] = '.';
            }       
        }
    }
    
}


/*Prints the grid in ASCII to the terminal*/
void printGrid(int occ[L][L]){
    int i,j;
    for(i=0;i<L;i++){
        printf("\n");
        for(j=0;j<L;j++){
            printf("%d ", occ[i][j]);
        }
    }
    printf("\n");
}

void printBond(char bond[2*L][L]){
    int i,j;
    for(i=0;i<2*L;i++){
        printf("\n");
        for(j=0;j<L;j++){
            if (i%2 == 0 && j==0){
                printf("    ");
            }
            if(i%2 == 1 && j == 0){
                printf(" ");
            }
            printf("%c    ", bond[i][j]);
            if (j==L-1){
                    printf("\n");
            }
        }
    }
    printf("\n");
}
