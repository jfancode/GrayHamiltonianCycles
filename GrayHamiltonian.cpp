/* GrayHamiltonianSearch.cpp

   Cilk program for finding a Hamiltonian cycle in a modular Gray graph.
   Created by jfan for CS49 Parallel Algorithms.

   usage: GrayHamiltonian [radices] -n [sequence length] -v [verbosity] 
     [radices]: a string of mixed or fixed integer radices encased in 
                quotation marks and delimited by spaces.
                e.g. GrayHamiltonian "2 3 4"

		In the case above, the mixed-radix tuple will be
		(2,3,4), listed from most significant to least
		significant radix. (However, the radices will be
		stored as an integer array from lsb to msb.)
  
     -n [sequence length]: (optional) an integer sequence length,
                required if we decide to perform a Gray Hamiltonian
                search only on that particular value of n. Otherwise,
                the program will perform the search on all possible
                values of n starting from n=2 and going to n=(product
                of the radices).  
		e.g. GrayHamiltonian "2 3 4" -n 13

		If n is NOT provided above, the program will search
		for a Hamiltonian cycle using the radices (2,3,4) and
		the sequence lengths n=3 to n=23.

		The program will alert the user if the input sequence
		length is invalid (most commonly, if n surpases the
		highest possible number representable with the
		specified radices)

     -v [verbosity]: (optional) changes the verbosity from default 2,
                which for each value of n prints a progress bar and,
                if no Hamiltonian cycle is found, an appropriate alert
                with the search status.
		e.g. GrayHamiltonian "2 3 4" -v 3

		-v 1 outputs only the values of n for which there is
                     no Hamiltonian cycle
		-v 2 outputs the above, plus progress bars and time 
		     lapse for each value of n tried
		-v 3 outputs the above, plus the Hamiltonian cycle if 
		     one is found
		-v 4 (not recommended) outputs the above, plus the 
		     stack trace as the recursive backtracking search 
		     is running

     --help: prints the help page from the local README file

   main function takes as an arg a set of mixed radices (input msb to
   lsb; stored lsb to msb) and optionally, the sequence length n,
   generates the Gray graph as an 2D adjacency array, and performs a
   modified DFS for a Hamiltonian cycle. It also performs
   preprocessing for those cases of radices and sequence lengths where
   it is known that a Hamiltonian cycle is impossible.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef NOCILK
#define cilk_spawn
#define cilk_sync
#define cilk_for for
#else
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#endif

#include "GrayHamiltonian.h"
#include "util.h"
#include "timer.h"

/* ===================== STRUCTS AND MACROS ======================= */

/* Verbosity argument can be reset in the command line using -v */
int verbosity = 2;

/* ========================= FUNCTIONS ============================ */

/* Creates the valueToRep map for easy conversion one way */
void populateConversionTable(int** valueToRep, int nHigh,
			     int k, int* radices) {			     
  /* Store the values and reps in the conversion tables in integer
     order */

  int* rep = (int *) calloc(nHigh, sizeof(int));
  valueToRep[0] = rep;
  //std::map<int*, int> repToValueLocal;
  //repToValueLocal[rep] = 0;
  //conversionTable->repToValue[rep] = 0;

  for (int v = 1; v < nHigh; v++) {
    rep = increment(rep, k, radices);
    valueToRep[v] = rep;
    //repToValueLocal[rep] = v;
    //printf("v = %d; rep = ", v); printv(k, rep); printf("\n");
  }
  //conversionTable->repToValue = repToValueLocal;
  /*  
  for (int v = 0; v < nHigh; v++) {
    printf("v = %d = ", v); printv(k, conversionTable->valueToRep[v]);
    printf(" = %d\n", conversionTable->repToValue[conversionTable->valueToRep[v]]);
  }
  */
}

/* Creates the Gray adjacency graph, a 2D integer array where graph[v]
   lists the integers that share the Gray-code property with v.
*/
void populateGraph(graph_t graph, int** valueToRep,
		   int n_, int k, int* radices) {
  /*
  for (int v = 0; v < n_; v++) {
    printf("v = %d = ", v); printv(k, valueToRep[v]); printf("\n");
    //printf(" = %d\n", conversionTable->repToValue[conversionTable->valueToRep[v]]);
  }
  */

  // build the adjacency list for each vertex v. The list is at most
  // size 2k
  for (int v = 0; v < n_; v++) {
    graph[v] = (int*) calloc(2*k, sizeof(int));
    int* V = valueToRep[v];
    
    int index = 0;
    //printf("v = %d = ", v); printv(k, V); printf(": ");

    // at most two neighbors for each digit
    for (int i = 0; i < k; i++) {

      // upper neighbor u
      int* U = copy(V, k);
      U[i] = (U[i]+1) % radices[i];
      int u = decimal(U, k, radices);
      if (u < n_) {
	graph[v][index++] = u;
	/*
	if (!(checkGray(U, V, k, radices)))
	  printf("ERROR");
	printf("u = %d = ", u); printv(k, U); printf(", ");
	*/
      }
      free(U);
      
      // lower neighbor w
      if (radices[i] != 2) {
	int* W = copy(V, k);
	W[i] = (radices[i]+W[i]-1) % radices[i];
	int w = decimal(W, k, radices);
	if (w < n_) {
	  graph[v][index++] = w;
	  /*
	  if (!(checkGray(W, V, k, radices)))
	    printf("ERROR");
	  printf("w = %d = ", w); printv(k, W); printf(", ");
	  */
	}
	free(W);
      }
    }
    //printf("\n");
  }
}

/* Increments the digit representation given the radices */
int* increment(int* rep, int k, int* radices) {
  int* repCopy = copy(rep, k);
  for (int i = 0; i < k; i++) {
    if (repCopy[i] + 1 == radices[i])
      repCopy[i] = 0;
    else {
      repCopy[i]++;
      break;
    }
  }
  return repCopy;
}

/* Checks whether two digit representations hold the modular Gray code
   property with the given radices */
int checkGray(int* rep1, int* rep2, int k, int* radices) {
  int differ = 0;
  // check every digit
  for (int i = 0; i < k; i++) {

    // if the digits differ, check that the Gray-code property holds
    if (rep1[i] != rep2[i]) {
      if (++differ > 1)
	return 0;
      
      bool Gray = (abs(rep1[i]-rep2[i]) ||
		   (rep1[i] == 0 && rep2[i] == radices[i]-1) ||
		   (rep2[i] == 0 && rep1[i] == radices[i]-1));	
      if (!(Gray))
	return 0;
    }
  }
  return 1;
}

/* Main algorithm: perform recursive backtracking Hamiltonian search
   on the given graph. Meanwhile, check off vertices in visited and
   build the cycle in solution */
search_t* HamiltonianSearch(graph_t graph, int** valueToRep,
			    int* visited, int* solution, int endIndex,
			    int n_, int k, int* radices) {

  /* are there any more vertices to visit?
  int recurse = 1;
  for (int i = 0; i < n_; i++) {
    recurse |= visited[i];
    if (recurse == 0)
      break;
  }
  */

  search_t* searchResults;
  if (endIndex == n_-1) {
    //printf("reached base case\n"); fflush(stdout);
    searchResults = (search_t *) calloc(1, sizeof(search_t));;
    // base case: no vertices left to visit; try to return to 0
    if (checkGray(valueToRep[0], valueToRep[solution[endIndex]],
		  k, radices)) {
      searchResults->found = 1;
      searchResults->cycle = solution;
    }
    // else, searchResults->found = 0 (default)
  } else {
    //printf("reached recursive case, endIndex = %d\n", endIndex);fflush(stdout);
    /* Recursively build the solution Hamiltonian cycle */

    int v = solution[endIndex];
      
    // iterate through all the neighbors u of v
    for (int i = 0; i < 2*k; i++) {
      int u = graph[v][i];
      
      if (visited[u])
	continue;
      else {
	int* visitedCopy = copy(visited, n_);
	int* solutionCopy = copy(solution, n_);

	// visit vertex u and add it to solution
	visitedCopy[u] = 1;
	solutionCopy[endIndex+1] = u;

	searchResults = HamiltonianSearch(graph, valueToRep,
					  visitedCopy, solutionCopy,
					  endIndex+1,
					  n_, k, radices);

	// free up memory; DO NOT TRY TO ACCESS SOLUTION
	free(visitedCopy);
	free(solutionCopy);

	// return right away if a solution has been found
	if (searchResults->found)
	  return searchResults;
	else free(searchResults);
      }
    }
    // otherwise, if we haven't iterated to base case and found a
    // solution, return newly allocated searchResult with no findings
    searchResults = (search_t *) calloc(1, sizeof(search_t));
  }
  return searchResults;
}

/* Part of main that can be executed in parallel
void mainParallel(int P, int** valueToRep,
		  int nLow, int nHigh, int k, int* radices) {

}
*/

/* Driver to perform Hamiltonian search on a Gray graph */
int main(int argc, char **argv) {
  cs49_timer_t timer;/* a timer */
  
  if (argc < 2 || askHelp(argc, argv)) {
    displayReadme();
    exit(1);
  }

  args_t* args = parseArguments(argc, argv);
  if (args->error) {
    printf("ERROR: %s", args->errMsg);
    freeArgs(args);
    exit(1);
  }

  // print the number of processors
  int P;
#ifdef NOCILK
  P = 1;
#else
  P = __cilkrts_get_nworkers();
#endif
  printf("P = %d\n", P);
  
  // make local copies of option arguments
  int n = args->n;
  int k = args->k;
  printf("r = "); printv(k, args->radices);
  int* radices = reverse(args->radices, k);
  verbosity = args->v;
  
  // find the range of n values to iterate through
  int nLow, nHigh;
  if (n != 0) {
    nLow = n;
    nHigh = n;
    printf("; n = %d\n\n", n);
  } else {
    nLow = 2;
    nHigh = radixProduct(k, radices);
    printf("; n = 2 to n = %d\n\n", nHigh);
  }

  // create map that allows easy conversion to digit representation
  int** valueToRep = (int**) calloc(nHigh, sizeof(int*));
  populateConversionTable(valueToRep, nHigh, k, radices);
  
  /* Build list of n values each proc will be responsible for. We
     assign n-value responsibilities in a snake-like fashion */

  int** procN = (int **) calloc(P, sizeof(int*));
  for (int p = 0; p < P; p++) {
    procN[p] = (int *) calloc(nHigh, sizeof(int));

    for (int n_ = 2+p; n_ <= nHigh; n_ += (2*(P-p)-1))
      procN[p][n_] = 1;
  }

  /* Start the timer. */
  TIMER_RESET(timer);
  TIMER_START(timer);
  
  /* Here, we parallelize */
  cilk_for (int p = 0; p < P; p++) {
    // compute Hamiltonian cycles for n_ = nLow to n_ = nHigh that
    // this proc is responsible for
    for (int n_ = nLow; n_ <= nHigh; n_++) {
      // skip computation if another process will handle it
      if (!(procN[p][n_]))
	continue;
      
      graph_t graph = (graph_t) calloc(n_, sizeof(int *));

#ifdef NOCILK
      printf("n_ = %d\n", n_); fflush(stdout);
#else
      printf("p = %d, n_ = %d\n", __cilkrts_get_worker_number(), n_); fflush(stdout);
#endif
      populateGraph(graph, valueToRep, n_, k, radices);

      int* visited = (int *) calloc(n_, sizeof(int));
      int* solution = (int *) calloc(n_, sizeof(int));
      visited[0] = 1;
      search_t* searchResults =  HamiltonianSearch(graph, valueToRep,
						 visited, solution, 0,
						 n_, k, radices);      
      if (!(searchResults->found))
	printf("not found\n");
      else if (args->v >= 3) {
	for (int x = 0; x < n_; x++) {
	  int* V = valueToRep[searchResults->cycle[x]];
	  /*
	    int* U = valueToRep[searchResults->cycle[(x+1) % n_]];
	    if (!(checkGray(V, U, k, radices))) {
	    printf("ERROR: not Gray\n");
	    exit(1);
	    }
	  */	   
	  printvReverse(k, V);
	  printf("\n");
	}
      }
      //printf("\n");

      free(searchResults);
      free(visited);
      free(solution);
      freeGraph(graph, n_);
    }
  }

  /* Stop the timer. */
  TIMER_STOP(timer);
  printf("Time = %f\n", TIMER_EVAL(timer));
  
  freeArgs(args);
  freeConversionTable(valueToRep, nHigh);
  freeProcN(procN, P);
  return 0;
}
