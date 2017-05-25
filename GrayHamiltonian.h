/* GrayHamiltonianSearch.h
   
   Types and macros for GrayHamiltonianSearch.
   Created by jfan for CS49 Parallel Algorithms.
*/

#ifndef GRAY_HAMILTONIAN
#define GRAY_HAMILTONIAN

#include <math.h>
#include <map>

/* ===================== STRUCTS AND MACROS ======================= */

/* A set of LUTs that allows lookup from integer value to digit
   representation and vice versa
typedef struct LUTs {
  int** valueToRep;    // array of int* representations. Index v into
                       // valueToRep to get the representation of v
  std::map<int*, int> repToValue;     // map of int* representations
                                      // to their int values
} LUTs; 
*/

/* A Gray adjacency graph on which to perform our modified DFS for a
   Hamiltonian cycle 
*/
typedef int** graph_t;

typedef struct search_s {
  int found;       // 1 means a cycle has been found, 0 means none
  int* cycle;      // if a cycle is found, we return the solution
} search_t;

/* ========================= FUNCTIONS ============================ */
  
/* Creates the conversion tables of type LUTs, which map each integer
   value between 0 and nHigh-1 to their digit representations in the
   given radices, and vice versa. */
void populateConversionTable(int** valueToRep,
			     int nHigh, int k, int* radices);

/* Frees the conversion table */
void freeConversionTable(int** valueToRep, int nHigh) {
  for (int v = 0; v < nHigh; v++)
    free(valueToRep[v]);
  free(valueToRep);
  //free(conversionTable);
}

/* Creates the Gray adjacency graph, a 2D integer array where graph[v]
   lists the integers that share the Gray code property with v.
*/
void populateGraph(graph_t graph, int** valueToRep,
		   int n_, int k, int* radices);

/* Frees the Gray adjaceny graph */
void freeGraph(graph_t graph, int n_) {
  for (int v = 0; v < n_; v++)
    free(graph[v]);
  free(graph);
}

/* Frees the processor responsibility array */
void freeProcN(int** procN, int P) {
  for (int p = 0; p < P; p++)
    free(procN[p]);
  free(procN);
}

/* Increments the digit representation given the radices */
int* increment(int* rep, int k, int* radices);

/* Checks whether two digit representations hold the modular Gray code
   property with the given radices */
int checkGray(int* rep1, int* rep2, int k, int* radices);

/* Main algorithm: perform recursive backtracking Hamiltonian search
   on the given graph. Meanwhile, check off vertices in visited and
   build the cycle in solution */
search_t* HamiltonianSearch(graph_t graph, int** valueToRep,
			    int* visited, int* solution, int endIndex,
			    int n_, int k, int* radices);

/* Part of main that can be executed in parallel 
void mainParallel(int P, int** valueToRep,
		  int nLow, int nHigh, int k, int* radices);
*/

#endif
