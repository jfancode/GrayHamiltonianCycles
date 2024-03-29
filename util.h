/* util.h
   
   Print, error, and housekeeping utilities for GrayHamiltonianSearch.
   Created by jfan for CS49 Parallel Algorithms.
*/

#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* ===================== STRUCTS AND MACROS ======================= */

/* Struct used to hold arguments received from the command line */
typedef struct args_s {
  int error;             // 0 means OK, 1 means a parsing problem
  const char *errMsg;    // error message to print of status not OK
  int *radices;          // array of mixed radices
  int k;                 // the number of radices in the radices array
  int n;                 // sequence length n
  int v;                 // verbosity argument
} args_t;

/* ========================= FUNCTIONS ============================ */

/* Usage report printed from this directory's README file. */
void displayReadme();

/* Loops through argv and determines whether the caller asked for
   help. 
*/
int askHelp(int argc, char **argv);

/* Populates the args_t struct with the specified error message */
args_t* setError(args_t* argStruct, const char* errMsg);

/* Evaluates whether the arguments parse correctly. Returns a struct
   with the argument results.
 */
args_t* parseArguments(int argc, char **argv);

/* Frees all dynamically allocated memory in argStruct */
void freeArgs(args_t* argStruct);

/* Prints the integer array */
void printv(int length, int* array);

/* Prints the integer array in reverse */
void printvReverse(int length, int* array);

/* Reverses the array so that the msb is located at the end of the
   array
 */
int* reverse(int* toReverse, int length);

/* Copies the array */
int* copy(int* toCopy, int length);

/* Computes the product of the radices */
int radixProduct(int k, int* radices);

/* Computes the decimal value from a digit representation */
int decimal(int* rep, int k, int* radices);

#endif
