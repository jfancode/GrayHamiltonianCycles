/* util.h
   
   Print, error, and housekeeping utilities for GrayHamiltonianSearch.
   Created by jfan for CS49 Parallel Algorithms.
*/

#include <math.h>

#include "util.h"

/* ===================== STRUCTS AND MACROS ======================= */

/* ========================= FUNCTIONS ============================ */

/* Usage report printed from this directory's README file. */
void displayReadme() {
  # define BUF_SIZE 100
  char buf[BUF_SIZE];
  size_t n;
  
  FILE *fp = fopen("README", "r");
  
  // README file cannot be read
  if (!fp) {
    fprintf(stderr, "Could not display help from README file");
    return;
  }
  
  while ((n = fread(buf, 1, sizeof buf, fp)) > 0)
    fwrite(buf, 1, n, stdout);

  // error while reading
  if (ferror(fp)) {
    fprintf(stderr, "An error occured while reading the README file");
    return;
  }

  fclose(fp);
}

/* Loops through argv and determines whether the caller asked for
   help. 
*/
int askHelp(int argc, char **argv) {
  while (argc-- != 0) {
    if (strcmp(*argv, "--help") == 0)
      return 1;
    argv++;
  }
  return 0;
}

/* Populates the args_t struct with the specified error message */
args_t* setError(args_t* argStruct, const char* errMsg) {
  argStruct->error = 1;
  argStruct->errMsg = errMsg;
  return argStruct;
}

/* Evaluates whether the arguments parse correctly. Returns a struct
   with the argument results.
 */
args_t* parseArguments(int argc, char **argv) {
  args_t* argStruct = (args_t *) calloc(1, sizeof(args_t));

  /* Process and store radices */

  // first argument must be radices
  char* radixStr = argv[1];

  // initialize argStruct's k and radices
  argStruct->k = 0;
  argStruct->radices = (int *) calloc(strlen(radixStr), sizeof(int));
  char *radixC;
  int radix;
  int i = 0;

  radixC = strtok(radixStr, " ");
  // parse each radix from the arguments
  while (radixC != NULL) {
    // check if radix can be converted into an int in valid range
    radix = atoi(radixC);
    if (radix < 2)
      return setError(argStruct, "Each radix must have value > 1.\n");
    argStruct->radices[i++] = radix;
    argStruct->k++;

    radixC = strtok(NULL, " ");
  }

  /* Process and store option arguments */

  int arg;
  while ((arg = getopt(argc, argv, ":n:v:")) != -1) {
    switch (arg) {
    // sequence length n
    case 'n':
      {
	int n = atoi(optarg);
	// check that n is in valid range
	if (n < 2 || n > radixProduct(argStruct->k, argStruct->radices))
	  return setError(argStruct, "n must be > 1 and less than the product of radices.\n");
	argStruct->n = n;
	break;
      }
      
    // verbosity argument v
    case 'v':
      {
	int v = atoi(optarg);
	// check that v is a valid verbosity
	if (v < 1 || v > 4)
	  return setError(argStruct, "Inoperable verbosity. Must be between 0 and 3.\n");
	argStruct->v = v;
	break;
      }
      
    case ':':
      return setError(argStruct, "Option missing an argument.\n");
      
    // unknown option
    case '?':
      return setError(argStruct, "Unknown option");
    }
  }
  return argStruct;
}

/* Frees all dynamically allocated memory in argStruct */
void freeArgs(args_t* argStruct) {
  if (argStruct->radices != NULL)
    free(argStruct->radices);
  free(argStruct);
}

/* Prints the integer array */
void printv(int length, int* array) {
  printf("[");  
  for (int i = 0; i < length; i++) {
    printf("%d", array[i]);
    if (i < length-1)
      printf(", ");
  }
  printf("]");
}

/* Prints the integer array in reverse */
void printvReverse(int length, int* array) {
  printf("[");  
  for (int i = length-1; i >= 0; i--) {
    printf("%d", array[i]);
    if (i > 0)
      printf(", ");
  }
  printf("]");  
}

/* Reverses the array so that the msb is located at the end of the
   array
 */
int* reverse(int* toReverse, int length) {
  int temp;
  for (int i = 0; i < floor(length/2); i++) {
    temp = toReverse[i];
    toReverse[i] = toReverse[length-1-i];
    toReverse[length-1-i] = temp;
  }
  return toReverse;
}

/* Copies the array */
int* copy(int* toCopy, int length) {
  int* copied = (int *) calloc(length, sizeof(int));
  for (int i = 0; i < length; i++)
    copied[i] = toCopy[i];
  return copied;
}

/* Computes the product of the radices */
int radixProduct(int k, int* radices) {		       
  int product = 1;
  for (int i = 0; i < k; i++)
    product *= radices[i];
  return product;
}

/* Computes the decimal value from a digit representation */
int decimal(int* rep, int k, int* radices) {
  int decimal = 0;
  int radixProduct = 1;
  for (int i = 0; i < k; i++) {
    decimal += (rep[i] * radixProduct);
    radixProduct *= radices[i];
  }
  return decimal;
}
