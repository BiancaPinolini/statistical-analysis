/*********************************************************
 *
 *  Pairwise_Sum.c
 *
 *  Compute the Euler-Mascheroni constant with pairwise
 *  sum method.
 *
 *  gamma = (sum_i^N)(1/i) - ln(N)
 *
 *
 ********************************************************/

#define MAIN_C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/******** Functions for the Pairwise sum algorithm ********/


long double *nums_ld;
float *nums_f;
long double sld;
float sf;

// FLOAT

void PWf(int start, int end){

    int len = end - start + 1;

    if (len<=2){
      for (int i = 0; i<len; i++) {
        sf += nums_f[start+i];
      }
    }

    else {
        int m = (int)(len/2);

        PWf(start,start+m-1);
        PWf(start+m,end);
    }

}

void Pairwise_Sum_f(int n){
    // Compute the sum of 1/i with i in [1, n] with
    // the Pairwise sum algorithm


    int j = 0;

    if (n == 10){

      nums_f = (float*)calloc(n, sizeof(float));

      printf("Ha creato nums_f FLOAT\n");

      if (nums_f == NULL) {
          printf("FALLITO!\n");
      }

      for(int i=1; i<=n; i++){
          nums_f[j] = 1/(float)i;
          j ++;
      }

      printf("Ha riempito nums_f FLOAT\n");

      PWf(0, n-1);

      printf("Ha fatto la somma FLOAT\n");

    } else{

      nums_f = (float*)calloc((int)(n/4), sizeof(float));

      printf("Ha creato nums_f FLOAT\n");

      if (nums_f == NULL) {
          printf("FALLITO!\n");
      }

      float a, b, k;


      for(int i=1; i<=n; i=i+4){
          k = (float) i;
          a = 1/k + 1/(k+1.);
          b = 1/(k+2.) + 1/(k+3.);

          nums_f[j] = a+b;
          j ++;
      }

      printf("Ha riempito nums_f FLOAT\n");

      PWf(0, (int)(n/4)-1);

      printf("Ha fatto la somma FLOAT\n");

    }



}

float gamma_P_f(int n){
    // Compute an approximation of the gamma function
    // up to n iterations using the Pairwise sum algotrithm

    float gamma_t;
    Pairwise_Sum_f(n);
    gamma_t = sf - logl((long double)n);
    return gamma_t;
}

// LONG DOUBLE

void PWld(int start, int end){

    int len = end - start + 1;

    if (len<=2){
      for (int i = 0; i<len; i++) {
        //printf("nums_ld[%d] = %.50Lf\t", start+i, nums_ld[start+i]);
        sld += nums_ld[start+i];
        //printf("sld = %Lf\n", sld);
      }
    }

    else {
        int m = (int)(len/2);

/*
        printf("start = %d\t", start);
        printf("end = %d\t", end);
        printf("m = %d\n", m);
*/

        PWld(start,start+m-1);
        PWld(start+m,end);
    }

}

void Pairwise_Sum_ld(int n){
    // Compute the sum of 1/i with i in [1, n] with
    // the Pairwise sum algorithm

    int j = 0;

    if (n == 10) {

      nums_ld = (long double*)calloc(n, sizeof(long double));

      printf("Ha creato nums_ld LONG DOUBLE\n");

      if (nums_ld == NULL) {
          printf("FALLITO!\n");
      }

      for(int i=1; i<=n; i++){
          nums_ld[j] = 1/(long double)i;
          j ++;
      }

      printf("Ha riempito nums_ld LONG DOUBLE\n");

      PWld(0, n-1);

      printf("Ha fatto la somma LONG DOUBLE\n");

    } else {

      nums_ld = (long double*)calloc((int)(n/4), sizeof(long double));

      printf("Ha creato nums_ld LONG DOUBLE\n");

      if (nums_ld == NULL) {
          printf("FALLITO!\n");
      }

      long double a, b, k;


      for(int i=1; i<=n; i=i+4){
          k = (long double) i;
          a = 1/k + 1/(k+1.);
          b = 1/(k+2.) + 1/(k+3.);

          nums_ld[j] = a+b;
          j ++;
      }

      printf("Ha riempito nums_ld LONG DOUBLE\n");

      PWld(0, (int)(n/4)-1);

      printf("Ha fatto la somma LONG DOUBLE\n");

    }



}

long double gamma_P_ld(int n){
    // Compute an approximation of the gamma function
    // up to n iterations using the Pairwise sum algotrithm

    long double gamma_t;
    Pairwise_Sum_ld(n);
    gamma_t = sld - logl((long double)n);
    return gamma_t;
}

//
////////////////////////////////////////////////////////////

int main (){

    // True values of gamma in the desired sizes
    long double gamma_true_ld = 0.57721566490153286;
    float       gamma_true_f  = 0.5772156;

    // Files
    FILE* errors;
    FILE* gammas;

    float       f, err_f;
    long double ld, err_ld;
    int         iter = 0;

    errors = fopen("Pairwise_Sum.txt", "w+");
    gammas = fopen("P_S_gammas.txt", "w+");

    for (int i=1; i<=9; i++){

        printf("\n \n %d\n", i);

        iter = (int)pow(10, i);

        sf = 0.;
        sld = 0.;
        free(nums_f);
        free(nums_ld);


        f     = gamma_P_f(iter);
        err_f = (f - gamma_true_f)/gamma_true_f;
        fprintf(errors, "%.70f\n", fabsf(err_f));
        fprintf(gammas, "%.70f\n", f);

        ld     = gamma_P_ld(iter);
        err_ld = (ld - gamma_true_ld)/gamma_true_ld;
        fprintf(errors, "%.70Lf\n", fabsl(err_ld));
        fprintf(gammas, "%.70Lf\n", ld);
    }

printf("THE END\n");



return 0;
}
