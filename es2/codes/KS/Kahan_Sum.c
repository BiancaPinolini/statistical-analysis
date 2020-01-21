/*********************************************************
 * 
 *  Kahan_Sum.c
 *
 *  Compute the Euler-Mascheroni constant with Kahan 
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


/******** Functions for the Kahan sum algorithm ********/

// FLOAT

float Kahan_Sum_f(int n){
    // Compute the sum of 1/i with i in [1, n] with
    // the Kahan sum algorithm

    float s = 0.;
    float c = 0.;
    float t = 0.;
    float y = 0.;

    for (int i = 1; i<=n; i+=1){
        y = 1./(float)i - c;
        t = s + y;
        c = (t - s) - y;
        s = t;
    }
    return s;
}

float gamma_K_f(int n){
    // Compute an approximation of the gamma function
    // up to n iterations using the Kahan sum algotrithm 

    float gamma_t = 0.;
    gamma_t = Kahan_Sum_f(n) - logf((float)n);
    return gamma_t;
}

// LONG DOUBLE

long double Kahan_Sum_ld(int n){
    // Compute the sum of 1/i with i in [1, n] with
    // the Kahan sum algorithm

    long double s = 0.;
    long double c = 0.;
    long double t = 0.;
    long double y = 0.;

    for (int i = 1; i<=n; i+=1){
        y = 1./(float)i - c;
        t = s + y;
        c = (t - s) - y;
        s = t;
    }

    return s;
}

long double gamma_K_ld(int n){
    // Compute an approximation of the gamma function
    // up to n iterations using the Kahan sum algotrithm 


    long double gamma_t = 0.;
    gamma_t = Kahan_Sum_ld(n) - logl((long double)n);
    return gamma_t;
}

//////////////////////////////////////////////////////////////

int main (){

    // True values of gamma in the desired sizes
    long double gamma_true_ld = 0.57721566490153286;
    float       gamma_true_f  = 0.5772156;


    FILE* errors;
    FILE* gammas;
    
    float       f, err_f;
    long double ld, err_ld;
    int         iter = 0;

    errors = fopen("Kahan_Sum.txt", "w+");
    gammas = fopen("K_S_gammas.txt", "w+");

    for (int i=1; i<=9; i++){

        iter = (int)pow(10, i);

        f     = gamma_K_f(iter);
        err_f = (f - gamma_true_f) / gamma_true_f;
        fprintf(errors, "%.70f\n", fabsf(err_f));
        fprintf(gammas, "%.70f\n", f);

        ld     = gamma_K_ld(iter);
        err_ld = (ld - gamma_true_ld) / gamma_true_ld;
        fprintf(errors, "%.70Lf\n", fabsl(err_ld));
        fprintf(gammas, "%.70Lf\n", ld);


        printf("%d\n", i);
    }

}