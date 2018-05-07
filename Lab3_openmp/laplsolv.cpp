//
// Created by David Tran on 2018-05-04.
//
#include <ctime>
#include <iostream>
#include <vector>
#include <omp.h>
int main(int argc, char* argv[]){
    /*Variable initialization*/
    int n{1000}, maxiter{1000}, i{0}, j{0}, k{0};
    double tol{1.0e-3}, error{0.0}, x{0.0};
    std::vector<std::vector<double>> T(n+1,std::vector<double>(n+1));
    std::vector<double> tmp1(n), tmp2(n);
    std::string str;
    std::clock_t t0, t1;

    /*Set boundaries and initial values for the unknowns*/
    for(int i = 0; i < n+1; i++) {
        T[i][0] = 1.0;
        T[i][n] = 1.0;
        T[n][i] = 2.0;
    }
    T[n][0] = 2.0;

    t0 = std::clock();

    /* Do something */

    t1 = std::clock();
    double time_elapsed_ms = 1000.0 * (t1-t0) / CLOCKS_PER_SEC;
    printf("Number of iterations: %d \n", k);
    printf("CPU time used: %g ms \n", time_elapsed_ms);

    return 0;
}