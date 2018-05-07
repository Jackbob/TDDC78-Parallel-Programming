//
// Created by David Tran on 2018-05-04.
//

#include <ctime>
#include <iostream>
#include <vector>
#include <omp.h>
#include <algorithm>
#include <cmath>
int main(int argc, char* argv[]){
    /*Variable initialization*/
    int n{5}, maxiter{1000}, i{0}, j{1}, k{1};
    double tol{1.0e-3}, error{0.0}, x{0.0};
    std::vector<std::vector<double>> T(n+2,std::vector<double>(n+2));
    std::vector<double> tmp1(n), tmp2(n), vec1(n), vec2(n), vec3(n);
    std::string str;
    std::clock_t t0, t1;

    /*Set boundaries and initial values for the unknowns*/
    for(int i = 0; i < n+2; i++) {
        T[0][i] = 1.0;
        T[n+1][i] = 1.0;
        T[i][n+1] = 2.0;
    }
    T[0][n+1] = 2.0;

    t0 = std::clock();

    //Heat conduction
    for(k = 1; k <= maxiter;k++){
        tmp1.assign(T[0].begin()+1, T[0].end()-1);
        error = 0.0;

        for(j = 1; j <= n; j++){
            tmp2.assign(T[j].begin()+1, T[j].end()-1);

            vec1.assign(T[j].begin(), T[j].end()-2);
            vec2.assign(T[j].begin()+2, T[j].end());
            vec3.assign(T[j+1].begin()+1, T[j+1].end()-1);

            std::transform(vec1.begin(), vec1.end(), vec2.begin(), vec1.begin(), std::plus<double>());
            std::transform(vec1.begin(), vec1.end(), vec3.begin(), vec1.begin(), std::plus<double>());
            std::transform(vec1.begin(), vec1.end(), tmp1.begin(), vec1.begin(), std::plus<double>());

            std::for_each(vec1.begin(), vec1.end(), [](auto& val){
                val/=4.0;
            });

            std::copy(vec1.begin(), vec1.end(), T[j].begin()+1);
            std::vector<double> temp(n);
            std::transform(tmp2.begin(), tmp2.end(), T[j].begin()+1, temp.begin(), [&temp](auto a, auto b){
                return std::abs(a-b);
            });

            error = std::max(error, *std::max_element(temp.begin(), temp.end()));
            tmp1.assign(tmp2.begin(), tmp2.end());
        }
        if(error < tol){
            break;
        }
    }

    t1 = std::clock();
    double time_elapsed_ms = 1000.0 * (t1-t0) / CLOCKS_PER_SEC;
    printf("Number of iterations: %d \n", k);
    printf("CPU time used: %g ms \n", time_elapsed_ms);

    for(auto r : T) {
        std::copy(r.begin(), r.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\n";
    }

    return 0;
}