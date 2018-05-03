//
// Created by David Tran on 2018-05-04.
//
#include <ctime>
#include <iostream>
#include <omp.h>
int main(int argc, char* argv[]){
    int target_thread_num = 4;
    omp_set_num_threads(target_thread_num);
    unsigned long times[target_thread_num];

// Initialize all the times
#pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        times[thread_id] = start_time();

        std::cout << "Thread number: " << omp_get_thread_num() << std::endl;

        times[thread_id] = end_time();
    }
}