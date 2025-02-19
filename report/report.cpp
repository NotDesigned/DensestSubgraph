#include "report.h"
void Report::add_graph_size(int n, int m){
    graph_size.push_back(std::make_pair(n,m));
}
void Report::add_density(int n, double d){
    if(!density.empty() && density.back().second == d){
        return;
    }
    density.push_back(std::make_pair(n,d));
}
void Report::add_total_run_time(double t){
    total_run_time = t;
}
void Report::add_final_density(double d){
    final_density = d;
}
void Report::add_total_iteration(int i){
    total_iteration = i;
}
void Report::print(){
    std::cout << "Graph size: " << std::endl;
    for(auto i: graph_size){
        std::cout << i.first << " " << i.second << std::endl;
    }
    std::cout << "Density: " << std::endl;
    for(auto i: density){
        std::cout << i.first << " " << i.second << std::endl;
    }
    std::cout << "Total run time: " << total_run_time << std::endl;
    std::cout << "Final density: " << final_density << std::endl;
    std::cout << "Total iteration: " << total_iteration << std::endl;
}