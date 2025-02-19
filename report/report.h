#include<iostream>
#include<vector>
class Report{
    private:
        std::vector<std::pair<int, int> > graph_size;
        std::vector<std::pair<int,double> > density; // iter, density
        double total_run_time;
        double final_density;
        int total_iteration;
    public:
        void add_graph_size(int n, int m);
        void add_density(int n, double d);
        void print(); 
        void add_total_run_time(double t);
        void add_final_density(double d);
        void add_total_iteration(int i);
};