#ifndef DENSESTSUBGRAPH_LP_H
#define DENSESTSUBGRAPH_LP_H

#include <vector>
#include <iostream>
#include "types.h"
#include "graph.h"
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include <cmath>
#include <queue>


struct Alpha{
    ui id_first;
    ui id_second;
    double weight_first;
    double weight_second;
};

class LinearProgramming{

private:
    bool is_directed_;

public:
    ui nodes_count_;
    ui edges_count_;
    ui type_;
    ui sort_type;
    ui cur_iter_num;
    double result, last_result;
    std::vector<int> perm;
    std::vector<std::pair<double,double> > w;
public:
    std::vector<std::vector<double>> r;
    std::vector<Alpha> alpha;
    std::vector<Alpha> beta;
    std::vector<double> weight;

    class adam{
    public:
        double beta1, beta2, gamma;
        double beta1_t, beta2_t;
        double alpha;
        double epsilon;
        double limit;
        double mn;
        int iter;
        std::vector<std::vector<double>> m;
        std::vector<std::vector<double>> v;
        adam(){
            m.resize(2);
            v.resize(2);
        }
        adam(double beta1, double beta2, double gamma, double alpha, double epsilon, double limit, double iter, ui edges_count){
            initadam(beta1,beta2,gamma,alpha,epsilon,limit,iter,edges_count);
        }
        void initadam(double beta1,double beta_2, double gamma, double alpha, double epsilon, double limit, double iter, ui edges_count){
            this->beta1 = beta1;
            this->beta1_t = 1;
            this->beta2 = beta_2;
            this->beta2_t = 1;
            this->gamma = gamma;
            this->alpha = alpha;
            this->epsilon = epsilon;
            this->limit = limit;
            this->iter = iter;
            this->mn = std::numeric_limits<double>::max();
            m.resize(2);
            v.resize(2);
            for(ui i = 0; i < 2; i++){
                m[i].assign(edges_count,0);
                v[i].assign(edges_count,0);
            }
        }
    public:
        void update(std::vector<std::vector<double>> &grad, std::vector<Alpha> &theta, ui cur_iter_num, std::pair<double,double> (*Proj) (double,double))
        {
            ui edges_count = theta.size();
            alpha *= gamma;
	    //if (alpha < limit) alpha = limit;
            beta1_t *= beta1;
            beta2_t *= beta2;
            for(ui j = 0; j < edges_count; j++){
                double _u = theta[j].id_first, _v=theta[j].id_second;
                double g0 = grad[0][_u],g1 = grad[1][_v];
                
                //m[0][j] = beta1 * m[0][j] + (1 - beta1) * g0;
                //m[1][j] = beta1 * m[1][j] + (1 - beta1) * g1;
                
                v[0][j] = beta2 * v[0][j] + (1 - beta2) * g0 * g0;
                v[1][j] = beta2 * v[1][j] + (1 - beta2) * g1 * g1;
                //double m_hat0 = m[0][j] / (1 - beta1_t),m_hat1 = m[1][j] / (1 - beta1_t);
                double v_hat0 = v[0][j] / (1 - beta2_t), v_hat1 = v[1][j] / (1 - beta2_t);

                double minimum0 = g0 * limit, minimum1 = g1 * limit;
                //theta[j].weight_first -= alpha * g0 / (sqrt(v_hat0) + epsilon);
                //theta[j].weight_second -= alpha * g1 / (sqrt(v_hat1) + epsilon);

                theta[j].weight_first -= std::max(alpha * g0 / (sqrt(v_hat0) + epsilon), minimum0);
                theta[j].weight_second -= std::max(alpha * g1 / (sqrt(v_hat1) + epsilon), minimum1);

                //theta[j].weight_first -= std::max(alpha * m_hat0 / (sqrt(v_hat0) + epsilon), minimum0);
                //theta[j].weight_second -= std::max(alpha * m_hat1 / (sqrt(v_hat1) + epsilon), minimum1);

                auto [a,b] = Proj(theta[j].weight_first,theta[j].weight_second);
                theta[j].weight_first = a;
                theta[j].weight_second = b;
            }
        }
        void clear(){
            m.clear();
            v.clear();
        }
    }Adam;


public:
    ~LinearProgramming();
    explicit LinearProgramming(bool is_directed, ui type = 0, ui vertices_count = 0, ui edge_count = 0, ui sort = 0);

    void Iterate(double learning_rate, double ratio = 0, bool is_synchronous = false);

    void FistaIterate(double learning_rate, double t, double ratio = 0, bool is_synchronous = false, bool is_random = false);

    void MWUIterate(ui t, bool is_synchronous = false);

    void Init(Graph &graph, double ratio = 0);

    void sort(Graph &graph);
};

#endif //DENSESTSUBGRAPH_FLOWNETWORK_H
