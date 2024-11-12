#ifndef DENSESTSUBGRAPH_REDUCTION_H
#define DENSESTSUBGRAPH_REDUCTION_H

#include "utility/graph.h"
#include "utility/types.h"
#include "utility/xycore.h"


#include <queue>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>
#include <cmath>
#include <climits>
#include "lp.h"
#include "wcore.h"

class Reduction {
public:
    void xyCoreReduction(Graph &graph, Graph &x_y_core, std::pair<double, double> ratios, double &l, double &r,
                         bool &is_init, bool is_dc, bool is_divide_by_number, bool is_exact, bool is_map,
                         bool is_res, ui res_width, bool is_copy);
    void kCoreReduction(Graph &graph, double &l, double &r);
    void stableSetReduction(Graph &graph, LinearProgramming &lp, std::vector<std::pair<VertexID, VertexID>> &edges, bool stable_set_reduction, bool is_map = false);
    void wCoreReduction(Graph &graph, Graph &subgraph, WCore &w_core);
    void UndirectedkCoreReduction(Graph &graph, LinearProgramming &lp, bool weighted = false);
    void UndirectedStableReduction(Graph &graph, LinearProgramming &lp);
    void InitXYCore(Graph &graph,bool is_exact);
public:
    std::vector<ui> core;
    //XYCore xycore,xycore_inherit;
    struct XYCoreBase{
        int x,y;
        XYCore xycore;
        Graph x_y_core;
        XYCoreBase():x_y_core(true,0){
            x=y=0;
        }
        XYCoreBase(const int &_x, const int &_y, const XYCore &_xycore,const Graph &_x_y_core)
        :x(_x),y(_y),xycore(_xycore),x_y_core(_x_y_core){}
    };
    std::vector<XYCoreBase> xycore_bases;
    XYCoreBase& NearestXYCore(int x2, int y2){// return XYCoreBase with x <= x2 and y <= y2 and minimum edge num 
        XYCoreBase *res = nullptr;
        for(auto &xycore_base: xycore_bases){
            if(xycore_base.x == 0 && xycore_base.y == 0) continue;
            if(xycore_base.x <= x2 && xycore_base.y <= y2){
                if(res == nullptr || 
                                xycore_base.x_y_core.getEdgesCount() 
                                < res->x_y_core.getEdgesCount())
                    res = &xycore_base;
            }
        }
        if(res == nullptr){
            return xycore_bases[0];
        }
        return *res;

    }

private:
    void coreDecomposition(const Graph &graph);

};

#endif //DENSESTSUBGRAPH_REDUCTION_H
