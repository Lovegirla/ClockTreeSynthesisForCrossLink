#include <iostream>
#include <limits>
#include <stdexcept>
#include <boost/math/tools/roots.hpp>
#include <boost/math/constants/constants.hpp>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
const double max_double_num = std::numeric_limits<unsigned int>::max();
const double r = 0.0001;
const double c = 0.0002;
const double Cap_Limit = 70;
struct  TreeNode
{
    /* data */
    int id;
    double x;
    double y;
    double slew;
    double Ceff;
    double cap;
    double delay;
    bool is_sink;
    bool is_buffer;
    TreeNode *lchild;
    TreeNode *rchild;

    TreeNode(){}
    TreeNode(int _id,double _x , double _y,double _cap){
        id=_id;
        x = _x;
        y = _y;
        cap = _cap;
    }
};
typedef struct Equation01 {
    double target_y;

    Equation01(double y) : target_y(y) {}

    std::tuple<double, double> operator()(double x) const {
            return std::make_tuple(
                x * (1 - x * (1 - exp(-1 / x))) - target_y,
                1 - x * (1 - exp(-1 / x))
            );
        }
} Slew_Equation;
/* 
* ==================================================================
* description:  获取两点之间得到曼哈顿距离
* function:   GetL1Distance  
* Author: WenJiazhi.
* ===================================================================
*/
double GetL1Distance(TreeNode *p1, TreeNode *p2){
    return abs(p1->x - p2->x) + abs(p1->y -p2->y);
}
/* 
* ==================================================================
* description: Input_Slew_Equation
* function: 给出输出slew根据反向传播计算输入slew
            根具子节点slew求父节点slew
* Author: WenJiazhi.
* ===================================================================
*/
double Input_Slew_Equation(double target_y){

    Slew_Equation equation(target_y);

    double initial_guess = 0.5; // 初始猜测值
    double min_x = 0.0; // x的最小搜索范围
    double max_x = 20.0; // x的最大搜索范围

    boost::uintmax_t max_iter = 100; // 最大迭代次数
    double tolerance = 1e-6; // 收敛容忍度
    double solution = boost::math::tools::newton_raphson_iterate(equation, initial_guess, min_x, max_x, tolerance, max_iter);

    return solution;
}
/* 
* ==================================================================
* description: 查找计算时钟各节点中的slew值
* function:  FindSlew
* Author: WenJiazhi.
* ===================================================================
*/
double FindSlew(TreeNode *p){
    double slew;
    if(p->is_sink){
        return p->slew;
    }else{
        if(p->lchild && p->rchild){
            double slew_l = p->lchild ? p->lchild->slew : max_double_num;
            double slew_r = p->rchild ? p->rchild->slew : max_double_num;
            double R1 = GetL1Distance(p,p->lchild)*r;
            double R2 = GetL1Distance(p,p->rchild)*r;
            double C1 = p->lchild->Ceff + 0.5*c*GetL1Distance(p,p->lchild);
            double C2 = p->rchild->Ceff + 0.5*c*GetL1Distance(p,p->rchild);
            double y1 = R1*C1/slew_l;
            double y2 = R2*C2/slew_r;
            double x1 = Input_Slew_Equation(y1);
            double x2 = Input_Slew_Equation(y2);
            double slew1 =  R1*C1/x1;
            double slew2 =  R2*C2/x2;
            return std::min(slew1,slew2);
        }else if(p->lchild){
            double slew_l = p->lchild ? p->lchild->slew : max_double_num;
            double R1 = GetL1Distance(p,p->lchild)*r;
            double C1 = p->lchild->Ceff + 0.5*c*GetL1Distance(p,p->lchild);
            double y1 = R1*C1/slew_l;
            double x1 = Input_Slew_Equation(y1);
            return R1*C1/x1;
        }else{
            double slew_r = p->rchild ? p->rchild->slew : max_double_num;
            double R2 = GetL1Distance(p,p->rchild)*r;
            double C2 = p->rchild->Ceff + 0.5*c*GetL1Distance(p,p->rchild);
            double y2 = R2*C2/slew_r;
            double x2 = Input_Slew_Equation(y2);
            return R2*C2/x2;
        }
    }

}
/* 
* ==================================================================
* description: 计算节点p的有效下游电容
* function: FindEffectiveCap
* Author: WenJiazhi.
* ===================================================================
*/
double FindEffectiveCap(TreeNode* p){
    double Ceff1;
    double Ceff2;
    if(p->is_sink){
        return p->cap;
    }else{
        double K1,C1,K2,C2,Ceff;
        if(p->lchild){
            Ceff1 = p->lchild->Ceff;
            double R1 = r*GetL1Distance(p,p->lchild);
            C1 = Ceff1 + 0.5*c*GetL1Distance(p,p->lchild);
            double K1 = R1*C1/p->slew; 
        }
        if(p->rchild){
            Ceff2 = p->rchild->Ceff;
            double R2 = r*GetL1Distance(p,p->rchild);
            C2 = Ceff2 + 0.5*c*GetL1Distance(p,p->rchild);
            double K2 = R2*C2/p->slew; 
        }
        Ceff = K1*C1 + K2*C2 + 0.05*c*(GetL1Distance(p,p->lchild)+GetL1Distance(p,p->rchild));
        return Ceff;
    }
}
/* 
* ==================================================================
* description:  DME中合并需要的线段长度计算，返回需要总线长
* function: DME_Merge_Cost
* Author: WenJiazhi.
* ===================================================================
*/
double DME_Merge_Cost(TreeNode *Ti,TreeNode *Tj){
    double distance = GetL1Distance(Ti,Tj);
    double __delay = abs(Ti->delay - Tj->delay);
    if(Ti->delay > Tj->delay){
        double this_delay = r*c*distance*distance/2 + r*distance*Tj->Ceff;
        if(this_delay >= __delay){
            return distance;
        }else{
            return max_double_num;
        }
    }else if(Ti->delay < Tj->delay){
        double this_delay = r*c*distance*distance/2 + r*distance*Ti->Ceff;
        if(this_delay >= __delay){
            return distance;
        }else{
            return max_double_num;
        }
    }else{
        return distance;
    }
}
/* 
* ==================================================================
* description: 使用zero_merge时计算合并点在连线的距离lp_i
* function: GetSegmentLength
* Author: WenJiazhi.
* ===================================================================
*/
double GetSegmentLength(TreeNode *Ti,TreeNode *Tj){
    // D1 = r*c*l1*l1/2 + r*l1*Ti->Ceff; 
    // D2 = r*c*l2*l2/2 + r*l2*Tj->Ceff; 
    // D1 + Ti->delay  == D2 + Tj->delay;
    double epsilon = 1e-6; // Desired precision
    double left = 0.0;     // Minimum possible value of l1
    double right = GetL1Distance(Ti, Tj); // Maximum possible value of l1

    while (right - left > epsilon) {
        double mid = (left + right) / 2.0;
        double l1 = mid;
        double l2 = GetL1Distance(Ti, Tj) - l1;

        double D1 = r * c * l1 * l1 / 2 + r * l1 * Ti->Ceff;
        double D2 = r * c * l2 * l2 / 2 + r * l2 * Tj->Ceff;

        if (D1 + Ti->delay < D2 + Tj->delay) {
            left = mid;
        } else {
            right = mid;
        }
    }

    return left;
}
/* 
* ==================================================================
* description: 在Ti和Tj线段上更根据l1（合并点到Ti的距离），
                找到合并点，并返回
* function: FindNodeInSegment
* Author: WenJiazhi.
* ===================================================================
*/
TreeNode* FindNodeInSegment(TreeNode *Ti,TreeNode *Tj,double l1){
    double minDistance = l1;
    TreeNode* bestNode = nullptr;
    
    for (double t = 0; t <= 1; t += 0.001) {
        double xMerge = Ti->x + t * (Tj->x - Ti->x);
        double yMerge = Ti->y + t * (Tj->y - Ti->y);
        
        double distance = abs(Ti->x - xMerge) + abs(Ti->y - yMerge);
        
        if (abs(distance - l1) < minDistance) {
            minDistance = abs(distance - l1);
            bestNode = new TreeNode();
            bestNode->x = xMerge;
            bestNode->y = yMerge;
        }
    }

    return bestNode;
}
/* 
* ==================================================================
* description: 找到两个点合并，并返回合并出来的新节点
* function: DME_Merge
* Author: WenJiazhi.
* ===================================================================
*/
TreeNode* DME_Merge(TreeNode *Ti,TreeNode *Tj){
    TreeNode* p = nullptr;
    double distance = GetL1Distance(Ti,Tj);
    double __delay = abs(Ti->delay - Tj->delay);
     if(Ti->delay > Tj->delay){
        double this_delay = r*c*distance*distance/2 + r*distance*Tj->Ceff;
        if(this_delay >= __delay){
            double l1 = GetSegmentLength(Ti,Tj);
            p = FindNodeInSegment(Ti,Tj,l1);
        }else{
            p->x = Ti->x;
            p->y = Ti->y;
        }
    }else if(Ti->delay < Tj->delay){
        double this_delay = r*c*distance*distance/2 + r*distance*Ti->Ceff;
        if(this_delay >= __delay){
           double l1 = GetSegmentLength(Ti,Tj);
           p = FindNodeInSegment(Ti,Tj,l1);
        }else{
            p->x = Tj->x;
            p->y = Tj->y;
        }
    }else{
        double l1 = GetSegmentLength(Ti,Tj);
        p = FindNodeInSegment(Ti,Tj,l1);
    }
    p->lchild = Ti;
    p->rchild = Tj;
    return p;
} 
/* 
* ==================================================================
* description: 计算两个点之间合并需要的代价并返回
* function: FindMergingCost
* Author: WenJiazhi.
* ===================================================================
*/
double FindMergingCost(TreeNode *Ti,TreeNode *Tj){
    double Cost = DME_Merge_Cost(Ti,Tj);
    TreeNode *p = DME_Merge(Ti,Tj);
    p->lchild = Ti;
    p->rchild = Tj;
    double EDSC = FindEffectiveCap(p);
    if(EDSC < Cap_Limit){
        return Cost;
    }
    return max_double_num;
}
/* 
* ==================================================================
* description: 从没有合并的集合中，找到一对点并返回
* function: GetSubTreesToBeMerged
* Author: WenJiazhi.
* ===================================================================
*/
std::pair<TreeNode *,TreeNode *> GetSubTreesToBeMerged(std::vector<TreeNode*> &U,std::vector<TreeNode*> &F){
    int PairsFound = 0;
    TreeNode *node_1 = nullptr;
    TreeNode *node_2 = nullptr;
    double MergingCost = max_double_num; 
    while(PairsFound != 1 and U.size()>1){
        double delay = max_double_num;
        for(auto node : U){
            if(node->delay <= delay){
                node_1 = node;
                delay = node->delay;
            }
        }
        for(auto node : U){
            if(node == node_1){
                continue;
            }else{
                double cost = FindMergingCost(node_1,node);
                if(cost < MergingCost){
                    node_2 = node;
                    MergingCost = cost;
                }
            }
        }
        if(MergingCost != max_double_num){
          PairsFound = 1;  
        }else{
            U.erase(std::remove(U.begin(), U.end(), node_1), U.end());
            F.push_back(node_1);
        }
    }
    if(MergingCost != max_double_num){
        return std::pair<TreeNode *,TreeNode *>(node_1,node_2);
    }else{
        return std::pair<TreeNode *,TreeNode *>(NULL,NULL);
    }
}
/* 
* ==================================================================
* description: 读取相应的输入文件，将sinks相关信息读出来
* function: ReadInputfile
* Author: WenJiazhi.
* ===================================================================
*/
std::vector<TreeNode*> ReadInputfile(std::string filename){
    std::vector<TreeNode*> sinks;
    std::ifstream infile(filename);
    std::string line;
    getline(infile, line);
    std::istringstream iss(line);

    while (getline(infile, line)) {
        if (line.find("sink") != std::string::npos) {
            break;
        }
    }
    while (getline(infile, line)) {
        if (line.find("num wirelib") != std::string::npos) {
            break;
        }
        std::istringstream iss(line);
        int id;
        double capacitance, x, y;

        iss >> id >> x >> y >> capacitance;
        TreeNode *sink = new TreeNode(id,x,y,capacitance);
        sink->is_sink = true;
        sinks.push_back(sink);
        
    }
    return sinks;
}
/* 
* ==================================================================
* description: 对点集插入buffer,并返回点集
* function: InsertBuffer
* Author: WenJiazhi.
* ===================================================================
*/
std::vector<TreeNode*> InsertBuffer(std::vector<TreeNode*> nodes){
    std::vector<TreeNode*> new_nodes;
    for (auto node : nodes)
    {
        node->is_buffer = true;
        new_nodes.push_back(node);
    }
    return new_nodes;
    
}
/* 
* ==================================================================
* description: 平衡树top_level框架，描述构建过程
* function: TopLevelFlow
* Author: WenJiazhi.
* ===================================================================
*/
void TopLevelFlow(std::vector<TreeNode*> sinks){
    
    std::vector<TreeNode*> F;
    std::vector<TreeNode*> U;//初始sinknode
    for(auto sink : sinks){
        U.push_back(sink);
    }
    while((U.size()+F.size())>1){
        auto pair = GetSubTreesToBeMerged(U,F);
        if(pair.first != NULL && pair.second != NULL){
            TreeNode* p =  DME_Merge(pair.first,pair.second);
            p->slew = FindSlew(p);
            p->Ceff = FindEffectiveCap(p);
            U.erase(std::remove(U.begin(), U.end(), pair.first), U.end());
            U.erase(std::remove(U.begin(), U.end(), pair.second), U.end());
            U.push_back(p);
        }else{
            F = InsertBuffer(F);
            //跟新点的信息；
            for(auto node : F){
                U.push_back(node);
            }
            F.clear();
        }
    }

}
int main(int  argc,  char*  argv[]){
    std::vector<TreeNode*> sinks = ReadInputfile(argv[1]);
    TopLevelFlow(sinks);
    return 0;
}