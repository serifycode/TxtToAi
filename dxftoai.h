#ifndef DXFTOAI_H
#define DXFTOAI_H

#include <vector>
#include <iostream>

using namespace std;

// #define LINE 1
// #define BESEIER 4

struct POINTF{
    double x;
    double y;
    //重载运算符
    POINTF operator-(const POINTF& pt)
    {
        POINTF rr;
        rr.x=this->x-pt.x;
        rr.y=this->y-pt.y;
        return rr;
    }
    POINTF operator+(const POINTF& pt)
    {
        POINTF rr;
        rr.x=this->x+pt.x;
        rr.y=this->y+pt.y;
        return rr;
    }
    POINTF operator*(const double& k)
    {
        POINTF rr;
        rr.x=this->x*k;
        rr.y=this->y*k;
        return rr;
    }
};



struct ai_besier
{
    int kind;//0：start,1:line,3:besier
    POINTF pt;
};

class dxftoai
{
    typedef vector<POINTF> dxf_identity;
    typedef vector<ai_besier> ai_identity;
public:
    void readdxf(const string& filename);//读入dxf

    void readtxt(const string& filename);//读入txt，生成开始的点

    void simplise();//简化点

    void handleper();//将点转成贝塞尔

    void writeai(string filename);//写出ai
private:
    vector<dxf_identity> original_pp;//初始读入的点集

    vector<dxf_identity> simple_pp;//简化后的点集

    vector<ai_identity> ai_pp;

private:
    vector<ai_besier> depart(vector<POINTF>& curve);//

    vector<ai_besier> triangle(vector<POINTF>& dots);//法1：将点集转成四点控制贝塞尔

    void adjust(vector<POINTF>& K,vector<POINTF>& A,POINTF& Q,
                double& ll,double& rr,int& pos,double& dis_all);//输出

    vector<ai_besier> matrix_algorithm(vector<POINTF>& dots);//法2：矩阵相乘求控制点

    int PointToLine(const POINTF& P1,const POINTF& P2,const POINTF& T);//判断三点为直线还是曲线
};





#endif