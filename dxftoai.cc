#include "dxftoai.h"
#include "Entity.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "assert.h"
#include <cmath>

using namespace std;

#define pi 3.1415926

//矩阵相乘
vector<vector<double>> MatrixMulti(vector<vector<double>>& A,vector<vector<double>>& B)
{
    int m1=A.size();
    int n1=A[0].size();
    int m2=B.size();
    int n2=B[0].size();
    assert(n1==m2);//判断是否满足矩阵乘法
    vector<vector<double>> res;
    vector<double> line;//单行
    for(int i=0;i<m1;i++)//m1行
    {
        
        for(int j=0;j<n2;j++)//求每个元素
        {
            double ele=0;
            for(int t=0;t<n1;t++)
            {
                ele+=A[i][t]*B[t][j];
            }
            line.push_back(ele);
        }
        res.push_back(line);
        line.clear();
    }
    return res;
}

//重载一下矩阵和点集相乘
vector<POINTF> MatrixMulti(vector<vector<double>>& A,vector<POINTF>& B)
{
    int m=A.size();
    int n=A[0].size();
    assert(B.size()==n);   
    vector<POINTF> res;
    POINTF pt;
    for(int i=0;i<m;i++)
    {
        pt.x=0;pt.y=0;
        for(int j=0;j<n;j++)
        {
            pt.x+=A[i][j]*B[j].x;
            pt.y+=A[i][j]*B[j].y;
        }
        res.push_back(pt);
    }
    return res;
}

//矩阵转置
vector<vector<double>> MatrixTrans(vector<vector<double>>& A)
{
    int m=A.size();
    int n=A[0].size();
    vector<vector<double>> res;
    vector<double> tmp;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            tmp.push_back(A[j][i]);
        }
        res.push_back(tmp);
        tmp.clear();
    }
    return res;
}


//矩阵求逆
bool MatrixInver(vector<vector<double>>& A, vector<vector<double>>& B)
{
    int n=A.size();
    int m=A[0].size();
    assert(m==n);
    int i, j, k;
    double max, temp;
    vector<vector<double>> t=A;  //临时矩阵
    for (i = 0; i < n; i++)       
    {
        for (j = 0; j < n; j++)
        {
            t[i][j] = A[i][j];
        }
    }
    //初始化B矩阵为单位阵
    for (i = 0; i < n; i++)       
    {
        for (j = 0; j < n; j++)
        {
            B[i][j] = (i == j) ? (double)1 : 0;
        }
    }
    for (i = 0; i < n; i++)
    {
        //寻找主元
        max = t[i][i];
        k = i;
        for (j = i+1; j < n; j++)
        {
            if (fabs(t[j][i]) > fabs(max))
            {
                max = t[j][i];
                k = j;
            }
        }
        //如果主元所在行不是第i行，进行行交换
        if (k != i)
        {
            for (j = 0; j < n; j++)
            {
                temp = t[i][j];
                t[i][j] = t[k][j];
                t[k][j] = temp;
                //B伴随交换
                temp = B[i][j];
                B[i][j] = B[k][j];
                B[k][j] = temp;
            }
        }
        //判断主元是否为0, 若是, 则矩阵A不是满秩矩阵,不存在逆矩阵
        if (t[i][i] == 0)
        {
            cout << "There is no inverse matrix!";
            return false;
        }
        //消去A的第i列除去i行以外的各行元素
        temp = t[i][i];
        for (j = 0; j < n; j++)
        {
            t[i][j] = t[i][j] / temp;        //主对角线上的元素变为1
            B[i][j] = B[i][j] / temp;        //伴随计算
        }
        for (j = 0; j < n; j++)        //第0行->第n行
        {
            if (j != i)                //不是第i行
            {
                temp = t[j][i];
                for (k = 0; k < n; k++)        //第j行元素 - i行元素*j列i行元素
                {
                    t[j][k] = t[j][k] - t[i][k]*temp;
                    B[j][k] = B[j][k] - B[i][k]*temp;
                }
            }
        }
    }
    return true;
}

void dxftoai::readdxf(const string& filename)
{
    ifstream f;
    f.open(filename);
    if(!f)
        perror("read dxf error");
    string s;
    double xmax,xmin,ymax,ymin;
    bool flag=false;
    while(getline(f,s))
    {
        if(!flag)
        {
            if(s=="$EXTMIN")
            {
                getline(f,s);
                getline(f,s);
                xmin=stof(s);
                getline(f,s);
                getline(f,s);
                ymin=stof(s);
            }
            if(s=="$EXTMAX")
            {
                getline(f,s);
                getline(f,s);
                xmax=stof(s);
                getline(f,s);
                getline(f,s);
                ymax=stof(s);
            }
            if(s=="ENTITIES")
            {
                flag=true;
            }
            else
            {
                continue;
            }           
        }
        else
        {
            if(s=="AcDbEntity")
            {
                getline(f,s);
                getline(f,s);
                getline(f,s);
                getline(f,s);
                if(s=="AcDbCircle")//圆
                {
                    static double cpara=0.552284749831;//圆转beseir系数，四个半圆
                    getline(f,s);
                    getline(f,s);
                    double x= stof(s);
                    getline(f,s);
                    getline(f,s);
                    double y= stof(s);
                    getline(f,s);
                    getline(f,s);
                    double z= stof(s);
                    getline(f,s);
                    getline(f,s);
                    double r= stof(s);
                    ai_identity circle;
                    circle.push_back({0,{x+r,y}});
                    circle.push_back({3,{x+r,y+r*cpara}});
                    circle.push_back({3,{x+r*cpara,y+r}});
                    circle.push_back({3,{x,y+r}});
                    circle.push_back({3,{x-r*cpara,y+r}});
                    circle.push_back({3,{x-r,y+r*cpara}});
                    circle.push_back({3,{x-r,y}});
                    circle.push_back({3,{x-r,y-r*cpara}});
                    circle.push_back({3,{x-r*cpara,y-r}});
                    circle.push_back({3,{x,y-r}});
                    circle.push_back({3,x+r*cpara,y-r});
                    circle.push_back({3,{x+r,y-r*cpara}});
                    circle.push_back({3,{x+r,y}});
                    ai_pp.push_back(circle);
                }
                if(s=="AcDbPolyline")//多线段
                {
                    Polyline_ poly;
                    getline(f,s);
                    getline(f,s);
                    poly.num=stoi(s);
                    getline(f,s);
                    getline(f,s);
                    if(stoi(s)==1)
                        poly.isclosed=true;
                    else
                        poly.isclosed=false;
                    getline(f,s);
                    getline(f,s);
                    for(int i=0;i<poly.num;i++)
                    {
                        POINTF pt;
                        getline(f,s);
                        if(stoi(s)==42)
                        {
                            getline(f,s);
                            poly.bulge.push_back(stof(s));
                            getline(f,s);
                        }
                        else
                        {
                            poly.bulge.push_back(0);
                        }
                        getline(f,s);
                        pt.x=stof(s);
                        getline(f,s);
                        getline(f,s);
                        pt.y=stof(s);
                        poly.center.push_back(pt);                        
                    }
                    getline(f,s);
                    if(stoi(s)==0)
                        poly.bulge.push_back(0);
                    else
                    {
                        getline(f,s);
                        poly.bulge.push_back(stof(s));
                    }
                    poly.bulge.erase(poly.bulge.begin());
                    
                    //将polylien点直接写入original_pp
                    dxf_identity dxf;
                    for(int i=0;i<poly.num;i++)
                    {
                        dxf.push_back(poly.center[i]);
                    }
                    if(poly.isclosed)
                    {
                        dxf.push_back(poly.center[0]);
                    }
                    original_pp.push_back(dxf);
                }
            }
        }
        
    }
    
}

void dxftoai::readtxt(const string& filename)
{
    ifstream f;
    f.open(filename);
    if(!f)
    {
        perror("read txt error");
    }
    dxf_identity tmp;
    string s;
    double x,y;
    while(getline(f,s))
    {
        stringstream ss(s);
        ss>>x;
        ss>>y;
        if(x==0&&y==0)
        {
            //tmp.push_back(tmp[0]);//封闭
            original_pp.push_back(tmp);
            tmp.clear();
        }
        else
        {
            tmp.push_back({x,y});
        }        
    }
    f.close();
}

void dxftoai::writeai(string filename)
{
    ofstream file;
    file.open(filename);
    file<<"%!PS-Adobe-3.0\n";
    file<<"%%BoundingBox:0 0 250 250\n";
    file<<"%AI5_FileFormat 2.0\n";
    file<<"%%BeginProlog\n";
    file<<"%%EndProlog\n";
    file<<"%%BeginSetup\n";
    file<<"%%EndSetup\n";
    file<<"%AI5_BeginLayer\n";
    file<<"1 1 1 1 0 0 1 0 0 0 Lb\n";
    file<<"(0) Ln\n";
    file<<"u\n";
    
    for(auto& ai_identity:ai_pp)
    {
        file<<"u\n";
        file<<"0 R\n";
        file<<"0 G\n";
        file<<"0.01 w\n";
        file<<"0 J 0 j 0 w []0 d\n";
        file<<ai_identity[0].pt.x<<" "<<ai_identity[0].pt.y<<" m\n";
        for(int i=1;i<ai_identity.size();i++)
        {
            if(ai_identity[i].kind==1)
            {
                file<<ai_identity[i].pt.x<<" "<<ai_identity[i].pt.y<<" L\n";
            }
            if(ai_identity[i].kind==3)
            {
                file<<ai_identity[i].pt.x<<" "<<ai_identity[i].pt.y<<" "
                <<ai_identity[i+1].pt.x<<" "<<ai_identity[i+1].pt.y<<" "
                <<ai_identity[i+2].pt.x<<" "<<ai_identity[i+2].pt.y
                <<" C\n";
                i+=2;
            }
        }
        file<<"s\n";
        file<<"U\n";
    }
    /*for(auto& ai_identity:ai_pp)//输出控制点
    {
        file<<"u\n";
        file<<"0 R\n";
        file<<"0 G\n";
        file<<"0.01 w\n";
        file<<"0 J 0 j 0 w []0 d\n";
        file<<ai_identity[0].pt.x<<" "<<ai_identity[0].pt.y<<" m\n";
        for(int i=1;i<ai_identity.size();i++)
        {
            if(ai_identity[i].kind==1)
            {
                file<<ai_identity[i].pt.x<<" "<<ai_identity[i].pt.y<<" L\n";
            }
            if(ai_identity[i].kind==3)
            {
                file<<ai_identity[i].pt.x<<" "<<ai_identity[i].pt.y<<" L\n";
                file<<ai_identity[i+1].pt.x<<" "<<ai_identity[i+1].pt.y<<" L\n";
                file<<ai_identity[i+2].pt.x<<" "<<ai_identity[i+2].pt.y<<" L\n";
                i+=2;
            }
        }
        file<<"s\n";
        file<<"U\n";
    }*/
    /*for(auto& ai_identity:simple_pp)
    {
        file<<"u\n";
        file<<"0 R\n";
        file<<"0 G\n";
        file<<"0.01 w\n";
        file<<"0 J 0 j 0 w []0 d\n";
        file<<ai_identity[0].x<<" "<<ai_identity[0].y<<" m\n";
        for(int i=1;i<ai_identity.size();i++)
        {
            file<<ai_identity[i].x<<" "<<ai_identity[i].y<<" L\n";
        }
        file<<"s\n";
        file<<"U\n";
    }*/
    file<<"U\n";
    file<<"LB\n";
    file<<"%AI5_EndLayer--\n";
    file<<"%%PageTrailer\n";
    file<<"gsave annotatepage grestore showpage\n";
    file<<"%%Trailer\n";
    file<<"%%EOF\n";
    file.close();
}

void dxftoai::simplise()
{    
    for(auto& identity:original_pp)//单个实体简化
    {
        dxf_identity tmp;
        POINTF P1=identity[0];
        tmp.push_back(P1);
        POINTF P2=identity[1];
        int i=2;
        POINTF T;
        while(i<identity.size())
        {
            T=identity[i];
            int res=PointToLine(P1,P2,T);
            if(res==1)
                P2=T;
            else
            {
                tmp.push_back(P2);
                P1=P2;
                P2=T;
            }
            i++;            
        }
        tmp.push_back(T);
        simple_pp.push_back(tmp);
    }
}


void dxftoai::handleper()
{
    for(auto& simidentity:simple_pp)
    {
        ai_identity ai;
        POINTF P1=simidentity[0];
        ai.push_back({0,P1});//start
        POINTF P2=simidentity[1];
        bool flag=false;
        int i=2;
        POINTF T;
        vector<POINTF> curvepp;
        while(i<simidentity.size())
        {
            T=simidentity[i];
            int res=PointToLine(P1,P2,T);
            i++;
            if(res==2)//曲线集，保存进行拟合
            {
                if(!flag)
                {
                    flag=true;
                    curvepp.push_back(P1);
                    curvepp.push_back(P2);
                    curvepp.push_back(T);
                }
                else
                {
                    curvepp.push_back(T);
                }                
            }
            else
            {
                if(flag)//前面点集为曲线，进行下一步
                {
                    if(curvepp.size()>3)//转成贝塞尔
                    {
                        vector<ai_besier> tmp=depart(curvepp);
                        //将tmp加到ai中

                        for(auto pt:tmp)
                        {
                            ai.push_back(pt);
                        }

                        curvepp.clear();//进行下一轮
                    }
                    else//直接作为点
                    {
                        for(int r=1;r<curvepp.size();r++)//第一个点已经在上一次保存，不用重复保存
                        {
                            ai.push_back({1,curvepp[r]});
                        }
                        curvepp.clear();
                    }
                    flag=false;
                }
                else
                {
                    //将P2直线点
                    ai.push_back({1,P2});
                }
                
            }
            P1=P2;
            P2=T;
            
        }
        if(!flag)//如果最后是直线，最后一点直接加
            ai.push_back({1,T});
        else//结尾是曲线点集，最后求一次控制点
        {
            vector<ai_besier> tmp=depart(curvepp);
            //将tmp加到ai中
            for(auto pt:tmp)
            {
                ai.push_back(pt);
            }
            curvepp.clear();//进行下一轮
        }
        ai_pp.push_back(ai);
    }

}

vector<ai_besier> dxftoai::depart(vector<POINTF>& curve)//分割点集，单独求曲线
{
    POINTF T1=curve[0];
    POINTF T2=curve[1];
    POINTF P1=curve[1];
    POINTF P2;
    vector<int> index;//分割点
    index.push_back(0);
    int i=2;
    double theta;
    double the=0;
    while(i<curve.size())
    {
        P2=curve[i];
        theta=acos(((P2.x-P1.x)*(T2.x-T1.x)+(P2.y-P1.y)*(T2.y-T1.y))/
        (sqrt(pow(P2.x-P1.x,2)+pow(P2.y-P1.y,2))*sqrt(pow(T2.x-T1.x,2)+pow(T2.y-T1.y,2))));
        if(theta>pi/2||(theta<the&&(i-*(index.end()-1)>5)))//拐角大于pi/2或者出现角度减小
        {
            index.push_back(i);
            the=0;
            T1=curve[i-1];
            T2=curve[i];
            P1=curve[i];
        }
        else
        {
            the=theta;
            P1=P2;
        }
        i++;        
    }
    int size=index.size();
    int last=index[size-1];
    //最后一段需要整理
    int h=curve.size()-1;
    int left=h-last;
    if(left>=0&&left<4)
    {
        index[size-1]=h;//直接加到最后一段
    }
    else
    {
        index.push_back(h);//直接自成一段
    }

    vector<ai_besier> res;
    for(int i=1;i<index.size();i++)
    {
        vector<POINTF> dots;
        for(int j=index[i-1];j<=index[i];j++)
            dots.push_back(curve[j]);
        vector<ai_besier> tmp=triangle(dots);//返回后面三个控制点
        for(auto pt:tmp)
            res.push_back(pt);
    }
    return res;

}

vector<ai_besier> dxftoai::triangle(vector<POINTF>& dots)
{
    int n=dots.size()-1;
    double theta1=atan((dots[1].y-dots[0].y)/(dots[1].x-dots[0].x));
    double theta2=atan((dots[2].y-dots[1].y)/(dots[2].x-dots[1].x));
    if(abs(theta1-theta2)>pi/2)
        theta1=pi+pi+theta1;
    double theta3=1.5*theta1-0.5*theta2;
    double theta4=atan((dots[n-1].y-dots[n-2].y)/(dots[n-1].x-dots[n-2].x));
    double theta5=atan((dots[n].y-dots[n-1].y)/(dots[n].x-dots[n-1].x));
    if(abs(theta4-theta5)>pi/2)
        theta4=pi+theta4;
    double theta6=1.5*theta5-0.5*theta4;
    double x1=dots[0].x; double y1=dots[0].y;
    double x2=dots[n].x; double y2=dots[n].y;
    double k1=tan(theta3); double k2=tan(theta6);
    //求两条直线交点
    double x=(y2-y1+k1*x1-k2*x2)/(k1-k2);
    double y=k1*(x-x1)+y1;
    POINTF S=dots[0];
    POINTF Q={x,y};
    POINTF E=dots[n];
    POINTF QS=Q-S;
    POINTF QE=Q-E;
    vector<POINTF> res(4);
    res[0]=S; res[3]=E;
    double left=0; double right=1;
    double ll,rr,dis_all; int pos;
    double LL,RR,MM;
    res[1]=S+QS*left;
    res[2]=E+QE*left;
    adjust(res,dots,Q,ll,rr,pos,dis_all);
    LL=dis_all;
    res[1]=S+QS*right;
    res[2]=E+QE*right;
    adjust(res,dots,Q,ll,rr,pos,dis_all);
    RR=dis_all;
    MM=RR;
    int times=0;
    double mid;
    while(abs(MM)>0.001&&times<20)//二分法，几何和为0
    {
        times++;
        mid=(left+right)/2;
        res[1]=S+QS*mid;
        res[2]=E+QE*mid;
        adjust(res,dots,Q,ll,rr,pos,dis_all);
        MM=dis_all;
        if(MM*LL<0)
        {
            right=mid;
            RR=MM;
        }
        else
        {
            left=mid;
            LL=MM;
        }
    }
    vector<ai_besier> result;
    if(times!=20)
    {
        if((abs(ll)>0.01||abs(rr)>0.01)&&dots.size()>=7&&pos!=0)//分两段再求
        {
            vector<POINTF> tmp;
            for(int i=0;i<=pos;i++)
            {
                tmp.push_back(dots[i]);
            }
            result=triangle(tmp);
            tmp.clear();
            //后部分再加到后面
            for(int i=pos;i<=n;i++)
            {
                tmp.push_back(dots[i]);
            }
            vector<ai_besier> next=triangle(tmp);
            for(int i=0;i<next.size();i++)
            {
                result.push_back(next[i]);
            }
        }
        else
        {
            for(int i=1;i<4;i++)
                result.push_back({3,res[i]});
        }
        
    }
    else
    {
        //直接矩阵求取
        result=matrix_algorithm(dots);
    }
    return result;
    
}

//K：控制点 A：点集 Q顶点 ll:左边误差 rr:右边误差 pos:中间位置 dis_all：误差和(+-)
void dxftoai::adjust(vector<POINTF>& K,vector<POINTF>& A,POINTF& Q,
                double& ll,double& rr,int& pos,double& dis_all)//输出
{
    ll=0; rr=0; pos=0; dis_all=0;
    double t=0;
    vector<POINTF> PP;
    vector<POINTF> P;
    double dis=0;
    double tmp;
    int j=0;//A的点
    vector<vector<double>> parameter(1);
    vector<double> distribute;
    for(int i=0;i<=1000;i++)
    {
        t=i*0.001;
        parameter[0]={pow(1-t,3),3*t*pow((1-t),2),3*pow(t,2)*(1-t),pow(t,3)};
        P=MatrixMulti(parameter,K);
        tmp=sqrt(pow(A[j].x-P[0].x,2)+pow(A[j].y-P[0].y,2));//离得最近的点
        if(dis<=tmp)
        {
            PP.push_back(P[0]);
            //直线公式，看Q和P在直线的位置
            double a=A[j+1].y-A[j].y;
            double b=A[j].x-A[j+1].x;
            double c=A[j+1].x*A[j].y-A[j].x*A[j+1].y;

            double q=a*Q.x+b*Q.y+c;
            double f=a*P[0].x+b*P[0].y+c;
            if(q*f>0)
            {
                distribute.push_back(dis);
                ll+=dis;
            }
            else
            {
                dis=-dis;
                distribute.push_back(dis);
                rr+=dis;
            }
            dis_all+=dis;
            j++;
            dis=sqrt(pow(A[j].x-P[0].x,2)+pow(A[j].y-P[0].y,2));//下一个点的距离了
            if(j==A.size()-1)
                break;         
        }
        else
        {
            dis=tmp;
        }   

    }
    PP.push_back(A[A.size()-1]);
    for(int i=3;i<A.size()-3;i++)//分割点必然是左右必须有四个，否则不进行分割
    {
        if(abs(distribute[i])<abs(distribute[i-1])&&abs(distribute[i])<abs(distribute[i+1]))
        {
            pos=i;
            break;
        }            
    }
    double corr=0;
    for(int i=0;i<A.size()/2;i++)
    {
        corr+=distribute[i];
    }
    if(corr<0)
    {
        tmp=ll;
        ll=rr;
        rr=tmp;
    }
}


vector<ai_besier> dxftoai::matrix_algorithm(vector<POINTF>& dots)
{
    int n=dots.size();
    double d_all=0;
    vector<double> d;
    for(int i=1;i<n;i++)
    {
        double tmp=sqrt(pow(dots[i].x-dots[i-1].x,2)+pow(dots[i].y-dots[i-1].y,2));
        d_all+=tmp;
        d.push_back(tmp);
    }
    double t=0;
    vector<vector<double>> CA;
    CA.push_back({1,0,0,0});
    for(int i=0;i<d.size();i++)
    {
        t+=d[i]/d_all;
        CA.push_back({pow(1-t,3),3*t*pow(1-t,2),3*pow(t,2)*(1-t),pow(t,3)});
    }
    vector<vector<double>> CA_T=MatrixTrans(CA);
    vector<vector<double>> CACA=MatrixMulti(CA_T,CA);
    vector<vector<double>> res=CACA;
    bool flag=MatrixInver(CACA,res);
    vector<ai_besier> result;
    if(!flag)
    {
        perror("not inverse");
        return result;
    }
    res=MatrixMulti(res,CA_T);
    vector<POINTF> dot_p=MatrixMulti(res,dots);
    dot_p[0]=dots[0];
    dot_p[3]=dots[n-1];//使其贴合端点
    
    for(int i=1;i<dot_p.size();i++)
    {
        result.push_back({3,dot_p[i]});//第一个点不要
    }
    return result;    
    
}


int dxftoai::PointToLine(const POINTF& P1,const POINTF& P2,const POINTF& T)
{
    double A=P2.y-P1.y;
    double B=P1.x-P2.x;
    double C=P2.x*P1.y-P1.x*P2.y;
    double d_TtoL=abs(A*T.x+B*T.y+C)/sqrt(A*A+B*B);//距离公式
    double d_P1P2=sqrt(pow(P1.x-P2.x,2)+pow(P1.y-P2.y,2));//P1P2 distance
    double d_P2T=sqrt(pow(T.x-P2.x,2)+pow(T.y-P2.y,2));
    //P1P2与P2T向量夹角cos
    double theta=((P2.x-P1.x)*(T.x-P2.x)+(P2.y-P1.y)*(T.y-P2.y))/(d_P1P2*d_P2T);

    if(d_TtoL<d_P1P2/50&&theta>0.98)//直线
        return 1;
    if(theta>0.9&&d_P1P2/d_P2T<2&&d_P1P2/d_P2T>0.5)//曲线
        return 2;
    return 0;//都不是
}
