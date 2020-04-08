#include "dxftoai.h"
#include <iostream>
#include <string>

using namespace std;

int main(int argc,char *argv[])//
{
    string readname="LUC-007.dxf";
    //string readname="facep.txt";
    string writename="LUC-007.ai";
    if(argc==3)
    {
        readname=argv[1];
        writename=argv[2];
    }
    dxftoai func;
    string::size_type pos;
    if((pos=readname.find("."))==string::npos)
    {
        perror("read error");
        return -1;
    }
        
    if(readname.substr(pos+1)=="txt")
        func.readtxt(readname);
    if(readname.substr(pos+1)=="dxf")
        func.readdxf(readname);
    func.simplise();
    func.handleper();
    func.writeai(writename);
    return 0;
}


/*
示例：
可以直接在里面修改参数：
readname：
输入的文件，可以使txt，每个实体以0 0"结束；或者是dxf，可包含圆，和多线段
writename：
输出文件为，为矢量图ai。

也可debug一次，后面直接调用exe文件完成输出，第一个参数为输入文件，第二个参数为输出文件
*/