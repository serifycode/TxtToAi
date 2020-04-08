#ifndef ENTITY_H
#define ENTITY_H

#include "dxftoai.h"
#include <vector>


using namespace std;

struct Line_
{
    POINTF start;
    POINTF end;
};

struct Circle_
{
    POINTF center;
    double radius;
};

struct Arc_
{
    POINTF center;
    double radius;
    double start_angle;
    double end_angle;
};

struct Polyline_
{
    int num;
    bool isclosed;
    vector<double> bulge;
    vector<POINTF> center;
};

struct Ellipse_
{
    POINTF center;
    double majoraxis;
    double minoraixs;
};

struct Text_
{
    char* text;
    char* content;
};


#endif