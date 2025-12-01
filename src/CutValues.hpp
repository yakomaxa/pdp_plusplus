#ifndef PDP_CUTVALUES_H
#define PDP_CUTVALUES_H

#include<string>

class CutValues {
public:
    double s_min;
    int site2;
    bool first_cut;

    double AD;

    CutValues();
    ~CutValues();
    std::string toString();

};

#endif
