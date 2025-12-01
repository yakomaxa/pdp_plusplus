#include "CutValues.hpp"

CutValues::CutValues() : s_min(100), site2(0), first_cut(true) {}

CutValues::~CutValues() {}

std::string CutValues::toString() {
    return "CutValues [s_min=" + std::to_string(s_min) + ", site2=" + std::to_string(site2) +
           ", AD=" + std::to_string(AD) + "]";
}

