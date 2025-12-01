#ifndef CUT_HPP
#define CUT_HPP

#include "Domain.hpp"
#include "CutDomain.hpp"
#include "PDPParameters.hpp"
#include "PDPDistanceMatrix.hpp"
#include "Segment.hpp"
#include "SegmentComparator.hpp"

const int MAXLEN = PDPParameters::MAXLEN;
const int MAXSIZE = PDPParameters::MAXSIZE;
const int RG = PDPParameters::RG;
const int TD = PDPParameters::TD;
const int ENDS =PDPParameters::ENDS;
const int ENDSEND = PDPParameters::ENDSEND;

class Cut {
private:
    
public:
    int cut(std::vector<Atom>& ca,
            Domain& dom,
            CutValues& val,
            std::vector<std::vector<int>>& dist,
            PDPDistanceMatrix& pdpMatrix);
    };
#endif
