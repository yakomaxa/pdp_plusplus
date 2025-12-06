#include "PDPParameters.hpp"


int PDPParameters::maxIndex=0;

const int PDPParameters::MAXSIZE = 350;
int PDPParameters::MAXLEN = 20000;
const int PDPParameters::MAXDOM = 1000;
const int PDPParameters::MAX_CUTS = 1000;
const int PDPParameters::MIN_DOMAIN_LENGTH = 35;
const int PDPParameters::ENDS = 12 ;
const int PDPParameters::ENDSEND = 9;
const float PDPParameters::RG1 = 1.0f;
const float PDPParameters::RG = 0.0f;
const float PDPParameters::TD1 = 40.f;
const float PDPParameters::TD = 25.f;
const float PDPParameters::DBL = .05f;
const float PDPParameters::CUT_OFF_VALUE = .50f;
const float PDPParameters::CUT_OFF_VALUE1 = .29f;
const float PDPParameters::CUT_OFF_VALUE2 = .44f;
const float PDPParameters::CUT_OFF_VALUE1S = .19f;
const float PDPParameters::CUT_OFF_VALUE1M = .21f;
bool PDPParameters::VERBOSE = false;

void PDPParameters::setMAXLEN(int MAXLEN){
  this->MAXLEN = MAXLEN;
}
