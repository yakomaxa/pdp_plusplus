#ifndef PDPPARAMETERS_H
#define PDPPARAMETERS_H

class PDPParameters {
public:

    static int maxIndex ;
  
    static const int MAXSIZE;
    static int MAXLEN ;
    static const int MAXDOM ;
    static const int MAX_CUTS ;
    static const int MIN_DOMAIN_LENGTH;
    static const int ENDS;
    static const int ENDSEND;
    static const float RG1;
    static const float RG;
    static const float TD1;
    static const float TD;
    static const float DBL;
    static const int MAXCONT;

    static const float CUT_OFF_VALUE; /* decide to cut */
    static float CUT_OFF_VALUE1; /* decide to combine */
    static const float CUT_OFF_VALUE2; /* decide to double cut */
    static float CUT_OFF_VALUE1S; /* decide to combine small domains */
    static float CUT_OFF_VALUE1M; /* decide to combine medium domains */

    static bool VERBOSE;


  void setMAXLEN(int maxlen);
};

#endif /* PDPPARAMETERS_H */

