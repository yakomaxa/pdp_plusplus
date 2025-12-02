#ifndef ATOM_HPP
#define ATOM_HPP

#include <stdio.h>
#include <string>

class Atom {
private:
  float x, y, z;
  std::string chain;
  std::string residue;
  int indexorg;
  int chainid;
    
public:
    Atom();
    ~Atom();
    void setX(float x);
    void setY(float y);
    void setZ(float z);
    float getX() ;
    float getY() ;
    float getZ() ;
  
    void setChain(std::string chain) ;
    std::string getChain() ;

    void setResidue(std::string residue) ;
    std::string getResidue() ;

    void setChainId(int chainid) ;
    int getChainId() ;

    void setIndexOrg(int resi) ;
    int getIndexOrg() ;
};

#endif

