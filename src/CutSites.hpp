#ifndef CUTSITES_HPP
#define CUTSITES_HPP

#include <vector>
#include <string>
#include <algorithm>

#include "PDPParameters.hpp"

class CutSites {
private:
    int ncuts;
    CutSites& operator=(const CutSites& other);
public:
  std::vector<int> cut_sites;
  CutSites();
  ~CutSites();
  CutSites(const CutSites& other);
  std::string toString();
  int getNcuts();
  void setNcuts(int ncuts);
  void addNcuts(int add); 
  std::vector<int> getCutSites();
  void setCutSites(std::vector<int> cut_sites);
  void pushbackCutSites(int site);
  void addCutSites(int index, int site);
  
};

#endif // CUTSITES_HPP

