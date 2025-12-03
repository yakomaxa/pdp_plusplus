#include "CutSites.hpp"
#include "PDPParameters.hpp"
#include <iostream>

// Default constructor
CutSites::CutSites() : cut_sites(PDPParameters::MAX_CUTS){
  this->ncuts=0;
}

CutSites::CutSites(CutSites const &cs): cut_sites(PDPParameters::MAX_CUTS){
  this->ncuts=0;
}

// Destructor
CutSites::~CutSites() {
    // Nothing to do here
}

// Getter for ncuts
int CutSites::getNcuts(){
    return ncuts;
}

// Setter for ncuts
void CutSites::setNcuts(int ncuts) {
    this->ncuts = ncuts;
}

void CutSites::addNcuts(int add) {
    this->ncuts += add;
}

// Getter for cut_sites
std::vector<int> CutSites::getCutSites() {
    return cut_sites;
}

void CutSites::addCutSites(int index, int site){
    this->cut_sites[index] = site;
};

void CutSites::pushbackCutSites(int site){
  this->cut_sites.push_back(site);
};


// Setter for cut_sites

void CutSites::setCutSites(std::vector<int> cut_sites) {
   this->cut_sites = cut_sites;    
}


// toString method
std::string CutSites::toString() {
    std::string result = "CutSites [ncuts=" + std::to_string(ncuts) + ", cut_sites=[";
    for (int i = 0; i < ncuts; i++) {
        result += std::to_string(cut_sites[i]);
        if (i < ncuts - 1) {
            result += ", ";
        }
    }
    result += "]]";
    return result;
}

