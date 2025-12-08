#ifndef CLUSTER_DOMAINS_HPP
#define CLUSTER_DOMAINS_HPP

#include <vector>
#include "Domain.hpp"
#include "PDPDistanceMatrix.hpp"

class ClusterDomains {
public:
  ClusterDomains();
  static int getTotalContacts(std::vector<Domain> &domains,PDPDistanceMatrix& pdpDistMatrix, Domain& d1,Domain& d2);
  static int isContacting(Domain& d1,Domain& d2, const std::vector<int>& iclose, const std::vector<int>& jclose, int nclose);
  //  static bool isContacting( Domain& d1,Domain& d2, std::vector<int> &iclose,std::vector<int> &jclose,int nclose);
  static std::vector<std::vector<Domain>> cluster(std::vector<Domain>& domains,PDPDistanceMatrix& pdpDistMatrix);
  static int ndom;
  static std::vector<int> visibleDomains;    
private:
  
};

#endif // CLUSTER_DOMAINS_HPP

