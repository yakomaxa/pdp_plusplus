#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>

#include "ClusterDomains.hpp"
#include "PDPDistanceMatrix.hpp"
#include "Domain.hpp"
#include "CutDomain.hpp"
#include "PDPParameters.hpp"

// K.S. 2025-Dec-06:
// compare the result of original PDP and PDP++ for AF-P14775-F1-model_v6.pdb
// the former returns 1–96, but the latter does not return anything, removing all the small fragments
// This clustering part of PDP++ might need debugging, or the original PDP is flawed; if the original is executed with -v, it shows many super large double (float) for the variable named avd, which suggests merging is superficial and misled by these weird scores: I think the insanely large value is something from NaN or overflow. But it's float, so it causes mis-scoring!



std::vector<int> ClusterDomains::visibleDomains;
int ClusterDomains::ndom;
ClusterDomains::ClusterDomains(){

}

static void listdomains(std::vector<Domain>& domains) {
  int i = -1;
  for (Domain& dom : domains) {
    i++;
    std::cout << "create DOMAIN" << i << ", resi ";
    std::vector<Segment>& segments = dom.getSegments();

    int flag=0;
    for (Segment& s : segments) {
      if (flag>0){
	std::cout << "+" ;
      }
      std::cout << s ;
      flag++;
    }
    std::cout << ";" << "\n";
  }
};

int calc_S(const int a1,
            const int b1,
            const int a2,
            const int b2,
            PDPDistanceMatrix & pdpDistMatrix) {

    int contacts = 0;
    std::vector<std::vector<int>> dist = pdpDistMatrix.getDist();

    for (int i = a1; i <= b1; i++) {
        for (int j = a2; j <= b2; j++) {
            contacts += dist[i][j];
        }
    }

    return contacts;
};

int ClusterDomains::getTotalContacts(std::vector<Domain>& domains,PDPDistanceMatrix& pdpDistMatrix,Domain& i,Domain& j) {
  int total_contacts = 0;
  std::vector<std::vector<int>> dist = pdpDistMatrix.getDist();
  for(int k=0; k<i.getNseg(); k++) {
    int a1=i.getSegmentAtPos(k).getFrom();
    int b1=i.getSegmentAtPos(k).getTo();      
    for(int l=0; l<j.getNseg(); l++) {
      int a2=j.getSegmentAtPos(l).getFrom();
      int b2=j.getSegmentAtPos(l).getTo();
      for (int ii = a1; ii <= b1; ii++) {
        for (int jj = a2; jj <= b2; jj++) {
	  total_contacts += dist[ii][jj];;
        }
      }
    }
  }
  return total_contacts;
};

int ClusterDomains::isContacting(Domain& i,Domain& j,const std::vector<int>& iclose, const std::vector<int>& jclose, int nclose) {
  for(int k=0; k<i.getNseg(); k++) {
    int fromi=i.getSegmentAtPos(k).getFrom();
    int toi=i.getSegmentAtPos(k).getTo();
    for(int l=0; l<j.getNseg(); l++) {
      int fromj=j.getSegmentAtPos(l).getFrom();
      int toj=j.getSegmentAtPos(l).getTo();     
      for (int n = 0 ; n < nclose ; n++){
	if (fromi <= iclose[n] && toi >= iclose[n] &&
	    fromj <= jclose[n] && toj >= jclose[n] 
	    ){
	  return true;
	}
	if (fromi <= jclose[n] && toi >= jclose[n] &&
	    fromj <= iclose[n] && toj >= iclose[n]
	    ){
	  return true;
	}
      }
    }
  }
  return false;
};
    
std::vector<Domain> combine(std::vector<Domain> &domains, int Si, int Sj, double maximum_value, std::vector<std::vector<int>> &contacts) {
  bool verbose_cluster = PDPParameters::VERBOSE;
  if (verbose_cluster){
        std::cout << "  +++  combining domains " << Si << " " << Sj << "\n";
  }
  for (int i = 0; i < domains[Sj].getNseg(); i++) {
    domains[Si].getSegmentAtPos(domains[Si].getNseg()).setFrom(domains[Sj].getSegmentAtPos(i).getFrom());
    domains[Si].getSegmentAtPos(domains[Si].getNseg()).setTo(domains[Sj].getSegmentAtPos(i).getTo());
    domains[Si].addNseg(1);  
  }

  std::erase(ClusterDomains::visibleDomains,Sj);
  
  for (int k : domains[Sj].getContacted()){
    domains[Si].removeContacted(k);
    domains[Si].pushbackContacted(k);
    domains[k].removeContacted(Si);
    domains[k].pushbackContacted(Si);

    domains[k].removeContacted(Sj);
    domains[Sj].removeContacted(k);
    
    contacts[Si][k] += contacts[Sj][k];
    contacts[k][Si] += contacts[k][Sj];
  }
  domains[Si].removeContacted(Sj);
  domains[Sj].removeContacted(Si);
  domains[Si].addSize(domains[Sj].getSize());

  ClusterDomains::ndom -=1;
  return std::move(domains);
};


std::vector<std::vector<Domain>> ClusterDomains::cluster(
                            std::vector<Domain>& domains,
                            PDPDistanceMatrix& pdpDistMatrix){

  bool verbose_cluster = PDPParameters::VERBOSE;
  ClusterDomains::ndom = (int)domains.size();
  int Si = -1;
  int Sj = -1;
  int Sis = -1;
  int Sjs = -1;
  int Sim = -1;
  int Sjm = -1;
  
  int total_max_contacts = 0;
  
  double maximum_values = PDPParameters::CUT_OFF_VALUE1S;
  double maximum_valuem = PDPParameters::CUT_OFF_VALUE1M;
  double maximum_value  = PDPParameters::CUT_OFF_VALUE1;

  std::vector<std::vector<Domain>> h_domain;
  if (ClusterDomains::ndom < 2){
    h_domain.push_back(domains);
    return h_domain;
  }
  
  Domain d1;
  Domain d2;
  std::vector<int> i_can_contact(ClusterDomains::ndom*ClusterDomains::ndom);
  std::vector<int> j_can_contact(ClusterDomains::ndom*ClusterDomains::ndom);
  int n_can=0;
  
  std::vector<int> icontacted(ClusterDomains::ndom*ClusterDomains::ndom);
  std::vector<int> jcontacted(ClusterDomains::ndom*ClusterDomains::ndom);
  int i,j;
  
  int nclose_raw = pdpDistMatrix.getNclose_raw();
  std::vector<int> iclose_raw = pdpDistMatrix.getIclose_raw();
  std::vector<int> jclose_raw = pdpDistMatrix.getJclose_raw();

  for(int i=0;i<ClusterDomains::ndom;i++) {
    d1 = domains.at(i);
    for(int j=i ;j<ClusterDomains::ndom;j++) {
      d2 = domains.at(j);
      if (ClusterDomains::isContacting(d1,d2,iclose_raw,jclose_raw,nclose_raw)){
	i_can_contact[n_can]=i;
	j_can_contact[n_can]=j;	    
	n_can++;	
      }
    }
  }


  std::vector<std::vector<int>> contacts_list((ClusterDomains::ndom), std::vector<int>(ClusterDomains::ndom));
  std::vector<std::vector<int>> dist = pdpDistMatrix.getDist();
  //int tmp=0;
  int total_contacts = 0;
  for(int n = 0 ; n < n_can; n++) {
    total_contacts = 0 ;
    i = i_can_contact[n];
    j = j_can_contact[n];
    
    //    if (i >= j){
    //      tmp = i;
    //      i = j;
    //      j = tmp;	
    //    }
    d1 = domains.at(i);    
    d2 = domains.at(j);


    for(int k=0; k<d1.getNseg(); k++) {
      int a1=d1.getSegmentAtPos(k).getFrom();
      int b1=d1.getSegmentAtPos(k).getTo();      
      for(int l=0; l<d2.getNseg(); l++) {
	int a2=d2.getSegmentAtPos(l).getFrom();
	int b2=d2.getSegmentAtPos(l).getTo();
	for (int ii = a1; ii <= b1; ii++) {
	  for (int jj = a2; jj <= b2; jj++) {
	    total_contacts += dist[ii][jj];;
	  }
	}
      }
    }

    if (verbose_cluster){
      std::cout << " pos: d1:" << i << " vs d2:" << j << " d1:" << d1.getSegmentAtPos(0).getFrom() << "-" << d1.getSegmentAtPos(0).getTo() << " " <<  d2.getSegmentAtPos(0).getFrom() << "-" << d2.getSegmentAtPos(0).getTo() << " " << total_contacts << "\n" ;
    }
    
    contacts_list[i][j]=contacts_list[j][i]=total_contacts;
     if (total_contacts > 0){
      domains[i].pushbackContacted(j);
      domains[j].pushbackContacted(i);
     }
  }
  
  for (int i = 0 ; i < (int)domains.size() ;i++){
    ClusterDomains::visibleDomains.push_back(i);
  }
  int prev_ndom = ClusterDomains::visibleDomains.size();


  //  std::vector<Domain> olddomains = domains; 
  do {
    //    std::vector<std::vector<int>> contacts;    
    for(int i : ClusterDomains::visibleDomains){      
      for (int j : domains[i].getContacted()){
       	if (j<=i){
	  continue;
	}
	int total_contacts = contacts_list[i][j];
	int size1dom1=domains[i].getSize();
	int size2dom2=domains[j].getSize();
	double minDomSize=std::min(size1dom1,size2dom2);
	double maxDomSize=std::max(size1dom1,size2dom2);
	// set some limits on how big the domains can get
	if(minDomSize>150&&maxDomSize>1.5*minDomSize){
	  maxDomSize=1.5*minDomSize;
	}else if(maxDomSize>2*minDomSize){
	  maxDomSize=2*minDomSize;
	}
	int size1= std::min(PDPParameters::MAXSIZE,(int)minDomSize);
	int size2= std::min(PDPParameters::MAXSIZE,(int)maxDomSize);
	minDomSize=std::min(pow(minDomSize,1.6/3)+PDPParameters::RG1,pow(minDomSize,1.4/3)+pow(PDPParameters::TD1,1.6/3)+PDPParameters::RG1);
	maxDomSize=std::min(pow(maxDomSize,1.6/3)+PDPParameters::RG1,pow(maxDomSize,1.4/3)+pow(PDPParameters::TD1,1.6/3)+PDPParameters::RG1);
	
	total_max_contacts=(int)(minDomSize*maxDomSize*10);
	if(size1>130){
	  total_max_contacts=(int)(minDomSize*maxDomSize*9);
	}
	
	double S_value= double(total_contacts)/(double)total_max_contacts;
		
      if(verbose_cluster) {
	printf("size1=%i size2=%i minDomSize=%f maxDomSize=%f total_contacts = %i \n", size1,size2,minDomSize,maxDomSize,total_contacts);
	printf(" total_contacts = %i total_max_contacts = %i\n", total_contacts, total_max_contacts);
	printf(" maximum_value = %f S_value = %f\n",maximum_value, S_value);
      }
      
      if (S_value  > maximum_value) {
	maximum_value = S_value;
	Si = i;
	Sj = j;
      };
      
      if (S_value > maximum_valuem && size1<70) {
	maximum_valuem = S_value;
	Sim = i;
	Sjm = j;
      };
      
      if (S_value > maximum_values && size1<52) {
	maximum_values = S_value;
	Sis = i;
	Sjs = j;
      };
      
      total_contacts = 0;
      
      total_max_contacts = 0;
      
      }
    }

    if (maximum_value > PDPParameters::CUT_OFF_VALUE1) {

      domains = combine(domains,Si, Sj, maximum_value,contacts_list);

      maximum_value = PDPParameters::CUT_OFF_VALUE1-.1;
      maximum_values = PDPParameters::CUT_OFF_VALUE1S-.1;
      maximum_valuem = PDPParameters::CUT_OFF_VALUE1M-.1;
      
    }else if (maximum_valuem > PDPParameters::CUT_OFF_VALUE1M) {

      domains = combine(domains, Sim, Sjm, maximum_valuem,contacts_list);

      maximum_value =  PDPParameters::CUT_OFF_VALUE1-.1;
      maximum_values = PDPParameters::CUT_OFF_VALUE1S-.1;
      maximum_valuem = PDPParameters::CUT_OFF_VALUE1M-.1;
    }else if (maximum_values > PDPParameters::CUT_OFF_VALUE1S) {

      domains = combine(domains, Sis, Sjs, maximum_values,contacts_list);

      maximum_value = PDPParameters::CUT_OFF_VALUE1-.1;
      maximum_values = PDPParameters::CUT_OFF_VALUE1S-.1;
      maximum_valuem = PDPParameters::CUT_OFF_VALUE1M-.1;
      
    }else {
      
      if (visibleDomains.size() < prev_ndom){
	std::vector<Domain> newdoms;
	for (auto i : ClusterDomains::visibleDomains){
	  newdoms.push_back(domains[i]);
	}
	h_domain.push_back(newdoms);
      }
      prev_ndom = visibleDomains.size();
      PDPParameters::CUT_OFF_VALUE1=PDPParameters::CUT_OFF_VALUE1*0.99;
      PDPParameters::CUT_OFF_VALUE1S=PDPParameters::CUT_OFF_VALUE1S*0.99;
      PDPParameters::CUT_OFF_VALUE1M=PDPParameters::CUT_OFF_VALUE1M*0.99;
      printf("PREV_NDOM %i %f %f %f\n",prev_ndom, PDPParameters::CUT_OFF_VALUE1, PDPParameters::CUT_OFF_VALUE1S, PDPParameters::CUT_OFF_VALUE1M);
      if (prev_ndom < 2){
	break;
      }
    }
    //    listdomains(domains);
    //  } while (maximum_value > 0.0 || maximum_values > 0.0 || maximum_valuem > 0.0);
  } while (visibleDomains.size()>=1);

  return h_domain;
};


