#include <fstream>

#include "Structure.hpp"
#include "Atom.hpp"
//#include "Chain.hpp"
#include "GetDistanceMatrix.hpp"
#include "PDPDistanceMatrix.hpp"
#include "CutSites.hpp"
#include "CutDomain.hpp"
#include "ClusterDomains.hpp"
#include "ShortSegmentRemover.hpp"
#include "Domain.hpp"
#include "PDPParameters.hpp"
//#include "LocalDomainParser.hpp"

// This function writes to standard output
static void listdomains(std::vector<Domain>& domains) {
  int i = -1;
  for (Domain& dom : domains) {
    i++;
    std::cout << "create DOMAIN" << i << ", ";
    //    std::vector<Segment>& segments = dom.getSegments();

    int flag=0;
    for (int i=0 ; i< dom.getNseg();i++) {
      if (flag>0){
	std::cout << "+" ;
      }
      std::cout << dom.getSegmentAtPos(i) ;
      flag++;
    }
    std::cout << ";" << "\n";
  }
};

// This function writes to a file
static void listdomains(std::vector<Domain>& domains, const std::string& filename) {
  std::ofstream output_file(filename);

  if (!output_file) {
    std::cerr << "Cannot open the output file: " << filename << "\n";
    return;
  }

  int i = -1;
  for (Domain& dom : domains) {
    i++;
    output_file << "create DOMAIN" << i << ", ";

    int flag=0;
    for (int i=0 ; i< dom.getNseg();i++) {
      if (flag>0){
	output_file << "+" ;
      }
      output_file << dom.getSegmentAtPos(i) ;
      flag++;
    }
    output_file << ";" << "\n";
  }
  output_file << "set grid_mode,1" << "\n";
  output_file << "set grid_slot,2,DOM*" << "\n";
  
  output_file.close();
};

int main(int argc, char *argv[]){
  std::string filename=argv[1];

  PDPParameters::VERBOSE=false;
  if (argc == 3) {
      std::string ch=argv[2];
      if (ch == "-v"){
	PDPParameters::VERBOSE=true;
      }
  }
  if (argc == 4){
    std::string verbose=argv[3];
    if (verbose == "-v"){
      PDPParameters::VERBOSE=true;
    }
  }
  
  if (PDPParameters::VERBOSE){
    printf("---------Reading structure\n");
  }
  Structure s = Structure(filename);
  if (PDPParameters::VERBOSE){
    printf("---------Reading structure Done\n");
  }
  std::vector<Domain> domains;
  PDPParameters param;
  param.setMAXLEN(s.numResidues);
  if (PDPParameters::VERBOSE){
    printf("---------Get Repr atoms\n");
  }
  std::vector<Atom> ca = s.getRepresentativeAtomArray();
  if (PDPParameters::VERBOSE){
    printf("---------Get Repr atoms Done\n");
  }
  GetDistanceMatrix distMtxCalculator;

  if (PDPParameters::VERBOSE){
    printf("---------distMat creation\n");
  }
  PDPDistanceMatrix pdpMatrix = distMtxCalculator.getDistanceMatrix(ca);
  if (PDPParameters::VERBOSE){
    printf("---------distMat creation Done\n");
  }

  Domain dom;
  //Chain c = ca[0].getGroup().getChain();
  //dom.setId("D"+c.getStructure().getPDBCode()+c.getId()+"1");
  if (PDPParameters::VERBOSE){
    printf("---------Setting domain info\n");
  }
  dom.setId("testDomain");
  dom.setSize((int)ca.size());
  dom.setNseg(1);

  dom.getSegmentAtPos(0).setFrom(0);
  dom.getSegmentAtPos(0).setTo(int(ca.size())-1);
  CutSites cutSites = CutSites();
  if (PDPParameters::VERBOSE){
    printf("---------Setting domain info done\n");
    printf("---------Initial splitting\n");
  }
  // Do the initial splitting
  

  std::vector<int> init_cutsites;
  // add the head residue of the chains except for the first one
  s.tailofchain.pop_back();;
  for (int i : s.tailofchain){
    if (PDPParameters::VERBOSE){
      printf("ADDING INITAIAL SITE%i\n",i);
    }
    init_cutsites.push_back(i+1);
  }
  CutValues val;
  CutDomain cutDomain(ca,pdpMatrix, init_cutsites);
  cutDomain.cutDomain(dom, cutSites,pdpMatrix, val);
  if (PDPParameters::VERBOSE){
    printf("---------Initial splitting done\n");
  }
  domains =  cutDomain.getDomains();  
  // Cluster domains
  if (PDPParameters::VERBOSE){
    printf("---------Clustering domains\n");
  }
  domains = ClusterDomains::cluster(domains, pdpMatrix);
  if (PDPParameters::VERBOSE){
    printf("---------Clustering domains Done\n");
  }

  for (int i= 0 ; i < (int)domains.size(); i++){
    for (int j = 0 ; j < domains[i].getNseg();j++){
      domains[i].getSegmentAtPos(j).setFromOrg(ca[domains[i].getSegmentAtPos(j).getFrom()].getIndexOrg());
      domains[i].getSegmentAtPos(j).setToOrg(ca[domains[i].getSegmentAtPos(j).getTo()].getIndexOrg());
      domains[i].getSegmentAtPos(j).setChain(ca[domains[i].getSegmentAtPos(j).getTo()].getChain());
      domains[i].getSegmentAtPos(j).setChainId(ca[domains[i].getSegmentAtPos(j).getTo()].getChainId());      
    }
  }

  listdomains(domains);

  for(int j=0;j<(int)domains.size();j++) {
    std::sort(domains[j].getSegments().begin(),
    	      domains[j].getSegments().end(),SegmentComparator());
    for (int i=0;i<domains[j].getNseg()-1;i++){
      std::cout << domains[j].getSegmentAtPos(i) << std::endl;
      
      if((domains[j].getSegmentAtPos(i).getToOrg())==(domains[j].getSegmentAtPos(i+1).getFromOrg()-1)
	 && (domains[j].getSegmentAtPos(i).getChain()==(domains[j].getSegmentAtPos(i+1).getChain()))){
	
       	domains[j].getSegmentAtPos(i).setToOrg(domains[j].getSegmentAtPos(i+1).getToOrg());
	domains[j].getSegmentAtPos(i).setTo(domains[j].getSegmentAtPos(i+1).getTo());
	domains[j].addNseg(-1);
	if (PDPParameters::VERBOSE){
	  printf("NSEG=%i\n",domains[j].getNseg());
	}
	for(int l=i+1;l<domains[j].getNseg();l++){
	  domains[j].getSegmentAtPos(l).setToOrg(domains[j].getSegmentAtPos(l+1).getToOrg());
	  domains[j].getSegmentAtPos(l).setTo(domains[j].getSegmentAtPos(l+1).getTo());	  
	  domains[j].getSegmentAtPos(l).setFromOrg(domains[j].getSegmentAtPos(l+1).getFromOrg());
	  domains[j].getSegmentAtPos(l).setFrom(domains[j].getSegmentAtPos(l+1).getFrom());
	  domains[j].getSegmentAtPos(l).setChain(domains[j].getSegmentAtPos(l+1).getChain());
	}
	if (PDPParameters::VERBOSE){
	  printf("%i\n",i);
	}
	i--;
      }
    }
  }


  if (PDPParameters::VERBOSE){
    listdomains(domains);
  }
  listdomains(domains,"naive.pml"); //should allow user to specify file name
  // Remove short segments
  if (PDPParameters::VERBOSE){
    printf("---------Cleanup \n");
  }
  ShortSegmentRemover::cleanup(domains);
  if (PDPParameters::VERBOSE){
    printf("---------Cleanup Done\n");  
    printf("FINAL!!\n");
  }
  if (PDPParameters::VERBOSE){
    listdomains(domains);
  }  
  listdomains(domains,"removed.pml"); //should allow user to specify file name

  // PDB file exporter should be here
  
  return 0;
}

