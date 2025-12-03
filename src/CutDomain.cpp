#include "CutDomain.hpp"

CutDomain::CutDomain(std::vector<Atom> &ca, PDPDistanceMatrix &pdpMatrix, std::vector<int> &init_cutsites){
    this->dist = pdpMatrix.getDist();
    this->ca = ca;
    this->ndom = 1;
    this->domains = std::vector<Domain>();
    this->init_cutsites = init_cutsites;
}

void CutDomain::cutDomain(Domain& dom, CutSites& cut_sites,PDPDistanceMatrix& pdpMatrix) {
  int site;
  CutValues val;
  val.s_min = 100;
  val.site2 = 0;
  val.first_cut = true;

  if (init_cutsites.size()<=0){
      printf("CUTTING DE NOVO\n");
      Cut cut;
      site = cut.cut(ca, dom, val, dist, pdpMatrix);
    }else{
      printf("CUTTING OF GIVEN\n");
      site = init_cutsites.back();
      init_cutsites.pop_back();
    }
    printf("site %i \n",site)   ;
    if (site < 0) {
      dom.setScore(val.s_min);
      domains.push_back(dom);
      ndom++;
      std::cout << "HOGE NDOM" << ndom << std::endl;
      return;
    }
    
    cut_sites.addNcuts(1);
    cut_sites.pushbackCutSites(site);
    std::cout << "HOGE SITE " << site << std::endl;
    std::cout << "HOGE SITE2 " << val.site2 << std::endl;    
    Domain dom1;
    dom1.setNseg(0);
    dom1.setSize(0);

    Domain dom2;
    dom2.setNseg(0);
    dom2.setSize(0);
    if (val.site2 == 0) {
        for (int i = 0; i < dom.getNseg(); i++) {
            if (site > dom.getSegmentAtPos(i).getTo()) {
                dom1.getSegmentAtPos(dom1.getNseg()).setTo(dom.getSegmentAtPos(i).getTo());
                dom1.getSegmentAtPos(dom1.getNseg()).setFrom(dom.getSegmentAtPos(i).getFrom());
                dom1.addNseg(1);
                dom1.addSize(dom.getSegmentAtPos(i).getTo() - dom.getSegmentAtPos(i).getFrom() + 1);
            } else if (site < dom.getSegmentAtPos(i).getFrom()) {
                dom2.getSegmentAtPos(dom2.getNseg()).setTo(dom.getSegmentAtPos(i).getTo());
                dom2.getSegmentAtPos(dom2.getNseg()).setFrom(dom.getSegmentAtPos(i).getFrom());
                dom2.addNseg(1);
                dom2.addSize(dom.getSegmentAtPos(i).getTo() - dom.getSegmentAtPos(i).getFrom() + 1);
            } else if (site > dom.getSegmentAtPos(i).getFrom() && site < dom.getSegmentAtPos(i).getTo()) {
	      dom1.getSegmentAtPos(dom1.getNseg()).setFrom(dom.getSegmentAtPos(i).getFrom());
	      dom1.getSegmentAtPos(dom1.getNseg()).setTo(site-1);
	      dom1.addNseg(1);
	      dom1.addSize(site - dom.getSegmentAtPos(i).getFrom());
	      dom2.getSegmentAtPos(dom2.getNseg()).setTo(dom.getSegmentAtPos(i).getTo());
	      dom2.getSegmentAtPos(dom2.getNseg()).setFrom(site);
	      dom2.addNseg(1);
	      dom2.addSize(dom.getSegmentAtPos(i).getTo() - site + 1);
            }

        }
    } else if(val.site2>0) { /* double cut */
      for(int i=0;i<dom.getNseg();i++) {
	if(site>dom.getSegmentAtPos(i).getTo()||val.site2<dom.getSegmentAtPos(i).getFrom()) {
	  dom1.getSegmentAtPos(dom1.getNseg()).setTo(dom.getSegmentAtPos(i).getTo());
	  dom1.getSegmentAtPos(dom1.getNseg()).setFrom(dom.getSegmentAtPos(i).getFrom());
	  dom1.addNseg(1);
	  dom1.addSize(dom.getSegmentAtPos(i).getTo() - dom.getSegmentAtPos(i).getFrom() + 1);
	}
	else if(site<dom.getSegmentAtPos(i).getFrom()&&val.site2>dom.getSegmentAtPos(i).getTo()) {
	  dom2.getSegmentAtPos(dom1.getNseg()).setTo(dom.getSegmentAtPos(i).getTo());
	  dom2.getSegmentAtPos(dom1.getNseg()).setFrom(dom.getSegmentAtPos(i).getFrom());
	  dom2.addNseg(1);
	  dom2.addSize(dom.getSegmentAtPos(i).getTo() - dom.getSegmentAtPos(i).getFrom() + 1);
	}
	else if(site>dom.getSegmentAtPos(i).getFrom() &&
		site<dom.getSegmentAtPos(i).getTo()) {
	  dom1.getSegmentAtPos(dom1.getNseg()).setTo(site);
	  dom1.getSegmentAtPos(dom1.getNseg()).setFrom(dom.getSegmentAtPos(i).getFrom());
	  dom1.addSize(dom1.getSegmentAtPos(dom1.getNseg()).getTo() - dom1.getSegmentAtPos(dom1.getNseg()).getFrom() + 1);
	  dom1.addNseg(1);
	  dom2.getSegmentAtPos(dom2.getNseg()).setFrom(site+1);
	  if(val.site2>dom.getSegmentAtPos(i).getFrom() &&
	     val.site2<dom.getSegmentAtPos(i).getTo()) {
	    dom2.getSegmentAtPos(dom2.getNseg()).setTo(val.site2-1);
	    dom2.addSize(dom2.getSegmentAtPos(dom2.getNseg()).getTo() - dom2.getSegmentAtPos(dom2.getNseg()).getFrom() + 1);
	    dom2.addNseg(1);
	    dom1.getSegmentAtPos(dom1.getNseg()).setFrom( val.site2);
	    dom1.getSegmentAtPos(dom1.getNseg()).setTo( dom.getSegmentAtPos(i).getTo());
	    dom1.addSize(dom1.getSegmentAtPos(dom1.getNseg()).getTo() - dom1.getSegmentAtPos(dom1.getNseg()).getFrom() + 1);
	    dom1.addNseg(1);
	  }
	  else {
	    dom2.getSegmentAtPos(dom2.getNseg()).setTo(dom.getSegmentAtPos(i).getTo());
	    dom2.addSize(dom2.getSegmentAtPos(dom2.getNseg()).getTo() - dom2.getSegmentAtPos(dom2.getNseg()).getFrom() + 1);
	    dom2.addNseg(1);
	  }
	}
	else if(val.site2>dom.getSegmentAtPos(i).getFrom() &&
		val.site2<dom.getSegmentAtPos(i).getTo()) {
	  dom2.getSegmentAtPos(dom2.getNseg()).setTo(val.site2-1);
	  dom2.getSegmentAtPos(dom2.getNseg()).setFrom(dom.getSegmentAtPos(i).getFrom());
	  dom2.addSize(dom2.getSegmentAtPos(dom2.getNseg()).getTo() - dom2.getSegmentAtPos(dom2.getNseg()).getFrom() + 1);
	  dom2.addNseg(1);
	  dom1.getSegmentAtPos(dom1.getNseg()).setFrom(val.site2);
	  dom1.getSegmentAtPos(dom1.getNseg()).setTo( dom.getSegmentAtPos(i).getTo());
	  dom1.addSize(dom1.getSegmentAtPos(dom1.getNseg()).getTo() - dom1.getSegmentAtPos(dom1.getNseg()).getFrom() + 1);
	  dom1.addNseg(1);
	}
      }
    }

    /**
    if(verbose){
      std::cout << "  CUTR dom1 ...  nse" << dom1.getNseg() << 
;
    }

    if ( verbose){
      for(i=0;i<dom1.getNseg();i++){
	std::cout << "F ... from %d to %d" << dom1.getSegmentAtPos(i).getFrom() << " "
		  << dom1.getSegmentAtPos(i).getTo() << 
;
      }
    }
    **/
    
    cutDomain(dom1, cut_sites, pdpMatrix);

    /**
    if(verbose){
      std::cout << "  CUTR dom2 ...  nse" << dom2.getNseg() << 
;
    }

    if ( verbose){
      for(i=0;i<dom2.getNseg();i++){
	std::cout << "F ... from %d to %d" << dom2.getSegmentAtPos(i).getFrom() << " " 
		  << dom2.getSegmentAtPos(i).getTo() << 
;
      }
    }
    **/

    cutDomain(dom2, cut_sites, pdpMatrix);

};

std::vector<Domain> CutDomain::getDomains(){
    return domains;
};



