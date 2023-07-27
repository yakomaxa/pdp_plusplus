#include "CutDomain.hpp"

bool CutDomain::verbose = true;
int ndom;

CutDomain::CutDomain(std::vector<Atom> &ca, PDPDistanceMatrix &pdpMatrix){
    dist = pdpMatrix.getDist();
    this->ca = ca;
    ndom = 0;
    domains = std::vector<Domain>();
}

void CutDomain::cutDomain(Domain& dom, CutSites& cut_sites, PDPDistanceMatrix& pdpMatrix) {
 
    int i, site;
    
    Domain dom1 = Domain();
    Domain dom2 = Domain();
    
    CutValues val = CutValues();
    val.s_min = 100;
    val.site2 = 0;
    val.first_cut = true;
    
    Cut cut = Cut();
    
    site = cut.cut(ca, dom, val, dist, pdpMatrix);
    printf("site %i \n",site)   ;
    if (site < 0) {
      //printf("site<0\n");
      dom.setScore(val.s_min);
      domains.push_back(dom);
      ndom++;
      return;
    }

    printf("CUT_SITE Ncuts= %i\n",    cut_sites.getNcuts());
    cut_sites.addNcuts(1);
    printf("CUT_SITE Ncuts= %i\n",    cut_sites.getNcuts());
    cut_sites.cut_sites[cut_sites.getNcuts()]=site;

    dom1.setSize(0);
    dom1.setNseg(0);
    dom2.setSize(0);
    dom2.setNseg(0);
    if (val.site2 == 0) {
        for (i = 0; i < dom.getNseg(); i++) {
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
      for(i=0;i<dom.getNseg();i++) {
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

    if(verbose){
      std::cout << "  CUTR dom1 ...  nse" << dom1.getNseg() << std::endl;
    }

    if ( verbose){
      for(i=0;i<dom1.getNseg();i++){
	std::cout << "F ... from %d to %d" << dom1.getSegmentAtPos(i).getFrom() << " "
		  << dom1.getSegmentAtPos(i).getTo() << std::endl;
      }
    }
       
    cutDomain(dom1, cut_sites, pdpMatrix);

    
    if(verbose){
      std::cout << "  CUTR dom2 ...  nse" << dom2.getNseg() << std::endl;
    }

    if ( verbose){
      for(i=0;i<dom2.getNseg();i++){
	std::cout << "F ... from %d to %d" << dom2.getSegmentAtPos(i).getFrom() << " " 
		  << dom2.getSegmentAtPos(i).getTo() << std::endl;
      }
    }


    cutDomain(dom2, cut_sites, pdpMatrix);
};

std::vector<Domain> CutDomain::getDomains(){
    return domains;
};



