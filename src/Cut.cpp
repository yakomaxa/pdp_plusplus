#include <math.h>
#include "Cut.hpp"
#include "Atom.hpp"
#include "CutValues.hpp"
#include "PDPDistanceMatrix.hpp"
#include "PDPParameters.hpp"

bool verbose = false;

int Cut::cut(std::vector<Atom>& ca,Domain& dom,CutValues& val,
             std::vector<std::vector<int>>& dist,
             PDPDistanceMatrix& pdpMatrix){
    
        int nclose = pdpMatrix.getNclose();
	////printf("nclose %i\n",nclose);
        std::vector<int> iclose = pdpMatrix.getIclose();
        std::vector<int> jclose = pdpMatrix.getJclose();


	int nclose_raw = pdpMatrix.getNclose_raw();
	std::vector<int> iclose_raw = pdpMatrix.getIclose_raw();
        std::vector<int> jclose_raw = pdpMatrix.getJclose_raw();

	
        std::vector<int> contacts(PDPParameters::MAXLEN);
        std::vector<double> max_contacts(PDPParameters::MAXLEN);
        std::vector<double> contact_density(PDPParameters::MAXLEN);
        double average_density,x,y;

        int endsf,endst;
        int k,l,nc, i,j;
        nc=0;
        int size1,size2;
        int size1t,size2t;
        int size11,size22,size0;
        int contactsd;
        int iseg,jseg,kseg;
        int from,to,from1,to1,from2,to2,lseg;


        int site_min = -1;


        // AP add sort here..
        // qsort(dom.segment,dom.getNseg(),sizeof(struct Segment),segcmp);
        // what is going on with the segments??

	//        std::vector<Segment> segments = dom.getSegments();
	std::sort(dom.getSegments().begin(), dom.getSegments().end(), SegmentComparator());
    
        average_density = 0.0;
        size0=0;
	//printf("Starting loop1\n");
        for(iseg=0;iseg<dom.getNseg();iseg++) {
	  contactsd=1;
	  size1t=0;
	  size2t=0;
	  //printf("-Starting loop2\n");
	  for(jseg=0;jseg<iseg;jseg++){	    
	    size1t+=(dom.getSegmentAtPos(jseg).getTo() - dom.getSegmentAtPos(jseg).getFrom() + 1);
	  }
	  
	  //printf("--Starting loop3\n");
	  for(jseg=iseg+1;jseg<dom.getNseg();jseg++){
	    size2t+=(dom.getSegmentAtPos(jseg).getTo() - dom.getSegmentAtPos(jseg).getFrom() + 1);
	  }
	  
	  //printf("---Starting loop4\n");

	  for (int n = 0; n < nclose_raw; n++){
	    i=iclose_raw[n];
	    j=jclose_raw[n];
	    for(jseg=0;jseg<iseg;jseg++) {
	      from1 = dom.getSegmentAtPos(jseg).getFrom();
	      to1 = dom.getSegmentAtPos(jseg).getTo();
	      //printf("----Starting loop5\n");
	      if(from1<= i && i<to1){
		//printf("-----Starting loop6\n");
		for(kseg=iseg+1;kseg<dom.getNseg();kseg++) {
		  from2 = dom.getSegmentAtPos(kseg).getFrom();
		  to2 = dom.getSegmentAtPos(kseg).getTo();
		  //printf("------Starting loop7\n");
		  if(from2<=j && j<to2){
		    if(abs(i-j)>4){
		      contactsd+=(dist[i][j]);
		    }
		  }
		}
	      }

	      if(from1<= j && j<to1){
		//printf("-----Starting loop6\n");
		for(kseg=iseg+1;kseg<dom.getNseg();kseg++) {
		  from2 = dom.getSegmentAtPos(kseg).getFrom();
		  to2 = dom.getSegmentAtPos(kseg).getTo();
		  //printf("------Starting loop7\n");
		  if(from2<=i && i<to2){
		    if(abs(i-j)>4){
		      contactsd+=(dist[i][j]);
		    }
		  }
		}
	      }	      
	    }
	  }
	  
	  from = dom.getSegmentAtPos(iseg).getFrom();
	  to = dom.getSegmentAtPos(iseg).getTo();
	  //printf("--------Starting loop8\n");
	  for(k=from;k<to;k++) {
	    contacts[k] = contactsd;
	    size11=size1t+(k-from+1);
	    size22=size2t+(to-k);
	    //printf("---------Starting loop9\n");
	    for (int n = 0; n < nclose_raw; n++){
	      i=iclose_raw[n];
	      j=jclose_raw[n];
	      if(from<=i && i<=k){
		//printf("----------Starting loop10\n");
		for(kseg=iseg+1;kseg<dom.getNseg();kseg++) {
		  from2 = dom.getSegmentAtPos(kseg).getFrom();
		  to2 = dom.getSegmentAtPos(kseg).getTo();
		  //printf("-----------Starting loop11\n");
		  if(from2<=j && j<=to2){
		    //printf("------------Starting loop12\n");
		    if(abs(i-j)>4){
		      contacts[k]+=(dist[i][j]);
		    }
		  }
		}
	      }
	      
	      if(from<=j && j<=k){
		//printf("----------Starting loop10\n");
		for(kseg=iseg+1;kseg<dom.getNseg();kseg++) {
		  from2 = dom.getSegmentAtPos(kseg).getFrom();
		  to2 = dom.getSegmentAtPos(kseg).getTo();
		  //printf("-----------Starting loop11\n");
		  if(from2<=i && i<=to2){
		    //printf("------------Starting loop12\n");
		    if(abs(i-j)>4){
		      contacts[k]+=(dist[i][j]);
		    }
		  }
		}
	      }
	      
	      if (from <= i && i <=k ){
		if (k+1 <= j && j <=to ){
		  if(abs(i-j)>4){
		    //printf("------------Starting loop13\n");
		    contacts[k]+=(dist[i][j]);
		  }
		}
	      }

	      if (from <= j && j <=k ){
		if (k+1 <= i && i <=to ){
		  if(abs(i-j)>4){
		    //printf("------------Starting loop13\n");
		    contacts[k]+=(dist[i][j]);
		  }
		}
	      }
	      
	      if ( k+1 <= i  && i<=to) {
		for(kseg=0;kseg<iseg;kseg++) {
		  from2 = dom.getSegmentAtPos(kseg).getFrom();
		  to2 = dom.getSegmentAtPos(kseg).getTo();
		  if(from2 <= j && j<to2){
		    if(abs(i-j)>4){
		      contacts[k]+=(dist[j][i]);
		    }
		  }
		}
	      }

	      if ( k+1 <= j  && j<=to) {
		for(kseg=0;kseg<iseg;kseg++) {
		  from2 = dom.getSegmentAtPos(kseg).getFrom();
		  to2 = dom.getSegmentAtPos(kseg).getTo();
		  if(from2 <= i && i<to2){
		    if(abs(i-j)>4){
		      contacts[k]+=(dist[j][i]);
		    }
		  }
		}
	      }
	    }
	    size1=std::min(size11,size22);
	    size2=std::max(size11,size22);
	    x=std::min(PDPParameters::MAXSIZE,size1);
	    y=std::min(PDPParameters::MAXSIZE,size2);
	    if(x>150&&y>1.5*x){
	      y=1.5*x;
	    }else if(y>2*x){
	      y=2*x;
	    };
	    x=std::min(pow(x,1.3/3)+PDPParameters::RG,pow(x,1.1/3)+pow(PDPParameters::TD,1.3/3)+PDPParameters::RG);
	    y=std::min(pow(y,1.3/3)+PDPParameters::RG,pow(y,1.1/3)+pow(PDPParameters::TD,1.3/3)+PDPParameters::RG);
	    
	    max_contacts[k] = 10*x*y;
	    if(size1>150){
	      max_contacts[k] = 9*x*y;
	      
	    };
	    contact_density[k]=contacts[k]/max_contacts[k];
	    
	    if(from==0){
	      endsf = PDPParameters::ENDSEND;
	    }else{
	      endsf = PDPParameters::ENDS;
	    };
	    
	    if(to==(int)ca.size()-1){
	      endst = PDPParameters::ENDSEND;
	    }else{
	      endst = PDPParameters::ENDS;
	    };
	    
	    if((contact_density[k])<val.s_min && k > from + endsf && k< to - endst) {
	      val.s_min = (contact_density[k]);
	      site_min=k+1;
	    }
	    if(k>from+endsf&&k<to-endst) {
	      average_density+=contact_density[k];
	      size0++;
	    }
	  };
	};
	if (size0>0){
	  average_density/=size0;
	}else{
	  return -1;	  
	}
	
	if(val.first_cut) {
	  val.AD = average_density;
	};
	val.AD = average_density;
	
	val.s_min/=val.AD;
	
	k=0;	
	nc=0;
	for(l=0;l<nclose;l++) {
	  /************ find iseg, jseg ****************/
	  iseg=jseg=-1;
	  for(kseg=0;kseg<dom.getNseg();kseg++) {
	    from=dom.getSegmentAtPos(kseg).getFrom();
	    to=dom.getSegmentAtPos(kseg).getTo();
	    if(from==0){
	      endsf = PDPParameters::ENDSEND;		      
	    }else{
	      endsf = PDPParameters::ENDS;
	    };
	    if(to==(int)ca.size()-1){
	      endst = PDPParameters::ENDSEND;
	    }else{
	      endst = PDPParameters::ENDS;
	    }
	    
	    if(iclose[l]>from+endsf&&iclose[l]<to-endst){
	      iseg=kseg;
	    }
	    if(jclose[l]>from+endsf&&jclose[l]<to-endst){
	      jseg=kseg;
	    }
	  }
	  if(iseg<0||jseg<0) continue;
	  /*********************************************/
	  
	  from=dom.getSegmentAtPos(iseg).getFrom();
	  to=dom.getSegmentAtPos(iseg).getTo();
	  from1=dom.getSegmentAtPos(jseg).getFrom();
	  to1=dom.getSegmentAtPos(jseg).getTo();
	  
	  /************ count contacts *****************/
	  contacts[nc] = 1;
	  
	  /******* contacts between [0,iseg[ and ]iseg,jseg[ ********/
	  for(kseg=0;kseg<iseg;kseg++)
	    for(lseg=iseg+1;lseg<jseg;lseg++)
	      for( int i=dom.getSegmentAtPos(kseg).getFrom();i<dom.getSegmentAtPos(kseg).getTo();i++)
		for(int j=dom.getSegmentAtPos(lseg).getFrom();j<dom.getSegmentAtPos(lseg).getTo();j++) {
		  contacts[nc]+=(dist[i][j]);
		}
	  
	  /******* contacts between ]jseg,nseg[ and ]iseg,jseg[ ********/
	  for(kseg=jseg+1;kseg<dom.getNseg();kseg++)
	    for(lseg=iseg+1;lseg<jseg;lseg++)
	      for(int i=dom.getSegmentAtPos(kseg).getFrom();i<dom.getSegmentAtPos(kseg).getTo();i++)
		for(int j=dom.getSegmentAtPos(lseg).getFrom();j<dom.getSegmentAtPos(lseg).getTo();j++) {
		  contacts[nc]+=(dist[j][i]);
		}
	  /*************************************************************/
	  /*************************************************************/
	  /**** contacts between [from,iclose] in iseg and ]iseg,jseg[ ****/
	  if(iseg==jseg) {
	    for(int i=from;i<=iclose[l];i++) {
	      for (int j=iclose[l]+1;j<=jclose[l];j++) {
		contacts[nc]+=(dist[i][j]);
	      }
	    }
	    for (int j=iclose[l]+1;j<jclose[l];j++) {
	      for(kseg=0;kseg<iseg;kseg++)
		for(int i=dom.getSegmentAtPos(kseg).getFrom();i<dom.getSegmentAtPos(kseg).getTo();i++) {
		  contacts[nc]+=(dist[i][j]);
		}
	      for(int i=jclose[l];i<to;i++) {
		contacts[nc]+=(dist[j][i]);
	      }
	      for(kseg=iseg+1;kseg<dom.getNseg();kseg++)
		for(int i=dom.getSegmentAtPos(kseg).getFrom();i<dom.getSegmentAtPos(kseg).getTo();i++) {
		  contacts[nc]+=(dist[j][i]);
		}
	    }
	  }
	  else {
	    for(int i=from;i<=iclose[l];i++) {
	      for(kseg=iseg+1;kseg<jseg;kseg++)
		for(int j=dom.getSegmentAtPos(kseg).getFrom();j<dom.getSegmentAtPos(kseg).getTo();j++) {
		  contacts[nc]+=(dist[i][j]);
		}
	      for(int j=from1;j<jclose[l];j++) {
		contacts[nc]+=(dist[i][j]);
	      }
	      for(int j=iclose[l]+1;j<to;j++) {
		contacts[nc]+=(dist[i][j]);
	      }
	    }
	    for(int i=iclose[l]+1;i<to;i++) {
	      for(kseg=0;kseg<iseg;kseg++)
		for(int j=dom.getSegmentAtPos(kseg).getFrom();j<dom.getSegmentAtPos(kseg).getTo();j++) {
		  contacts[nc]+=(dist[j][i]);
		}
	      for(kseg=jseg+1;kseg<dom.getNseg();kseg++)
		for(int j=dom.getSegmentAtPos(kseg).getFrom();j<dom.getSegmentAtPos(kseg).getTo();j++) {
		  contacts[nc]+=(dist[i][j]);
		}
	      for(int j=jclose[l];j<=to1;j++) {
		contacts[nc]+=(dist[i][j]);
	      }
	    }
	    for (int i=from1;i<jclose[l];i++) {
	      for(kseg=0;kseg<iseg;kseg++)
		for(int j=dom.getSegmentAtPos(kseg).getFrom();j<dom.getSegmentAtPos(kseg).getTo();j++) {
		  contacts[nc]+=(dist[j][i]);
		}
	      for(kseg=jseg+1;kseg<dom.getNseg();kseg++)  {
		for(int j=dom.getSegmentAtPos(kseg).getFrom();j<dom.getSegmentAtPos(kseg).getTo();j++)
		  contacts[nc]+=(dist[i][j]);
	      }
	      for(int j=jclose[l];j<to1;j++) {
		contacts[nc]+=(dist[i][j]);
	      }
	    }
	    for(int i=jclose[l];i<to1;i++)
	      for(kseg=iseg+1;kseg<jseg;kseg++)
		for(int j=dom.getSegmentAtPos(kseg).getFrom();j<dom.getSegmentAtPos(kseg).getTo();j++) {
		  contacts[nc]+=(dist[j][i]);
		}
	  }
	  /*******************************************************************/
	  /*******************************************************************/
	  /*******************************************************************/
	  size11=0;
	  size22=0;
	  for(kseg=0;kseg<iseg;kseg++){
	    size11+=(dom.getSegmentAtPos(kseg).getTo()-dom.getSegmentAtPos(kseg).getFrom()+1);
	  }
	  for(kseg=jseg+1;kseg<dom.getNseg();kseg++){
	    size11+=(dom.getSegmentAtPos(kseg).getTo()-dom.getSegmentAtPos(kseg).getFrom()+1);
	  }
	  size11+=(iclose[l]-from+1);
	  size11+=(to1-jclose[l]+1);
	  
	  for(kseg=iseg+1;kseg<jseg;kseg++){
	    size22+=(dom.getSegmentAtPos(kseg).getTo()-dom.getSegmentAtPos(kseg).getFrom()+1);
	  }
	  
	  if(iseg==jseg){
	    size22+=(jclose[l]-iclose[l]);
	  }else{
	    size22+=(jclose[l]-from1);
	    size22+=(to-iclose[l]);
	  };
	  
	  size1=std::min(size11,size22);
	  size2=std::max(size11,size22);
	  x=std::min(PDPParameters::MAXSIZE,size1);
	  y=std::max(PDPParameters::MAXSIZE,size2);
	  if(y>2*x){
	    y=2*x;
	  };
	  
	  x=std::min(pow(x,1.3/3)+PDPParameters::RG,pow(x,1.1/3)+pow(PDPParameters::TD,1.3/3)+PDPParameters::RG);
	  y=std::min(pow(y,1.3/3)+PDPParameters::RG,pow(y,1.1/3)+pow(PDPParameters::TD,1.3/3)+PDPParameters::RG);
	  
	  max_contacts[nc] = x*y*10;
	  if(size1>150){
	    max_contacts[k] = 9*x*y;
	  }
	  contact_density[nc]=contacts[nc]/max_contacts[nc];
	  if((contact_density[nc]/val.AD+PDPParameters::DBL)<val.s_min&&contact_density[nc]/val.AD+PDPParameters::DBL<PDPParameters::CUT_OFF_VALUE2) {
	    val.s_min = (contact_density[nc]/val.AD)+PDPParameters::DBL;
	    site_min=iclose[l];
	    val.site2=jclose[l];
	  }
	  nc++;
	  if ( nc >= PDPParameters::MAXSIZE)
	    nc = PDPParameters::MAXSIZE-1;
	}
	val.first_cut=false;
	if(val.s_min> PDPParameters::CUT_OFF_VALUE){
	  return -1;
	}
	
	return(site_min);
}


