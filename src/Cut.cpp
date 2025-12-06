#include <math.h>
#include "Cut.hpp"
#include "Atom.hpp"
#include "CutValues.hpp"
#include "PDPDistanceMatrix.hpp"
#include "PDPParameters.hpp"

bool verbose = PDPParameters::VERBOSE;;
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

  std::vector<int> contacts(nclose_raw);
  std::vector<int> max_contacts(nclose_raw);
  // (K.S. 2025-Dec-06) I used to define this max_contacts as vector<double>, which caused slight but persistent inconsistency from the original because of floating point in 10*x*y or 9*x*y to be injected into max_contacts[nc]. Now it's made int, and 10*x*y and 9*x*y are explicitly cast to int as (int)(10*x*y), realizing perfect consistency to the original PDP.

  std::vector<double> contact_density(nclose_raw);
  double average_density,x,y;
  
  int endsf,endst;
  int nc;
  int size1,size2;
  int size1t,size2t;
  int size11,size22,size0;
  int contactsd;
  int from,to,from1,to1,from2,to2;
  int i, j;
  
  int site_min = -1;
  val.s_min = 100;
  
  // AP add sort here..
  // qsort(dom.segment,dom.getNseg(),sizeof(struct Segment),segcmp);
  // what is going on with the segments??
  
  //        std::vector<Segment> segments = dom.getSegments();
  std::sort(dom.getSegments().begin(), dom.getSegments().end(), SegmentComparator());
  
  average_density = 0.0;
  size0=0;
  //printf("Starting loop1\n");
  for(int iseg=0;iseg<dom.getNseg();iseg++) {
    contactsd=1;
    size1t=0;
    size2t=0;
    //printf("-Starting loop2\n");
    for(int jseg=0;jseg<iseg;jseg++){	    
      size1t+=(dom.getSegmentAtPos(jseg).getTo() - dom.getSegmentAtPos(jseg).getFrom() + 1);
    }
    
    //printf("--Starting loop3\n");
    for(int jseg=iseg+1;jseg<dom.getNseg();jseg++){
      size2t+=(dom.getSegmentAtPos(jseg).getTo() - dom.getSegmentAtPos(jseg).getFrom() + 1);
    }
    
    //printf("---Starting loop4\n");
    for (int n = 0; n < nclose_raw; n++){
      i=iclose_raw[n];
      j=jclose_raw[n];
      if(abs(i-j)>4){
	for(int jseg=0;jseg<iseg;jseg++) {
	  from1 = dom.getSegmentAtPos(jseg).getFrom();
	  to1 = dom.getSegmentAtPos(jseg).getTo();
	  if(from1<= i && i<to1){
	    for(int kseg=iseg+1;kseg<dom.getNseg();kseg++) {
	      from2 = dom.getSegmentAtPos(kseg).getFrom();
	      to2 = dom.getSegmentAtPos(kseg).getTo();
	      if(from2<=j && j<to2){
		contactsd+=(dist[i][j]);
	      }
	    }
	  }else if(from1<= j && j<to1){
	    for(int kseg=iseg+1;kseg<dom.getNseg();kseg++) {
	      from2 = dom.getSegmentAtPos(kseg).getFrom();
	      to2 = dom.getSegmentAtPos(kseg).getTo();
	      if(from2<=i && i<to2){
		contactsd+=(dist[j][i]);
	      }
	    }
	  }    
	}
      }
    }
    
    from = dom.getSegmentAtPos(iseg).getFrom();
    to = dom.getSegmentAtPos(iseg).getTo();	  
    //printf("--------Starting loop8\n");
    for(int k=from;k<to;k++) {
      contacts[k] = contactsd;
      if (verbose){
	printf("init contacts = %d\n",contacts[k]);
      }
      size11=size1t+(k-from+1);
      size22=size2t+(to-k);
      //printf("---------Starting loop9\n");
      for (int n = 0; n < nclose_raw; n++){
	i=iclose_raw[n];
	j=jclose_raw[n];
	if(abs(i-j)>4){
	  if(from<=i && i<=k){ // condition 1
	    for(int kseg=iseg+1;kseg<dom.getNseg();kseg++) {
	      from2 = dom.getSegmentAtPos(kseg).getFrom();
	      if(from2<=j){ // made nest rather than AND like if (from2 <= j && j<= to)
		to2 = dom.getSegmentAtPos(kseg).getTo();
		if(j<=to2){
		  contacts[k]+=(dist[i][j]);
		}
	      }
	    }
	    if (k+1 <= j && j <=to ){
	      contacts[k]+=(dist[i][j]);
	    }
	  }else if(from<=j && j<=k){ // i-j swap of condition 1
	    for(int kseg=iseg+1;kseg<dom.getNseg();kseg++) {
	      from2 = dom.getSegmentAtPos(kseg).getFrom();
	      if(from2<=i){
		to2 = dom.getSegmentAtPos(kseg).getTo();
		if (i<=to2){
		  contacts[k]+=(dist[i][j]);
		}
	      }
	    }	    
	    if (k+1 <= i && i <=to ){
	      contacts[k]+=(dist[i][j]);
	    }
	  }else if ( k+1 <= i  && i<=to) { // condition 2
	    for(int kseg=0;kseg<iseg;kseg++) {
	      from2 = dom.getSegmentAtPos(kseg).getFrom();
	      if(from2 <= j){
		to2 = dom.getSegmentAtPos(kseg).getTo();
		if (j<to2){
		  contacts[k]+=(dist[j][i]);
		}
	      }
	    }
	  }else if ( k+1 <= j  && j<=to) { // i-j swap of condition
	    for(int kseg=0;kseg<iseg;kseg++) {
	      from2 = dom.getSegmentAtPos(kseg).getFrom();
	      if(from2 <= i){
		to2 = dom.getSegmentAtPos(kseg).getTo();
		if( i<to2){
		  contacts[k]+=(dist[j][i]);
		}
	      }
	    }
	  }
	}
      }
      size1=std::min(size11,size22);
      size2=std::max(size11,size22);
      x=std::min(PDPParameters::MAXSIZE,size1);
      y=std::min(PDPParameters::MAXSIZE,size2);
      if(verbose){
	printf("init contacts = %d\n",contacts[k]);
      }
      
      if(x>150&&y>1.5*x){
	y=1.5*x;
      }else if(y>2*x){
	y=2*x;
      };
      x=std::min(pow(x,1.3/3)+PDPParameters::RG,pow(x,1.1/3)+pow(PDPParameters::TD,1.3/3)+PDPParameters::RG);
      y=std::min(pow(y,1.3/3)+PDPParameters::RG,pow(y,1.1/3)+pow(PDPParameters::TD,1.3/3)+PDPParameters::RG);
      
      max_contacts[k] = (int)(10*x*y);
      if(size1>150){
	max_contacts[k] = (int)(9*x*y);
	
      };
      contact_density[k]=(double)contacts[k]/(double)max_contacts[k];
      
      if(verbose){
	printf("data %i  %i      %i      %f      %f      %i      %i      %f\n",k,size1,size2,x,y,max_contacts[k],contacts[k],contact_density[k]);
      }
      
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
      
      //      if((contact_density[k])<val.s_min && k > from + endsf && k< to - endst) {
      //	val.s_min = (contact_density[k]);
      //	site_min=k+1;
      //      }
      if(k>from+endsf&&k<to-endst) {
	if((contact_density[k])<val.s_min){
	  val.s_min = (contact_density[k]);
	  site_min=k+1;
	}
	average_density+=contact_density[k];
	size0++;
      }
    }
  }
  
  if (!size0==0){ // K.S.: Added this to avoid 0-division.
    average_density/=size0;
    if(verbose){
      printf("Trying to cut domain of size %d having %d segments and  average cont_density %f\n",dom.getSize(),dom.getNseg(),average_density);
      for(int kseg=0;kseg<dom.getNseg();kseg++){
	//to=dom.getSegmentAtPos(iseg).getTo();
	printf("Trying segment %d from %d to %d\n",kseg,dom.getSegmentAtPos(kseg).getFrom(),dom.getSegmentAtPos(kseg).getTo());
      }
    }
  }else{
    val.AD = 123456;
    if(verbose){
      printf("could have had NaN because size0 == 0, forced to exit!");
      printf("at the end of cut: s_min %f CUTOFF %f site_min %d *site2 %d\n",val.s_min,PDPParameters::CUT_OFF_VALUE,site_min,val.site2);
    }
    return -1;
  }
    
  if(val.first_cut) { // K.S.: What does this if sentence do...? I think the flag first_cut is doing nothing.
    val.AD = average_density;
  };
  val.AD = average_density;
  
  if(verbose){
    printf("AD=%f\n", average_density);
  }

  if (!val.AD == 0){ // K.S.: Added this to avoid 0-division.
    val.s_min/=val.AD;
  }else{
    if(verbose){
      printf("could have had NaN because val.AD  == 0, forced to exit!");
      printf("at the end of cut: s_min %f CUTOFF %f site_min %d *site2 %d\n",val.s_min,PDPParameters::CUT_OFF_VALUE,site_min,val.site2);
    }
    return -1;
  }
  
  if(verbose){
    printf("after single cut: s_min = %f site_min = %d\n",val.s_min,site_min);
  }
  
  nc=0;
  for(int l=0;l<nclose;l++) {
    /************ find iseg, jseg ****************/
    int iseg=-1;
    int jseg=-1;
    for(int kseg=0;kseg<dom.getNseg();kseg++) {
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
    
    if(iseg<0||jseg<0){
      continue;
    }
    /*********************************************/
    
    from=dom.getSegmentAtPos(iseg).getFrom();
    to=dom.getSegmentAtPos(iseg).getTo();
    from1=dom.getSegmentAtPos(jseg).getFrom();
    to1=dom.getSegmentAtPos(jseg).getTo();
    
    /************ count contacts *****************/
    contacts[nc] = 1;
    
    /******* contacts between [0,iseg[ and ]iseg,jseg[ ********/
    for(int kseg=0;kseg<iseg;kseg++){
      for(int lseg=iseg+1;lseg<jseg;lseg++){
	for( int ii=dom.getSegmentAtPos(kseg).getFrom();ii<dom.getSegmentAtPos(kseg).getTo();ii++){
	  for(int jj=dom.getSegmentAtPos(lseg).getFrom();jj<dom.getSegmentAtPos(lseg).getTo();jj++) {
	    contacts[nc]+=(dist[ii][jj]);
	  }
	}
      }
    }
    
    /******* contacts between ]jseg,nseg[ and ]iseg,jseg[ ********/
    for(int kseg=jseg+1;kseg<dom.getNseg();kseg++){
      for(int lseg=iseg+1;lseg<jseg;lseg++){
	for(int ii=dom.getSegmentAtPos(kseg).getFrom();ii<dom.getSegmentAtPos(kseg).getTo();ii++){
	  for(int jj=dom.getSegmentAtPos(lseg).getFrom();jj<dom.getSegmentAtPos(lseg).getTo();jj++) {
	    contacts[nc]+=(dist[jj][ii]);
	  }
	}
      }
    }
    /*************************************************************/
    /*************************************************************/
    /**** contacts between [from,iclose] in iseg and ]iseg,jseg[ ****/
    if(iseg==jseg) {
      for(int ii=from;ii<=iclose[l];ii++) {
	for (int jj=iclose[l]+1;jj<=jclose[l];jj++) {
	  contacts[nc]+=(dist[ii][jj]);
	}
      }
      for (int jj=iclose[l]+1;jj<jclose[l];jj++) {
	for(int kseg=0;kseg<iseg;kseg++){
	  for(int ii=dom.getSegmentAtPos(kseg).getFrom();ii<dom.getSegmentAtPos(kseg).getTo();ii++) {
	    contacts[nc]+=(dist[ii][jj]);
	  }
	}
	for(int ii=jclose[l];ii<to;ii++) {
	  contacts[nc]+=(dist[jj][ii]);
	}
	for(int kseg=iseg+1;kseg<dom.getNseg();kseg++){
	  for(int ii=dom.getSegmentAtPos(kseg).getFrom();ii<dom.getSegmentAtPos(kseg).getTo();ii++) {
	    contacts[nc]+=(dist[jj][ii]);
	  }
	}
      }
    }
    else {
      for(int ii=from;ii<=iclose[l];ii++) {
	for(int kseg=iseg+1;kseg<jseg;kseg++){
	  for(int jj=dom.getSegmentAtPos(kseg).getFrom();jj<dom.getSegmentAtPos(kseg).getTo();jj++) {
	    contacts[nc]+=(dist[ii][jj]);
	  }
	}
	for(int jj=from1;jj<jclose[l];jj++) {
	  contacts[nc]+=(dist[ii][jj]);
	}
	for(int jj=iclose[l]+1;jj<to;jj++) {
	  contacts[nc]+=(dist[ii][jj]);
	}
      }
      for(int ii=iclose[l]+1;ii<to;ii++) {
	for(int kseg=0;kseg<iseg;kseg++){
	  for(int jj=dom.getSegmentAtPos(kseg).getFrom();jj<dom.getSegmentAtPos(kseg).getTo();jj++) {
	    contacts[nc]+=(dist[jj][ii]);
	  }
	}
	for(int kseg=jseg+1;kseg<dom.getNseg();kseg++){
	  for(int jj=dom.getSegmentAtPos(kseg).getFrom();jj<dom.getSegmentAtPos(kseg).getTo();jj++) {
	    contacts[nc]+=(dist[ii][jj]);
	  }
	}
	for(int jj=jclose[l];jj<=to1;jj++) {
	  contacts[nc]+=(dist[ii][jj]);
	}
      }
      for (int ii=from1;ii<jclose[l];ii++) {
	for(int kseg=0;kseg<iseg;kseg++){
	  for(int jj=dom.getSegmentAtPos(kseg).getFrom();jj<dom.getSegmentAtPos(kseg).getTo();jj++) {
	    contacts[nc]+=(dist[jj][ii]);
	  }
	}
	for(int kseg=jseg+1;kseg<dom.getNseg();kseg++)  {
	  for(int jj=dom.getSegmentAtPos(kseg).getFrom();jj<dom.getSegmentAtPos(kseg).getTo();jj++)
	    contacts[nc]+=(dist[ii][jj]);
	}
	for(int jj=jclose[l];jj<to1;jj++) {
	  contacts[nc]+=(dist[ii][jj]);
	}
      }
      for(int ii=jclose[l];ii<to1;ii++){
	for(int kseg=iseg+1;kseg<jseg;kseg++){
	  for(int jj=dom.getSegmentAtPos(kseg).getFrom();jj<dom.getSegmentAtPos(kseg).getTo();jj++) {
	    contacts[nc]+=(dist[jj][ii]);
	  }
	}
      }
    }
    /*******************************************************************/
    /*******************************************************************/
    /*******************************************************************/
    size11=0;
    size22=0;
    for(int kseg=0;kseg<iseg;kseg++){
      size11+=(dom.getSegmentAtPos(kseg).getTo()-dom.getSegmentAtPos(kseg).getFrom()+1);
    }
    for(int kseg=jseg+1;kseg<dom.getNseg();kseg++){
      size11+=(dom.getSegmentAtPos(kseg).getTo()-dom.getSegmentAtPos(kseg).getFrom()+1);
    }
    size11+=(iclose[l]-from+1);
    size11+=(to1-jclose[l]+1);
    
    for(int kseg=iseg+1;kseg<jseg;kseg++){
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
    
    max_contacts[nc] = (int)(x*y*10);
    if(size1>150){
      max_contacts[nc] = (int)(9*x*y); // the original code seems to have wrong index max_contacts[k] instead of [nc] here. 
    }
    contact_density[nc]=(double)contacts[nc]/(double)max_contacts[nc]; 
    
    if(verbose){
      std::cout << max_contacts[nc] << std::endl;
      printf(" double cut: %i %s %i %i c=%d mc=%f x=%f y=%f s1=%i s2=%i cd=%f cd/ad=%f\n",l,ca[iclose[l]].getResidue(),iclose[l],jclose[l],contacts[nc],max_contacts[nc],x,y,size11,size22,contact_density[nc],contact_density[nc]/val.AD);
    }
    
    if((contact_density[nc]/val.AD+PDPParameters::DBL)<val.s_min&&contact_density[nc]/val.AD+PDPParameters::DBL<PDPParameters::CUT_OFF_VALUE2) {
      val.s_min = (contact_density[nc]/val.AD)+PDPParameters::DBL;
      site_min=iclose[l];
      val.site2=jclose[l];
    }
    nc++;
    if (verbose){
      std::cout << "NC++" << nc << std::endl;
    }
    //if ( nc >= PDPParameters::MAXSIZE)
    //	    nc = PDPParameters::MAXSIZE-1;
  }
  val.first_cut=false;
  if(verbose){
    printf("at the end of cut: s_min %f CUTOFF %f site_min %d *site2 %d\n",val.s_min,PDPParameters::CUT_OFF_VALUE,site_min,val.site2);
  }
  if(val.s_min > PDPParameters::CUT_OFF_VALUE){
    return -1;
  }
  
  return(site_min);
}


