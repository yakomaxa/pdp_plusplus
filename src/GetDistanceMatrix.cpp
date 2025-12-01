#include <iostream>
#include <math.h>
#include "Group.hpp"
#include "Atom.hpp"
#include "PDPDistanceMatrix.hpp"
#include "PDPParameters.hpp"
#include "GetDistanceMatrix.hpp"

float getDistance(Atom atom1,Atom atom2);
float getDistance(Atom atom1,Atom atom2){
    float dist=0;
    float dx=0;
    float dy=0;
    float dz=0;
    
    dx = atom1.getX() - atom2.getX();
    dy = atom1.getY() - atom2.getY();
    dz = atom1.getZ() - atom2.getZ();
    
    dist = sqrt(dx*dx + dy*dy + dz*dz);
    
    return dist;
}


PDPDistanceMatrix GetDistanceMatrix::getDistanceMatrix(std::vector<Atom>& protein)
{
    
  float dx=0;
  float dy=0;
  float dz=0;
  
  std::vector<std::vector<int>> dist(PDPParameters::MAXLEN+3,
				     std::vector<int>(PDPParameters::MAXLEN+3));
  int i,j;
  float d,dt1,dt2,dt3,dt4;
  int nc=0;  

  std::cerr << protein.size() << " protein.len < MAXLEN " << PDPParameters::MAXLEN << "\n";
  if((int)protein.size() >= PDPParameters::MAXLEN) {
    std::cerr << protein.size() << " protein.len > MAXLEN " << PDPParameters::MAXLEN << "\n";
    //return 0;
    exit;
  }

  dt1=81;
  dt2=64;
  dt3=49;
  dt4=36;

  // counter loop to count loosest contacts 
  for(i=0; i < (int)protein.size()-1; i++) {
    Atom ca1 = protein.at(i);
    for(j=i+1; j < (int)protein.size(); j++) {
      Atom ca2 = protein.at(j);          
      dx = ca1.getX() - ca2.getX();
      dy = ca1.getY() - ca2.getY();
      dz = ca1.getZ() - ca2.getZ();	
      d = dx*dx + dy*dy + dz*dz;
      if(d<dt1) {
	nc++;
      }
    }
  }

  // allocate vector based on the number of contacts
  std::vector<int> iclose(nc);
  std::vector<int> jclose(nc);  
  std::vector<int> iclose_raw(nc);
  std::vector<int> jclose_raw(nc);

  int nclose_raw=0;
  int nclose=0;

  // make dist (not distance but contacts although...) matrix
  for(i=0; i < (int)protein.size(); i++) {
    dist[i][i]=0;       
  }

  for(i=0; i < (int)protein.size()-1; i++) {
    Atom ca1 = protein.at(i);
    std::string  chain1 = protein.at(i).getChain();
    for(j=i+1; j < (int)protein.size(); j++) {
      std::string  chain2 = protein.at(j).getChain();
      Atom ca2 = protein.at(j);          
      dx = ca1.getX() - ca2.getX();
      dy = ca1.getY() - ca2.getY();
      dz = ca1.getZ() - ca2.getZ();	
      d = dx*dx + dy*dy + dz*dz;
      
      if(d<dt1){
	dist[i][j]=1;
	dist[j][i]=1;

	iclose_raw[nclose_raw]=i;
	jclose_raw[nclose_raw]=j;
	nclose_raw++;
			
	if(d<dt2) {
	  dist[i][j]=2;
	  dist[j][i]=2;	  
	  if(j-i>35 || (chain1!=chain2)){
	    iclose[nclose]=i;
	    jclose[nclose]=j;
	    nclose++;
	  }
	  if(d<dt3) {
	    dist[i][j]=4;
	    dist[j][i]=4;
	    if(d<dt4) {
	      dist[i][j]=6;
	      dist[j][i]=6;
	    }
	  }
	}
      }
    }
  }

  for(int i=1;i<(int)protein.size()-1;i++) {
    for(int j=i;j<(int)protein.size()-1;j++) {	
      if(dist[i][j]>=2&&j-i>5) {
	if((dist[i-1][j-1]>=2&&dist[i+1][j+1]>=2)||(dist[i-1][j+1]>=2&&dist [i+1][j-1]>=2))  {
	  dist[i][j]+=4;
	  dist[j][i]+=4;
	  //	    printf("COND1 I=%i J=%i\n",i,j);
	}
	else if(i>2&&j<(int)protein.size()-2) {
	  if((dist[i-3][j-3]>=1&&dist[i+3][j+3]>=1)||(dist[i-3][j+3]>=1&&dist[i+3][j-3]>=1)) {
	    dist[i][j]+=4;
	    dist[j][i]+=4;
	    //	      printf("COND2 I=%i J=%i\n",i,j);
	  }
	  else if(i>3&&j<(int)protein.size()-3) {
	    if(((dist[i-3][j-3]>=1||dist[i-3][j-4]>=1||dist[i-4][j-3]>=1||dist[i-4][j-4]>=1)&&
		(dist[i+4][j+4]>=1||dist[i+4][j+3]>=1||dist[i+3][j+3]>=1||dist[i+3] [j+4]>=1))
	       ||((dist[i-4][j+4]>=1||dist[i-4][j+3]>=1||dist[i-3][j+4]>=1||dist[i-3][j+3]>=1)&&
		  (dist[i+4][j-4]>=1||dist[i+4][j-3]>=1||dist[i+3][j-4]>=1||dist[i+3][j-3]>=1))) {
	      dist[i][j]+=4;
	      dist[j][i]+=4;
	      //		printf("COND3 I=%i J=%i\n",i,j);
	    }	     
	  }
	}
      }
    }
  }
  PDPDistanceMatrix matrix;
  matrix.setNclose(nclose);
  matrix.setIclose(iclose);
  matrix.setJclose(jclose);
  
  matrix.setNclose_raw(nclose_raw);
  matrix.setIclose_raw(iclose_raw);
  matrix.setJclose_raw(jclose_raw);
  
  matrix.setDist(dist);
  return matrix;
}


