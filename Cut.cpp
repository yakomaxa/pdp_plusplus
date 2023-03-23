#include <math.h>
#include "Cut.hpp"
#include "Atom.hpp"
#include "CutValues.hpp"
#include "PDPDistanceMatrix.hpp"
#include "PDPParameters.hpp"

bool verbose = false;

int Cut::cut(std::vector<Atom> ca,Domain dom,CutValues& val,
             std::vector<std::vector<int>> dist,
             PDPDistanceMatrix pdpMatrix){
    
        int nclose = pdpMatrix.getNclose();
	printf("nclose %i\n",nclose);
        std::vector<int> iclose = pdpMatrix.getIclose();
        std::vector<int> jclose = pdpMatrix.getJclose();
    
        std::vector<int> contacts(PDPParameters::MAXLEN);
        std::vector<double> max_contacts(PDPParameters::MAXLEN);
        std::vector<double> contact_density(PDPParameters::MAXLEN);
        double average_density,x,y;

        int endsf,endst;
        int k,l,nc;
        nc=0;
        int size1,size2;
        int size1t,size2t;
        int size11,size22,size0;
        int contactsd;
        int iseg,jseg,kseg;
        int from,to,from1,to1,from2,to2,lseg;


        int site_min = -1;


        // AP add sort here..
        //qsort(dom.segment,dom.getNseg(),sizeof(struct Segment),segcmp);
        // what is going on with the segments??

        std::vector<Segment> segments = dom.getSegments();
        std::sort(segments.begin(), segments.end(), SegmentComparator());
    
        average_density = 0.0;
        size0=0;
        for(iseg=0;iseg<dom.getNseg();iseg++) {
                contactsd=1;
                size1t=0;
                size2t=0;
                for(jseg=0;jseg<iseg;jseg++)
                    size1t+=(dom.getSegmentAtPos(jseg).getFrom() - dom.getSegmentAtPos(jseg).getFrom() + 1);
                for(jseg=iseg+1;jseg<dom.getNseg();jseg++)
                    size2t+=(dom.getSegmentAtPos(jseg).getTo() - dom.getSegmentAtPos(jseg).getFrom() + 1);
                for(jseg=0;jseg<iseg;jseg++) {
                    from1 = dom.getSegmentAtPos(jseg).getFrom();
                    to1 = dom.getSegmentAtPos(jseg).getTo();
                    for(int i=from1;i<to1;i++) {
                        for(kseg=iseg+1;kseg<dom.getNseg();kseg++) {
                            from2 = dom.getSegmentAtPos(kseg).getFrom();
                            to2 = dom.getSegmentAtPos(kseg).getFrom();
                            for(int j=from2;j<to2;j++)
                                if(abs(i-j)>4) contactsd+=(dist[i][j]);
                        }
                    }
                }
                from = dom.getSegmentAtPos(iseg).getFrom();
                to = dom.getSegmentAtPos(iseg).getTo();
                for(k=from;k<to;k++) {
                    contacts[k] = contactsd;
                    /*
        if(k==392) printf("init contacts = %d\n",contacts[k]);
                     */
                    size11=size1t+(k-from+1);
                    size22=size2t+(to-k);
                    for(int i=from;i<=k;i++) {
                        for(kseg=iseg+1;kseg<dom.getNseg();kseg++) {
                            from2 = dom.getSegmentAtPos(kseg).getFrom();
                            to2 = dom.getSegmentAtPos(kseg).getTo();
                            for(int j=from2;j<=to2;j++)
                                if(abs(i-j)>4) contacts[k]+=(dist[i][j]);
                        }
                    }
                    /*
        if(k==392) printf("[from,k]x]iseg,nseg[ = %d\n",contacts[k]);
                     */
                    for(int i=from;i<=k;i++) {
                        for (int j=k+1;j<=to;j++)
                            if(abs(i-j)>4) contacts[k]+=(dist[i][j]);
                    }
                    /*
        if(k==392) printf("[from,k]x]k,to[ = %d\n",contacts[k]);
                     */
                    for(int i=k+1;i<=to;i++) {
                        for(kseg=0;kseg<iseg;kseg++) {
                            from2 = dom.getSegmentAtPos(kseg).getFrom();
                            to2 = dom.getSegmentAtPos(kseg).getTo();
                            for(int j=from2;j<to2;j++)
                                if(abs(i-j)>4) contacts[k]+=(dist[j][i]);
                        }
                    }
                    /*
        if(k==392) printf("]k,to]x]0,iseg[ = %d\n",contacts[k]);
                     */
                    size1=std::min(size11,size22);
                    size2=std::max(size11,size22);
                    x=std::min(PDPParameters::MAXSIZE,size1);
                    y=std::min(PDPParameters::MAXSIZE,size2);
                    if(x>150&&y>1.5*x) y=1.5*x;
                    else if(y>2*x) y=2*x;
                    /*
                     */
                    x=std::min(pow(x,1.3/3)+PDPParameters::RG,pow(x,1.1/3)+pow(PDPParameters::TD,1.3/3)+PDPParameters::RG);
                    y=std::min(pow(y,1.3/3)+PDPParameters::RG,pow(y,1.1/3)+pow(PDPParameters::TD,1.3/3)+PDPParameters::RG);
                    /* max_ contacts depend on the size of domains */
                    /* stella wanted comments at this point */
                    /*
                max_contacts[k] = min(MAXCONT,x*y);
                     */
                    max_contacts[k] = 10*x*y;
                    if(size1>150) max_contacts[k] = 9*x*y;
                    contact_density[k]=contacts[k]/max_contacts[k];
                    /*
                if(contact_density[k]>2.5) contact_density[k]=2.5;
                     */
                    /*
        if(first_cut)
            if(k==277 && !first_cut)
            printf("k=%d s1=%d s2=%d x=%f y=%f mc=%d c=%d cd=%f\n",k,size1,size2,x,y,max_contacts[k],contacts[k],contact_density[k]);
                     */
                    //if(verbose) System.out.println(String.format("%d    %d    %d    %f    %f    %f    %d    %f",k,size1,size2,x,y,max_contacts[k],contacts[k],contact_density[k]));
                    /*
                     */
                    if(from==0) endsf = PDPParameters::ENDSEND;
                    else endsf = PDPParameters::ENDS;
                    if(to==ca.size()-1) endst = PDPParameters::ENDSEND;
                    else endst = PDPParameters::ENDS;
                    if((contact_density[k])<val.s_min&&k>from+endsf&&k<to-endst) {
                        val.s_min = (contact_density[k]);
                        site_min=k+1;
                    }
                    if(k>from+endsf&&k<to-endst) {
                        average_density+=contact_density[k];
                        size0++;
                    }
                }
            }
            average_density/=size0;
          //  if(verbose) System.out.printf("  --- Trying to cut domain of size %d having %d segments and  average cont_density %f%n%n",dom.size,dom.getNseg(),average_density);

           // if ( verbose )
             //   for(kseg=0;kseg<dom.getNseg();kseg++)
               //     System.out.printf("  --- segment %d from %d to %d av density: %f%n",kseg,dom.getSegmentAtPos(kseg).getFrom(),dom.getSegmentAtPos(kseg).getTo(), average_density);


            if(val.first_cut) {
                /*
            for(k=from+ENDS;k<to-ENDS;k++)
                printf("%d    %s    %d    %d    %d    %d    %f\n",k,protein.res[k].type,k-from,to-k,contacts[k],max_contacts[k],contact_density[k]);
                 */
                val.AD = average_density;
//                if(verbose) System.out.printf("  --- AD=%f%n", average_density);
                /*
                 */
            }
            val.AD = average_density;

            val.s_min/=val.AD;

           // if(verbose) System.out.printf("  --- after single cut: s_min = %f site_min = %d%n",val.s_min,site_min);

            ///
            k=0;

            /* check double cuts */
            if ( verbose )
            //System.out.println("  --- checking double cuts up to: " + nclose);
            nc=0;
            for(l=0;l<nclose;l++) {

                /************ find iseg, jseg ****************/
                iseg=jseg=-1;
                for(kseg=0;kseg<dom.getNseg();kseg++) {
                    from=dom.getSegmentAtPos(kseg).getFrom();
                    to=dom.getSegmentAtPos(kseg).getTo();
                    if(from==0) endsf = PDPParameters::ENDSEND;
                    else endsf = PDPParameters::ENDS;
                    if(to==ca.size()-1) endst = PDPParameters::ENDSEND;
                    else endst = PDPParameters::ENDS;
                    if(iclose[l]>from+endsf&&iclose[l]<to-endst)
                        iseg=kseg;
                    if(jclose[l]>from+endsf&&jclose[l]<to-endst)
                        jseg=kseg;
                }
                if(iseg<0||jseg<0) continue;
                /*
        printf("l = %d iclose[l] = %d jclose[l] = %d\n",l,iclose[l],jclose[l]);
                 */
                /*********************************************/

                from=dom.getSegmentAtPos(iseg).getFrom();
                to=dom.getSegmentAtPos(iseg).getTo();
                from1=dom.getSegmentAtPos(jseg).getFrom();
                to1=dom.getSegmentAtPos(jseg).getTo();

                /************ count contacts *****************/
                contacts[nc] = 1;
                //no=0;

                /******* contacts between [0,iseg[ and ]iseg,jseg[ ********/
                for(kseg=0;kseg<iseg;kseg++)
                    for(lseg=iseg+1;lseg<jseg;lseg++)
                        for( int i=dom.getSegmentAtPos(kseg).getFrom();i<dom.getSegmentAtPos(kseg).getTo();i++)
                            for(int j=dom.getSegmentAtPos(lseg).getFrom();j<dom.getSegmentAtPos(lseg).getTo();j++) {
                                contacts[nc]+=(dist[i][j]);
                            }

                //System.out.println(String.format("[0,iseg[ - ]iseg,jseg[ : %d\n",contacts[nc]-no));

                //no=contacts[nc];
                /**********************************************************/
                /**********************************************************/
                /******* contacts between ]jseg,nseg[ and ]iseg,jseg[ ********/
                for(kseg=jseg+1;kseg<dom.getNseg();kseg++)
                    for(lseg=iseg+1;lseg<jseg;lseg++)
                        for(int i=dom.getSegmentAtPos(kseg).getFrom();i<dom.getSegmentAtPos(kseg).getTo();i++)
                            for(int j=dom.getSegmentAtPos(lseg).getFrom();j<dom.getSegmentAtPos(lseg).getTo();j++) {
                                contacts[nc]+=(dist[j][i]);
                            }
                /*
            printf("]jseg,nseg] - ]iseg,jseg[ : %d\n",contacts[nc]-no);
                 */
                //no=contacts[nc];
                /*************************************************************/
                /*************************************************************/
                /**** contacts between [from,iclose] in iseg and ]iseg,jseg[ ****/
                if(iseg==jseg) {
                    //System.out.println(" CONTACT:  " + from + " " + iclose[l] + " " + iseg + " " + jseg);
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
                    /*
            printf("iclose==jclose : %d\n",contacts[nc]-no);
                     */
                    //no=contacts[nc];
                }
                else {
                    //System.out.println(" ISEG!=JSEG " + " " + from + " " + iclose[l]);
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
                                /*
                            if(iclose[l]==33&&jclose[l]==69&&dist[i][j]) printf("%d %s %d %s %d\n",i,protein.res[i].type,j,protein.res[j].type,dist[i][j]);
                                 */
                                contacts[nc]+=(dist[j][i]);
                            }
                }
                /*******************************************************************/
                /*******************************************************************/
                /*******************************************************************/
                size11=0;
                size22=0;
                for(kseg=0;kseg<iseg;kseg++)
                    size11+=(dom.getSegmentAtPos(kseg).getTo()-dom.getSegmentAtPos(kseg).getFrom()+1);
                for(kseg=jseg+1;kseg<dom.getNseg();kseg++)
                    size11+=(dom.getSegmentAtPos(kseg).getTo()-dom.getSegmentAtPos(kseg).getFrom()+1);
                size11+=(iclose[l]-from+1);
                size11+=(to1-jclose[l]+1);
                /*
        printf("size11 = %d from = %d to1 = %d \n",size11,from,to1);
                 */
                for(kseg=iseg+1;kseg<jseg;kseg++)
                    size22+=(dom.getSegmentAtPos(kseg).getTo()-dom.getSegmentAtPos(kseg).getFrom()+1);
                if(iseg==jseg)
                    size22+=(jclose[l]-iclose[l]);
                else {
                    size22+=(jclose[l]-from1);
                    size22+=(to-iclose[l]);
                }
                /*
        printf("size22 = %d from1 = %d jclose %d iclose %d \n",size22,from1,jclose[l],iclose[l]);
                 */
                size1=std::min(size11,size22);
                size2=std::max(size11,size22);
                x=std::min(PDPParameters::MAXSIZE,size1);
                y=std::max(PDPParameters::MAXSIZE,size2);
                if(y>2*x) y=2*x;
                /*
                 */
                x=std::min(pow(x,1.3/3)+PDPParameters::RG,pow(x,1.1/3)+pow(PDPParameters::TD,1.3/3)+PDPParameters::RG);
                y=std::min(pow(y,1.3/3)+PDPParameters::RG,pow(y,1.1/3)+pow(PDPParameters::TD,1.3/3)+PDPParameters::RG);
                /*
            max_contacts[nc] = min(MAXCONT,x*y);
                 */
                max_contacts[nc] = x*y*10;
                if(size1>150) max_contacts[k] = 9*x*y;
                contact_density[nc]=contacts[nc]/max_contacts[nc];
                /*
        if(first_cut)
                 */
                //if ( verbose)
                //    System.out.println("double cut" + l + " " + ca[iclose[l]].getGroup().getType() + " " + iclose[l] + " " + jclose[l] + " " + contacts[nc]+
                //            " " + max_contacts[nc] + " " + x+ " " + y + " " + size11 + " " + size22 + " " + contact_density[nc] + " " + contact_density[nc]/val.AD );
                //if(verbose) System.out.println(String.format(" double cut: %d    %s %d %d c=%d mc=%d x=%f y=%f s1=%d s2=%d cd=%f cd/ad=%f\n",l,ca[iclose[l]].getGroup().getType(),iclose[l],jclose[l],contacts[nc],max_contacts[nc],x,y,size11,size22,contact_density[nc],contact_density[nc]/val.AD));
                if((contact_density[nc]/val.AD+PDPParameters::DBL)<val.s_min&&contact_density[nc]/val.AD+PDPParameters::DBL<PDPParameters::CUT_OFF_VALUE2) {
                    /*
        printf("________________\n");
                     */
                    val.s_min = (contact_density[nc]/val.AD)+PDPParameters::DBL;
                    site_min=iclose[l];
                    val.site2=jclose[l];
                }

                nc++;
                if ( nc >= PDPParameters::MAXSIZE)
                    nc = PDPParameters::MAXSIZE-1;
            }
            val.first_cut=false;
    if(verbose){
        //System.out.printf("  --- E ... at the end of cut: s_min %f CUTOFF %f site_min %d *site2 %d%n",val.s_min,PDPParameters::CUT_OFF_VALUE,site_min,val.site2);
    }
            if(val.s_min> PDPParameters::CUT_OFF_VALUE) return -1;

            return(site_min);
}


