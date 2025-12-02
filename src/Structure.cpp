//
//  Structure.cpp
//  PDP_translate
//
//  Created by Koya Sakuma on 2023/03/23.
//
#include <iostream>
#include "Structure.hpp"
#include "PDPParameters.hpp"

gemmi::Structure openCoordinatefile(const std::string & infilename){
  gemmi::MaybeGzipped inputfile(infilename);
  gemmi::CoorFormat infileformat = gemmi::coor_format_from_ext(inputfile.basepath());
  if(infileformat != gemmi::CoorFormat::Unknown){
    return gemmi::read_structure(inputfile, infileformat);
  }else{
    return gemmi::read_structure(inputfile, gemmi::CoorFormat::Pdb);
  }
};

Structure::Structure(std::string filename){
    this->structure=openCoordinatefile(filename);
    this->numResidues = 0;
    for (gemmi::Model& model : this->structure.models){
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                for (gemmi::Atom &atom : residue.atoms) {
                    std::string elementname = atom.element.name();
                    if (atom.name == "CA" && elementname == "C"){
                        this->numResidues += 1;
                    }
                }
            }
        }
    }
};

std::vector<Atom> Structure::getRepresentativeAtomArray(){
  std::vector<Atom> Atoms(this->numResidues);
  int index=-1;
  int resi=0;
  int CA_flag=0;
  int chainid=0;
  int maxindex=0;
  for (gemmi::Model& model : this->structure.models){
    for (gemmi::Chain& chain : model.chains) {
      CA_flag=0;
      for (gemmi::Residue& residue : chain.residues) {
	for (gemmi::Atom &atom : residue.atoms) {
	  std::string elementname = atom.element.name();
	  if (atom.name == "CA" && elementname == "C"){
	    index += 1;
	    Atoms[index].setX(atom.pos.x);
	    Atoms[index].setY(atom.pos.y);
	    Atoms[index].setZ(atom.pos.z);
	    Atoms[index].setChain(chain.name);
	    resi=stoi(residue.seqid.str());
	    Atoms[index].setIndexOrg(resi);
	    Atoms[index].setChainId(chainid);
	    Atoms[index].setResidue(residue.name);
	    if (maxindex < index){
	      maxindex = index;
	    }	      
	    CA_flag=1;
	  }
	  if (atom.name == "CB" && elementname == "C"){
	    Atoms[index].setX(atom.pos.x);
	    Atoms[index].setY(atom.pos.y);
	    Atoms[index].setZ(atom.pos.z);
	  }
	}
      }
      if(CA_flag==1){
	this->tailofchain.push_back(index);
      }
      chainid++;
    }
  }
  PDPParameters::maxIndex = maxindex;
  return Atoms;
};
