//
//  Structure.cpp
//  PDP_translate
//
//  Created by Koya Sakuma on 2023/03/23.
//
#include <iostream>
#include "Structure.hpp"
#include "PDPParameters.hpp"

static std::vector<std::vector<std::string>> extractPdpRows(
    const gemmi::cif::Document& doc) {
  static const int PDP_TAG_COUNT = 9;
  std::vector<std::vector<std::string>> rows;
  if (doc.blocks.empty()) return rows;
  gemmi::cif::Column col =
      const_cast<gemmi::cif::Block&>(doc.blocks[0]).find_loop("_pdp.ordinal");
  if (!col) return rows;
  const gemmi::cif::Loop* loop = col.get_loop();
  if (!loop || (int)loop->tags.size() != PDP_TAG_COUNT) return rows;
  int w = PDP_TAG_COUNT;
  for (int i = 0; i + w <= (int)loop->values.size(); i += w)
    rows.emplace_back(loop->values.begin() + i, loop->values.begin() + i + w);
  return rows;
}

Structure::Structure(std::string filename){
  gemmi::MaybeGzipped inputfile(filename);
  gemmi::CoorFormat fmt = gemmi::coor_format_from_ext(inputfile.basepath());
  if (fmt == gemmi::CoorFormat::Mmcif) {
    gemmi::cif::Document doc = gemmi::cif::read(std::move(inputfile));
    this->pdp_rows = extractPdpRows(doc);
    this->structure = gemmi::make_structure(doc);
  } else {
    gemmi::MaybeGzipped inputfile2(filename);
    this->structure = gemmi::read_structure(
        inputfile2, fmt == gemmi::CoorFormat::Unknown ? gemmi::CoorFormat::Pdb : fmt);
  }
    this->numResidues = 0;
    int n_model = 0;
    for (gemmi::Model& model : this->structure.models){
      n_model += 1;
      if (n_model > 1){
	break;
      }
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
  int n_model = 0;
  for (gemmi::Model& model : this->structure.models){
    n_model += 1;
    if (n_model > 1){
      break;
    }
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
