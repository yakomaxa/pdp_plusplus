//
//  Structure.cpp
//  PDP_translate
//
//  Created by Koya Sakuma on 2023/03/23.
//
#include <iostream>
#include "Structure.hpp"
#include "PDPParameters.hpp"

bool is_integer(const std::string& s) {
    if (s.empty()) return false;
    size_t start = (s[0] == '-') ? 1 : 0;    
    if (start == 1 && s.size() == 1) return false;
    return std::all_of(s.begin() + start, s.end(), [](unsigned char c) {
        return std::isdigit(c);
    });
}

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
  gemmi::CoorFormat infileformat = gemmi::coor_format_from_ext(inputfile.basepath());
  if(infileformat == gemmi::CoorFormat::Mmcif || infileformat == gemmi::CoorFormat::Mmjson){
    this->structure = gemmi::read_structure(inputfile, infileformat);
    gemmi::cif::Document doc = gemmi::make_mmcif_document(this->structure);
    this->pdp_rows = extractPdpRows(doc);
    PDPParameters::INPUT_FILETYPE="pdbx";
  }else if(infileformat == gemmi::CoorFormat::Pdb){
    this->structure = gemmi::read_structure(inputfile, gemmi::CoorFormat::Pdb);
    PDPParameters::INPUT_FILETYPE="pdb";    
  }
  if (PDPParameters::VERBOSE){
    std::cout << "[Structure.cpp] The input file type was recognized as " << PDPParameters::INPUT_FILETYPE << std::endl;
  }
 
  this->numResidues = 0;
  int n_model = 0;
  for (gemmi::Model& model : this->structure.models){
    n_model += 1;
    if (n_model > 1){
      break;
    }
    for (gemmi::Chain& chain : model.chains) {
      for (gemmi::ResidueSpan& sub : chain.subchains()){
	std::string sub_id = sub.subchain_id();
	if ( sub.length() <= PDPParameters::MIN_CHAIN_LENGTH){
	  if (PDPParameters::VERBOSE){
	    std::cout << "[Structure.cpp]: Subchain " << sub_id << " was skipped as it's shorter than or equal to " << PDPParameters::MIN_CHAIN_LENGTH << std::endl;
	  }	  
	  continue;
	}

	//int resi = 1;
	int num_calpha_of_chain = 0;
	for (gemmi::Residue& residue : sub) {
	  if (!is_integer(residue.label_seq.str())){
	    if (PDPParameters::VERBOSE){
	      std::cout << "[Structure.cpp]: As label_seq =" <<  residue.label_seq.str() << " for the residue with residue.name=" << residue.name << " label_asym_id=" << sub_id << " and auth_asym_id=" << chain.name  <<  " seems not to be integer, residue.seqid = " << residue.seqid.str() << " is copied to label_seq (i.e. label_seq_id for subchain)" << std::endl;
	    }	  	    
	    residue.label_seq = stoi(residue.seqid.str()) ;
	    //residue.label_seq = resi;
	    //resi += 1;
	  }
	  for (const gemmi::Atom &atom : residue.atoms) {
	    std::string elementname = atom.element.name();
	    if (atom.name == "CA" && elementname == "C"){
	      num_calpha_of_chain  += 1;
	    }
	  }
	}
	if (num_calpha_of_chain <= PDPParameters::MIN_CHAIN_LENGTH){
	  if (PDPParameters::VERBOSE){
	    std::cout << "[Structure.cpp]: Subchain " << sub_id << " was not considered in residue counter as it's shorter than or equal to " << PDPParameters::MIN_CHAIN_LENGTH << std::endl;
	  }
	}else{
	  this->numResidues += num_calpha_of_chain;
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
      for (gemmi::ResidueSpan& sub : chain.subchains()){
	if (sub.length() <= PDPParameters::MIN_CHAIN_LENGTH){
	  continue;
	}
	std::string sub_id = sub.subchain_id();
	int num_calpha_of_chain = 0;
	for (gemmi::Residue& residue : sub) {
	  for (const gemmi::Atom &atom : residue.atoms) {
	    std::string elementname = atom.element.name();
	    if (atom.name == "CA" && elementname == "C"){
	      num_calpha_of_chain += 1;
	    }
	  }
	}
	if (num_calpha_of_chain <= PDPParameters::MIN_CHAIN_LENGTH){
	  continue;
	}
	int chain_has_calpha = 0;
	for (gemmi::Residue& residue : sub) {
	  CA_flag=0;
	  for (const gemmi::Atom &atom : residue.atoms) {
	    std::string elementname = atom.element.name();
	    if (atom.name == "CA" && elementname == "C"){
	      index += 1;
	      Atoms[index].setX(atom.pos.x);
	      Atoms[index].setY(atom.pos.y);
	      Atoms[index].setZ(atom.pos.z);
	      if (PDPParameters::INPUT_FILETYPE == "pdb"){
		Atoms[index].setChain(chain.name);
		resi=stoi(residue.seqid.str());
	      }else{
		Atoms[index].setChain(sub_id);
		resi=stoi(residue.label_seq.str());
	      }

	      Atoms[index].setIndexOrg(resi);
	      Atoms[index].setChainId(chainid);
	      Atoms[index].setResidue(residue.name);
	      if (PDPParameters::VERBOSE){
		std::cout  << "[Structure.cpp] " << index << "-th Representative atom coordinate was taken from: " <<  Atoms[index].getChain() << " " << Atoms[index].getIndexOrg() << " "<< Atoms[index].getChainId()  << " " << Atoms[index].getResidue() <<  std::endl;
	      }
	      if (maxindex < index){
		maxindex = index;
	      }	      
	      CA_flag=1;
	      chain_has_calpha = 1;
	    }
	    if (atom.name == "CB" && elementname == "C" && CA_flag == 1){
	      if (PDPParameters::VERBOSE){
		std::cout  << "[Structure.cpp] " << index << "-th Representative atom coordinate was updated to C-beta from: "<<  Atoms[index].getChain() << " " << Atoms[index].getIndexOrg() << " "<< Atoms[index].getChainId()  << " " << Atoms[index].getResidue() <<  std::endl;
	      }
	      Atoms[index].setX(atom.pos.x);
	      Atoms[index].setY(atom.pos.y);
	      Atoms[index].setZ(atom.pos.z);
	    }
	  }
	}
	if (chain_has_calpha==1){
	  this->tailofchain.push_back(index);
	  chainid++;
	}
      }
    }
  }
  PDPParameters::maxIndex = maxindex;
  return Atoms;
};
