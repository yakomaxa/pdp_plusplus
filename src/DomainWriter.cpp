#include <fstream>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <iostream>

#define GEMMI_WRITE_IMPLEMENTATION
#include "./gemmi/to_pdb.hpp"
#include "./gemmi/to_mmcif.hpp"
#include "./gemmi/to_cif.hpp"

#include "PDPParameters.hpp"
#include "DomainWriter.hpp"
#include "Segment.hpp"
#include "SegmentComparator.hpp"

// TODO: set from git tag at build time
static const char* SOFTWARE_NAME    = "PDP++";
static const char* SOFTWARE_VERSION = "v1.0";
//static const char* REFERCENCE = "Under Preparation";
static const char* REFERCENCE = "?";

static std::string getSourceId(const gemmi::Structure& structure) {
  const std::string& entry_id = structure.get_info("_entry.id");
  return entry_id.empty() ? structure.name : entry_id;
}

static std::string makeTimestamp() {
  auto now = std::chrono::system_clock::now();
  auto tt  = std::chrono::system_clock::to_time_t(now);
  std::tm tm_utc{};
  gmtime_r(&tt, &tm_utc);
  char buf[32];
  std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm_utc);
  return std::string(buf);
}

static DomainSeq filterDomainStructure(Domain& dom,
                                              gemmi::Structure& structure) {

  gemmi::Structure out_struct = structure;
  std::vector<std::vector<std::string>> data_for_json;
  out_struct.models.clear();
  DomainSeq out_domain;
  if (structure.models.empty()){
    out_domain.setValues(out_struct,"",data_for_json);
    return out_domain;
  }
  std::string fasta;
  int n_model = 0;  
  for (gemmi::Model& model : structure.models){
    n_model += 1;
    if (n_model > 1){
      break;
    }
    gemmi::Model out_model = model;
    out_model.chains.clear();
    for (int si = 0; si < dom.getNseg(); si++) {
      Segment& seg = dom.getSegmentAtPos(si);      
      for (gemmi::Chain& chain : model.chains) {
	gemmi::Chain out_chain = chain;
	if (PDPParameters::INPUT_FILETYPE == "pdb"){
	  out_chain.residues.clear();
	  out_chain.name = chain.name;
	  if (seg.getChain() == chain.name){
	    for (gemmi::Residue& res : chain.residues) {
	      int seqid = stoi(res.seqid.str());
	      if (seqid >= seg.getFromOrg() && seqid <= seg.getToOrg()){
		out_chain.residues.push_back(res);
	      }
	    }
	  }
	  gemmi::ResidueSpan polymer=out_chain.get_polymer();
	  if (!polymer.empty()){
	    std::string seq = gemmi::make_one_letter_sequence(polymer);
	    std::vector<std::string> data_tmp  =  {std::to_string(si),
						   "auth_asym_id", chain.name,
						   "auth_seq_id", std::to_string(seg.getFromOrg()), std::to_string(seg.getToOrg()), seq};
	    data_for_json.push_back(data_tmp);
	    fasta += ">auth_asym_id=" + chain.name + " segment=" + std::to_string(si) + " start:auth_seq_id=" + std::to_string(seg.getFromOrg()) + " end:auth_seq_id=" + std::to_string(seg.getToOrg())  + "\n" + seq + "\n";
	  if (!out_chain.residues.empty()){
	    out_model.chains.push_back(std::move(out_chain));
	  }

	  }
	}else{
	  for (gemmi::ResidueSpan& sub : chain.subchains()){
	    out_chain.residues.clear();
	    std::string sub_id = sub.subchain_id();
	    if (seg.getChain() == sub_id){
	      for (gemmi::Residue& res : sub) {
		int seqid = std::stoi(res.label_seq.str());
		if (seqid >= seg.getFromOrg() && seqid <= seg.getToOrg()){
		  out_chain.residues.push_back(res);
		}
	      }
	    }
	    gemmi::ResidueSpan polymer=out_chain.get_polymer();
	    if (!polymer.empty()){
	      std::string seq = gemmi::make_one_letter_sequence(out_chain.get_polymer());
	      std::vector<std::string> data_tmp  =  {std::to_string(si),
						     "label_asym_id", sub_id,
						     "label_seq_id", std::to_string(seg.getFromOrg()), std::to_string(seg.getToOrg()), seq};
	      data_for_json.push_back(data_tmp);
	      fasta += ">label_asym_id=" + sub_id + " segment=" + std::to_string(si) + " start:label_seq_id=" + std::to_string(seg.getFromOrg()) + " end:label_seq_id=" + std::to_string(seg.getToOrg())  + "\n" + seq + "\n";
	    }
	    if (!out_chain.residues.empty()){
	      out_model.chains.push_back(std::move(out_chain));
	    }
	  }
	}
      }
    }
    out_struct.models.push_back(std::move(out_model));
  }
  // set domain container
  out_domain.setValues(out_struct, fasta, data_for_json);
  return out_domain;
}


void listdomains(std::vector<Domain>& domains) {
  int i = -1;
  for (Domain& dom : domains) {
    i++;
    std::cout << "create DOMAIN" << i << ", ";

    int flag = 0;
    for (int j = 0; j < dom.getNseg(); j++) {
      if (flag > 0) std::cout << "+";
      std::cout << dom.getSegmentAtPos(j);
      flag++;
    }
    std::cout << ";\n";
  }
}

void listdomains(std::vector<Domain>& domains, const std::string& filename) {
  std::ofstream output_file(filename);
  if (!output_file) {
    std::cerr << "Cannot open the output file: " << filename << "\n";
    return;
  }

  int i = -1;
  for (Domain& dom : domains) {
    i++;
    output_file << "create DOMAIN" << i << ", ";

    int flag = 0;
    for (int j = 0; j < dom.getNseg(); j++) {
      if (flag > 0) output_file << "+";
      output_file << dom.getSegmentAtPos(j);
      flag++;
    }
    output_file << ";\n";
  }
  output_file << "set grid_mode,1\n";
  output_file << "set grid_slot,2,DOM*\n";
}

static const std::vector<std::string> PDP_TAGS = {
  "ordinal", "source_id", "domain_number", "domain_count",
  "execution_time", "software_name", "software_version", "reference",
  "source_path"
};

void writeDomainFiles(std::vector<Domain>& domains,
		      std::vector<Domain>& naive_domains,
                      Structure& s,
                      const std::string& prefix,
                      std::vector<std::string>& OutFormats,
                      const std::string& input_path,
		      bool is_naive,
                      bool include_path
		      ) {
  const std::string timestamp = makeTimestamp();
  const std::string source_id = getSourceId(s.structure);
  const int total = (int)naive_domains.size();

  auto survivedCleanup = [&](Domain& naive) {
    Segment& s0 = naive.getSegmentAtPos(0);
    for (Domain& d : domains)
      if (d.getSegmentAtPos(0).getChain() == s0.getChain() &&
          d.getSegmentAtPos(0).getFromOrg() == s0.getFromOrg())
        return true;
    return false;
  };

  std::vector<std::vector<std::vector<std::string>>> data_for_json;
  std::vector<std::string> tf_table;
  
  int json_flag = 0; 
  for (int di = 0; di < total; di++) {
    DomainSeq output = filterDomainStructure(naive_domains[di], s.structure);
    std::string survived = survivedCleanup(naive_domains[di]) ? "True" : "False";
    gemmi::Structure out_struct =  output.structure;
    std::string out_sequence =  output.sequence;
    data_for_json.push_back(output.data);
    tf_table.push_back(survived);
    //std::cout << out_sequence << std::endl;
    
    const int domain_num = di + 1;
    if (survived == "False" && !is_naive){
          continue;
    }
    for (std::string format: OutFormats){
      std::string ext = (format == "CIF") ? ".cif" :
	(format == "PDB") ? ".pdb" :
	(format == "JSON") ? ".json" :
	(format == "FASTA") ? ".fasta" :
	".txt";
      std::string outname = prefix + "_domain" + std::to_string(domain_num) + ext;      
      if (format == "CIF") {
	std::ofstream out(outname);
	if (!out) {
	  std::cerr << "Cannot open output file: " << outname << "\n";
	  continue;
	}
	gemmi::cif::Document doc;
	doc.blocks.resize(1);
	gemmi::cif::Block& block = doc.blocks[0];
	block.name = source_id + "_domain" + std::to_string(domain_num);

	const auto& existing_rows = s.pdp_rows;
	int ordinal = (int)existing_rows.size() + 1;
	
	std::vector<std::string> new_row = {
	  std::to_string(ordinal),
	  source_id,
	  std::to_string(domain_num),
	  std::to_string(total),
	  gemmi::cif::quote(timestamp),
	  SOFTWARE_NAME,
	  SOFTWARE_VERSION,
	  REFERCENCE,
	  include_path ? gemmi::cif::quote(input_path) : "?"
	};
	
	gemmi::cif::Loop& pdp_loop = block.init_loop("_pdp.", PDP_TAGS);
	for (auto& row : existing_rows){
	  pdp_loop.add_row(row);
	}
	pdp_loop.add_row(new_row);
	
	gemmi::MmcifOutputGroups groups(false);
	groups.atoms     = true;
	groups.group_pdb = true;
	gemmi::update_mmcif_block(out_struct, block, groups);
	
	gemmi::cif::write_cif_to_stream(out, doc);
      }
      
      if (format == "PDB") {
	std::ofstream out(outname);
	if (!out) {
	  std::cerr << "Cannot open output file: " << outname << "\n";
	  continue;
	}
	gemmi::write_minimal_pdb(out_struct, out);
      }
      
      if (format == "FASTA") {
	std::ofstream out(outname);
	if (!out) {
	  std::cerr << "Cannot open output file: " << outname << "\n";
	  continue;
	}
	out << out_sequence;
      }

      if (format == "JSON") {
	json_flag = 1;
      }
    }
  }
  if (json_flag == 1) {
    writeDomainJson(naive_domains, data_for_json, tf_table, prefix, OutFormats, is_naive);
  }
}

void writeDomainJson(std::vector<Domain>& naive_domains,
		     const std::vector<std::vector<std::vector<std::string>>>& data_for_json,
		     const std::vector<std::string>& tf_table,
                     const std::string& prefix,
		     const std::vector<std::string>& OutFormats,
		     bool is_naive) {
  std::string filename = prefix + ".json";
  std::ofstream out(filename);
  if (!out) {
    std::cerr << "Cannot open output file: " << filename << "\n";
    return;
  }

  std::vector<std::string> exts;
  for (std::string format : OutFormats){
    if (format == "CIF"){
      exts.push_back("cif");
    }
    if (format == "PDB"){
      exts.push_back("pdb");
    }
    if (format == "FASTA"){
      exts.push_back("fasta");
    }  
  }
    
  out << "{\"domain_data\": [\n";
  for (int di = 0; di < (int)naive_domains.size(); di++) {
    Domain& dom = naive_domains[di];
    std::string survived = tf_table[di];
    
    out << "  {\n";
    out << "    \"domain_index\": " << di+1 << ",\n";
    if ( exts.size() >= 1 && (survived == "True" || is_naive)){
      out << "    \"files\": {";
      for (int ei = 0 ; ei < (int)exts.size(); ei++){
	  out << "\"" << exts[ei] << "\": " << "\"" << prefix << "_domain" << std::to_string(di+1) << "." << exts[ei] << "\"";
	  if (ei < (int)exts.size() - 1){
	    out << ",";
	  }
	}
	out << "},\n";
    }
    out << "    \"segments\": [\n";
    int domain_size = 0;
    std::vector<std::vector<std::string>> data_dom = data_for_json[di];
    for (int si = 0; si < dom.getNseg(); si++) {
      std::vector<std::string>  data = data_dom[si];
      if (PDPParameters::VERBOSE){
	std::cout << "[DomainWriter.cpp] domainID="<< di << " segmentID=" << data[0] << " asym_id=("<< data[1] << " "<< data[2] << ") seq_id=("<< data[3] << " from="<< data[4] << " to="<< data[5]  << ") sequence="<< data[6] << std::endl;
      }
      Segment& seg = dom.getSegmentAtPos(si);
      domain_size +=  seg.getToOrg() - seg.getFromOrg() + 1;

      out << "      {\"segment_id\": " << data[0] 
	  << ", \"chain\": \"" << seg.getChain() << "\""
          << ", \"from\": " << seg.getFromOrg()
          << ", \"to\": " << seg.getToOrg()
	  << ", \"sequence\": \"" << data[6] << "\""
	  << ", \"seq_id\": {\"type\": \"" << data[3] << "\", \"from\": " << data[4] << ", \"to\": " << data[5] << "}"
	  << ", \"asym_id\": {\"type\": \"" << data[1] << "\", \"id\": \"" << data[2] << "\"}"
	  << "}";
      if (si < dom.getNseg() - 1) out << ",";
      out << "\n";
    }
    out << "    ],\n";
    out << "    \"n_residue\": " << domain_size << ",\n";
    out << "    \"survived_cleanup\": \"" << survived << "\"\n";
    out << "  }";
    if (di < (int)naive_domains.size() - 1) out << ",";
    out << "\n";
  }
  out << "]}\n";
}
