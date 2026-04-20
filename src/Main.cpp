#include <fstream>
#include <algorithm>

#include "Structure.hpp"
#define GEMMI_WRITE_IMPLEMENTATION
#include "./gemmi/to_pdb.hpp"
#include "./gemmi/to_mmcif.hpp"
#include "./gemmi/to_cif.hpp"
#include "Atom.hpp"
//#include "Chain.hpp"
#include "GetDistanceMatrix.hpp"
#include "PDPDistanceMatrix.hpp"
#include "CutSites.hpp"
#include "CutDomain.hpp"
#include "ClusterDomains.hpp"
#include "ShortSegmentRemover.hpp"
#include "Domain.hpp"
#include "PDPParameters.hpp"
//#include "LocalDomainParser.hpp"

// This function writes to standard output
static void listdomains(std::vector<Domain>& domains) {
  int i = -1;
  for (Domain& dom : domains) {
    i++;
    std::cout << "create DOMAIN" << i << ", ";
    //    std::vector<Segment>& segments = dom.getSegments();

    int flag=0;
    for (int i=0 ; i< dom.getNseg();i++) {
      if (flag>0){
	std::cout << "+" ;
      }
      std::cout << dom.getSegmentAtPos(i) ;
      flag++;
    }
    std::cout << ";" << "\n";
  }
};

// This function writes to a file
static void listdomains(std::vector<Domain>& domains, const std::string& filename) {
  std::ofstream output_file(filename);

  if (!output_file) {
    std::cerr << "Cannot open the output file: " << filename << "\n";
    return;
  }

  int i = -1;
  for (Domain& dom : domains) {
    i++;
    output_file << "create DOMAIN" << i << ", ";

    int flag=0;
    for (int i=0 ; i< dom.getNseg();i++) {
      if (flag>0){
	output_file << "+" ;
      }
      output_file << dom.getSegmentAtPos(i) ;
      flag++;
    }
    output_file << ";" << "\n";
  }
  output_file << "set grid_mode,1" << "\n";
  output_file << "set grid_slot,2,DOM*" << "\n";

  output_file.close();
};

enum class OutFormat { PDB, CIF };
enum class DomainStage { NAIVE, REMOVED };

struct OutFlags {
  bool pdb  = false;
  bool cif  = false;
  bool pml  = false;
  bool json = false;
};

static gemmi::Structure filterDomainStructure(Domain& dom,
                                              const gemmi::Structure& structure) {
  gemmi::Structure out_struct = structure;
  if (out_struct.models.size() > 1)
    out_struct.models.erase(out_struct.models.begin() + 1, out_struct.models.end());

  for (gemmi::Chain& chain : out_struct.models[0].chains) {
    std::vector<std::pair<int,int>> ranges;
    for (int si = 0; si < dom.getNseg(); si++) {
      Segment& seg = dom.getSegmentAtPos(si);
      if (seg.getChain() == chain.name)
        ranges.emplace_back(seg.getFromOrg(), seg.getToOrg());
    }
    chain.residues.erase(
      std::remove_if(chain.residues.begin(), chain.residues.end(),
        [&](const gemmi::Residue& res) {
          int seqid = std::stoi(res.seqid.str());
          for (auto& [from, to] : ranges)
            if (seqid >= from && seqid <= to) return false;
          return true;
        }),
      chain.residues.end());
  }

  out_struct.models[0].chains.erase(
    std::remove_if(out_struct.models[0].chains.begin(),
                   out_struct.models[0].chains.end(),
                   [](const gemmi::Chain& c) { return c.residues.empty(); }),
    out_struct.models[0].chains.end());

  return out_struct;
}

static void writeDomainFiles(std::vector<Domain>& domains,
                             const gemmi::Structure& structure,
                             const std::string& prefix,
                             OutFormat format) {
  for (int di = 0; di < (int)domains.size(); di++) {
    gemmi::Structure out_struct = filterDomainStructure(domains[di], structure);

    std::string ext = (format == OutFormat::CIF) ? ".cif" : ".pdb";
    std::string outname = prefix + std::to_string(di) + ext;
    std::ofstream out(outname);
    if (!out) {
      std::cerr << "Cannot open output file: " << outname << "\n";
      continue;
    }

    if (format == OutFormat::CIF) {
      gemmi::cif::Document doc = gemmi::make_mmcif_document(out_struct);
      gemmi::cif::write_cif_to_stream(out, doc);
    } else {
      gemmi::write_pdb(out_struct, out);
    }
  }
}

static bool survivedCleanup(Domain& naive, std::vector<Domain>& cleaned) {
  Segment& s0 = naive.getSegmentAtPos(0);
  for (Domain& d : cleaned)
    if (d.getSegmentAtPos(0).getChain() == s0.getChain() &&
        d.getSegmentAtPos(0).getFromOrg() == s0.getFromOrg())
      return true;
  return false;
}

static void writeDomainJson(std::vector<Domain>& naive_domains,
                            std::vector<Domain>& cleaned_domains,
                            const std::string& prefix) {
  std::string filename = prefix + ".json";
  std::ofstream out(filename);
  if (!out) {
    std::cerr << "Cannot open output file: " << filename << "\n";
    return;
  }

  out << "[\n";
  for (int di = 0; di < (int)naive_domains.size(); di++) {
    Domain& dom = naive_domains[di];
    int survived = survivedCleanup(dom, cleaned_domains) ? 1 : 0;

    out << "  {\n";
    out << "    \"domain_index\": " << di << ",\n";
    out << "    \"segments\": [\n";
    for (int si = 0; si < dom.getNseg(); si++) {
      Segment& seg = dom.getSegmentAtPos(si);
      out << "      {\"chain\": \"" << seg.getChain()
          << "\", \"from\": " << seg.getFromOrg()
          << ", \"to\": " << seg.getToOrg() << "}";
      if (si < dom.getNseg() - 1) out << ",";
      out << "\n";
    }
    out << "    ],\n";
    out << "    \"survived_cleanup\": " << survived << "\n";
    out << "  }";
    if (di < (int)naive_domains.size() - 1) out << ",";
    out << "\n";
  }
  out << "]\n";
}

static void printUsage(const char* prog) {
  std::cerr << "Usage: " << prog << " <input> [-o prefix] [-f pdb|cif|pml|json|all ...] [-s naive|removed] [-v|-vv]\n";
  std::cerr << "  -v   print domain listing to stdout\n";
  std::cerr << "  -vv  print computation details to stdout (exclusive with -v)\n";
}

int main(int argc, char *argv[]){
  if (argc < 2 || argv[1][0] == '-') {
    printUsage(argv[0]);
    return 1;
  }

  std::string filename = argv[1];
  std::string outprefix = "DOMAIN";
  OutFlags outflags;
  DomainStage outstage = DomainStage::REMOVED;
  int verbosity = 0;

  for (int i = 2; i < argc; i++) {
    std::string arg = argv[i];
    if (arg == "-v" || arg == "-vv") {
      int level = (arg == "-vv") ? 2 : 1;
      if (verbosity != 0 && verbosity != level) {
        std::cerr << "-v and -vv are exclusive.\n";
        return 1;
      }
      verbosity = level;
    } else if (arg == "-o" && i + 1 < argc) {
      outprefix = argv[++i];
    } else if (arg == "-f") {
      int consumed = 0;
      while (i + 1 < argc && argv[i+1][0] != '-') {
        std::string fmt = argv[++i];
        consumed++;
        if      (fmt == "pdb")  outflags.pdb  = true;
        else if (fmt == "cif")  outflags.cif  = true;
        else if (fmt == "pml")  outflags.pml  = true;
        else if (fmt == "json") outflags.json = true;
        else if (fmt == "all")  outflags = {true, true, true, true};
        else {
          std::cerr << "Unknown format '" << fmt << "'. Use pdb, cif, pml, json, or all.\n";
          return 1;
        }
      }
      if (consumed == 0) {
        std::cerr << "-f requires at least one format.\n";
        printUsage(argv[0]);
        return 1;
      }
    } else if (arg == "-s" && i + 1 < argc) {
      std::string stage = argv[++i];
      if (stage == "naive") {
        outstage = DomainStage::NAIVE;
      } else if (stage == "removed") {
        outstage = DomainStage::REMOVED;
      } else {
        std::cerr << "Unknown stage '" << stage << "'. Use naive or removed.\n";
        return 1;
      }
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      printUsage(argv[0]);
      return 1;
    }
  }

  PDPParameters::VERBOSE = (verbosity >= 2);

  if (PDPParameters::VERBOSE){
    printf("---------Reading structure\n");
  }
  Structure s = Structure(filename);
  if (PDPParameters::VERBOSE){
    printf("---------Reading structure Done\n");
  }
  std::vector<Domain> domains;
  PDPParameters param;
  param.setMAXLEN(s.numResidues);
  if (PDPParameters::VERBOSE){
    printf("---------Get Repr atoms\n");
  }
  std::vector<Atom> ca = s.getRepresentativeAtomArray();
  if (PDPParameters::VERBOSE){
    printf("---------Get Repr atoms Done\n");
  }
  GetDistanceMatrix distMtxCalculator;

  if (PDPParameters::VERBOSE){
    printf("---------distMat creation\n");
  }
  PDPDistanceMatrix pdpMatrix = distMtxCalculator.getDistanceMatrix(ca);
  if (PDPParameters::VERBOSE){
    printf("---------distMat creation Done\n");
  }

  Domain dom;
  //Chain c = ca[0].getGroup().getChain();
  //dom.setId("D"+c.getStructure().getPDBCode()+c.getId()+"1");
  if (PDPParameters::VERBOSE){
    printf("---------Setting domain info\n");
  }
  dom.setId("testDomain");
  dom.setSize((int)ca.size());
  dom.setNseg(1);

  dom.getSegmentAtPos(0).setFrom(0);
  dom.getSegmentAtPos(0).setTo(int(ca.size())-1);
  CutSites cutSites = CutSites();
  if (PDPParameters::VERBOSE){
    printf("---------Setting domain info done\n");
    printf("---------Initial splitting\n");
  }
  // Do the initial splitting


  std::vector<int> init_cutsites;
  // add the head residue of the chains except for the first one
  s.tailofchain.pop_back();;
  for (int i : s.tailofchain){
    if (PDPParameters::VERBOSE){
      printf("ADDING INITAIAL SITE%i\n",i);
    }
    init_cutsites.push_back(i+1);
  }
  CutValues val;
  CutDomain cutDomain(ca,pdpMatrix, init_cutsites);
  cutDomain.cutDomain(dom, cutSites,pdpMatrix, val);
  if (PDPParameters::VERBOSE){
    printf("---------Initial splitting done\n");
  }
  domains =  cutDomain.getDomains();
  // Cluster domains
  if (PDPParameters::VERBOSE){
    printf("---------Clustering domains\n");
  }
  domains = ClusterDomains::cluster(domains, pdpMatrix);
  if (PDPParameters::VERBOSE){
    printf("---------Clustering domains Done\n");
  }

  for (int i= 0 ; i < (int)domains.size(); i++){
    for (int j = 0 ; j < domains[i].getNseg();j++){
      domains[i].getSegmentAtPos(j).setFromOrg(ca[domains[i].getSegmentAtPos(j).getFrom()].getIndexOrg());
      domains[i].getSegmentAtPos(j).setToOrg(ca[domains[i].getSegmentAtPos(j).getTo()].getIndexOrg());
      domains[i].getSegmentAtPos(j).setChain(ca[domains[i].getSegmentAtPos(j).getTo()].getChain());
      domains[i].getSegmentAtPos(j).setChainId(ca[domains[i].getSegmentAtPos(j).getTo()].getChainId());
    }
  }

  if (verbosity >= 1) listdomains(domains);

  for(int j=0;j<(int)domains.size();j++) {
    std::sort(domains[j].getSegments().begin(),
    	      domains[j].getSegments().end(),SegmentComparator());
    for (int i=0;i<domains[j].getNseg()-1;i++){
      if (PDPParameters::VERBOSE) std::cout << domains[j].getSegmentAtPos(i) << std::endl;

      if((domains[j].getSegmentAtPos(i).getToOrg())==(domains[j].getSegmentAtPos(i+1).getFromOrg()-1)
	 && (domains[j].getSegmentAtPos(i).getChain()==(domains[j].getSegmentAtPos(i+1).getChain()))){

       	domains[j].getSegmentAtPos(i).setToOrg(domains[j].getSegmentAtPos(i+1).getToOrg());
	domains[j].getSegmentAtPos(i).setTo(domains[j].getSegmentAtPos(i+1).getTo());
	domains[j].addNseg(-1);
	if (PDPParameters::VERBOSE){
	  printf("NSEG=%i\n",domains[j].getNseg());
	}
	for(int l=i+1;l<domains[j].getNseg();l++){
	  domains[j].getSegmentAtPos(l).setToOrg(domains[j].getSegmentAtPos(l+1).getToOrg());
	  domains[j].getSegmentAtPos(l).setTo(domains[j].getSegmentAtPos(l+1).getTo());
	  domains[j].getSegmentAtPos(l).setFromOrg(domains[j].getSegmentAtPos(l+1).getFromOrg());
	  domains[j].getSegmentAtPos(l).setFrom(domains[j].getSegmentAtPos(l+1).getFrom());
	  domains[j].getSegmentAtPos(l).setChain(domains[j].getSegmentAtPos(l+1).getChain());
	}
	if (PDPParameters::VERBOSE){
	  printf("%i\n",i);
	}
	i--;
      }
    }
  }


  if (PDPParameters::VERBOSE){
    listdomains(domains);
  }
  if (outflags.pml) listdomains(domains, outprefix + "_naive.pml");
  std::vector<Domain> naive_domains = domains;
  // Remove short segments
  if (PDPParameters::VERBOSE){
    printf("---------Cleanup \n");
  }
  ShortSegmentRemover::cleanup(domains);
  if (PDPParameters::VERBOSE){
    printf("---------Cleanup Done\n");
    printf("FINAL!!\n");
  }
  if (PDPParameters::VERBOSE){
    listdomains(domains);
  }
  if (outflags.pml)  listdomains(domains, outprefix + "_removed.pml");

  std::vector<Domain>& dump_domains = (outstage == DomainStage::NAIVE) ? naive_domains : domains;
  if (outflags.pdb)  writeDomainFiles(dump_domains, s.structure, outprefix, OutFormat::PDB);
  if (outflags.cif)  writeDomainFiles(dump_domains, s.structure, outprefix, OutFormat::CIF);
  if (outflags.json) writeDomainJson(naive_domains, domains, outprefix);

  return 0;
}
