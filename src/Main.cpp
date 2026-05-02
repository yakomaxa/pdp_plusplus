#include <algorithm>

#include "Structure.hpp"
#include "Atom.hpp"
#include "GetDistanceMatrix.hpp"
#include "PDPDistanceMatrix.hpp"
#include "CutSites.hpp"
#include "CutDomain.hpp"
#include "ClusterDomains.hpp"
#include "ShortSegmentRemover.hpp"
#include "Domain.hpp"
#include "PDPParameters.hpp"
#include "SegmentComparator.hpp"
#include "DomainWriter.hpp"

enum class DomainStage { NAIVE, REMOVED };

struct OutFlags {
  bool pdb          = false;
  bool cif          = false;
  bool pml          = false;
  bool json         = false;
  bool fasta         = false;
  bool include_path = false;
};

static void resolveOriginalCoords(std::vector<Domain>& domains,
                                  std::vector<Atom>& ca) {
  for (Domain& dom : domains) {
    for (int j = 0; j < dom.getNseg(); j++) {
      Segment& seg = dom.getSegmentAtPos(j);
      seg.setFromOrg(ca[seg.getFrom()].getIndexOrg());
      seg.setToOrg(ca[seg.getTo()].getIndexOrg());
      seg.setChain(ca[seg.getTo()].getChain());
      seg.setChainId(ca[seg.getTo()].getChainId());
    }
  }
}

static void mergeAdjacentSegments(std::vector<Domain>& domains) {
  for (Domain& dom : domains) {
    std::sort(dom.getSegments().begin(), dom.getSegments().end(), SegmentComparator());
    for (int i = 0; i < dom.getNseg() - 1; i++) {
      if (PDPParameters::VERBOSE) std::cout << dom.getSegmentAtPos(i) << "\n";
      Segment& cur  = dom.getSegmentAtPos(i);
      Segment& next = dom.getSegmentAtPos(i + 1);
      if (cur.getToOrg() == next.getFromOrg() - 1 && cur.getChain() == next.getChain()) {
        cur.setToOrg(next.getToOrg());
        cur.setTo(next.getTo());
        dom.addNseg(-1);
        if (PDPParameters::VERBOSE) printf("NSEG=%i\n", dom.getNseg());
        for (int l = i + 1; l < dom.getNseg(); l++) {
          Segment& s = dom.getSegmentAtPos(l);
          Segment& t = dom.getSegmentAtPos(l + 1);
          s.setFromOrg(t.getFromOrg());
          s.setFrom(t.getFrom());
          s.setToOrg(t.getToOrg());
          s.setTo(t.getTo());
          s.setChain(t.getChain());
        }
        i--;
      }
    }
  }
}

static void printUsage(const char* prog) {
  std::cerr << "Usage: " << prog << " <input> [-o prefix] [-f pdb|cif|pml|json|fasta|all ...] [-s naive|removed] [--include-path] [-v|-vv]\n";
  std::cerr << "  -v              print domain listing to stdout\n";
  std::cerr << "  -vv             print computation details to stdout (exclusive with -v)\n";
  std::cerr << "  --include-path  include source file path in _pdp.source_path (CIF output)\n";
}

int main(int argc, char* argv[]) {
  if (argc < 2 || argv[1][0] == '-') {
    printUsage(argv[0]);
    return 1;
  }

  const std::string filename = argv[1];
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
	else if (fmt == "fasta") outflags.fasta = true;
        else if (fmt == "all")  outflags.pdb = outflags.cif = outflags.pml = outflags.json = true;
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
      if (stage == "naive")        outstage = DomainStage::NAIVE;
      else if (stage == "removed") outstage = DomainStage::REMOVED;
      else {
        std::cerr << "Unknown stage '" << stage << "'. Use naive or removed.\n";
        return 1;
      }
    } else if (arg == "--include-path") {
      outflags.include_path = true;
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      printUsage(argv[0]);
      return 1;
    }
  }

  PDPParameters::VERBOSE = (verbosity >= 2);

  if (PDPParameters::VERBOSE) printf("---------Reading structure\n");
  Structure s(filename);
  
  if (s.numResidues < 10){
    if (verbosity >= 1){
      std::cout << "PDP skipped " <<  argv[1] << " as this structure has C-alpha atoms less than 10." << std::endl;
    }
    return 1;
  }
  
  if (PDPParameters::VERBOSE) printf("---------Reading structure Done\n");

  PDPParameters param;
  param.setMAXLEN(s.numResidues);
  
  if (PDPParameters::VERBOSE) printf("---------Get Repr atoms\n");
  std::vector<Atom> ca = s.getRepresentativeAtomArray();

  if (PDPParameters::VERBOSE) printf("---------Get Repr atoms Done\n");
  if (PDPParameters::VERBOSE) printf("---------distMat creation\n");
  PDPDistanceMatrix pdpMatrix = GetDistanceMatrix().getDistanceMatrix(ca);
  if (PDPParameters::VERBOSE) printf("---------distMat creation Done\n");

  if (PDPParameters::VERBOSE) printf("---------Setting domain info\n");
  Domain dom;
  dom.setId("testDomain");
  dom.setSize((int)ca.size());
  dom.setNseg(1);
  dom.getSegmentAtPos(0).setFrom(0);
  dom.getSegmentAtPos(0).setTo((int)ca.size() - 1);
  if (PDPParameters::VERBOSE) printf("---------Setting domain info done\n");

  s.tailofchain.pop_back();
  std::vector<int> init_cutsites;
  for (int i : s.tailofchain) {
    if (PDPParameters::VERBOSE) printf("ADDING INITAIAL SITE%i\n", i);
    init_cutsites.push_back(i + 1);
  }

  if (PDPParameters::VERBOSE) printf("---------Initial splitting\n");
  CutSites cutSites;
  CutValues val;
  CutDomain cutDomain(ca, pdpMatrix, init_cutsites);
  cutDomain.cutDomain(dom, cutSites, pdpMatrix, val);
  if (PDPParameters::VERBOSE) printf("---------Initial splitting done\n");

  if (PDPParameters::VERBOSE) printf("---------Clustering domains\n");
  std::vector<Domain> domains = cutDomain.getDomains();
  domains = ClusterDomains::cluster(domains, pdpMatrix);
  if (PDPParameters::VERBOSE) printf("---------Clustering domains Done\n");

  resolveOriginalCoords(domains, ca);



  if (verbosity >= 1) listdomains(domains);

  mergeAdjacentSegments(domains);

  if (PDPParameters::VERBOSE) listdomains(domains);
  if (outflags.pml) listdomains(domains, outprefix + "_naive.pml");

  std::vector<Domain> naive_domains = domains;

  if (PDPParameters::VERBOSE) printf("---------Cleanup\n");
  ShortSegmentRemover::cleanup(domains);
  if (PDPParameters::VERBOSE) printf("---------Cleanup Done\nFINAL!!\n");
  if (PDPParameters::VERBOSE) listdomains(domains);
  if (outflags.pml) listdomains(domains, outprefix + "_removed.pml");

  std::vector<Domain>& dump_domains = (outstage == DomainStage::NAIVE) ? naive_domains : domains;
  std::vector<std::string> formats;
  if (outflags.pdb){
    formats.push_back("PDB");
  }  
  if (outflags.cif){
    formats.push_back("CIF");
      }
  if (outflags.fasta){ // define fasta later
    formats.push_back("FASTA");
  }
  if (outflags.json){
    formats.push_back("JSON");
  }
  if (formats.size()>0){
    writeDomainFiles(dump_domains, naive_domains, s, outprefix, formats, filename, outflags.include_path);
  }

  return 0;
}
