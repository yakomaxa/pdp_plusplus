#ifndef DomainWriter_hpp
#define DomainWriter_hpp

#include <string>
#include <vector>

#include "Structure.hpp"
#include "Domain.hpp"

void listdomains(std::vector<Domain>& domains);
void listdomains(std::vector<Domain>& domains, const std::string& filename);

void writeDomainFiles(std::vector<Domain>& domains,
		      std::vector<Domain>& naive_domains,
                      Structure& s,
                      const std::string& prefix,
                      std::vector<std::string> formats,
                      const std::string& input_path = "",
                      bool include_path = false);

void writeDomainJson(std::vector<Domain>& naive_domains,
                     std::vector<Domain>& cleaned_domains,
                     const std::string& prefix);

class DomainSeq {
public:
  gemmi::Structure structure;
  std::string sequence;
  void setValues(gemmi::Structure s, std::string n){
    structure = s;
    sequence = n;
  }
};

#endif /* DomainWriter_hpp */
