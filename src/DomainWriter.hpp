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
                      std::vector<std::string>&formats,
                      const std::string& input_path = "",
		      bool is_naive = false,
                      bool include_path = false);

void writeDomainJson(std::vector<Domain>& naive_domains,
		     const std::vector<std::vector<std::vector<std::string>>>& data_for_json,
		     const std::vector<std::string>& tf_table,
                     const std::string& prefix,
		     const std::vector<std::string>& formats,
		     bool is_naive)  ;

class DomainSeq {
public:
  gemmi::Structure structure;
  std::string sequence;
  std::vector<std::vector<std::string>> data;
  void setValues(gemmi::Structure s, std::string n, std::vector<std::vector<std::string>> v){
    structure = s;
    sequence = n;
    data = v;
  }
};

#endif /* DomainWriter_hpp */
