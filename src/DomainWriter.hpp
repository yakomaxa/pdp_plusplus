#ifndef DomainWriter_hpp
#define DomainWriter_hpp

#include <string>
#include <vector>

#include "Structure.hpp"
#include "Domain.hpp"

enum class OutFormat { PDB, CIF };

void listdomains(std::vector<Domain>& domains);
void listdomains(std::vector<Domain>& domains, const std::string& filename);

void writeDomainFiles(std::vector<Domain>& domains,
                      const Structure& s,
                      const std::string& prefix,
                      OutFormat format,
                      const std::string& input_path = "",
                      bool include_path = false);

void writeDomainJson(std::vector<Domain>& naive_domains,
                     std::vector<Domain>& cleaned_domains,
                     const std::string& prefix);

#endif /* DomainWriter_hpp */
