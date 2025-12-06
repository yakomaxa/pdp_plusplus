#include "Segment.hpp"
#include "PDPParameters.hpp"

// Constructors
Segment::Segment() : from_(0), to_(0), score_(0.0), chainid_(0) {}
Segment::Segment(int from, int to, double score) : from_(from), to_(to), score_(score) {}

// Getters and setters
int Segment::getFrom() const { return from_; }
void Segment::setFrom(int from) { from_ = from; }
int Segment::getTo() const { return to_; }
void Segment::setTo(int to) { to_ = to; }

int Segment::getFromOrg() const { return from_org_; }
void Segment::setFromOrg(int from) { from_org_ = from; }
int Segment::getToOrg() const { return to_org_; }
void Segment::setToOrg(int to) { to_org_ = to; }

double Segment::getScore() const { return score_; }
void Segment::setScore(double score) { score_ = score; }

std::string Segment::getChain() const { return chain_; }
void Segment::setChain(std::string chain) { chain_ = chain; }

int Segment::getChainId() const { return chainid_; }
void Segment::setChainId(int chainid) { chainid_ = chainid; }


// Overloaded operators
bool Segment::operator==(const Segment& other) const {
  return (from_ == other.from_) && (to_ == other.to_) && (score_ == other.score_) && (chainid_ == other.chainid_);
}
bool Segment::operator!=(const Segment& other) const { return !(*this == other); }

bool Segment::operator<(const Segment& other) const {
  int tmp1 = chainid_ * PDPParameters::maxIndex * 10;
  int tmp2 = other.chainid_ * PDPParameters::maxIndex * 10;
  if ( from_ + tmp1 != other.from_ + tmp2 ) {
    if (PDPParameters::VERBOSE){
      printf("A = %i, B = %i\n", from_ + tmp1,other.from_ + tmp2);
    }
    return (from_ + tmp1) < (other.from_ + tmp2);
  }
  if (PDPParameters::VERBOSE){
    printf("A = %i, B = %i\n", to_ + tmp1, other.to_ + tmp2);
  }
  return (to_ + tmp1) < (other.to_ + tmp2);
}
bool Segment::operator>(const Segment& other) const { return other < *this; }
bool Segment::operator<=(const Segment& other) const { return !(other < *this); }
bool Segment::operator>=(const Segment& other) const { return !(*this < other); }

std::ostream& operator<<(std::ostream& os, const Segment& segment) {
  //os << "Segment [from=" << segment.from_ << ", to=" << segment.to_ << ", score=" << segment.score_ << "]";
  //  os <<  segment.from_ + 1 << "-" << segment.to_ + 1;
  os << "(resi " << segment.from_org_ << "-" << segment.to_org_ << " and chain " << segment.chain_ + ")";
  return os;
}
