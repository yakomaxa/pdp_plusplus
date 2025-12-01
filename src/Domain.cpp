#include "Domain.hpp"

Domain::Domain() : size(0), nseg(0), score(0.0) {}

std::string Domain::getId() const {
    return id;
}

void Domain::setId(const std::string& id) {
    this->id = id;
}

std::vector<Segment>& Domain::getSegments() {
    return segments;
}

Segment& Domain::getSegmentAtPos(int pos) {
    int size = (int)segments.size();
    while (pos >= size) {
        segments.emplace_back();
        size++;
    }
    return segments[pos];
}

int Domain::getNSeg() const {
    return nseg;
}

bool Domain::operator<(const Domain& other) const {
    if (id.empty()) {
        return true;
    } else if (other.id.empty()) {
        return false;
    } else {
        return id < other.id;
    }
}

int Domain::getSize() const {
    return size;
}

void Domain::setSize(int size) {
    this->size = size;
}

void Domain::addSize(int size) {
    this->size += size ;
}


int Domain::getNseg() const {
    return nseg;
}

void Domain::setNseg(int nseg) {
    this->nseg = nseg;
}

void Domain::addNseg(int nseg) {
    this->nseg += nseg ;
}

double Domain::getScore() const {
    return score;
}

void Domain::setScore(double score) {
    this->score = score;
}

std::vector<int> Domain::getContacted() const{
  return contacted;
}

void Domain::pushbackContacted(int i){
  contacted.push_back(i);
};

void Domain::removeContacted(int i){
  std::erase(contacted,i);
}
