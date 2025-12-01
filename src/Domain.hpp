#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <vector>
#include "Segment.hpp"

class Domain {
public:
    Domain();

    std::string getId() const;
    void setId(const std::string& id);

    std::vector<Segment>& getSegments();
    Segment& getSegmentAtPos(int pos);
    int getNSeg() const;

    bool operator<(const Domain& other) const;

    int getSize() const;
    void setSize(int size);
    void addSize(int size);

    int getNseg() const;
    void setNseg(int nseg);
    void addNseg(int size);

    double getScore() const;
    void setScore(double score);

    std::vector<int> getContacted() const;
    void pushbackContacted(int i);
    void removeContacted(int i);

  //    std::vector<int> contacts;
private:
    std::string id;
    int size;
    int nseg;
    double score;
    std::vector<Segment> segments;
    std::vector<int> contacted;
};

#endif // DOMAIN_HPP

