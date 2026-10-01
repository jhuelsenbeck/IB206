#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>



class Node {

    public:
                            Node(void);
        void                addDescendant(Node* p) { descendants.insert(p); }
        Node*               getAncestor(void) { return ancestor; }
        double              getBranchLength(void) { return branchLength; }
        std::set<Node*>&    getDescendants(void) { return descendants; }
        int                 getIndex(void) { return index; }
        bool                getIsTip(void) { return isTip; }
        std::string         getName(void) { return name; }
        double              getTime(void) { return time; }
        int                 getNucleotide(void) { return nucleotide; }
        void                setAncestor(Node* p) { ancestor = p; }
        void                setBranchLength(double x) { branchLength = x; }
        void                setIndex(int x) { index = x; }
        void                setIsTip(bool tf) { isTip = tf; }
        void                setName(std::string n) { name = n; }
        void                setTime(double x) { time = x; }
        void                setNucleotide(int x) { nucleotide = x; }
    
    private:
        std::set<Node*>     descendants;
        Node*               ancestor;
        std::string         name;
        int                 index;
        bool                isTip;
        double              branchLength;
        double              time;
        int                 nucleotide;
};

#endif
