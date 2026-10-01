#ifndef Tree_hpp
#define Tree_hpp

#include <string>
#include <vector>
class Node;
class RandomVariable;



class Tree {

    public:
                                    Tree(void) = delete;
                                    Tree(std::string newickString);
                                    Tree(RandomVariable* rng, double lambda, double mu, double duration);
                                   ~Tree(void);
        const std::vector<Node*>&   getDownPassSequence(void) { return downPassSequence; }
        std::string                 getNewickString(void);
        int                         getNumTaxa(void) { return numTaxa; }
        Node*                       getRoot(void) { return root; }
        void                        print(void);
        
    private:
        Node*                       addNode(void);
        void                        initializeDownPassSequence(void);
        void                        newickNode(Node* p, std::stringstream& strm);
        void                        passDown(Node* p);
        void                        showNode(Node* p, int indent);
        std::vector<std::string>    tokenizeNewickString(std::string newickString);
        Node*                       root;
        std::vector<Node*>          nodes;
        std::vector<Node*>          downPassSequence;
        int                         numTaxa;
};

#endif
