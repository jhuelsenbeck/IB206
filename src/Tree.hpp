#ifndef Tree_hpp
#define Tree_hpp

class Node;



class Tree {

    public:
                                    Tree(void) = delete;
                                    Tree(std::string newickString);
        void                        print(void);
        
    private:
        Node*                       addNode(void);
        void                        initializeDownPassSequence(void);
        void                        passDown(Node* p);
        void                        showNode(Node* p, int indent);
        std::vector<std::string>    tokenizeNewickString(std::string newickString);
        Node*                       root;
        std::vector<Node*>          nodes;
        std::vector<Node*>          downPassSequence;
};

#endif
