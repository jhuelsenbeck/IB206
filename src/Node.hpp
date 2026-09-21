#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>



class Node {

    public:
                            Node(void);
        void                addDescendant(Node* p) { descendants.insert(p); }
        Node*               getAncestor(void) { return ancestor; }
        std::set<Node*>&    getDescendants(void) { return descendants; }
        int                 getIndex(void) { return index; }
        bool                getIsTip(void) { return isTip; }
        std::string         getName(void) { return name; }
        void                setAncestor(Node* p) { ancestor = p; }
        void                setIndex(int x) { index = x; }
        void                setIsTip(bool tf) { isTip = tf; }
        void                setName(std::string n) { name = n; }
    
    private:
        std::set<Node*>     descendants;
        Node*               ancestor;
        std::string         name;
        int                 index;
        bool                isTip;
};

#endif
