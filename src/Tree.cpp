#include <iostream>
#include "Msg.hpp"
#include "Node.hpp"
#include "Tree.hpp"



Tree::Tree(std::string newickString) {

    // break the Newick string into tokens
    std::vector<std::string> newickTokens = tokenizeNewickString(newickString);
    
    Node* p = nullptr;
    for (int i=0; i<newickTokens.size(); i++)
        {
        std::string token = newickTokens[i];
        
        if (token == "(")
            {
            // add a new node as a descendant of p
            Node* newNode = addNode();
            if (p == nullptr)
                {
                root = newNode;
                }
            else 
                {
                p->addDescendant(newNode);
                newNode->setAncestor(p);
                }
            p = newNode;
            }
        else if (token == ")" || token == ",")
            {
            // move to the ancestor of p
            if (p == nullptr)
                Msg::error("No node to move to");
            p = p->getAncestor();
            }
        else if (token == ";")
            {
            // we should end at the root of the tree
            if (p != root)
                Msg::error("Expected reading of Newick string to end at the root");
            }
        else if (token == ":")
            {
            // the next token is a branch length
            //double brlen = std::stod(token);
            }
        else 
            {
            // add a tip as a descendant of node p
            Node* newNode = addNode();
            newNode->setName(token);
            newNode->setIsTip(true);
            p->addDescendant(newNode);
            newNode->setAncestor(p);
            p = newNode;
            }
        }
        
    // initialzie down pass sequence
    initializeDownPassSequence();
    
    // index tips first
    int tipIdx = 0;
    for (Node* p : downPassSequence)
        {
        if (p->getIsTip() == true)
            p->setIndex(tipIdx++);
        }
    
    for (Node* p : downPassSequence)
        {
        if (p->getIsTip() == false)
            p->setIndex(tipIdx++);
        }
}

Node* Tree::addNode(void) {

    Node* newNode = new Node;
    nodes.push_back(newNode);
    return newNode;
}

void Tree::initializeDownPassSequence(void) {

    downPassSequence.clear();
    passDown(root);
}

void Tree::passDown(Node* p) {

    if (p == nullptr)
        return;
        
    std::set<Node*>& pDescendants = p->getDescendants();
    for (Node* d : pDescendants)
        passDown(d);
        
    downPassSequence.push_back(p);
}

void Tree::print(void) {

    showNode(root, 0);
    
    std::cout << "Postorder: ";
    for (int i=0; i<downPassSequence.size(); i++)
        std::cout << downPassSequence[i]->getIndex() << " ";
    std::cout << std::endl;
}

void Tree::showNode(Node* p, int indent) {

    if (p == nullptr)
        return;
        
    std::set<Node*>& pDescendants = p->getDescendants();
    
    for (int i=0; i<indent; i++)
        std::cout << " ";
    std::cout << p->getIndex() << " ( ";
    if (p->getAncestor() == nullptr)
        std::cout << "a_N ";
    else 
        std::cout << "a_" << p->getAncestor()->getIndex() << " ";
    for (Node* d : pDescendants)
        std::cout << d->getIndex() << " ";
    std::cout << ") ";
    std::cout << p->getName() << " ";
    std::cout << std::endl;
    
    for (Node* d : pDescendants)
        showNode(d, indent + 3);
}

std::vector<std::string> Tree::tokenizeNewickString(std::string newickString) {

    std::vector<std::string> tokens;
    
    std::string token = "";
    for (int i=0; i<newickString.length(); i++)
        {
        char c = newickString[i];
        if (c == '(' || c == ')' || c == ',' || c == ';' || c == ':')
            {
            if (token != "")
                {
                tokens.push_back(token);
                token = "";
                }
            tokens.push_back( std::string(1,c) );
            }
        else 
            {
            token += std::string(1,c);
            }
        }
        
    return tokens;
}
