#include <iomanip>
#include <iostream>
#include <sstream>
#include "Msg.hpp"
#include "Node.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"



Tree::Tree(std::string newickString) {

    // break the Newick string into tokens
    std::vector<std::string> newickTokens = tokenizeNewickString(newickString);
    
    Node* p = nullptr;
    numTaxa = 0;
    bool readingBranchLength = false;
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
            readingBranchLength = true;
            }
        else 
            {
            if (readingBranchLength == false)
                {
                // add a tip as a descendant of node p
                Node* newNode = addNode();
                newNode->setName(token);
                newNode->setIsTip(true);
                p->addDescendant(newNode);
                newNode->setAncestor(p);
                p = newNode;
                numTaxa++;
                }
            else 
                {
                // reading the branch length and assigning to p
                double brlen = std::stod(token);
                p->setBranchLength(brlen);
                readingBranchLength = false;
                }
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

Tree::Tree(RandomVariable* rng, double lambda, double mu, double duration) {

    double massExtinctionRate = 0.5;
    double survivalProbability = 0.1;
    
    // initialize tree
    root = addNode();
    Node* p = addNode();
    root->addDescendant(p);
    p->setAncestor(root);
    
    // initialize the list of active nodes
    std::vector<Node*> activeNodes;
    activeNodes.push_back(p);
    
    double t = 0.0;
    int fossilNumber = 0;
    while (t < duration && activeNodes.size() > 0)
        {
        double rate = (lambda + mu) * activeNodes.size() + massExtinctionRate;
        t += -log(rng->uniformRv()) / rate;
        if (t < duration)
            {
            if (rng->uniformRv() < massExtinctionRate / rate)
                {
                // mass extinction event
                std::cout << "MASS EXTINCTION!!!!!" << std::endl;
                getchar();
                std::vector<Node*> survivors;
                for (int i=0; i<activeNodes.size(); i++)
                    {
                    if (rng->uniformRv() > survivalProbability)
                        {
                        // this lineage goes extinct!
                        fossilNumber++;
                        activeNodes[i]->setTime(t);
                        activeNodes[i]->setName("Fossil_" + std::to_string(fossilNumber));
                        activeNodes[i]->setIsTip(true);
                        }
                    else 
                        survivors.push_back(activeNodes[i]);
                    }
                activeNodes = survivors;
                }
            else 
                {
                // speciation or extinction event affecting one of the active nodes
                int whichNode = (int)(rng->uniformRv() * activeNodes.size());
                p = activeNodes[whichNode];
                p->setTime(t);
                if (rng->uniformRv() < lambda / (lambda + mu))
                    {
                    // speciation event
                    Node* new1 = addNode();
                    Node* new2 = addNode();
                    p->addDescendant(new1);
                    p->addDescendant(new2);
                    new1->setAncestor(p);
                    new2->setAncestor(p);
                    activeNodes[whichNode] = new1;
                    activeNodes.push_back(new2);
                    }
                else    
                    {
                    // extinction event
                    fossilNumber++;
                    p->setName("Fossil_" + std::to_string(fossilNumber));
                    p->setIsTip(true);
                    activeNodes[whichNode] = activeNodes[activeNodes.size()-1];
                    activeNodes.pop_back();
                    }
                }
            }
        }
        
    // clean up
    numTaxa = (int)activeNodes.size();
    for (int i=0; i<activeNodes.size(); i++)
        {
        p = activeNodes[i];
        p->setTime(duration);
        p->setIsTip(true);
        p->setName("Taxon_" + std::to_string(i+1));
        p->setIndex(i);
        }
    
    // initialzie down pass sequence
    initializeDownPassSequence();
    
    // set the branch lengths from the node times
    for (Node* p : downPassSequence)
        {
        if (p != root)
            p->setBranchLength(p->getTime() - p->getAncestor()->getTime());
        }
    
    // index tips first
    int tipIdx = (int)activeNodes.size();
    
    for (Node* p : downPassSequence)
        {
        if (p->getIsTip() == false)
            p->setIndex(tipIdx++);
        }
}

Tree::~Tree(void) {

    for (int i=0; i<nodes.size(); i++)
        delete nodes[i];
}

Node* Tree::addNode(void) {

    Node* newNode = new Node;
    nodes.push_back(newNode);
    return newNode;
}

std::string Tree::getNewickString(void) {

    std::stringstream strm;
    newickNode(root, strm);
    strm << ";";
    return strm.str();
}

void Tree::initializeDownPassSequence(void) {

    downPassSequence.clear();
    passDown(root);
}

void Tree::newickNode(Node* p, std::stringstream& strm) {

    if (p == nullptr)
        return;
        
    std::set<Node*>& pDescendants = p->getDescendants();
    
    if (p->getIsTip() == false)
        strm << "(";
    else 
        strm << p->getName() << ":" << p->getBranchLength();
        
    bool isFirst = true;
    for (Node* d : pDescendants)
        {
        if (isFirst == false)
            strm << ",";
        isFirst = false;
        newickNode(d, strm);
        }
        
    if (p->getIsTip() == false)
        strm << "):" << p->getBranchLength();
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
    std::cout << std::fixed << std::setprecision(6) << p->getBranchLength() << " ";
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
