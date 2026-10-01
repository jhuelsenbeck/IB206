#include <iostream>
#include "Alignment.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"



Alignment::Alignment(RandomVariable* rng, Tree* tree, int sequenceLength) {

    numTaxa  = tree->getNumTaxa();
    numSites = sequenceLength;
    
    matrix = new int*[numTaxa];
    matrix[0] = new int[numTaxa*numSites];
    for (int i=1; i<numTaxa; i++)
        matrix[i] = matrix[i-1] + numSites;
    for (int i=0; i<numTaxa; i++)
        for (int j=0; j<numSites; j++)
            matrix[i][j] = 0;
            
    // parameters of GTR
    double alpha = 0.5;
    std::vector<double> baseFrequencies = { 0.4, 0.3, 0.2, 0.1 };
    std::vector<double> exchangeabilities = { 1.0, 5.0, 1.0, 1.0 , 5.0, 1.0 }; // AC, (AG), AT, CG, (CT), GT
            
    // initialize a rate matrix
    double Q[4][4];
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            Q[i][j] = exchangeabilities[k] * baseFrequencies[j];
            Q[j][i] = exchangeabilities[k] * baseFrequencies[i];
            k++;
            }
        }
    for (int i=0; i<4; i++)
        {
        double sum = 0.0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += Q[i][j];
            }
        Q[i][i] = -sum;
        }
        
    double averageRate = 0.0;
    for (int i=0; i<4; i++)
        {
        for (int j=0; j<4; j++)
            {
            if (i != j)
                averageRate += baseFrequencies[i] * Q[i][j];
            }
        }
    double factor = 1.0 / averageRate;
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            Q[i][j] *= factor;
        
    for (int i=0; i<4; i++)
        {
        for (int j=0; j<4; j++)
            std::cout << Q[i][j] << " ";
        std::cout << std::endl;
        }
        
    // simulate the alignment site-by-site
    taxonNames.resize(numTaxa);
    const std::vector<Node*>& dpSeq = tree->getDownPassSequence();
    for (Node* p : dpSeq)
        {
        if (p->getIsTip() == true)
            taxonNames[p->getIndex()] = p->getName();
        }
    for (int i=0; i<taxonNames.size(); i++)
        std::cout << i << " " << taxonNames[i] << std::endl;
    
    for (int site = 0; site<numSites; site++)
        {
        double r = Probability::Gamma::rv(rng, alpha, alpha);
        for (int n=(int)dpSeq.size()-1; n>=0; n--)
            {
            Node* p = dpSeq[n];
            if (p == tree->getRoot())
                {
                // draw the nucleotide state from the stationary distribution
                double u = rng->uniformRv();
                double sum = 0.0;
                int nuc = 0;
                for (int i=0; i<4; i++)
                    {
                    sum += baseFrequencies[i];
                    if (u < sum)
                        {
                        nuc = i;
                        break;
                        }
                    }
                p->setNucleotide(nuc);
                }
            else 
                {
                // simulate from the ancestor of the branch to the other end
                int currentNuc = p->getAncestor()->getNucleotide();
                double t = 0.0;
                double v = p->getBranchLength() * r;
                while (t < v)
                    {
                    double rate = -Q[currentNuc][currentNuc];
                    t += -log(rng->uniformRv()) / rate;
                    if (t < v)
                        {
                        double u = rng->uniformRv();
                        double sum = 0.0;
                        for (int j=0; j<4; j++)
                            {
                            if (j != currentNuc)
                                {
                                sum += Q[currentNuc][j] / rate;
                                if (u < sum)
                                    {
                                    currentNuc = j;
                                    break;
                                    }
                                }
                            }
                        }
                    }
                p->setNucleotide(currentNuc);
                if (p->getIsTip() == true)
                    matrix[p->getIndex()][site] = currentNuc;
                }
            }
        }

}

Alignment::~Alignment(void) {

    delete [] matrix[0];
    delete [] matrix;
}

char Alignment::convertNuc(int x) {

    if (x == 0)
        return 'A';
    else if (x == 1)
        return 'C';
    else if (x == 2)
        return 'G';
    return 'T';
}

void Alignment::print(void) {

    int longestNameLen = 0;
    for (int i=0; i<numTaxa; i++)
        {
        if (taxonNames[i].length() > longestNameLen)
            longestNameLen = (int)taxonNames[i].length();
        }
        
    std::cout << "#NEXUS" << std::endl << std::endl;
    std::cout << "begin data;" << std::endl;
    std::cout << "   dimensions ntax=" << numTaxa << " nchar=" << numSites << ";" << std::endl;
    std::cout << "   format datatype=dna;" << std::endl;
    std::cout << "   matrix" << std::endl;
    for (int i=0; i<numTaxa; i++)
        {
        std::cout << "   " << taxonNames[i] << " ";
        for (int j=0; j<longestNameLen-taxonNames[i].length(); j++)
            std::cout << " ";
        for (int j=0; j<numSites; j++)
            std::cout << convertNuc(matrix[i][j]);
        std::cout << std::endl;
        }
    std::cout << "   ;" << std::endl;
    std::cout << "end;" << std::endl;
}
