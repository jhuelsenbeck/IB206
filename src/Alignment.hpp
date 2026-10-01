#ifndef Alignment_hpp
#define Alignment_hpp

#include <string>
#include <vector>
class RandomVariable;
class Tree;



class Alignment {

    public:
                    Alignment(void) = delete;
                    Alignment(RandomVariable* rng, Tree* tree, int sequenceLength);
                   ~Alignment(void);
        void        print(void);

    private:
        char        convertNuc(int x);
        int         numTaxa;
        int         numSites;
        int**       matrix;
        std::vector<std::string>    taxonNames;
};

#endif 
