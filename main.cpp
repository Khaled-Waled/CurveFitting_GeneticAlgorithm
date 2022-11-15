#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <set>

#define pop_size 1000
//probability of crossover as a fraction of a thousand
#define prob_CO 587
//probability of crossover as a fraction of a thousand
#define prob_M 190
#define K_Tournament 2
#define DEP_FACTOR 1.5

using namespace std;

int degree;
int numGenes;
vector<pair<double, double>> dataset;
double minVal = -10, maxVal = 10;
std::mt19937 RNG(time(nullptr));

//Randomisation helper functions
long long gen(long long mn = 1, long long mx = 1000)
{
    return RNG()%(mx - mn + 1) + mn;
}

double genF(double mn, double mx)
{
    if(mn > mx)
        swap(mn, mx);
    double f = (double)RNG() /RNG.max();
    return mn + f * (mx - mn);
}

bool test(int p)
{
    int num = gen();
    return num <= p;
}

struct Chromosome
{
    vector<double> genes;
    double fitness=0;
    void initialize();
    Chromosome()
    {
        initialize();
        for (int i=0; i<numGenes; i++)
            genes.push_back(genF(minVal, maxVal));
        updateFitness();
    }
    Chromosome(vector<double> nums)
    {
        initialize();
        genes = nums;
        updateFitness();
    }
    Chromosome(vector<double> a, vector<double> b, vector<double> c)
    {
        initialize();
        for(auto x : a)
            genes.push_back(x);
        for(auto x : b)
            genes.push_back(x);
        for(auto x : c)
            genes.push_back(x);
        updateFitness();
    }

    void updateFitness()
    {
        double mse = 0;

        for(auto point : dataset)
        {
            double x = point.first;
            double polynomial = 0;
            for (int i=0; i<degree+1; i++)
                polynomial += (genes[i] * pow(x,i));
            mse+= pow((polynomial - point.second),2);
        }

        fitness = -mse;
    }

    bool operator<(const Chromosome& ch1)const
    {
        return this->fitness < ch1.fitness;
    }

//    bool operator==(const Chromosome& ch1)const
//    {
//        return this->genes == ch1.genes;
//    }

    void mutate(double t, double T)
    {
        //non uniform mutation
        for (int i=0; i<genes.size(); i++)
        {
            if(!test(prob_M))
                continue;

            double delta;
            bool low = test(500);

            if(low)
            {
                delta = genF(0, genes[i]-minVal);
                delta *= -1;
            }
            else
                delta = genF(0, maxVal - genes[i]);

            double r = genF(0, 1);
            genes[i] += delta * (1.0- pow(r, pow((1.0- t/T),DEP_FACTOR)));
        }
        updateFitness();
    }

    void print()
    {
        for (int i=0; i<genes.size()-1; i++)
        {
            cout<<genes[i]<<',';
        }
        cout<<genes[genes.size()-1]<<endl;
        cout<<"Fitness = "<<fitness<<endl<<endl;
    }

    void print2()
    {
        cout<<'[';
        for (int i=0; i<genes.size()-1; i++)
        {
            cout<<genes[i]<<',';
        }
        cout<<genes[genes.size()-1]<<"]  || Fitness = "<<fitness<<endl;

    }

    void printDesmos ()const
    {
        for (int i=0; i<genes.size(); i++)
        {
            if(i)
                cout << " + ";
            cout<<genes[i]<<"x^"<<i;
        }
        cout<<endl;
    }
};

void printDataset()
{
    for(auto point: dataset)
    {
        printf("%lf, %lf\n", point.first, point.second);
    }
}

void Chromosome::initialize()
{
    numGenes = degree+1;
}


vector<Chromosome> population;
vector<Chromosome> newGeneration;

vector<double> GetSubChromosome(const Chromosome& ch, int i, int j)
{
    vector<double> result;
    for (int k=i; k<=j; k++)
        result.push_back(ch.genes[k]);
    return result;
}


void runTestCase()
{
    population.clear();
    newGeneration.clear();

    vector<Chromosome> matingPool;
    std::mt19937 shuffleRNG(time(nullptr));

    //Initialize Population
    population.resize(pop_size);
    sort(population.begin(), population.end());

    int iterations = 20;
    for (int currIter=1; currIter<=iterations; currIter++)
    {
        newGeneration.clear();

        //Tournament selection
        matingPool.clear();
        vector<int> availableIndices(pop_size);
        for(int i=0; i<pop_size; i++)
            availableIndices[i] = i;
        shuffle(availableIndices.begin(), availableIndices.end(), shuffleRNG);
        for (int i=0; i<pop_size/2; i++) {
            matingPool.push_back(max(population[availableIndices[i << 1]],
                                     population[availableIndices[i << 1 | 1]]));
        }



        //Crossover
        for (int i=0; i<pop_size; i++)
        {
            int cutOffPoint1 = gen(0,numGenes-1);
            int cutOffPoint2 = gen(0,numGenes-1);

            if(cutOffPoint1>cutOffPoint2)
                swap(cutOffPoint1,cutOffPoint2);


            int firstInd = gen(0, (int) matingPool.size() - 1);
            int secondInd = gen(0, (int) matingPool.size() - 1);
            if(!test(prob_CO))
            {
                continue;
            }
            Chromosome c1 = Chromosome(
                    GetSubChromosome(matingPool[firstInd],0,cutOffPoint1),
                    GetSubChromosome(matingPool[secondInd],cutOffPoint1+1,cutOffPoint2),
                    GetSubChromosome(matingPool[firstInd],cutOffPoint2+1,numGenes-1));

            Chromosome c2 = Chromosome(
                    GetSubChromosome(matingPool[secondInd],0,cutOffPoint1),
                    GetSubChromosome(matingPool[firstInd],cutOffPoint1+1,cutOffPoint2),
                    GetSubChromosome(matingPool[secondInd],cutOffPoint2+1,numGenes-1));

            newGeneration.push_back(c1);
            newGeneration.push_back(c2);
        }

//        int desIter = 45;
        //mutation
        for(auto c: newGeneration)
            c.mutate(currIter, iterations);

        //Elitist replacement
        population.insert(population.end(),newGeneration.begin(), newGeneration.end());
        std::sort(population.rbegin(), population.rend());
        population[0].printDesmos();
//        for(auto& ch:population)
//        {
//            ch.printDesmos();
//        }
        population.resize(pop_size);
    }
}

void takeInput()
{
    int dataSize;
    scanf("%d %d", &dataSize,&degree);
    for (int i=0; i<dataSize; i++)
    {
        double x,y;
        scanf("%lf %lf", &x,&y);
        dataset.push_back({x,y});
    }
}


int main() {

    freopen("../input.txt", "r", stdin);
    freopen("../output.txt", "w", stdout);
    int T;
    scanf("%d", &T);
    for(int tc = 1; tc <= T; tc++)
    {

        takeInput();
        printDataset();
        runTestCase();

        //print best solution
        cout<<"Case# "<<tc<<":\n";
        population[0].printDesmos();
        cout<<endl;
    }
}

