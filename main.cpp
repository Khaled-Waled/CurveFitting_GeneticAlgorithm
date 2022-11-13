#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <cstring>
#include <sstream>

#define pop_size 10
//probability of crossover as a fraction of a thousand
#define prob_CO 550
//probability of crossover as a fraction of a thousand
#define prob_M 60
#define K_Tournament 2
#define DEP_FACTOR 1

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
            genes.push_back(gen(minVal, maxVal));
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
                polynomial+=(genes[i] * pow(x,i));
            mse+= pow((polynomial - point.second),2);
        }

        fitness = 1.0/mse;
    }

    bool operator<(const Chromosome& ch1)const
    {
        return this->fitness < ch1.fitness;
    }

    void mutate(double t, double T)
    {
        if(!test(prob_M)) return;
        //non uniform mutation
       for (int i=0; i<genes.size(); i++)
       {
           double delta;
           bool low = test(5000);

           if(low)
           {
               delta = gen(0, genes[i]-minVal);
               delta *= -1;
           }
           else
               delta = gen(0, maxVal - genes[i]);

           double r = gen(0, 1000)/1000.0;
           genes[i] += delta * (1- pow(r, pow((1.0-t/T),DEP_FACTOR)));
       }
       updateFitness();
    }

    void print()
    {
        cout<<'[';
        for (int i=0; i<genes.size()-1; i++)
        {
            cout<<genes[i]<<',';
        }
        cout<<genes[genes.size()-1]<<"]  || Fitness = "<<fitness<<endl;

    }

    void printDesmos()
    {
        for (int i=0; i<genes.size(); i++)
        {
            cout<<genes[i]<<"x^"<<i<<" +";
        }
        cout<<endl;
    }
};

void printDataset()
{
    for(auto point: dataset)
    {
        printf("(%lf, %lf)\n", point.first, point.second);
    }
}

void Chromosome::initialize()
{
        numGenes = degree+1;
        //genes.resize(numGenes);
}


vector<Chromosome> population;
vector<Chromosome> newGeneration;

vector<double> GetSubChromosome(const Chromosome& ch, int i, int j)
{
    vector<double> result;
    for (int k=i; k<j+1; k++)
        result.push_back(ch.genes[k]);
    return result;
}


void runTestCase()
{
    population.clear();
    newGeneration.clear();

    //Initialize Population
    population.resize(pop_size);

    int iterations = 100;
    for (int currIter=1; currIter<=iterations; currIter++)
    {
        newGeneration.clear();
        //Tournament selection
        for (int i=0; i<pop_size/2; i++)
        {
            vector<int> groupA;
            vector<int> groupB;
            for (int j=0; j<K_Tournament; j++)
            {
                groupA.push_back(gen(0, population.size()-1));
                groupB.push_back(gen(0, population.size()-1));
            }
            //Get best of the two groups and store them
            std::sort(groupA.rbegin(), groupA.rend());
            std::sort(groupB.rbegin(), groupB.rend());
            newGeneration.push_back(population[groupA[0]]);
            newGeneration.push_back(population[groupB[0]]);
        }

        //Crossover
        for (int i=0; i<pop_size; i+=2)
        {
            if(!test(prob_CO)) return;

            int cutOffPoint1 = gen(0,numGenes-2);
            int cutOffPoint2 = gen(0,numGenes-2);

            if(cutOffPoint1>cutOffPoint2)
                swap(cutOffPoint1,cutOffPoint2);

            Chromosome c1 = Chromosome(
                    GetSubChromosome(newGeneration[i],0,cutOffPoint1),
                    GetSubChromosome(newGeneration[i+1],cutOffPoint1+1,cutOffPoint2),
                    GetSubChromosome(newGeneration[i],cutOffPoint2+1,numGenes));

            Chromosome c2 = Chromosome(
                    GetSubChromosome(newGeneration[i+1],0,cutOffPoint1),
                    GetSubChromosome(newGeneration[i],cutOffPoint1+1,cutOffPoint2),
                    GetSubChromosome(newGeneration[i+1],cutOffPoint2+1,numGenes));

            newGeneration[i]    = c1;
            newGeneration[i+1]  = c2;
        }

        //mutation
        for(auto c: newGeneration)
            c.mutate(currIter, iterations);

        //Elitist replacement
        population.insert(population.end(),newGeneration.begin(), newGeneration.end());
        std::sort(population.rbegin(), population.rend());
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

vector<string> splitString(const string& s, char delim)
{
    stringstream raw(s);
    string temp;
    vector<string> arr;
    while(getline(raw, temp, delim))
        arr.push_back(temp);
    return arr;
}


int main() {

    freopen("../input.txt", "r", stdin);
    freopen("../output.txt", "w", stdout);
    int T;
    bool flg = true;
    scanf("%d", &T);

    while (T--)
    {
        takeInput();
        if(flg)
        {
            flg = false;
            printDataset();
        }
        runTestCase();

        //print best solution
        cout<<"Case# "<<T<<":\n";
        population[0].printDesmos();
        cout<<endl;
    }
}

