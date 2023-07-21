#include <iostream>
#include <vector>
#include <array>

#include "preprocess.h"
#include "dijkstra.h"

int main() {
    //simple, lets consider 3 customers
    //distances: 1,2 = 2 , 1,3 = 1, 1,depot = 4, 2,3 = 3, 2,depot = 1, 3, depot = 6
    std::vector<std::vector<double>> dist { {0,2,1,4},{2,0,3,1}, {1,3,0,6}, {4,1,6,0} };
    //time taken: 1 = 1, 2 = 2, 3 = 1, depot = 0
    std::vector<double> timeTaken {1,2,1,0};
    //earliest Times: 1 = 10, 2 = 12, 3 = 8, 4 = 20
    std::vector<double> earliestTimes { 10,12,8,20 };
    //latest Times: 1 = 8, 2= 3, 3 = 5, 4 = 0
    std::vector<double> latestTimes { 8,3,5,0 };
    const int numCust = 3;

    std::cout << "Time Matrix" << std::endl;
    std::vector<std::vector<double>> timeMat = constructTimeTaken(dist, numCust, earliestTimes, latestTimes, timeTaken);
    for (int i = 0; i < numCust + 1; i++)
    {
        for (int j = 0; j < numCust + 1; j++)
        {
            std::cout << timeMat[i][j] << " ";
        }
        std::cout << std::endl;
    }

    //const int numNeigh = 1;
    const int numNeigh = 2;

    std::cout << "\nLA Neighbor Matrix" << std::endl;
    std::vector<std::vector<int>> LANeigh = constructLAneigh(timeMat, numCust, numNeigh);
    for (int i = 0; i < numCust; i++)
    {
        for (int j = 0; j < numNeigh; j++)
        {
            std::cout << LANeigh[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nPa Matrix" << std::endl;
    //Pa should be Neighbors, u, v
    std::vector<std::vector<int>> Pa = generatePa(numCust, numNeigh, LANeigh);
    for (int i = 0; i < Pa.size(); i++)
    {
        for (int j = 0; j < Pa[i].size(); j++)
        {
            std::cout << Pa[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nRp vals" << std::endl;
    std::pair<std::unordered_map<std::vector<int>, std::unordered_set<std::vector<int>, VectorHasher>, VectorHasher>, std::unordered_map<std::vector<int>, std::array<double, 3>, VectorHasher>> Rptauc =  generateRp(Pa, timeMat, earliestTimes, latestTimes);

    for (auto i : Rptauc.first)
    {
        std::cout << "\nNew Key" << std::endl;
        for (int j = 0; j < i.first.size(); j++)
        {
            std::cout << i.first[j] << " ";
        }
        std::cout << "\nKey's values" << std::endl;
        for (auto k : i.second)
        {
            std::cout << "\nNew path" << std::endl;
            for (auto l : k)
            {
                std::cout << l << " ";
            }
            std::cout << std::endl;
        }
    }

    int V = 9;
    Graph g(V);

    // making above shown graph
    g.addEdge(0, 1, 4);
    g.addEdge(0, 7, 8);
    g.addEdge(1, 2, 8);
    g.addEdge(1, 7, 11);
    g.addEdge(2, 3, 7);
    g.addEdge(2, 8, 2);
    g.addEdge(2, 5, 4);
    g.addEdge(3, 4, 9);
    g.addEdge(3, 5, 14);
    g.addEdge(4, 5, 10);
    g.addEdge(5, 6, 2);
    g.addEdge(6, 7, 1);
    g.addEdge(6, 8, 6);
    g.addEdge(7, 8, 7);

    //test dijstra
    std::vector<int> path = g.shortestPath(0,8);

    for (int i: path)
    {
        std::cout << i << " ";
    }
    std::cout << std::endl;

    return 0;
}