//establishing some norms
//undefined values in the array are not used
//dist is an NdxNd matrix with distances of [i][j] for all i \neq j
//earliestTimes, latestTimes, timeTaken is an Nd matrix each 
//numCust and numNeigh are just numbers set by user
//all matrices are 0-indexed, so customer 2 is at index 1

//IF THIS CODE TAKES TOO LONG THERE ARE AREAS TO OPTIMIZE
//MOST NOTABLY WE CAN OPTIMIZE THE PATH STORAGE THROUGH HASHMAPS -- O(1) time
//AND THE SORT TO USING A MINHEAP -- though I doubt this is needed

#include "helpers.h"
#include "preprocess.h"

#include <iostream>
#include <array>
#include <vector>
#include <tuple>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>

#define INF 0x3f3f3f3f

//returns a matrix of NdxNd size for time costs (with distance) -- Nd is num customers + 1 depot
std::vector<std::vector<double>> constructTimeTaken(std::vector<std::vector<double>> dist, const int numCust, std::vector<double> earliestTimes, std::vector<double> latestTimes, std::vector<double> timeTaken)
{
	std::vector<std::vector<double>> timeMat( numCust + 1, std::vector<double> (numCust + 1,0));

	for (int i = 0; i < (numCust + 1); i++)
	{
		timeMat[i][i] = INF;
		for (int j = i + 1; j < (numCust + 1); j++)
		{
			double earliestReachj = earliestTimes[i] - timeTaken[i] - dist[i][j];
			if (earliestReachj < latestTimes[j])
			{
				timeMat[i][j] = INF;
			}
			else
			{
				timeMat[i][j] = dist[i][j] + timeTaken[i];
			}
			double earliestReachi = earliestTimes[j] - timeTaken[j] - dist[j][i];
			if (earliestReachi < latestTimes[i])
			{
				timeMat[j][i] = INF;
			}
			else
			{
				timeMat[j][i] = dist[j][i] + timeTaken[j];
			}
		}
	}

	return timeMat;
}

std::vector<std::vector<int>> constructLAneigh(std::vector<std::vector<double>> timeMat, const int numCust, const int numNeigh)
{
	std::vector<std::vector<int>> neighMat(numCust, std::vector<int>(numNeigh, 0));

	int* indices = new int[numCust];
	for (int i = 0; i < numCust; i++)
	{
		std::vector<double> cur_time = timeMat[i];
		//sort will be used, but time can be improved using a MinHeap implementation if needed
		std::iota(indices, indices + numCust, 0);
		std::sort(indices, indices + numCust,
			[cur_time](size_t i, size_t j) { return cur_time[i] < cur_time[j]; });

		for (int j = 0; j < numNeigh; j++)
		{
			neighMat[i][j] = indices[j];
		}

	}
	delete[] indices;

	return neighMat;
}

//generate P^a
std::vector<std::vector<int>> generatePa(const int numCust, const int numNeigh, std::vector<std::vector<int>> neighMat)
{
	std::unordered_set<std::vector<int>, VectorHasher> Pa_help;
	std::vector<std::vector<int>> Pa;

	//how will we format Pa's elements? First element path, second element up, third element vp

	//now loop and generate Pa
	for (int u = 0; u < numCust; u++)
	{
		//create hashmap for v
		std::unordered_map<int, int> possible_v;
		for (int v = 0; v < numCust + 1; v++)
		{
			possible_v[v] = 1;
		}
		possible_v[u] = 0;
		for (int i = 0; i < numNeigh; i++)
		{
			possible_v[neighMat[u][i]] = 0;
		}

		//generate subsets for N_u
		std::vector<std::vector<int>> neighborSetsU = subsets(neighMat[u]);

		//first loop over all subsets for N_u
		for (int i = 0; i < neighborSetsU.size(); i++)
		{
			//loop over all v possibilities
			for (int v = 0; v < numCust + 1; v++)
			{
				if (possible_v[v] == 1)
				{
					//do u 
					std::vector<int> curP = neighborSetsU[i];
					curP.push_back(u);
					curP.push_back(v);
					if (Pa_help.find(curP) == Pa_help.end())
					{
						Pa_help.insert(curP);
						Pa.push_back(curP);
					}

					//do other neighbors in neighMat
					for (int j = 0; j < numNeigh; j++)
					{
						std::vector<int> curP = neighborSetsU[i];
						curP.push_back(neighMat[u][j]);
						curP.push_back(v);
						if (Pa_help.find(curP) == Pa_help.end())
						{
							Pa_help.insert(curP);
							Pa.push_back(curP);
						}
					}
				}

			}


		}
	}
	//now lets sort Pa by its size of the neighSubset
	std::sort(Pa.begin(), Pa.end(), [](const std::vector<int>& a, const std::vector<int>& b) { return a.size() < b.size(); });

	//return P^a
	return Pa;
}


//generate efficient frontier Rp, currently we do not return tauc but this can be easily changed
std::pair<std::unordered_map<std::vector<int>, std::unordered_set<std::vector<int>, VectorHasher>, VectorHasher>, std::unordered_map<std::vector<int>, std::array<double, 3>, VectorHasher>> generateRp(std::vector<std::vector<int>> Pa, std::vector<std::vector<double>> timeMat, std::vector<double> earliestTimes, std::vector<double> latestTimes)
{
	//Rp dictionary as follows
	std::unordered_map<std::vector<int>, std::unordered_set<std::vector<int>, VectorHasher>, VectorHasher> Rp;
	//we can ignore u,v vals in Rp safely

	//first element t1,second t2, third cr
	std::unordered_map<std::vector<int>, std::array<double, 3>, VectorHasher> tauc;

	for (std::vector<int> p : Pa)
	{
		int pLen = p.size();
		int u = p[pLen - 2];
		int v = p[pLen - 1];
		std::unordered_set<std::vector<int>, VectorHasher> paths = {};
		Rp[p] = paths;

		//this LA arc is empty
		if (pLen - 2 == 0)
		{
			std::vector<int> path = { u,v };
			
			if (tauc.find(path) == tauc.end())
			{
				double cArc = timeMat[u][v];
				std::array<double, 3> vals = { std::min(earliestTimes[u],earliestTimes[v] + cArc), std::max(latestTimes[u],latestTimes[v] + cArc), cArc };
				tauc[path] = vals;

				if (vals[2] <= earliestTimes[u])
				{
					Rp[p].insert(path);
				}
			}
		}
		//this LA arc is not empty
		else
		{
			//iterating over possible w...
			for (int i = 0; i < pLen - 2; i++)
			{
				int w = p[i];
				std::vector<int> pH = p;
				pH.erase(pH.begin() + i);

				double cArc = timeMat[u][w];

				//note index must exist in hashtable
				for (std::vector<int> rMinus : Rp[pH])
				{
					std::array<double, 3> taurMinus = tauc[rMinus];
					std::vector<int> path = { u };
					path.insert(path.end(), pH.begin(), pH.end());

					//compute values if not already computed
					if (tauc.find(path) == tauc.end())
					{
						std::array<double, 3> vals = { std::min(earliestTimes[u],taurMinus[0] + cArc), std::max(latestTimes[u],taurMinus[1] + cArc), taurMinus[2] + cArc };
						tauc[path] = vals;

						if (vals[2] <= earliestTimes[u])
						{
							Rp[p].insert(path);
						}
					}
				}
			}

			//filtering possible vals...
			for (std::vector<int> r : Rp[p])
			{
				for (std::vector<int> rMinus : Rp[p])
				{
					if (r != rMinus)
					{
						std::array<double, 3> rVals = tauc[r];
						std::array<double, 3> rMinusVals = tauc[rMinus];

						if ((rVals[2] >= rMinusVals[2]) && (rVals[1] >= rMinusVals[1]) && ((rVals[0] - rVals[2]) <= (rMinusVals[0] - rMinusVals[2])))
						{
							Rp[p].erase(r);
						}
					}
				}
			}
		}
	}

	std::pair<std::unordered_map<std::vector<int>, std::unordered_set<std::vector<int>, VectorHasher>, VectorHasher>, std::unordered_map<std::vector<int>, std::array<double, 3>, VectorHasher>> Rptauc{ Rp, tauc };

	return Rptauc;

}

std::vector<std::vector<int>> generateP(std::vector<std::vector<int>> neighMat, int numCust)
{
	std::vector<std::vector<int>> P;

	for (int i = 0; i < numCust + 1; i++)
	{
		for (int j = 0; j < numCust + 2; j++)
		{
			if (i != numCust)
			{
				std::vector<int> neighSet = neighMat[i];
				if (std::find(neighSet.begin(), neighSet.end(), j) != neighSet.end()) {
					continue;
				}
				else
				{
					for (std::vector<int> neighSubset : subsets(neighMat[i]))
					{
						std::vector<int> curP = neighSubset;
						curP.push_back(i);
						curP.push_back(j);
						P.push_back(curP);
					}
				}
			}
			else
			{
				P.push_back({ i,j });
			}
		}
	}

	return P;
}