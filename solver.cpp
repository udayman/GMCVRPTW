#include <tuple>
#include <vector>
#include <algorithm>

void solver(const int num_cust, const int capacity, const int start_time, std::vector<double> earliestTimes, std::vector<std::vector<double>> timeMat)
{
	std::vector<std::tuple<int, int, double>> nodes;
	std::vector<std::tuple<int, int, std::vector<int>>> arcs;

	//initial nodes
	nodes.push_back(std::tuple<int, int, double>{num_cust, capacity, start_time});
	nodes.push_back(std::tuple<int, int, double>{num_cust + 1, 0, 0});

	for (int i = 0; i < num_cust; i++)
	{
		nodes.push_back(std::tuple<int, int, double>{i,capacity,std::min(start_time - timeMat[num_cust][i],earliestTimes[i])});
		//initial arcs
		arcs.push_back(std::tuple<int, int, std::vector<int>>{num_cust, i, {});
	}

	while (true)
	{
		while (true)
		{
			//solver (need the implementation here)

		}
	}


}