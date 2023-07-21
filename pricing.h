#pragma once
#include <utility>
#include <unordered_set>
#include <vector>
#include <unordered_map>

#include "helpers.h"

std::tuple<std::vector<int>, std::unordered_map<std::pair<int, int>, int, PairHasher>, std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher>> solve_pricing(std::vector<double> duals, std::vector<std::vector<double>> timeMat, bool is_family, std::vector<double> beta, const int num_cust, const int capacity, const int time_start, std::vector<int> demands, std::vector<double> earliestTimes, std::vector<double> latestTimes, std::unordered_map<std::vector<int>, std::array<double, 3>, VectorHasher> tauc, std::unordered_map<std::vector<int>, std::unordered_set<std::vector<int>>, VectorHasher> Rp, std::vector<std::vector<int>> P, std::vector<std::vector<int>> neighMat);

std::pair<std::vector<std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>>>, std::unordered_map<int, std::unordered_set<int>>> create_initial_nodes(const int num_cust, const int capacity, const int time_start, std::vector<int> demands, std::vector<double> earliestTimes, std::vector<double> latestTimes, std::vector<std::vector<int>> neighMat);

double compute_eta(std::vector<std::vector<int>> P, std::vector<double> duals, std::vector<double> earliestTimes, std::vector<double> latestTimes, std::unordered_map<std::vector<int>, std::array<double, 3>, VectorHasher> tauc, std::unordered_map<std::vector<int>, std::unordered_set<std::vector<int>>, VectorHasher> Rp, std::vector<int> demands);

double compute_pi(std::vector<int> p, std::vector<double> duals);

int compute_demand(std::vector<int> p, std::vector<int> demands);

std::pair<double, std::vector<int>> compute_cpt(std::vector<int> p, const int t_1, const int t_2, std::unordered_map<std::vector<int>, std::array<double, 3>, VectorHasher> tauc, std::unordered_map<std::vector<int>, std::unordered_set<std::vector<int>>, VectorHasher> Rp);