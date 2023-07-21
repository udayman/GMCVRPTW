#pragma once

#include "helpers.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>


std::vector<std::vector<double>> constructTimeTaken(std::vector<std::vector<double>> dist, const int numCust, std::vector<double> earliestTimes, std::vector<double> latestTimes, std::vector<double> timeTaken);

std::vector<std::vector<int>> constructLAneigh(std::vector<std::vector<double>> timeMat, const int numCust, const int numNeigh);

std::vector<std::vector<int>> generatePa(const int numCust, const int numNeigh, std::vector<std::vector<int>> neighMat);

std::pair<std::unordered_map<std::vector<int>, std::unordered_set<std::vector<int>, VectorHasher>, VectorHasher>, std::unordered_map<std::vector<int>, std::array<double, 3>, VectorHasher>> generateRp(std::vector<std::vector<int>> Pa, std::vector<std::vector<double>> timeMat, std::vector<double> earliestTimes, std::vector<double> latestTimes);

std::vector<std::vector<int>> generateP(std::vector<std::vector<int>> neighMat, int numCust);