#include "helpers.h"

std::vector<std::vector<int>> subsets(std::vector<int> nums) {
    std::vector<std::vector<int>> subs = { {} };
    for (int num : nums) {
        int n = subs.size();
        for (int i = 0; i < n; i++) {
            subs.push_back(subs[i]);
            subs.back().push_back(num);
        }
    }
    return subs;
}