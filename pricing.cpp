#include "helpers.h"
#include "dijkstra.h"
#include "pricing.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <tuple>
#include <algorithm>
#include <iterator>

#define INF 0x3f3f3f3f

std::tuple<std::vector<int>, std::unordered_map<std::pair<int, int>, int, PairHasher>, std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher>> solve_pricing(std::vector<double> duals, std::vector<std::vector<double>> timeMat, bool is_family, std::vector<double> beta, const int num_cust, const int capacity, const int time_start, std::vector<int> demands, std::vector<double> earliestTimes, std::vector<double> latestTimes, std::unordered_map<std::vector<int>, std::array<double, 3>, VectorHasher> tauc, std::unordered_map<std::vector<int>, std::unordered_set<std::vector<int>>, VectorHasher> Rp, std::vector<std::vector<int>> P, std::vector<std::vector<int>> neighMat)
{
    //define initial nodes and edges
    std::pair<std::vector<std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>>>, std::unordered_map<int, std::unordered_set<int>>> nodesnedges = create_initial_nodes(num_cust, capacity, time_start, demands, earliestTimes, latestTimes, neighMat);
    std::vector<std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>>> nodes = nodesnedges.first;
    std::unordered_map<int, std::unordered_set<int>> edges = nodesnedges.second;

    //compute eta
    double eta = compute_eta(P, duals, earliestTimes, latestTimes, tauc, Rp, demands);

    //define chat_fg, d_fg, and p_fg used for edges
    std::unordered_map<std::pair<int, int>, int, PairHasher> demand_fg;
    std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher> path_fg;
    std::unordered_map<std::pair<int, int>, double, PairHasher> chat_fg;

    while (true)
    {
        //define graph
        Graph g(nodes.size());

        //compute associated edge terms. skip if already computed
        for (auto keyvalue : edges)
        {
            int src = keyvalue.first;
            std::unordered_set<int> recipients = keyvalue.second;

            int u_p = std::get<0>(nodes[src]);

            for (int end : recipients)
            {
                if (chat_fg.find(std::pair<int, int>{src, end}) != chat_fg.end())
                {
                    g.addEdge(src, end, chat_fg[std::pair<int, int>{src, end}]);
                }

                //first defining all paths
                int v_p = std::get<0>(nodes[end]);

                std::vector<std::vector<int>> P_fg;

                std::vector<std::vector<int>> neighsubsets;

                if (u_p == num_cust)
                {
                    neighsubsets.push_back({});
                }
                else
                {
                    neighsubsets = subsets(neighMat[u_p]);
                }

                for (std::vector<int> neighbor_set : neighsubsets)
                {
                    if (is_family == true)
                    {
                        bool skip = false;
                        if (!(beta[u_p] < beta[v_p]))
                        {
                            continue;
                        }
                        for (int neigh_cust : neighbor_set)
                        {
                            if (!(beta[u_p] < beta[neigh_cust]) || !(beta[v_p] > beta[neigh_cust]))
                            {
                                skip = true;
                                continue;
                            }
                        }
                        if (skip == true)
                        {
                            continue;
                        }
                    }

                    neighbor_set.push_back(u_p);
                    neighbor_set.push_back(v_p);

                    int d_p = compute_demand(neighbor_set, demands);
                    if ((std::get<2>(nodes[src]) - d_p >= std::get<1>(nodes[end])) && (std::get<1>(nodes[src]) - d_p <= std::get<2>(nodes[end])))
                    {
                        //first check (Mfminus intersect Np+ and vp must be zero)
                        std::vector<int> Np_plus_v = neighbor_set;

                        std::vector<int> Np_plus = Np_plus_v;

                        Np_plus.pop_back();

                        std::unordered_set<int> M_f_ms = std::get<5>(nodes[src]);

                        std::vector<int> M_f_m(M_f_ms.begin(), M_f_ms.end());

                        std::sort(Np_plus_v.begin(), Np_plus_v.end());

                        std::sort(M_f_m.begin(), M_f_m.end());

                        std::vector<int> overlap;

                        std::set_intersection(Np_plus_v.begin(), Np_plus_v.end(), M_f_m.begin(), M_f_m.end(), std::back_inserter(overlap));

                        //no overlap -- continue to second check that Mg minus must be a subset
                        if (overlap.size() == 0)
                        {
                            std::unordered_set<int> M_f_ps = std::get<6>(nodes[src]);

                            std::vector<int> M_f_p(M_f_ps.begin(), M_f_ps.end());

                            std::sort(M_f_p.begin(), M_f_p.end());

                            std::sort(Np_plus.begin(), Np_plus.end());

                            std::vector<int> overlap;

                            std::set_union(Np_plus.begin(), Np_plus.end(), M_f_p.begin(), M_f_p.end(), std::back_inserter(overlap));

                            std::unordered_set<int> M_g_ms = std::get<5>(nodes[end]);

                            std::vector<int> M_g_m(M_g_ms.begin(), M_g_ms.end());

                            std::sort(M_g_m.begin(), M_g_m.end());

                            //final check!! (is Mfm and Np+'s union a subset?)
                            if (std::includes(M_g_m.begin(), M_g_m.end(), overlap.begin(), overlap.end()))
                            {
                                std::unordered_set<int> M_g_ps = std::get<6>(nodes[end]);

                                std::vector<int> M_g_p(M_g_ps.begin(), M_g_ps.end());

                                std::sort(M_g_p.begin(), M_g_p.end());

                                std::vector<int> overlap;

                                std::set_union(M_f_m.begin(), M_f_m.end(), Np_plus.begin(), Np_plus.end(), overlap);

                                if (std::includes(overlap.begin(), overlap.end(), M_g_p.begin(), M_g_p.end()))
                                {
                                    std::vector<int> p = neighbor_set;
                                    p.push_back(u_p);
                                    p.push_back(v_p);
                                    P_fg.push_back(p);
                                }
                            }
                        }
                    }
                }

                //can't be negative for edge
                double best_chat_fg = INF;

                for (std::vector<int> p : P_fg)
                {
                    int d_p = compute_demand(p, demands);
                    int d_fgp;

                    //if the ending customer is sink
                    if (std::get<0>(nodes[end]) == num_cust + 1)
                    {

                        int d_fgp = std::max(std::get<1>(nodes[src]), d_p);

                    }
                    //otherwise default definition
                    else
                    {
                        int d_fgp = std::max(std::get<1>(nodes[src]) - std::get<2>(nodes[end]), d_p);

                    }
                    std::pair<double, std::vector<int>> cptnpath = compute_cpt(p, std::get<4>(nodes[src]), std::get<3>(nodes[end]), tauc, Rp);
                    double cpt = cptnpath.first;
                    double chat_fgp = eta * d_fgp - compute_pi(p, duals) + cpt;

                    if (chat_fgp < best_chat_fg)
                    {
                        std::pair<int, int> key{ src,end };
                        demand_fg[key] = d_fgp;
                        path_fg[key] = cptnpath.second;
                        chat_fg[key] = chat_fgp;

                        best_chat_fg = chat_fgp;
                    }
                }

                if (best_chat_fg != INF)
                {
                    g.addEdge(src, end, best_chat_fg);
                }
                else
                {
                    edges[src].erase(end);
                }
            }
        }

        //getting shortest path (source and sink will never change in terms of split)
        std::vector<int> path = g.shortestPath(0, num_cust + 1);

        //check cases and redefine edges and nodes + repeat

        //case 1: eta correctness violated
        int total_demand_path = 0;
        for (int i = 0; i < path.size() - 1; i++) {
            total_demand_path += demand_fg[std::pair<int, int>{path[i], path[i + 1]}];
        }

        if (total_demand_path < capacity)
        {
            int demand_so_far = capacity - demand_fg[std::pair<int, int>{0, 1}];
            for (int i = 1; i < path.size() - 1; i++)
            {
                if (demand_so_far <= std::get<2>(nodes[path[i]]) && demand_so_far > std::get<1>(nodes[path[i]]))
                {
                    std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>> new_node_one = nodes[path[i]];
                    std::get<2>(new_node_one) = demand_so_far - 1;
                    std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>> new_node_two = nodes[path[i]];
                    std::get<1>(new_node_two) = demand_so_far;

                    //add nodes
                    nodes[path[i]] = new_node_one;
                    nodes.push_back(new_node_two);

                    //add edges (edges for 1 added by default) and remove from chat
                    edges[nodes.size() - 1] = edges[path[i]];

                    for (int potend : edges[path[i]])
                    {
                        chat_fg.erase(std::pair<int, int>{path[i], potend});
                    }

                    for (auto edge : edges)
                    {
                        std::unordered_set<int> endings = edge.second;
                        if (endings.find(path[i]) != endings.end())
                        {
                            edges[edge.first].insert(nodes.size() - 1);
                            chat_fg.erase(std::pair<int, int>{edge.first, path[i]});
                        }
                    }
                }
                demand_so_far -= demand_fg[std::pair<int, int>{i, i + 1}];
            }
            continue;
        }

        //case 2: demand exceeds capacity
        if (total_demand_path > capacity)
        {
            int demand_so_far = capacity;
            for (int i = 1; i < path.size(); i++)
            {
                if (demand_so_far < std::get<2>(nodes[path[i]]) && demand_so_far >= std::get<1>(nodes[path[i]]))
                {
                    std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>> new_node_one = nodes[path[i]];
                    std::get<1>(new_node_one) = demand_so_far + 1;
                    std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>> new_node_two = nodes[path[i]];
                    std::get<2>(new_node_two) = demand_so_far;

                    //add nodes
                    nodes[path[i]] = new_node_one;
                    nodes.push_back(new_node_two);

                    //add edges (edges for 1 added by default) and remove from chat
                    edges[nodes.size() - 1] = edges[path[i]];

                    for (int potend : edges[path[i]])
                    {
                        chat_fg.erase(std::pair<int, int>{path[i], potend});
                    }

                    for (auto edge : edges)
                    {
                        std::unordered_set<int> endings = edge.second;
                        if (endings.find(path[i]) != endings.end())
                        {
                            edges[edge.first].insert(nodes.size() - 1);
                            chat_fg.erase(std::pair<int, int>{edge.first, path[i]});
                        }
                    }
                }
                demand_so_far -= demand_fg[std::pair<int, int>{i, i + 1}];
            }
            continue;
        }

        //case 3: time window infeasible
        bool time_infeasible = false;
        double time_remaining = time_start;
        for (int i = 0; i < path.size() - 1; i++)
        {
            std::vector<int> cur_customer_path = path_fg[std::pair<int, int>{path[i], path[i + 1]}];
            for (int cust_idx = 1; cust_idx < cur_customer_path.size(); cust_idx++)
            {
                time_remaining -= time_remaining - timeMat[cur_customer_path[cust_idx - 1]][cur_customer_path[cust_idx]];
                if (time_remaining < latestTimes[cur_customer_path[cust_idx]])
                {
                    time_infeasible = true;
                    break;
                }
                time_remaining = std::min(earliestTimes[cur_customer_path[cust_idx]], time_remaining);
            }
            if (time_infeasible == true)
            {
                break;
            }
        }


        //doing the check inside the splitting condition
        if (time_infeasible == true)
        {
            double arrival_time = time_start;
            bool neg_inf = false;
            for (int i = 0; i < path.size() - 1; i++)
            {
                std::vector<int> cur_customer_path = path_fg[std::pair<int, int>{path[i], path[i + 1]}];
                for (int cust_idx = 1; cust_idx < cur_customer_path.size(); cust_idx++)
                {
                    arrival_time -= arrival_time - timeMat[cur_customer_path[cust_idx - 1]][cur_customer_path[cust_idx]];
                    if (arrival_time < latestTimes[cur_customer_path[cust_idx]])
                    {
                        neg_inf = true;
                        break;
                    }
                    arrival_time = std::min(earliestTimes[cur_customer_path[cust_idx]], arrival_time);
                }
                if (neg_inf == true)
                {
                    break;
                }
                if ((arrival_time < std::get<4>(nodes[path[i]])) && (arrival_time >= std::get<3>(nodes[path[i]])))
                {
                    std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>> new_node_one = nodes[path[i]];
                    std::get<3>(new_node_one) = arrival_time + 1;
                    std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>> new_node_two = nodes[path[i]];
                    std::get<4>(new_node_two) = arrival_time;

                    //add nodes
                    nodes[path[i]] = new_node_one;
                    nodes.push_back(new_node_two);

                    //add edges (edges for 1 added by default) and remove from chat
                    edges[nodes.size() - 1] = edges[path[i]];

                    for (int potend : edges[path[i]])
                    {
                        chat_fg.erase(std::pair<int, int>{path[i], potend});
                    }

                    for (auto edge : edges)
                    {
                        std::unordered_set<int> endings = edge.second;
                        if (endings.find(path[i]) != endings.end())
                        {
                            edges[edge.first].insert(nodes.size() - 1);
                            chat_fg.erase(std::pair<int, int>{edge.first, path[i]});
                        }
                    }
                }
            }
            continue;
        }

        //case 4: cycle detected
        std::vector<int> complete_path; //no need to add depot at the end 
        for (int i = 0; i < path.size() - 1; i++)
        {
            std::vector<int> cur_customer_path = path_fg[std::pair<int, int>{path[i], path[i + 1]}];
            for (int cust_idx = 0; cust_idx < cur_customer_path.size() - 1; cust_idx++)
            {
                complete_path.push_back(cur_customer_path[cust_idx]);
            }
        }

        int shortest_path_length = INF;
        std::pair<int, int> shortest_path_pos;
        int shortest_path_customer;
        for (int i = 0; i < complete_path.size(); i++)
        {
            for (int j = i + 1; j < complete_path.size(); j++)
            {
                if (complete_path[i] == complete_path[j])
                {
                    if ((j - i + 1) < shortest_path_length)
                    {
                        shortest_path_length = j - i + 1;
                        shortest_path_customer = complete_path[i];
                        shortest_path_pos = std::pair<int, int>{ i,j };
                    }
                }
            }
        }

        //if there is a cycle
        if (shortest_path_length < INF)
        {
            int num_cust_crossed = 0;
            int prev_num_cust_crossed = num_cust_crossed;
            int node_i;
            int node_j;

            //finding node positions
            for (int i = 0; i < path.size() - 1; i++)
            {
                if (num_cust_crossed >= std::get<1>(shortest_path_pos))
                {
                    node_j = path[i];
                    break;
                }
                std::vector<int> cur_customer_path = path_fg[std::pair<int, int>{path[i], path[i + 1]}];
                num_cust_crossed += cur_customer_path.size() - 1;
                if ((prev_num_cust_crossed <= std::get<1>(shortest_path_pos)) && (num_cust_crossed > std::get<0>(shortest_path_pos)))
                {
                    node_i = path[i];
                }
            }

            //go through the nodes and split
            for (int i = node_i + 1; i < node_j; i++)
            {
                std::unordered_set<int> mg_minus = std::get<5>(nodes[node_i]);
                std::unordered_set<int> mg_plus = std::get<6>(nodes[node_i]);

                    if ((mg_minus.find(shortest_path_customer) == mg_minus.end()) && (mg_plus.find(shortest_path_customer) != mg_plus.end()))
                    {
                        //we split this node
                        std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>> new_node_one = nodes[node_i];
                        std::get<5>(new_node_one).insert(shortest_path_customer);
                        std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>> new_node_two = nodes[node_j];
                        std::get<6>(new_node_two).erase(shortest_path_customer);

                        //add nodes
                        nodes[node_i] = new_node_one;
                        nodes.push_back(new_node_two);

                        //add edges (edges for 1 added by default) and remove from chat
                        edges[nodes.size() - 1] = edges[node_i];

                        for (int potend : edges[node_i])
                        {
                            chat_fg.erase(std::pair<int, int>{node_i, potend});
                        }

                        for (auto edge : edges)
                        {
                            std::unordered_set<int> endings = edge.second;
                            if (endings.find(node_i) != endings.end())
                            {
                                edges[edge.first].insert(nodes.size() - 1);
                                chat_fg.erase(std::pair<int, int>{edge.first, node_i});
                            }
                        }
                    }
            }

            continue;
        }


        //otherwise return
        return std::tuple<std::vector<int>, std::unordered_map<std::pair<int, int>, int, PairHasher>, std::unordered_map<std::pair<int, int>, std::vector<int>, PairHasher>>{path, demand_fg, path_fg};
    }
}

std::pair<std::vector<std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>>>, std::unordered_map<int, std::unordered_set<int>>> create_initial_nodes(const int num_cust, const int capacity, const int time_start, std::vector<int> demands, std::vector<double> earliestTimes, std::vector<double> latestTimes, std::vector<std::vector<int>> neighMat)
{
    //creating all the nodes

    std::vector<std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>>> nodes = {};

    std::tuple<int, double, double, double, double, std::unordered_set<int>, std::unordered_set<int>> source(num_cust, capacity, capacity, time_start, time_start, {}, {});

    std::unordered_set<int> all_cust = {};

    for (int i = 0; i < num_cust; i++)

    {
        all_cust.insert(i);
    }

    std::tuple<int, double, double, double, double, std::unordered_set<int>, std::unordered_set<int>> sink(num_cust + 1, 0, capacity, 0, time_start, {}, all_cust);

    for (int i = 0; i < num_cust; i++)

    {

        std::unordered_set<int> cur_possible = all_cust;

        cur_possible.erase(i);

        std::tuple<int, double, double, double, double, std::unordered_set<int>, std::unordered_set<int>> customer(i, demands[i], capacity, earliestTimes[i], latestTimes[i], {}, cur_possible);

        nodes.push_back(customer);

    }

    nodes.push_back(source);

    nodes.push_back(sink);

    //creating the initial edges
    std::unordered_map<int, std::unordered_set<int>> edges;

    for (int i = 0; i < num_cust + 1; i++)
    {
        edges[i] = {};
        for (int j = 0; j < num_cust + 2; j++)
        {
            if (i != num_cust + 1)
            {
                std::vector<int> neighSet = neighMat[i];
                if (std::find(neighSet.begin(), neighSet.end(), j) != neighSet.end()) {
                    continue;
                }
            }

            if (((i != num_cust) || (j != num_cust + 1)) && (i != j))
            {
                edges[i].insert(j);
            }
        }

    }

    std::pair<std::vector<std::tuple<int, int, int, int, int, std::unordered_set<int>, std::unordered_set<int>>>, std::unordered_map<int, std::unordered_set<int>>> nodesnedges{ nodes,edges };
    return nodesnedges;
}

double compute_eta(std::vector<std::vector<int>> P, std::vector<double> duals, std::vector<double> earliestTimes, std::vector<double> latestTimes, std::unordered_map<std::vector<int>, std::array<double, 3>, VectorHasher> tauc, std::unordered_map<std::vector<int>, std::unordered_set<std::vector<int>>, VectorHasher> Rp, std::vector<int> demands)
{

    double min_so_far = INF;

    for (std::vector<int> p : P)
    {

        int d_p = compute_demand(p, demands);

        int pLen = p.size();

        int u = p[pLen - 2];

        int v = p[pLen - 1];

        double cpt = compute_cpt(p, earliestTimes[u], latestTimes[v], tauc, Rp).first;

        double pi = compute_pi(p, duals);

        double curValue = (cpt - pi) / d_p;

        if (curValue < min_so_far)
        {

            min_so_far = curValue;

        }

    }

    return (-1 * min_so_far);

}

double compute_pi(std::vector<int> p, std::vector<double> duals)

{

    int pLen = p.size();

    int u = p[pLen - 2];

    double totalPi = duals[u];



    for (int i = 0; i < pLen - 2; i++)

    {

        totalPi += duals[p[i]];

    }



    return totalPi;

}



//to compute d_p, note p is of the form N_p, u, v
int compute_demand(std::vector<int> p, std::vector<int> demands)

{

    int pLen = p.size();

    int u = p[pLen - 2];

    int totalDemand = demands[u];

    for (int i = 0; i < pLen - 2; i++)
    {

        totalDemand += demands[p[i]];

    }

    return totalDemand;

}

std::pair<double, std::vector<int>> compute_cpt(std::vector<int> p, const int t_1, const int t_2, std::unordered_map<std::vector<int>, std::array<double, 3>, VectorHasher> tauc, std::unordered_map<std::vector<int>, std::unordered_set<std::vector<int>>, VectorHasher> Rp)

{

    double cpt = INF;
    std::vector<int> best_path = {};

    for (std::vector<int> r : Rp[p])

    {

        if ((tauc[r][1] <= t_1) && (t_2 <= -tauc[r][2] + std::min((double)t_1, tauc[r][0])))

        {

            if (tauc[r][2] <= cpt)

            {

                cpt = tauc[r][2];
                best_path = r;
            }

        }

    }

    std::pair<double, std::vector<int>> cptnpath{ cpt,best_path };

    return cptnpath;
}

