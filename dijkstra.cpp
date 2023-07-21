// C++ Program to find Dijkstra's shortest path using
// priority_queue in STL
#include <utility>
#include "dijkstra.h"

#include <list>
#include <vector>
#include <queue>

#define INF 0x3f3f3f3f

// iPair ==> Integer Pair
typedef std::pair<int, int> iPair;

// Allocates memory for adjacency list
Graph::Graph(int V)
{
    this->V = V;
    adj = new std::list<iPair>[V];
}

void Graph::addEdge(int u, int v, int w)
{
    adj[u].push_back(std::make_pair(v, w));
    adj[v].push_back(std::make_pair(u, w));
}

// Give shortest path to ending vertex
std::vector<int> Graph::shortestPath(int src, int end)
{
    // Create a priority queue to store vertices that
    // are being preprocessed. This is weird syntax in C++.
    // Refer below link for details of this syntax
    // https://www.geeksforgeeks.org/implement-min-heap-using-stl/
    std::priority_queue<iPair, std::vector<iPair>, std::greater<iPair> >
        pq;

    // Create a vector for distances and initialize all
    // distances as infinite (INF)
    std::vector<int> dist(V, INF);

    //predecessor vector
    std::vector<int> pred(V, 0);

    for (int i = 0; i < V; i++)
    {
        pred[i] = i;
    }

    // Insert source itself in priority queue and initialize
    // its distance as 0.
    pq.push(std::make_pair(0, src));
    dist[src] = 0;

    /* Looping till priority queue becomes empty (or all
    distances are not finalized) */
    while (!pq.empty()) {
        // The first vertex in pair is the minimum distance
        // vertex, extract it from priority queue.
        // vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices
        // sorted distance (distance must be first item
        // in pair)
        int u = pq.top().second;
        pq.pop();

        // 'i' is used to get all adjacent vertices of a
        // vertex
        std::list<std::pair<int, int> >::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i) {
            // Get vertex label and weight of current
            // adjacent of u.
            int v = (*i).first;
            int weight = (*i).second;

            // If there is shorted path to v through u.
            if (dist[v] > dist[u] + weight) {
                // Updating distance of v
                dist[v] = dist[u] + weight;
                pred[v] = u;
                pq.push(std::make_pair(dist[v], v));
            }
        }
    }

    std::vector<int> shortest_path;
    int cur_node = end;
    while (true)
    {
        shortest_path.push_back(cur_node);
        if (pred[cur_node] == cur_node)
        {
            break;
        }
        cur_node = pred[cur_node];
    }
    reverse(shortest_path.begin(), shortest_path.end());

    return shortest_path;
}
