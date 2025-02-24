#include <bits/stdc++.h>
using namespace std;
// A very large weight to represent 'infinity'

const double INF = 1e15;

// Structure to store an edge
struct Edge
{
    int u, v;
    double w;
    // replacedId indicates which changed edge (by index) might replace this edge
    // Initialize with -1 meaning no replacement
    int replacedId = -1;

    Edge() : u(-1), v(-1), w(0.0), replacedId(-1) {}
    Edge(int u, int v, double w) : u(u), v(v), w(w), replacedId(-1) {}

    // Overloading the '==' operator to compare two edges
    bool operator==(const Edge &other) const
    {
        return (u == other.u && v == other.v && w == other.w && replacedId == other.replacedId);
    }
};

// Structure to store Rooted Tree information per vertex
struct RootedVertex
{
    int parent;
    int root;
    // maxE: The maximum weighted edge on the path from this vertex to the root
    Edge maxE;
    RootedVertex()
    {
        parent = -1;
        root = -1;
        maxE = Edge(-1, -1, -1);
    }
};

// Status of changed edges
enum EdgeStatus
{
    NONE, // added to remainder
    DEL,  // deleted from key edges
    INS,  // inserted to key edges
    RPL   // replaced an existing key edge
};

// Global variables or encapsulated data
// Assuming we have a graph with V vertices and initial sets of edges
int V;

// Rooted tree structure
vector<RootedVertex> RootedT;

// We will need adjacency lists for the key edges (Ex) to form the current spanning structure
// Adjacency list for the current tree (Ex)
vector<vector<pair<int, int>>> adjEx; // (neighbor, index_of_edge_in_Ex)

// Build adjacency from Ex
void buildExAdj(vector<Edge> &Ex)
{
    adjEx.assign(V, {});
    for (int i = 0; i < (int)Ex.size(); i++)
    {
        Edge &e = Ex[i];
        if (e.u < 0 || e.v < 0)
            continue;
        adjEx[e.u].push_back({e.v, i});
        adjEx[e.v].push_back({e.u, i});
    }
}

// Step 1: Create the Rooted Tree (Algorithm 3)
// We find a vertex with Minimu Height vertex as root and run BFS from it.
// If multiple components exist, we repeat.

int findMinHeightVertex()
{
    // Step 1: Initialize the graph and degree array
    vector<int> degree(V, 0); // Degree of each vertex
    queue<int> q;             // Queue to store leaves

    // Step 2: Build the graph and find initial leaves
    for (int i = 0; i < V; i++)
    {
        degree[i] = adjEx[i].size(); // Degree is the size of the adjacency list
        if (degree[i] == 1)
        {
            q.push(i); // Leaves are vertices with degree 1
        }
    }

    // Step 3: Trim leaves level by level
    int remainingVertices = V;
    while (remainingVertices > 2)
    {
        int sz = q.size();
        remainingVertices -= sz; // Remove all leaves from the graph

        // Remove leaves one by one
        for (int i = 0; i < sz; i++)
        {
            int leaf = q.front();
            q.pop();

            // Reduce the degree of each neighbor of the leaf
            for (int neighbor : adjEx[leaf])
            {
                degree[neighbor]--;

                // If the neighbor becomes a leaf, add it to the queue
                if (degree[neighbor] == 1)
                {
                    q.push(neighbor);
                }
            }
        }
    }

    // Step 4: After trimming, the remaining vertices in the queue are the centers (min height)
    // We need to decide which one to return based on the adjacency list size (degree)
    int vertex = -1;

    if (q.size() == 1)
    {
        // If only one vertex is remaining, return that vertex
        vertex = q.front();
    }
    else if (q.size() == 2)
    {
        // If two vertices remain, return the one with the larger adjacency list
        int first = q.front();
        q.pop();
        int second = q.front();

        // Compare their degrees (size of adjacency list)
        if (adjEx[first].size() > adjEx[second].size())
        {
            vertex = first;
        }
        else
        {
            vertex = second;
        }
    }

    return vertex; // Return the vertex with the maximum degree if there are two remaining, or the single remaining vertex
}

void bfsFromRoot(int root)
{
    // BFS from root
    RootedT[root].parent = root;
    RootedT[root].root = root;
    queue<int> q;
    // Push all neighbors of root
    for (auto &nbr : adjEx[root])
    {
        int v = nbr.first;
        int eidx = nbr.second;
        if (RootedT[v].parent == -1)
        {
            RootedT[v].parent = root;
            RootedT[v].root = root;
            RootedT[v].maxE = Ex[eidx]; // edge (root->v)
            q.push(v);
        }
    }

    while (!q.empty())
    {
        int cur = q.front();
        q.pop();
        for (auto &nbr : adjEx[cur])
        {
            int nxt = nbr.first;
            int eidx = nbr.second;
            if (RootedT[nxt].parent == -1)
            {
                RootedT[nxt].parent = cur;
                RootedT[nxt].root = RootedT[cur].root;
                // Check if this edge is max on path
                if (Ex[eidx].w > RootedT[cur].maxE.w)
                {
                    RootedT[nxt].maxE = Ex[eidx];
                }
                else
                {
                    RootedT[nxt].maxE = RootedT[cur].maxE;
                }
                q.push(nxt);
            }
        }
    }
}

void Create_Tree()
{
    // Initialize RootedT
    for (int i = 0; i < V; i++)
    {
        RootedT[i].parent = -1;
        RootedT[i].root = -1;
        RootedT[i].maxE = Edge(-1, -1, -1);
    }

    // Run BFS until all vertices have parents assigned
    bool allAssigned = false;
    while (!allAssigned)
    {
        allAssigned = true;
        for (int i = 0; i < V; i++)
        {
            if (RootedT[i].parent == -1)
            {
                // This vertex not assigned yet, become a root
                int r = findMinHeightVertex();
                RootedT[r].parent = r;
                RootedT[r].root = r;
                // BFS from r
                bfsFromRoot(r);
                break;
            }
        }

        for (int i = 0; i < V; i++)
        {
            if (RootedT[i].parent == -1)
            {
                allAssigned = false;
                break;
            }
        }
    }
}

// Utility: find root of a vertex using RootedT
int findRoot(int v)
{
    return RootedT[v].root;
}

// Find max weighted edge in path u->v

Edge findMaxEdgeOnPath(int u, int v)
{
    // If in different components
    if (RootedT[u].root != RootedT[v].root)
    {
        // Path passes through roots or is disconnected
        // If they are different roots => no direct path in the tree
        //  a dummy edge with w = -1
        return Edge(-1, -1, -1);
    }

    int root = RootedT[u].root;
    // If the maximum edges on u->root and v->root differ, O(1)
    // else O(h) by traversing up.

    Edge umx = RootedT[u].maxE;
    Edge vmx = RootedT[v].maxE;

    if (umx == vmx)
    {
    }
    else
    {
        return umx.w > vmx.w ? umx : vmx;
    }

    auto getPathToRoot = [&](int x)
    {
        vector<int> path;
        while (x != RootedT[x].parent)
        {
            path.push_back(x);
            x = RootedT[x].parent;
        }
        path.push_back(x);
        return path;
    };
}

// Step 2: Classify Edges
// CE: changed edges
// Status, Marked arrays
void Classify_Edges(vector<Edge> &CE, vector<EdgeStatus> &Status, vector<Edge> &Marked)
{
    // Initialize Status and Marked
    for (int i = 0; i < (int)CE.size(); i++)
    {
        Status[i] = NONE;
        Marked[i] = Edge(-1, -1, -1);
    }

    for (int i = 0; i < (int)CE.size(); i++)
    {
        Edge &E = CE[i];
        if (E.w > 1e14)
        {
            // Marked for deletion
            Status[i] = DEL;
        }
        else
        {
            // Marked for insertion
            int u = E.u, v = E.v;
            int ru = findRoot(u), rv = findRoot(v);
            if (ru != rv)
            {
                // connecting disconnected components
                Status[i] = INS;
            }
            else
            {
                // same component, check if can replace max weighted edge in path
                Edge MaxW = findMaxEdgeOnPath(u, v);
                if (MaxW.w > E.w && MaxW.w >= 0)
                {
                    int x = MaxW.replacedId;
                    // If x = -1 or CE[x].w > E.w
                    if (x == -1 || (x >= 0 && x < (int)CE.size() && CE[x].w > E.w))
                    {
                        Marked[i] = MaxW;
                        Status[i] = RPL;
                        // Mark MaxW replacedId = i
                        MaxW.replacedId = i;
                        break;
                    }
                }
                // else Status is NONE by default (added to remainder)
            }
        }
    }
}

// Step 3: Process edges by Status
void Process_Status(vector<Edge> &Ex, vector<Edge> &Er, vector<Edge> &CE, vector<EdgeStatus> &Status, vector<Edge> &Marked)
{
    // We might produce a new CE set if some replacements fail
    vector<Edge> newCE;
    vector<EdgeStatus> newStatus;
    vector<Edge> newMarked;

    for (int i = 0; i < (int)CE.size(); i++)
    {
        Edge &E = CE[i];
        EdgeStatus S = Status[i];
        if (S == DEL)
        {
            // Delete E from appropriate edge set
            // If E is in Ex, remove it:
            // We'll mark weight INF

            E.w = INF; // Mark deleted
            break;
        }
        else if (S == NONE)
        {
            Er.push_back(E);
        }
        else if (S == INS)
        {
            // Add to key edges
            Ex.push_back(E);
        }
        else if (S == RPL)
        {
            // Find edge to be replaced from Marked[i]
            Edge &Erpl = Marked[i];
            // Check if Erpl's replacedId matches i
            // We must find Erpl in Ex:
            int rpl_idx = Erpl.replacedId;

            if (rpl_idx != -1 && Ex[rpl_idx].replacedId == i)
            {
                // Perfect match
                // Add E to key edges
                Ex.push_back(E);
                // Mark Erpl as replaced
                Ex[rpl_idx].w = -1; // replaced mark
            }
            else
            {
                // Not match or replaced by another
                // Add E to new changed edges
                newCE.push_back(E);
            }
        }
    }

    // // After processing all, we run BFS again to fix parents and maxE
    // buildExAdj();
    // // Reassign parents by BFS from root
    // // To do this, reset RootedT and call Create_Tree again or a BFS update
    // for (int i = 0; i < V; i++)
    // {
    //     RootedT[i].parent = -1;
    //     RootedT[i].root = -1;
    //     RootedT[i].maxE = Edge(-1, -1, -1);
    // }
    // // A single BFS from all roots (Create_Tree)
    // Create_Tree();

    // // If we have newCE (rare case), we run steps 2 and 3 again
    // if (!newCE.empty())
    // {
    //     newStatus.resize(newCE.size(), NONE);
    //     newMarked.resize(newCE.size());
    //     Classify_Edges(newCE, newStatus, newMarked);
    //     Process_Status(Ex, Er, newCE, newStatus, newMarked);
    // }
}

// Step for repairing the tree
// Similar to Classify_Edges but only for insertion from Er if it can reconnect disconnected parts
void Repair_Tree(vector<Edge> &Er, vector<EdgeStatus> &Status, vector<Edge> &Marked)
{
    // Check only edges that can be inserted to reconnect parts
    for (int i = 0; i < (int)Er.size(); i++)
    {
        Edge &E = Er[i];
        // If E can reconnect edges marked as deleted (with INF weight)
        int ru = findRoot(E.u), rv = findRoot(E.v);
        if (ru != rv)
        {
            // Potentially reconnect
            // find max weighted edge in path:
            Edge MaxW = findMaxEdgeOnPath(E.u, E.v);
            // If MaxW has INF weight => means it was deleted
            if (MaxW.w == INF)
            {
                // Mark as RPL
                Status[i] = RPL;
                Marked[i] = MaxW;

                MaxW.replacedId = i;
                break;
            }
        }
    }
}

// Main algorithm (Algorithm 2)
// Input: Ex, Er, CE
// Output: Updated Ex

void ProcessAllEdges(vector<Edge> &Ex, vector<Edge> &Er, vector<Edge> &CE)
{
    // RootedT structure
    RootedT.resize(V);

    // Step 1: Create the Rooted Tree
    buildExAdj(Ex);
    Create_Tree();

    // Prepare Status and Marked arrays
    vector<EdgeStatus> Status(CE.size(), NONE);
    vector<Edge> Marked(CE.size());

    // While CE not empty
    while (!CE.empty())
    {
        // Step 2: Classify Edges
        Classify_Edges(CE, Status, Marked);

        // Step 3: Process by Status
        Process_Status(Ex, Er, CE, Status, Marked);
    }

    // Step 4: Repair Tree with remainder edges Er
    vector<EdgeStatus> StatusR(Er.size(), NONE);
    vector<Edge> MarkedR(Er.size());
    // Similar to Classify_Edges but only RPL case is relevant
    Repair_Tree(Er, StatusR, MarkedR);

    // Process Remainder edges with RPL status
    Process_Status(Ex, Er, Er, StatusR, MarkedR);
}

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    cin >> V;
    int V;
    vector<Edge> Ex;
    vector<Edge> Er;
    vector<Edge> CE;

    ProcessAllEdges(Ex, Er, CE);

      

    return 0;
}
