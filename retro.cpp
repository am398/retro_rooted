#include <bits/stdc++.h>
using namespace std;

#define mp make_pair

typedef pair<int, int> i;

const double INF = 1e15;

// Structure to store an edge
struct Edge
{
    int u, v, id;
    double w;
    // replacedId indicates which changed edge (by index) might replace this edge
    // -1 meaning no replacement
    int replacedId = -1;
    Edge() : u(-1), v(-1), w(0.0), replacedId(-1) {}
    Edge(int u, int v, double w, int id_) : u(u), v(v), w(w), id(id_), replacedId(-1) {}

    bool operator==(const Edge &other) const
    {
        return (u == other.u && v == other.v && w == other.w && replacedId == other.replacedId);
    }

    bool operator<(const Edge &other) const
    {
        return w < other.w;
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
        maxE = Edge();
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

vector<vector<pair<int, int>>> buildExAdj(vector<Edge> &Ex, int V)
{

    vector<vector<pair<int, int>>> adjEx(V); // (neighbor, index_of_edge_in_Ex)
    for (int i = 0; i < Ex.size(); i++)
    {
        Edge &e = Ex[i];
        if (e.u < 0 || e.v < 0)
            continue;
        adjEx[e.u].push_back({e.v, i});
        adjEx[e.v].push_back({e.u, i});
    }

    return adjEx;
}

int findMinHeightVertex(vector<vector<pair<int, int>>> &adjEx)
{
    int V = adjEx.size();
    vector<int> degree(V, 0); // Degree of each vertex
    queue<int> q;             // Queue to store leaves

    // Initialize the degree array and find initial leaves
    for (int i = 0; i < V; i++)
    {
        degree[i] = adjEx[i].size();
        if (degree[i] == 1)
        {
            q.push(i); // Add leaves to the queue
        }
    }

    // Trim leaves level by level
    int remainingVertices = V;
    while (remainingVertices > 2)
    {
        int sz = q.size();
        remainingVertices -= sz; // Remove current leaves

        for (int i = 0; i < sz; i++)
        {
            int leaf = q.front();
            q.pop();

            // Reduce the degree of neighbors
            for (const auto &[neighbor, w] : adjEx[leaf])
            {
                degree[neighbor]--;
                if (degree[neighbor] == 1)
                {
                    q.push(neighbor); // Add new leaf
                }
            }
        }
    }

    // The remaining vertices are the centers
    if (q.size() == 1)
    {
        return q.front();
    }
    else if (q.size() == 2)
    {
        int first = q.front();
        q.pop();
        int second = q.front();

        // Return either vertex (degree comparison is unnecessary as both are balanced)
        return first;
    }

    return -1; // Should not reach here for a valid tree
}

class RootedTree
{
public:
    int V;
    vector<RootedVertex> RootedT;
    vector<Edge> Ex;
    vector<Edge> Er;

    RootedTree() : V(0) {}

    RootedTree(int V)
    {
        this->V = V;
        RootedT.resize(V);
    }

    // Step 1: Create the Rooted Tree (Algorithm 3)
    // We find a vertex with Minimum Height vertex as root and run BFS from it.
    // If multiple components exist, we repeat.

    void bfsFromRoot(int root, vector<vector<pair<int, int>>> &adjEx)
    {

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

    void Create_Tree(vector<vector<pair<int, int>>> &adjEx)
    {
        V = adjEx.size();

        for (int i = 0; i < V; i++)
        {
            if (RootedT[i].parent == -1)
            {
                int r = findMinHeightVertex(adjEx);
                RootedT[r].parent = r;
                RootedT[r].root = r;
                // BFS from r
                bfsFromRoot(r, adjEx);
                break;
            }
        }
    }

    // Utility: find root of a vertex using Rooted Vertex
    int findRoot(int v)
    {
        return RootedT[v].root;
    }

    // Utility function to get path to root
    std::vector<int> getPathToRoot(int x)
    {
        std::vector<int> path;
        while (x != this->RootedT[x].parent)
        {
            path.push_back(x);
            x = RootedT[x].parent;
        }
        path.push_back(x);
        return path;
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
            return Edge();
        }

        int root = RootedT[u].root;
        // If the maximum edges on u->root and v->root differ, O(1)
        // else O(h) by traversing up.

        Edge umx = RootedT[u].maxE;
        Edge vmx = RootedT[v].maxE;

        if (umx.id == vmx.id)
        {
            std::vector<int> upath = getPathToRoot(u);
            std::vector<int> vpath = getPathToRoot(v);
            std::reverse(upath.begin(), upath.end());
            std::reverse(vpath.begin(), vpath.end());
            vector<vector<pair<int, int>>> adjEx = buildExAdj(Ex, V);

            int len = (int)std::min(upath.size(), vpath.size());
            int lca = -1;
            for (int i = 0; i < len; i++)
            {
                if (upath[i] == vpath[i])
                    lca = upath[i];
                else
                    break;
            }

            auto maxOnPath = [&](int start, int end)
            {
                double mxw = -1;
                Edge mxE(-1, -1, -1, -1);
                int cur = start;
                while (cur != end)
                {
                    int p = RootedT[cur].parent;
                    if (p == cur)
                        break;
                    // Find edge cur-p
                    bool found = false;
                    for (auto &nb : adjEx[cur])
                    {
                        if (nb.first == p)
                        {
                            Edge &ed = Ex[nb.second];
                            if (ed.w > mxw)
                            {
                                mxw = ed.w;
                                mxE = ed;
                            }
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                        break;
                    cur = p;
                }
                return mxE.w >= 0 ? mxE : Edge();
            };
            Edge mxE1 = maxOnPath(u, lca);
            Edge mxE2 = maxOnPath(v, lca);

            if (mxE1.w > mxE2.w)
                return mxE1;
            return mxE2;
        }
        return umx.w > vmx.w ? umx : vmx;
    }

    // Step 2: Classify Edges
    // CE: changed edges
    // Status, Marked arrays
    void Classify_Edges(vector<Edge> &CE, vector<EdgeStatus> &Status, vector<Edge> &Marked, vector<bool> &operation)
    {

        for (int i = 0; i < (int)CE.size(); i++)
        {
            Status[i] = NONE;
            Marked[i] = Edge();

            Edge &E = CE[i];
            bool op = operation[i];
            if (op == 1)
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
                            MaxW.replacedId = i;
                        }
                    }
                    // else Status is NONE by default (added to remainder)
                }
            }
        }
    }

    // Step 3: Process edges by Status
    void Process_Status(vector<Edge> &CE, vector<EdgeStatus> &Status, vector<Edge> &Marked)
    {
        vector<pair<Edge, bool>> newCE;
        for (int i = 0; i < (int)CE.size(); i++)
        {
            Edge &E = CE[i];
            EdgeStatus S = Status[i];
            if (S == DEL)
            {
                E.w = INF; // Mark deleted
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

                int rpl_idx = Erpl.replacedId;

                if (rpl_idx != -1 && rpl_idx == i)
                {
                    // Perfect match
                    // Add E to key edges
                    Ex.push_back(E);
                    // Mark Erpl as replaced
                    Erpl.w = -1;
                }
                else
                {
                    newCE.push_back({E, 0});
                }
            }
        }

        CE.clear();

        // // // // After processing all, we run BFS again to fix parents and maxE
        // vector<vector<pair<int, int>>> adjEx = buildExAdj(Ex);
        // Create_Tree(adjEx);

        // CE = std::move(newCE); // Swap to minimize copy overhead
        // newCE.clear();
    }

    // Step for repairing the tree
    void Repair_Tree(vector<Edge> &Er, vector<EdgeStatus> &Status, vector<Edge> &Marked)
    {
        // Check only edges that can be inserted to reconnect parts
        for (int i = 0; i < (int)Er.size(); i++)
        {
            Status[i] = NONE;
            Marked[i] = Edge();
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
                }
            }
        }
    }

    // Main algorithm (Algorithm 2)

    vector<Edge> ProcessAllEdges(vector<Edge> &CE, vector<bool> &operation)
    {
        vector<EdgeStatus> Status(CE.size(), NONE);
        vector<Edge> Marked(CE.size());

        while (!CE.empty())
        {
            Classify_Edges(CE, Status, Marked, operation);
            Process_Status(CE, Status, Marked);
        }

        // Step 4: Repair Tree with remainder edges Er
        vector<EdgeStatus> StatusR(Er.size(), NONE);
        vector<Edge> MarkedR(Er.size());
        Repair_Tree(Er, StatusR, MarkedR);

        // Process Remainder edges with RPL status
        Process_Status(Er, StatusR, MarkedR);
        return Ex;
    }

    void updateNode(int id, int x, int y, int w, bool op)
    {
        vector<Edge> CE(1, Edge(x, y, w, id));
        vector<bool> operation(1, op);
        ProcessAllEdges(CE, operation);
    }
};

// Disjoint Set Union-Find with path compression and union by rank
class DSU
{
public:
    vector<int> parent, rank;

    DSU(int n)
    {
        parent.resize(n);
        rank.resize(n, 0);
        for (int i = 0; i < n; ++i)
        {
            parent[i] = i;
        }
    }

    int find(int x)
    {
        if (parent[x] != x)
            parent[x] = find(parent[x]);
        return parent[x];
    }

    bool union_sets(int x, int y)
    {
        int rootX = find(x);
        int rootY = find(y);

        if (rootX != rootY)
        {
            // Union by rank
            if (rank[rootX] > rank[rootY])
            {
                parent[rootY] = rootX;
            }
            else if (rank[rootX] < rank[rootY])
            {
                parent[rootX] = rootY;
            }
            else
            {
                parent[rootY] = rootX;
                rank[rootX]++;
            }
            return true; // Edge was added to the MST
        }
        return false; // Edge would form a cycle, so not added to the MST
    }
};

// Function to find the MST and divide edges into Ex and Er
pair<vector<Edge>, vector<Edge>> findMST(vector<Edge> &edges, int n)
{
    sort(edges.begin(), edges.end());

    DSU dsu(n);
    vector<Edge> Ex; // Edges that are part of the MST
    vector<Edge> Er; // Edges that are not part of the MST

    for (Edge edge : edges)
    {
        if (dsu.union_sets(edge.u, edge.v))
        {
            edge.id = Ex.size() + 1;
            Ex.push_back(edge);
        }
        else
        {
            edge.id = Er.size() + 1;
            Er.push_back(edge);
        }
    }
    return {Ex, Er};
}

// Constants for maximum number of nodes and a block size based on square root decomposition
// root decomposition to be chnaged later on with ML Algorithm

const int N = 50011;
const int M = sqrt(N) + 10;

// Global variables

int tp;                          // Total number of blocks
int n, m, q;                     // Number of nodes, edges, and queries
vector<pair<Edge, bool>> que[N]; // Queries organized by time , 1 Added , 0 Deleted

// Function to compute the MST up to a certain time t

vector<Edge> getMST(int t, vector<RootedTree> &rt)
{
    vector<Edge> CE; // List of edges to be added or deleted
    vector<bool> operation;
    int bkt = t / M; // Determine the current block
    for (int i = bkt * M + 1; i <= t; ++i)
    {
        // I am allowing addition and deletion at the same time
        if (!que[i].empty())
        {
            for (int j = 0; j < que[i].size(); j++)
            {
                CE.push_back(que[i][j].first); // Add/delete edges scheduled at this time
                operation.push_back(que[i][j].second);
            }
        }
    }
    return rt[bkt].ProcessAllEdges(CE, operation);
}

// Function to add or delete an edge with specific parameters

void addAndDeleteEdge(int x, int y, int w, int t, int id, vector<RootedTree> &rt, bool op)
{
    que[t].push_back({Edge(x, y, w, id), op}); // Schedule the edge to be added/deleted at time t
    for (int i = (t + M - 1) / M; i < tp; ++i)
    {
        rt[i].updateNode(id, x, y, w, op);
    }
}

int main()
{
    clock_t tStart = clock();
    cout << tStart << endl;
    if (freopen("input/in000", "r", stdin) == nullptr)
    {
        perror("Error opening file");
        return 1;
    }
    tp = N / M + 2;                // Calculate the total number of blocks
    scanf("%d %d %d", &n, &m, &q); // Read number of nodes, edges, and queries

    vector<RootedTree> rt(tp, RootedTree(n + q)); // Array of Rooted Trees for different blocks
    vector<Edge> edges(m);

    for (int i = 0; i < m; ++i)
    {
        int x, y, w;
        scanf("%d %d %d", &x, &y, &w);
        edges[i] = Edge(--x, --y, w, i);
    }

    auto [Ex, Er] = findMST(edges, n);

    vector<vector<pair<int, int>>> adjEx = buildExAdj(Ex, n);

    // Initialize Rooted Trees for each block
    for (int i = 0; i < tp; ++i)
    {
        rt[i].Ex = Ex;
        rt[i].Er = Er;
        rt[i].Create_Tree(adjEx);
    }

    // Process queries
    for (int i = n + m; i < n + m + q; ++i)
    {
        int op;
        scanf("%d", &op);
        if (op == 0)
        {
            // If op == 0, it's an MST query
            int t;
            scanf("%d", &t);
            // vector<Edge> mst = getMST(t, rt);
            // double sum = 0.0;
            // for (auto &e : mst)
            // {
            //     sum += e.w;
            // }
            // cout << "MST weight: " << sum << endl;
        }
        else
        {
            // If op == 1, it's an edge addition
            // If op == 2, it's an edge deletion
            int x, y, w, t;
            scanf("%d %d %d %d", &x, &y, &w, &t);
            x--;
            y--;                                       // Read edge details and adjust to 0-based index
            // addAndDeleteEdge(x, y, w, t, i, rt, --op); // Schedule the edge to be added at time t
        }
    }

    fclose(stdin);

    printf("RootedTree + sqrt(n): %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    return 0;
}
