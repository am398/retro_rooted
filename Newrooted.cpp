#include <bits/stdc++.h>
using namespace std;

#define mp make_pair

typedef pair<int, int> i;

const double INF = 1e15;

struct Edge {
    int u, v, id;
    double w;
    int replacedId = -1;
    Edge() : u(-1), v(-1), w(0.0), id(-1), replacedId(-1) {}
    Edge(int u, int v, double w, int id_) : u(u), v(v), w(w), id(id_), replacedId(-1) {}

    bool operator==(const Edge &other) const {
        return (id == other.id) && (u == other.u) && (v == other.v) && (w == other.w) && (replacedId == other.replacedId);
    }

    bool operator<(const Edge &other) const {
        return w < other.w;
    }
};

struct RootedVertex {
    int parent;
    int root;
    Edge maxE;
    RootedVertex() : parent(-1), root(-1), maxE() {}
};

enum EdgeStatus {
    NONE, DEL, INS, RPL
};

vector<vector<pair<int, int>>> buildExAdj(vector<Edge> &Ex, int V) {
    vector<vector<pair<int, int>>> adjEx(V);
    for (int i = 0; i < Ex.size(); i++) {
        Edge &e = Ex[i];
        if (e.u < 0 || e.v < 0) continue;
        adjEx[e.u].push_back({e.v, i});
        adjEx[e.v].push_back({e.u, i});
    }
    return adjEx;
}

int findMinHeightVertex(vector<vector<pair<int, int>>> &adjEx) {
    int V = adjEx.size();
    vector<int> degree(V, 0);
    queue<int> q;

    for (int i = 0; i < V; i++) {
        degree[i] = adjEx[i].size();
        if (degree[i] == 1) q.push(i);
    }

    int remainingVertices = V;
    while (remainingVertices > 2) {
        int sz = q.size();
        remainingVertices -= sz;

        for (int i = 0; i < sz; i++) {
            int leaf = q.front();
            q.pop();

            for (const auto &[neighbor, eidx] : adjEx[leaf]) {
                if (--degree[neighbor] == 1) q.push(neighbor);
            }
        }
    }

    if (q.size() == 1) return q.front();
    else if (q.size() == 2) {
        int first = q.front(); q.pop();
        return first;
    }
    return -1;
}

class RootedTree {
public:
    int V;
    vector<RootedVertex> RootedT;
    vector<Edge> Ex;
    vector<Edge> Er;

    RootedTree() : V(0) {}
    RootedTree(int V) : V(V) { RootedT.resize(V); }

    void bfsFromRoot(int root, vector<vector<pair<int, int>>> &adjEx) {
        queue<int> q;
        for (auto &nbr : adjEx[root]) {
            int v = nbr.first, eidx = nbr.second;
            if (RootedT[v].parent == -1) {
                RootedT[v] = {root, root, Ex[eidx]};
                q.push(v);
            }
        }

        while (!q.empty()) {
            int cur = q.front(); q.pop();
            for (auto &nbr : adjEx[cur]) {
                int nxt = nbr.first, eidx = nbr.second;
                if (RootedT[nxt].parent == -1) {
                    RootedT[nxt].parent = cur;
                    RootedT[nxt].root = RootedT[cur].root;
                    RootedT[nxt].maxE = (Ex[eidx].w > RootedT[cur].maxE.w) ? Ex[eidx] : RootedT[cur].maxE;
                    q.push(nxt);
                }
            }
        }
    }

    void Create_Tree(vector<vector<pair<int, int>>> &adjEx) {
        RootedT.assign(V, RootedVertex());
        for (int i = 0; i < V; i++) {
            if (RootedT[i].parent == -1) {
                int r = findMinHeightVertex(adjEx);
                RootedT[r] = {r, r, Edge()};
                bfsFromRoot(r, adjEx);
            }
        }
    }

    int findRoot(int v) { return RootedT[v].root; }

    vector<int> getPathToRoot(int x) {
        vector<int> path;
        while (x != RootedT[x].parent) {
            path.push_back(x);
            x = RootedT[x].parent;
        }
        path.push_back(x);
        return path;
    }

    Edge findMaxEdgeOnPath(int u, int v) {
        if (findRoot(u) != findRoot(v)) return Edge();

        Edge umx = RootedT[u].maxE, vmx = RootedT[v].maxE;
        if (umx.id != vmx.id) return umx.w > vmx.w ? umx : vmx;

        vector<int> upath = getPathToRoot(u), vpath = getPathToRoot(v);
        reverse(upath.begin(), upath.end());
        reverse(vpath.begin(), vpath.end());
        int lca = -1, len = min(upath.size(), vpath.size());
        for (int i = 0; i < len; i++) {
            if (upath[i] != vpath[i]) break;
            lca = upath[i];
        }

        auto maxOnPath = [&](int start, int end) {
            Edge mxE;
            for (int cur = start; cur != end; cur = RootedT[cur].parent) {
                for (auto &nb : buildExAdj(Ex, V)[cur]) {
                    if (nb.first == RootedT[cur].parent) {
                        Edge &ed = Ex[nb.second];
                        if (ed.w > mxE.w) mxE = ed;
                        break;
                    }
                }
            }
            return mxE;
        };

        Edge mxE1 = maxOnPath(u, lca), mxE2 = maxOnPath(v, lca);
        return mxE1.w > mxE2.w ? mxE1 : mxE2;
    }

    void Classify_Edges(vector<Edge> &CE, vector<EdgeStatus> &Status, vector<Edge> &Marked, vector<bool> &operation) {
        for (int i = 0; i < CE.size(); i++) {
            Status[i] = NONE;
            Marked[i] = Edge();

            Edge &E = CE[i];
            if (operation[i]) { // Deletion
                Status[i] = DEL;
            } else { // Insertion
                if (findRoot(E.u) != findRoot(E.v)) {
                    Status[i] = INS;
                } else {
                    Edge MaxW = findMaxEdgeOnPath(E.u, E.v);
                    if (MaxW.w > E.w && MaxW.id != -1) {
                        Marked[i] = MaxW;
                        Status[i] = RPL;
                        MaxW.replacedId = i;
                    }
                }
            }
        }
    }

    void Process_Status(vector<Edge> &CE, vector<EdgeStatus> &Status, vector<Edge> &Marked) {
        vector<Edge> newCE;

        for (int i = 0; i < CE.size(); i++) {
            Edge &E = CE[i];
            EdgeStatus S = Status[i];

            if (S == DEL) {
                auto it = find(Ex.begin(), Ex.end(), E);
                if (it != Ex.end()) Ex.erase(it);
            } else if (S == NONE) {
                Er.push_back(E);
            } else if (S == INS) {
                Ex.push_back(E);
            } else if (S == RPL) {
                Edge &oldE = Marked[i];
                auto it = find(Ex.begin(), Ex.end(), oldE);
                if (it != Ex.end()) {
                    Ex.erase(it);
                    Ex.push_back(E);
                } else {
                    newCE.push_back(E);
                }
            }
        }

        CE = newCE;
        vector<vector<pair<int, int>>> adjEx = buildExAdj(Ex, V);
        Create_Tree(adjEx);
    }

    void Repair_Tree(vector<Edge> &Er, vector<EdgeStatus> &Status, vector<Edge> &Marked) {
        for (int i = 0; i < Er.size(); i++) {
            Edge &E = Er[i];
            if (findRoot(E.u) != findRoot(E.v)) {
                Status[i] = INS;
            } else {
                Edge MaxW = findMaxEdgeOnPath(E.u, E.v);
                if (MaxW.w == INF) {
                    Status[i] = RPL;
                    Marked[i] = MaxW;
                }
            }
        }
    }

    vector<Edge> ProcessAllEdges(vector<Edge> &CE, vector<bool> &operation) {
        vector<EdgeStatus> Status(CE.size(), NONE);
        vector<Edge> Marked(CE.size());

        while (!CE.empty()) {
            Classify_Edges(CE, Status, Marked, operation);
            Process_Status(CE, Status, Marked);
        }

        vector<EdgeStatus> StatusR(Er.size(), NONE);
        vector<Edge> MarkedR(Er.size());
        Repair_Tree(Er, StatusR, MarkedR);
        Process_Status(Er, StatusR, MarkedR);

        return Ex;
    }

    void updateNode(int id, int x, int y, int w, bool op) {
        vector<Edge> CE{Edge(x, y, w, id)};
        vector<bool> operation{op};
        ProcessAllEdges(CE, operation);
    }
};

class DSU {
public:
    vector<int> parent, rank;

    DSU(int n) {
        parent.resize(n);
        rank.assign(n, 0);
        for (int i = 0; i < n; ++i) parent[i] = i;
    }

    int find(int x) {
        return parent[x] == x ? x : parent[x] = find(parent[x]);
    }

    bool union_sets(int x, int y) {
        x = find(x), y = find(y);
        if (x == y) return false;
        if (rank[x] < rank[y]) swap(x, y);
        parent[y] = x;
        if (rank[x] == rank[y]) rank[x]++;
        return true;
    }
};

pair<vector<Edge>, vector<Edge>> findMST(vector<Edge> &edges, int n) {
    sort(edges.begin(), edges.end());
    DSU dsu(n);
    vector<Edge> Ex, Er;

    for (Edge &e : edges) {
        if (dsu.union_sets(e.u, e.v)) {
            e.id = Ex.size();
            Ex.push_back(e);
        } else {
            e.id = Er.size() + Ex.size();
            Er.push_back(e);
        }
    }
    return {Ex, Er};
}

const int N = 50011;
const int M = sqrt(N) + 10;

int tp, n, m, q;
vector<pair<Edge, bool>> que[N];

vector<Edge> getMST(int t, vector<RootedTree> &rt) {
    int bkt = t / M;
    vector<Edge> CE;
    vector<bool> operation;

    for (int i = bkt*M; i <= t; ++i)
        for (auto &[e, op] : que[i])
            CE.push_back(e), operation.push_back(op);

    return rt[bkt].ProcessAllEdges(CE, operation);
}

void addAndDeleteEdge(int x, int y, int w, int t, int id, vector<RootedTree> &rt, bool op) {
    que[t].push_back({Edge(x, y, w, id), op});
    for (int i = (t + M - 1)/M; i < tp; ++i)
        rt[i].updateNode(id, x, y, w, op);
}

int main() {
    clock_t tStart = clock();
    freopen("input.txt", "r", stdin);

    tp = N / M + 2;
    scanf("%d %d %d", &n, &m, &q);

    vector<Edge> edges(m);
    for (int i = 0; i < m; ++i) {
        int x, y, w;
        scanf("%d %d %d", &x, &y, &w);
        edges[i] = Edge(--x, --y, w, i);
    }

    auto [Ex, Er] = findMST(edges, n);
    vector<vector<pair<int, int>>> adjEx = buildExAdj(Ex, n);
    vector<RootedTree> rt(tp, RootedTree(n));

    for (int i = 0; i < tp; ++i) {
        rt[i].Ex = Ex;
        rt[i].Er = Er;
        rt[i].Create_Tree(adjEx);
    }

    for (int i = 0; i < q; ++i) {
        int op; scanf("%d", &op);
        if (op == 0) {
            int t; scanf("%d", &t);
            // Process MST query
        } else {
            int x, y, w, t;
            scanf("%d %d %d %d", &x, &y, &w, &t);
            x--; y--;
            addAndDeleteEdge(x, y, w, t, m + i, rt, op == 1);
        }
    }

    printf("Time: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
    return 0;
}