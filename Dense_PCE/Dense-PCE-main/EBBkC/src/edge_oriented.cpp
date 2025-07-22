#include <set>
#include <algorithm>
#include <unordered_map>
#include <omp.h>
#include <cassert>
#include <chrono>

#include "set_operation.h"
#include "edge_oriented.h"
#include "truss/util/graph/graph.h"
#include "truss/util/log/log.h"
#include "truss/util/timer.h"
#include "truss/decompose/parallel_all_edge_cnc.h"
#include "truss/util/reordering/reorder_utils.h"
#include "truss/decompose/iter_helper.h"

extern const int K, L;
extern unsigned long long N;

// EBBkC_Graph_t::EBBkC_Graph_t() = default;
EBBkC_Graph_t::EBBkC_Graph_t() {
    omp_init_lock(&clique_lock); // Initialize the lock
}

EBBkC_Graph_t::~EBBkC_Graph_t() {
    omp_destroy_lock(&clique_lock); // Destroy the lock
    int i, k = K, node_size = v_size, link_size = e_size;

    if (is_sub) {
        k = K - 2;
        node_size = truss_num;
        link_size = truss_num * (truss_num - 1) / 2;
    }

    delete [] edges;

    if (T) {
        for (i = 0; i < link_size; i++) delete [] T[i];
        delete [] T;
    }

    if (C) {
        for (i = 0; i < link_size; i++) delete [] C[i];
        delete [] C;
    }

    delete [] T_size;

    delete [] C_size;

    if (sub_v) {
        for (i = 0; i <= k; i++) delete [] sub_v[i];
        delete [] sub_v;
    }

    if (sub_e) {
        for (i = 0; i <= k; i++) delete [] sub_e[i];
        delete [] sub_e;
    }

    delete [] sub_v_size;

    delete [] sub_e_size;

    delete [] lab;

    if (DAG_deg) {
        for (i = 0; i <= k; i++) delete [] DAG_deg[i];
        delete [] DAG_deg;

    }

    if (G_deg) {
        for (i = 0; i <= k; i++) delete [] G_deg[i];
        delete [] G_deg;
    }


    delete [] col;

    if (DAG_adj) {
        for (i = 0; i < node_size; i++) delete [] DAG_adj[i];
        delete [] DAG_adj;
    }

    if (G_adj) {
        for (i = 0; i < node_size; i++) delete [] G_adj[i];
        delete [] G_adj;
    }

    if (used) {
        for (i = 0; i <= k; i++) delete [] used[i];
        delete [] used;
    }


    delete [] v_lab;

    if (out_v_size) {
        for (i = 0; i <= k; i++) delete [] out_v_size[i];
        delete [] out_v_size;
    }

    if (out_e_size) {
        for (i = 0; i <= k; i++) delete [] out_e_size[i];
        delete [] out_e_size;
    }

    delete [] F;

    delete [] P;

    delete [] lack_size;

    if (lack) {
        for (i = 0; i < node_size; i++) delete [] lack[i];
        delete [] lack;
    }

    delete [] lev;

    delete [] loc;
}

void EBBkC_Graph_t::build(bool sub) {
    int i, k = K, node_size = v_size, link_size = e_size;

    is_sub = sub;

    if (sub) {
        k = K - 2;
        node_size = truss_num;
        link_size = (truss_num) * (truss_num - 1) / 2;
    }

    sub_v = new int* [k + 1];

    sub_e = new int* [k + 1];

    sub_e_size = new int [k + 1];

    sub_v_size = new int [k + 1];

    for (i = 0; i < k; i++) sub_v[i] = new int [truss_num + 1];
    sub_v[k] = new int [node_size];

    for (i = 0; i < k; i++) sub_e[i] = new int [truss_num * (truss_num - 1) / 2];
    sub_e[k] = new int [link_size];

    sub_v_size[k] = 0;

    sub_e_size[k] = 0;

    lab = new int [node_size];
    for (i = 0; i < node_size; i++) lab[i] = k;

    DAG_deg = new int* [k + 1];
    for (i = 0; i <= k; i++) DAG_deg[i] = new int [node_size];

    G_deg = new int* [k + 1];
    for (i = 0; i <= k; i++) G_deg[i] = new int [node_size];

    col = new int [node_size];

    DAG_adj = new int* [node_size];
    for (i = 0; i < node_size; i++) DAG_adj[i] = new int [truss_num + 1];

    G_adj = new int* [node_size];
    for (i = 0; i < node_size; i++) G_adj[i] = new int  [truss_num + 1];

    used = new bool* [k + 1];
    for (i = 0; i <= k; i++) used[i] = new bool [node_size + 1]();

    v_lab = new int [node_size];
    for (i = 0; i < node_size; i++) v_lab[i] = k;

    out_v_size = new int* [k + 1];
    for (i = 0; i <= k; i++) out_v_size[i] = new int [link_size];

    out_e_size = new int* [k + 1];
    for (i = 0; i <= k; i++) out_e_size[i] = new int [link_size];

    F = new int [truss_num + 1];

    P = new int [truss_num + 1];

    lack_size = new int [node_size];

    lack = new int* [node_size];
    for (i = 0; i < node_size; i++) lack[i] = new int [L + 1];

    lev = new int [node_size]();

    loc = new int [node_size];
}

void EBBkC_Graph_t::branch(int e, EBBkC_Graph_t* g) {
    int c, i, j, k, p, e_, u, v, w, s, t, end, dist, l = K;
    int *old2new = new int[v_size];

    g->v_size = 0;
    g->e_size = 0;
    g->edges = new Edge_t[C_size[e] + 1];
    g->new2old = vector<int>(T_size[e]);

    for (i = 0; i < T_size[e]; i++) {
        u = T[e][i];
        old2new[u] = g->v_size;

        // --- Start of Modification ---
        // This is the fix. It ensures the subproblem's map ('g->new2old') stores
        // the final original vertex ID by looking it up from the parent's map ('this->new2old').
        g->new2old[g->v_size++] = this->new2old[u];
        // --- End of Modification ---
    }

    for (i = 0; i < C_size[e]; i++) {
        e_ = C[e][i];
        s = edges[e_].s;
        t = edges[e_].t;
        g->edges[g->e_size].s = old2new[s];
        g->edges[g->e_size++].t = old2new[t];
    }

    delete [] old2new;

    g->clear_stored_cliques();
    g->add_to_current_clique(this->new2old[edges[e].s]);
    g->add_to_current_clique(this->new2old[edges[e].t]);


    if (l == 3 || l == 4) {
        g->sub_v_size[l - 2] = g->v_size;
        g->sub_e_size[l - 2] = g->e_size;
        return;
    }

    if (g->v_size < l - 2 || g->e_size < (l - 2) * (l - 3) / 2) {
        g->sub_v_size[l - 2] = g->v_size;
        g->sub_e_size[l - 2] = g->e_size;
        return;
    }

    for (j = 0; j < g->v_size; j++) {
        g->col[j] = 0;
        g->DAG_deg[l - 2][j] = 0;
        g->G_deg[l - 2][j] = 0;
    }

    for (j = 0; j < g->e_size; j++) {
        s = g->edges[j].s;
        t = g->edges[j].t;
        g->G_adj[s][g->G_deg[l - 2][s]++] = t;
        g->G_adj[t][g->G_deg[l - 2][t]++] = s;
    }

    auto *list = new KeyVal_t[truss_num + 1];
    for (j = 0; j < g->v_size; j++) {
        list[j].key = j;
        list[j].val = g->G_deg[l - 2][j];
    }
    sort(list, list + g->v_size);

    for (j = 0; j < g->v_size; j++) {
        u = list[j].key;
        for (k = 0; k < g->G_deg[l - 2][u]; k++) {
            v = g->G_adj[u][k];
            g->used[l - 2][g->col[v]] = true;
        }
        for (c = 1; g->used[l - 2][c]; c++);
        g->col[u] = c;
        for (k = 0; k < g->G_deg[l - 2][u]; k++) {
            v = g->G_adj[u][k];
            g->used[l - 2][g->col[v]] = false;
        }
    }
    delete[] list;

    g->sub_v_size[l - 2] = 0;

    for (j = 0; j < g->v_size; j++) {
        g->sub_v[l - 2][g->sub_v_size[l - 2]++] = j;
    }

    g->sub_e_size[l - 2] = 0;
    for (j = 0; j < g->e_size; j++) {
        g->sub_e[l - 2][g->sub_e_size[l - 2]++] = j;
        s = g->edges[j].s;
        t = g->edges[j].t;
        g->edges[j].s = (g->col[s] > g->col[t]) ? s : t;
        g->edges[j].t = (g->col[s] > g->col[t]) ? t : s;
        s = g->edges[j].s;
        t = g->edges[j].t;

        g->DAG_adj[s][g->DAG_deg[l - 2][s]++] = t;
    }

    return;
}

void EBBkC_Graph_t::EBBkC_plus_plus(int l, unsigned long long *cliques) {
    int c, i, j, k, p, e, e_, u, v, w, s, t, end, dist;

    if (sub_v_size[l] < l || (l > 1 && sub_e_size[l] < l * (l - 1) / 2)) return;

    // --- Start of Base Case Modifications ---

    // Top-level call for finding 3-cliques (triangles)
    if (K == 3) {
        // The 'current_clique' already contains the initial edge's two vertices from branch().
        // We just need to add the third vertex to complete the triangle.
        for(i = 0; i < sub_v_size[l]; ++i) {
            u = sub_v[l][i];
            add_to_current_clique(new2old[u]);
            store_current_clique();
            remove_from_current_clique();
            (*cliques)++;
        }
        return;
    }

    // Top-level call for finding 4-cliques
    if (K == 4) {
        // The 'current_clique' has the initial edge. We need to find another edge
        // among the common neighbors to complete the 4-clique.
        for(i = 0; i < sub_e_size[l]; ++i) {
            e = sub_e[l][i];
            u = new2old[edges[e].s];
            v = new2old[edges[e].t];
            add_to_current_clique(u);
            add_to_current_clique(v);
            store_current_clique();
            remove_from_current_clique();
            remove_from_current_clique();
            (*cliques)++;
        }
        return;
    }

    // Base case for the recursion: finding the final edge (l=2)
    if (l == 2) {
        for (i = 0; i < sub_v_size[l]; i++) {
            u = sub_v[l][i];
            for (j = 0; j < DAG_deg[l][u]; j++) {
                v = DAG_adj[u][j];

                // Add the final two vertices, store the complete clique, and backtrack
                add_to_current_clique(new2old[u]);
                add_to_current_clique(new2old[v]);
                store_current_clique();
                remove_from_current_clique();
                remove_from_current_clique();
                
                (*cliques)++;
            }
        }
        return;
    }

    // Base case for the recursion: finding the final triangle (l=3)
    if (l == 3) {
        for (i = 0; i < sub_v_size[l]; i++) {
            u = sub_v[l][i];
            if (col[u] < l) continue;

            add_to_current_clique(new2old[u]); // Add u to the current clique path

            for (j = 0; j < DAG_deg[l][u]; j++) lab[DAG_adj[u][j]] = l - 1;

            for (j = 0; j < DAG_deg[l][u]; j++) {
                v = DAG_adj[u][j];
                if (col[v] < l - 1) continue;
                
                add_to_current_clique(new2old[v]); // Add v to path

                for (k = 0; k < DAG_deg[l][v]; k++) {
                    w = DAG_adj[v][k];
                    if (lab[w] == l - 1) {
                        add_to_current_clique(new2old[w]); // Final vertex found
                        store_current_clique();
                        remove_from_current_clique(); // Backtrack w
                        (*cliques)++;
                    }
                }
                remove_from_current_clique(); // Backtrack v
            }

            for (j = 0; j < DAG_deg[l][u]; j++) lab[DAG_adj[u][j]] = l;
            
            remove_from_current_clique(); // Backtrack u
        }
        return;
    }
    
    // --- End of Base Case Modifications ---

    if (can_terminate(l, cliques)) {
        return;
    }
    
    // --- Main Recursive Step with Backtracking ---
    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (col[u] < l) continue;

        // Add the current vertex 'u' to the clique path before recursing
        add_to_current_clique(new2old[u]);
        
        sub_v_size[l - 1] = 0;
        dist = 0;

        for (j = 0; j < DAG_deg[l][u]; j++) {
            v = DAG_adj[u][j];
            lab[v] = l - 1;
            sub_v[l - 1][sub_v_size[l - 1]++] = v;
            DAG_deg[l - 1][v] = 0;
            G_deg[l - 1][v] = 0;

            if (!used[l][col[v]]) {
                used[l][col[v]] = true;
                dist++;
            }
        }

        if (dist >= l - 1) {
            sub_e_size[l - 1] = 0;
            for (j = 0; j < sub_v_size[l - 1]; j++) {
                v = sub_v[l - 1][j];

                end = DAG_deg[l][v];
                for (k = 0; k < end; k++) {
                    w = DAG_adj[v][k];
                    if (lab[w] == l - 1) {
                        DAG_deg[l - 1][v]++;
                        sub_e_size[l - 1]++;
                        G_deg[l - 1][v]++;
                        G_deg[l - 1][w]++;
                    } else {
                        DAG_adj[v][k--] = DAG_adj[v][--end];
                        DAG_adj[v][end] = w;
                    }
                }
            }
            
            EBBkC_plus_plus(l - 1, cliques);
        }

        // Backtrack: Remove 'u' from the clique path after its recursive branch is fully explored
        remove_from_current_clique();

        for (j = 0; j < sub_v_size[l - 1]; j++) {
            v = sub_v[l - 1][j];
            lab[v] = l;
            used[l][col[v]] = false;
        }
    }
}

void EBBkC_Comb_list(int *list, int list_size, int start, int picked, int k, unsigned long long *cliques, EBBkC_Graph_t* graph, int* current_selection) {
    if (picked == k) {
        (*cliques)++;
        if (graph != nullptr) {
            // --- Start of Modification ---
            // Take the existing partial clique and add the new vertices to it.
            std::vector<int> full_clique = graph->current_clique;
            for (int i = 0; i < k; i++) {
                if (current_selection[i] < graph->new2old.size()) {
                    full_clique.push_back(graph->new2old[current_selection[i]]);
                }
            }

            // Store the complete clique.
            if (!full_clique.empty()) {
                // This function might be called from a single thread, but we use the lock
                // for safety and consistency with the main recursive path.
                omp_set_lock(&graph->clique_lock);
                graph->stored_cliques.push_back(full_clique);
                omp_unset_lock(&graph->clique_lock);
            }
            // --- End of Modification ---
        }
        return;
    }

    for (int i = start; i < list_size; i++) {
        current_selection[picked] = list[i];
        EBBkC_Comb_list(list, list_size, i + 1, picked + 1, k, cliques, graph, current_selection);
    }
}

void EBBkC_Graph_t::list_in_plex(int start, int p, int q, unsigned long long *cliques) {
    if (F_size < q) return;

    if (p == 0) {
        // if (q > F_size - q)
        //     EBBkC_Comb_list(F, F_size, 0, 0, F_size - q, cliques);
        // else
        //     EBBkC_Comb_list(F, F_size, 0, 0, q, cliques);

        // modified to store cliques
        int* temp_selection = new int[F_size];
        if (q > F_size - q)
            EBBkC_Comb_list(F, F_size, 0, 0, F_size - q, cliques, this, temp_selection);
        else
            EBBkC_Comb_list(F, F_size, 0, 0, q, cliques, this, temp_selection);
        delete[] temp_selection;
        return;
    }

    int i, j, u, v, vis = 0;

    for (i = start; i < P_size && P_act >= p; i++) {
        u = P[i];

        if (lev[u]) continue;

        for (j = 0; j < lack_size[u]; j++) {
            v = lack[u][j];
            if (loc[v] >= i && lev[v] == 0) {
                lev[v] = p;
                P_act--;
            }
        }

        list_in_plex(i + 1, p - 1, q, cliques);

        for (j = 0; j < lack_size[u]; j++) {
            v = lack[u][j];
            if (loc[v] >= i && lev[v] == p) {
                lev[v] = 0;
                P_act++;
            }
        }

        P_act--;
        vis++;
    }
    P_act += vis;
}

bool EBBkC_Graph_t::can_terminate(int l, unsigned long long *cliques) {
    int i, j, k, u, v, end, p_;

    if (sub_e_size[l] < sub_v_size[l] * (sub_v_size[l] - L) / 2) return false;

    if (sub_e_size[l] == sub_v_size[l] * (sub_v_size[l] - 1) / 2) {
        // if (l > sub_v_size[l] - l)
        //     EBBkC_Comb_list(sub_v[l], sub_v_size[l], 0, 0, sub_v_size[l] - l, cliques);
        // else
        //     EBBkC_Comb_list(sub_v[l], sub_v_size[l], 0, 0, l, cliques);

        // modified to store cliques
        int* temp_selection = new int[sub_v_size[l]];
        if (l > sub_v_size[l] - l)
            EBBkC_Comb_list(sub_v[l], sub_v_size[l], 0, 0, sub_v_size[l] - l, cliques, this, temp_selection);
        else
            EBBkC_Comb_list(sub_v[l], sub_v_size[l], 0, 0, l, cliques, this, temp_selection);
        delete[] temp_selection;

        return true;
    }

    if (L == 1) return false;

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (sub_v_size[l] - G_deg[l][u] > L) {
            return false;
        }
    }

    F_size = 0;
    P_size = 0;

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (G_deg[l][u] == sub_v_size[l] - 1) {
            loc[u] = -1;
            F[F_size++] = u;
            continue;
        }

        loc[u] = P_size;
        P[P_size++] = u;
    }

    int* e = new int [P_size * P_size]();

    for (i = 0; i < P_size; i++) {
        u = P[i];
        lack_size[u] = 0;

        end = DAG_deg[l][u];
        for (j = 0; j < end; j++) {
            v = DAG_adj[u][j];

            if (loc[v] != -1) {
                e[loc[u] * P_size + loc[v]] = 1;
                e[loc[v] * P_size + loc[u]] = 1;
            }
        }
    }

    for (i = 0; i < P_size * P_size; i++) {
        if (!e[i]) {
            j = i / P_size, k = i % P_size;
            u = P[j], v = P[k];
            lack[u][lack_size[u]++] = v;
        }
    }

    delete [] e;

    for(i = 0; i <= l; i++) {
        P_act = P_size;
        list_in_plex(0, i, l - i, cliques);
    }

    return true;
}

void EBBkC_Graph_t::truss_decompose(const char *dir) {
    graph_t g;

    //load the graph from file
    Graph G(dir);
    g.adj = G.edge_dst;
    g.num_edges = G.node_off;
    g.n = G.nodemax;
    g.m = G.edgemax;

    string reorder_method("core");

    vector <int32_t> new_vid_dict;
    vector <int32_t> old_vid_dict;
    ReorderWrapper(g, dir, reorder_method, new_vid_dict, old_vid_dict);

    //edge list array
    Timer get_eid_timer;

    Edge *edgeIdToEdge = (Edge *) malloc((g.m / 2) * sizeof(Edge));
    assert(edgeIdToEdge != nullptr);
    log_info("Malloc Time: %.9lf s", get_eid_timer.elapsed());
    get_eid_timer.reset();

    //Populate the edge list array
    getEidAndEdgeList(&g, edgeIdToEdge);
    log_info("Init Eid Time: %.9lf s", get_eid_timer.elapsed());
    get_eid_timer.reset();

    int *EdgeSupport = (int *) malloc(g.m / 2 * sizeof(int));
    assert(EdgeSupport != nullptr);

    auto max_omp_threads = omp_get_max_threads();
    omp_set_num_threads(max_omp_threads);
    log_info("Max Threads: %d", max_omp_threads);
    #pragma omp parallel for
    for (auto i = 0; i < max_omp_threads; i++) {
        auto avg = g.m / 2 / max_omp_threads;
        auto iter_beg = avg * i;
        auto iter_end = (i == max_omp_threads - 1) ? g.m / 2 : avg * (i + 1);
        memset(EdgeSupport + iter_beg, 0, (iter_end - iter_beg) * sizeof(int));
    }
    log_info("Init EdgeSupport Time: %.9lf s", get_eid_timer.elapsed());
    get_eid_timer.reset();

    Timer global_timer;
    truss_num = PKT_intersection(&g, EdgeSupport, edgeIdToEdge);

    #pragma omp single
    {
        int u, v, w;
        int *old2new = new int[N_NODES];
        for (int i = 0; i < N_NODES; i++) old2new[i] = -1;

        e_size = 0;
        edges = new Edge_t[g.m / 2];

        T = new int *[g.m / 2];
        T_size = new int[g.m / 2];

        for (int i = 0; i < g.m / 2; i++) {
            if (g.edge_truss[i] <= K) continue;

            Edge e = edgeIdToEdge[i];
            u = e.u;
            v = e.v;

            if (old2new[u] == -1) {
                old2new[u] = (int) new2old.size();
                new2old.push_back(u);
            }

            if (old2new[v] == -1) {
                old2new[v] = (int) new2old.size();
                new2old.push_back(v);
            }

            edges[e_size] = Edge_t(old2new[u], old2new[v], false);
            edge2id.insert(edges[e_size], e_size);
            rank.push_back(g.edge_rank[i]);

            int sz = g.v_set[i].size();
            T_size[e_size] = sz;

            T[e_size] = new int[truss_num + 1];
            for (int j = 0; j < sz; j++) {
                int w = g.v_set[i][j];

                if (old2new[w] == -1) {
                    old2new[w] = (int) new2old.size();
                    new2old.push_back(w);
                }

                T[e_size][j] = old2new[w];
            }

            e_size++;
        }

        C = new int *[e_size];
        C_size = new int[e_size];

        for (int i = 0; i < e_size; i++) {

            C_size[i] = 0;
            int sz = T_size[i];
            C[i] = new int[sz * (sz - 1) / 2];

            for (int j = 0; j < T_size[i]; j++) {
                for (int k = j + 1; k < T_size[i]; k++) {

                    Edge_t e = Edge_t(T[i][j], T[i][k], false);
                    int idx = edge2id.exist(e);

                    if (idx != -1 && rank[idx] > rank[i]) {
                        C[i][C_size[i]++] = idx;
                    }
                }
            }
        }
    };

    v_size = new2old.size();

    printf("|V| = %d, |E| = %d\n", v_size, e_size);
    printf("Truss number = %d\n", truss_num - 2);

    //Free memory
    free_graph(&g);
    free(edgeIdToEdge);
    free(EdgeSupport);

}

double EBBkC_t::list_k_clique(const char *file_name) {
    double runtime;
    struct rusage start, end;
    EBBkC_Graph_t G, g;

    printf("Reading edges from %s ...\n", file_name);
    GetCurTime(&start);
    G.truss_decompose(file_name);

    printf("Building necessary data structure ...\n");
    G.build(false);
    G.clear_stored_cliques(); // Clear main storage before starting

    printf("Iterate over all cliques\n");

    g.truss_num = G.truss_num;
    g.build(true);
    for (int i = 0; i < G.e_size; i++) {
        G.branch(i, &g);
        g.EBBkC_plus_plus(K - 2, &N);

        // --- Start of Modification ---
        // Collect cliques from the subproblem 'g' after each branch.
        if (!g.stored_cliques.empty()) {
            G.stored_cliques.insert(G.stored_cliques.end(), g.stored_cliques.begin(), g.stored_cliques.end());
        }
        // --- End of Modification ---
    }

    GetCurTime(&end);
    runtime = GetTime(&start, &end);

    // After the loop, copy the fully collected cliques to the global vector.
    stored_cliques = G.get_stored_cliques();

    return runtime;
}

double EBBkC_t::list_k_clique_parallel(const char *file_name) {
    int i;
    double start, end;
    EBBkC_Graph_t G;

    printf("Reading edges from %s ...\n", file_name);
    start = omp_get_wtime();
    G.truss_decompose(file_name);

    printf("Building necessary data structure ...\n");
    G.build(false);
    G.clear_stored_cliques(); // Clear main storage before starting

    printf("Iterate over all cliques\n");

    // MODIFICATION: Removed 'g' and 'i' from the private clause.
    #pragma omp parallel reduction(+:N)
    {
        EBBkC_Graph_t g; // 'g' is declared here, making it thread-private.
        int n_edges = 0; // for thread-local counting
        g.truss_num = G.truss_num;
        g.build(true);
        #pragma omp for schedule(dynamic, 1) nowait
        for (i = 0; i < G.e_size; i++) {
            n_edges++;
            G.branch(i, &g);
            g.EBBkC_plus_plus(K - 2, &N);

            if(!g.stored_cliques.empty()) {
                omp_set_lock(&G.clique_lock);
                G.stored_cliques.insert(G.stored_cliques.end(), g.stored_cliques.begin(), g.stored_cliques.end());
                omp_unset_lock(&G.clique_lock);
                g.clear_stored_cliques();
            }
        }
        printf("Thread: %d, handled %d edges\n", omp_get_thread_num(), n_edges);
    }

    end = omp_get_wtime();

    // After parallel execution, copy the results to the global vector for printing
    stored_cliques = G.get_stored_cliques();

    return (double) (end - start) * 1e3;
}

// clique-storing functions

void EBBkC_Graph_t::add_to_current_clique(int vertex) {
    current_clique.push_back(vertex);
}

void EBBkC_Graph_t::remove_from_current_clique() {
    if (!current_clique.empty()) {
        current_clique.pop_back();
    }
}

void EBBkC_Graph_t::store_current_clique() {
    if (!current_clique.empty()) {
        omp_set_lock(&clique_lock); // Set lock before modifying shared vector
        stored_cliques.push_back(current_clique);
        omp_unset_lock(&clique_lock); // Unset the lock
    }
}

void EBBkC_Graph_t::store_clique_from_list(int* list, int list_size, int k) {
    std::vector<int> clique;
    for (int i = 0; i < k; i++) {
        if (list[i] < new2old.size()) {
            clique.push_back(new2old[list[i]]);
        }
    }
    if (!clique.empty()) {
        stored_cliques.push_back(clique);
    }
}

std::vector<std::vector<int>> EBBkC_Graph_t::get_stored_cliques() const {
    return stored_cliques;
}

void EBBkC_Graph_t::clear_stored_cliques() {
    stored_cliques.clear();
    current_clique.clear();
}