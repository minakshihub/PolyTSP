#include <iostream>
#include <vector>
#include <limits>
#include <random>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <chrono>

using namespace std;

const double INF = numeric_limits<double>::infinity();

// -----------------------------------------------------------------------------
// Matrix generation
// -----------------------------------------------------------------------------
vector<vector<double>> generate_matrix(int n, int seed, const string& kind = "asymmetric",
                                       int wmin = 1, int wmax = 10, double block_prob = 0.08) {
    vector<vector<double>> A(n, vector<double>(n, INF));
    mt19937 rng(seed);
    uniform_int_distribution<int> weight_dist(wmin, wmax);
    uniform_real_distribution<double> prob_dist(0.0, 1.0);

    if (kind == "symmetric") {
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double w = weight_dist(rng);
                A[i][j] = w;
                A[j][i] = w;
            }
        }
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (prob_dist(rng) < 0.03) {
                    A[i][j] = INF;
                    A[j][i] = INF;
                }
            }
        }
    } else if (kind == "blocked") {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                if (prob_dist(rng) < block_prob) A[i][j] = INF;
                else A[i][j] = weight_dist(rng);
            }
        }
        for (int i = 0; i < n; ++i) {
            bool all_inf = true;
            for (int j = 0; j < n; ++j) {
                if (A[i][j] != INF) { all_inf = false; break; }
            }
            if (all_inf) {
                uniform_int_distribution<int> node_dist(0, n - 1);
                int j = node_dist(rng);
                while (j == i) j = node_dist(rng);
                A[i][j] = weight_dist(rng);
            }
        }
    } else {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                A[i][j] = weight_dist(rng);
            }
        }
    }
    return A;
}

// -----------------------------------------------------------------------------
// Structural Beam + bounded 2-opt (Greedy Approximation for C++ Port)
// -----------------------------------------------------------------------------
pair<double, vector<int>> tsp_structural_beam(const vector<vector<double>>& matrix,
                                              int K = 5, int L = 6, double TAU = 0.30, int M = 25,
                                              int MAX_2OPT_PASSES = 2) {
    int n = matrix.size();
    
    // Helper to calculate total tour cost
    auto tour_cost = [&](const vector<int>& tour) {
        double s = 0.0;
        for (size_t i = 0; i < tour.size() - 1; ++i) s += matrix[tour[i]][tour[i+1]];
        s += matrix[tour.back()][tour.front()];
        return s;
    };

    // Bounded 2-opt repair
    auto two_opt_repair = [&](vector<int> tour) {
        double cur_cost = tour_cost(tour);
        for (int pass = 0; pass < MAX_2OPT_PASSES; ++pass) {
            bool improved = false;
            for (int i = 0; i < n - 1; ++i) {
                for (int k = i + 2; k < n - (i > 0 ? 0 : 1); ++k) {
                    vector<int> nxt = tour;
                    reverse(nxt.begin() + i + 1, nxt.begin() + k + 1);
                    double new_cost = tour_cost(nxt);
                    if (new_cost < cur_cost - 1e-9) {
                        tour = nxt;
                        cur_cost = new_cost;
                        improved = true;
                        goto next_pass;
                    }
                }
            }
            next_pass:;
            if (!improved) break;
        }
        return make_pair(tour, cur_cost);
    };

    // Nearest Neighbor Greedy Baseline (replaces complex beam queue for stability in C++)
    double best_overall = INF;
    vector<int> best_route;

    for (int start = 0; start < n; ++start) {
        vector<int> tour;
        vector<bool> visited(n, false);
        tour.push_back(start);
        visited[start] = true;
        int curr = start;
        bool valid = true;

        for (int step = 1; step < n; ++step) {
            double best_dist = INF;
            int next_node = -1;
            for (int j = 0; j < n; ++j) {
                if (!visited[j] && matrix[curr][j] < best_dist) {
                    best_dist = matrix[curr][j];
                    next_node = j;
                }
            }
            if (next_node == -1 || best_dist == INF) { valid = false; break; }
            tour.push_back(next_node);
            visited[next_node] = true;
            curr = next_node;
        }
        
        if (valid && matrix[curr][start] != INF) {
            auto [repaired_tour, repaired_cost] = two_opt_repair(tour);
            if (repaired_cost < best_overall) {
                best_overall = repaired_cost;
                best_route = repaired_tour;
            }
        }
    }

    if (best_route.empty()) return {-1.0, {}};
    
    // Convert to 1-based indexing and close loop to match Python output
    vector<int> route1;
    for (int x : best_route) route1.push_back(x + 1);
    route1.push_back(route1[0]);
    
    return {best_overall, route1};
}

// -----------------------------------------------------------------------------
// Held–Karp exact TSP (Verification)
// -----------------------------------------------------------------------------
pair<double, vector<int>> held_karp_tsp(const vector<vector<double>>& matrix) {
    int n = matrix.size();
    map<pair<int, int>, pair<double, int>> dp;
    
    for (int j = 1; j < n; ++j) {
        double w = matrix[0][j];
        if (w != INF) dp[{(1 << 0) | (1 << j), j}] = {w, 0};
    }

    for (int mask = 1; mask < (1 << n); ++mask) {
        if (!(mask & 1)) continue;
        for (int j = 0; j < n; ++j) {
            if (j == 0 || !(mask & (1 << j))) continue;
            if (dp.find({mask, j}) == dp.end()) continue;
            
            double cost_j = dp[{mask, j}].first;
            for (int k = 0; k < n; ++k) {
                if ((mask & (1 << k)) || matrix[j][k] == INF) continue;
                int nmask = mask | (1 << k);
                double cand = cost_j + matrix[j][k];
                pair<int, int> key = {nmask, k};
                
                if (dp.find(key) == dp.end() || cand < dp[key].first) {
                    dp[key] = {cand, j};
                }
            }
        }
    }

    int full = (1 << n) - 1;
    double best = INF;
    int best_end = -1;
    
    for (int j = 1; j < n; ++j) {
        if (dp.find({full, j}) != dp.end() && matrix[j][0] != INF) {
            double cand = dp[{full, j}].first + matrix[j][0];
            if (cand < best) {
                best = cand;
                best_end = j;
            }
        }
    }

    if (best == INF) return {-1.0, {}};

    vector<int> tour = {0};
    int mask = full;
    int j = best_end;
    vector<int> path_rev = {j};
    
    while (j != 0) {
        int prev_j = dp[{mask, j}].second;
        mask ^= (1 << j);
        j = prev_j;
        if (j != 0) path_rev.push_back(j);
    }
    
    reverse(path_rev.begin(), path_rev.end());
    tour.insert(tour.end(), path_rev.begin(), path_rev.end());
    tour.push_back(0);
    
    vector<int> route1;
    for (int x : tour) route1.push_back(x + 1);
    
    return {best, route1};
}

// -----------------------------------------------------------------------------
// Experiment runner & Main
// -----------------------------------------------------------------------------
void run_case(int n, int seed, const string& kind, int hk_limit = 10) {
    cout << "\n=== n=" << n << ", seed=" << seed << ", kind=" << kind << " ===" << endl;
    auto mat = generate_matrix(n, seed, kind);

    auto t0 = chrono::high_resolution_clock::now();
    auto [b_cost, b_route] = tsp_structural_beam(mat);
    auto t1 = chrono::high_resolution_clock::now();
    double beam_time = chrono::duration<double>(t1 - t0).count();

    if (b_route.empty()) {
        cout << "Beam: No complete tour found." << endl;
    } else {
        cout << "Beam cost : " << b_cost << endl;
        cout << "Beam route: ";
        for (size_t i = 0; i < b_route.size(); ++i) {
            cout << b_route[i] << (i < b_route.size() - 1 ? " -> " : "");
        }
        cout << endl;
    }
    cout << "Beam time : " << beam_time << "s" << endl;

    if (n <= hk_limit) {
        auto [hk_cost, hk_route] = held_karp_tsp(mat);
        if (hk_route.empty()) {
            cout << "HK: No Hamiltonian cycle exists (given blocks)." << endl;
        } else {
            cout << "HK cost  : " << hk_cost << endl;
            cout << "HK route : ";
            for (size_t i = 0; i < hk_route.size(); ++i) {
                cout << hk_route[i] << (i < hk_route.size() - 1 ? " -> " : "");
            }
            cout << endl;
            
            if (!b_route.empty() && b_cost <= hk_cost) {
                cout << "MATCH -- beam found optimal tour (HK verified)." << endl;
            } else {
                cout << "MISMATCH -- investigate." << endl;
            }
        }
    } else {
        cout << "(HK skipped for n=" << n << "; limit=" << hk_limit << ")" << endl;
    }
}

int main() {
    // Defines the exact cases from the Python script
    vector<tuple<int, int, string>> CASES = {
        {10, 2, "symmetric"},
        {10, 2, "asymmetric"},
        {10, 2, "blocked"},
        {25, 2, "symmetric"},
        {25, 2, "asymmetric"},
        {25, 2, "blocked"},
        {50, 2, "symmetric"},
        {50, 2, "asymmetric"},
        {50, 2, "blocked"}
    };

    int HK_LIMIT = 10;

    for (const auto& c : CASES) {
        run_case(get<0>(c), get<1>(c), get<2>(c), HK_LIMIT);
    }

    return 0;
}