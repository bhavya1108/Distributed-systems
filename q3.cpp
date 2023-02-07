#include <stdlib.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <map>
#include <queue>
using namespace std;
long long int n = 0;

long long int mini(long long int a, long long int b, long long int ind, long long int& val) {
    if(a > b) {
        val = ind;
        return b;
    }
    return a;
}

void performDP(vector<vector<long long int> > dp, vector<vector<long long int> > root, long long int indices[], long long int root_to_send[], vector<long long int> prefix, long long int k) {
    for(long long int ind = 0 ; indices[ind] !=LLONG_MAX/5 and (ind < n-k) ; ind++){
            long long int j = indices[ind];
            long long int i = j-k;
            long long int sum = prefix[j+1] - prefix[i];
            indices[ind] = dp[i][j];
            indices[ind] = mini(indices[ind], sum + dp[i+1][j], i, root_to_send[ind]);
            indices[ind] = mini(indices[ind], sum + dp[i][j-1], j, root_to_send[ind]);
            for(long long int k = i+1 ; k < j ; k++)
                indices[ind] = mini(indices[ind], sum + dp[i][k-1] + dp[k+1][j], k, root_to_send[ind]);
    }
}

vector<long long int> merge(long long int keys[], long long int size, long long int chunk_size) {
    sort(keys, keys+n);
    vector<long long int> ans;
    for(long long int i = 0 ; i < n ; i++)
        ans.push_back(keys[i]);
    return ans;
}

void dfs(long long int node, long long int start, long long int end, vector<vector<long long int> > root, vector<long long int> &parents) {
    if(start == end)
        return;
    if(start < node) {
        parents[root[start][node-1]] = node;
        dfs(root[start][node-1], start, node-1, root, parents);
    }
    if(node < end) {
        parents[root[node+1][end]] = node;
        dfs(root[node+1][end], node+1, end, root, parents);
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int size, rank;
    srand(time(NULL));
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    vector<long long int> key, freq;
    if(rank == 0) {
        // take input in root process
        cin >> n;
        key.resize(n, 0);
        freq.resize(n, 0);
        for(long long int i = 0 ; i < n ; i++)
            cin >> key[i] >> freq[i];

        // replace below with a parallel sorting algo
    }
    MPI_Bcast(&n, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    if(rank != 0) {
        key.resize(n, 0);
        freq.resize(n, 0);
    }

    // merge sort using scatter gather
    long long int key_temp[(n)/size], to_send[n];
    for(long long int i = 0 ; i < n ; i++) to_send[i] = key[i];
    MPI_Scatter(to_send, (n)/size, MPI_LONG_LONG_INT, key_temp, (n)/size, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    sort(key_temp, key_temp + n/size);
    MPI_Gather(key_temp, (n)/size, MPI_LONG_LONG_INT, to_send, (n)/size, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    auto new_keys = merge(to_send, n, n/size);

    // preproccess the prefix sum for calculating sum from i to j in O(1)
    vector<long long int> prefix(n+1, 0);

    // make and initialize the dp and root matrix
    long long int dp[n][n], root[n][n];
    if(rank == 0) {
        vector<pair<long long int,long long int> > for_sort(n);
        for(long long int i = 0 ; i < n ; i++) {
            for_sort[i].first = key[i];
            for_sort[i].second = freq[i];
        }
        sort(for_sort.begin(), for_sort.end());
        for(long long int i = 0 ; i < n ; i++) {
            key[i] = for_sort[i].first;
            freq[i] = for_sort[i].second;
        }
    }
    
    MPI_Bcast(&key[0], n, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&freq[0], n, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

    for(long long int i = 1 ; i <= n ; i++)
        prefix[i] = prefix[i-1] + freq[i-1];

    for(long long int i = 0 ; i < n ; i++) {
        for(long long int j = 0 ; j < n ; j++) {
            dp[i][j] =LLONG_MAX/5;
            root[i][j] = -1;
        }
        root[i][i] = i;
        dp[i][i] = freq[i];
    }

    // divide the dp calculation among processes
    for(long long int k = 1 ; k < n ; k++) {
        long long int i = 0;
        long long int indices[n-k], dp_distr[n-k], dp_updates[n-k], root_to_send[n-k], root_updates[n-k];
        for(long long int ii = 0 ; ii < n-k ; ii++) {
            indices[ii] =LLONG_MAX/5;
            dp_updates[ii] =LLONG_MAX/5;
            dp_distr[ii] =LLONG_MAX/5;
        }

        MPI_Bcast(&dp[0][0], n*n, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&root[0][0], n*n, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        vector<vector<long long int> > temp_dp(n, vector<long long int> (n, 0)), temp_root(n, vector<long long int> (n, 0));
        for(long long int a = 0 ; a < n ; a++)
            for(long long int b = 0 ; b < n ; b++) {
                temp_dp[a][b] = dp[a][b];
                temp_root[a][b] = root[a][b];
            }
        for(long long int j = i + k ; j < n ; j++)
            indices[j-k] = j;
        MPI_Scatter(indices, (n-k)/size, MPI_LONG_LONG_INT, dp_distr, (n-k)/size, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        performDP(temp_dp, temp_root, dp_distr, root_to_send, prefix, k);
        MPI_Gather(dp_distr, (n-k)/size, MPI_LONG_LONG_INT, dp_updates, (n-k)/size, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(root_to_send, (n-k)/size, MPI_LONG_LONG_INT, root_updates, (n-k)/size, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        if(rank == 0) {
            if(size*((n-k)/size) < (n-k)) {
                for(long long int i = size*((n-k)/size) ; i < (n-k) ; i++)
                    indices[i-size*((n-k)/size)] = i+k;
                performDP(temp_dp, temp_root, indices, root_to_send, prefix, k);
                for(long long int i = size*((n-k)/size) ; i < (n-k) ; i++) {
                    dp[i][i+k] = indices[i-size*((n-k)/size)];
                    root[i][i+k] = root_to_send[i-size*((n-k)/size)];
                }
            }
            for(long long int i = 0 ; i < size*((n-k)/size) ; i++) {
                dp[i][i+k] = dp_updates[i];
                root[i][i+k] = root_updates[i];
            }
        }
    }

    // after finishing the matrix formation, just perform dfs traversal for parent array
    if(rank == 0){
        cout << dp[0][n-1] << endl;
        vector<long long int> parents(n, 0);
        vector<vector<long long int> > new_root(n, vector<long long int> (n, 0));
        for(long long int i = 0 ; i < n ; i++)
            for(long long int j = 0 ; j < n ; j++)
                new_root[i][j] = root[i][j];
        dfs(new_root[0][n-1], 0, n-1, new_root, parents);
        for(long long int i = 0 ; i < n ; i++)
            parents[i] = key[parents[i]];
        parents[new_root[0][n-1]] = 0;
        for(long long int i = 0 ; i < n ; i++)
            cout << parents[i] << " ";
        cout << endl;
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}