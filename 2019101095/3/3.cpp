#include<bits/stdc++.h>
using namespace std;
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

void dfs(vector<vector<long long int> > root, vector<long long int> &parents,long long int node, long long int start, long long int end ) {
    int flag=0;
    if(start == end){
        start = end;
        return;
    }
    if(flag==1){
        return;
    }
    // cout<<flag<<endl;
    if(start < node)
    {
        long long int x = root[start][node-1];
        parents[x] = node;
        dfs( root, parents ,x, start, node-1);
    }
    if(node < end) {
        long long int x = root[node+1][end];
        parents[x] = node;
        dfs( root, parents,x, node+1, end);
    }
}

long long int minimuni(long long int a, long long int b, long long int ind, long long int& val) {
    int fl=0;
    if(a > b) 
    {
        val = ind;
        if(fl==0)
            return b;
    }
    return a;
}
long long int n = 0;

void dp_computation( long long int k, vector<vector<long long int> > dp, vector<vector<long long int> > root, long long int ind_arr[], long long int root_send[], vector<long long int> prefix) {
    long long int ind=0;
    long long int i;
    long long int j;
    long long int sum;
    while(ind<(n-k) && ind_arr[ind] !=LLONG_MAX/5){
        j = ind_arr[ind];
        i = j-k;
        sum = prefix[j+1] - prefix[i];
        ind_arr[ind] = dp[i][j];
        long long int v1 = sum + dp[i+1][j];
        ind_arr[ind] = minimuni(ind_arr[ind], v1, i, root_send[ind]);
        long long int v2 = sum + dp[i][j-1];
        ind_arr[ind] = minimuni(ind_arr[ind], v2, j, root_send[ind]);
        long long int k = i+1;
        while(k<j){
            long long int x = sum + dp[i][k-1] + dp[k+1][j];
            ind_arr[ind] = minimuni(ind_arr[ind], x, k, root_send[ind]);
            k++;
        }
        ind++;
    }
}

vector<long long int> merge(long long int keys[], long long int size, long long int chunk_size)
{
    long long int i=0;
    vector<long long int> res(n,0);
    sort(keys, keys+n);
    while(i<n){
        res[i] = keys[i];
        i++;
    }
    return res;
}



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int process_count, process_id;
    srand(time(NULL));
    vector<long long int> key, freq;
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);

    if(process_id == 0) {
        cin >> n;
        for(long long int i = 0 ; i < n ; i++){
            long long int v1;
            long long int v2;
            cin >> v1 >> v2;
            key.push_back(v1);
            freq.push_back(v2);
        }
    }
    MPI_Bcast(&n, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    if(process_id != 0) {
        key.resize(n, 0);
        freq.resize(n, 0);
    }
    vector<long long int> prefix(n+1, 0);
    // merge sort using scatter gather
    long long int temp_keys[(n)/process_count];
    long long int to_send[n];
    long long int i=0;
    while(i<n){
        to_send[i] = key[i];
        i++;
    }
    long long int dp[n][n], root[n][n];

    MPI_Scatter(to_send, (n)/process_count, MPI_LONG_LONG_INT, temp_keys, (n)/process_count, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    long long int vae = n / process_count;
    sort(temp_keys, temp_keys + vae);
    long long int fl1;
    MPI_Gather(temp_keys, (n)/process_count, MPI_LONG_LONG_INT, to_send, (n)/process_count, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    fl1=0;
    auto new_keys = merge(to_send, n, n/process_count);
    if(process_id == 0)
    {
        long long int i=0;
        vector<pair<long long int,long long int> > for_sort;
        while(i<n){
            for_sort.push_back({key[i],freq[i]});
            i++;
        }
        sort(for_sort.begin(), for_sort.end());
        i=0;
        while(i<n){
            key[i] = for_sort[i].first;
            freq[i] = for_sort[i].second;
            i++;
        }
    }
    MPI_Bcast(&key[0], n, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    fl1=1;
    MPI_Bcast(&freq[0], n, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    i=1;
    while(i<=n){
        prefix[i] = prefix[i-1] + freq[i-1];
        i++; 
    }
    i=0;
    while(i<n){
        long long int j=0;
        while(j<n){
            dp[i][j] =LLONG_MAX/5;
            root[i][j] = -1;
            j++;
        }
        i++;
    }
    i=0;
    while(i<n){
        root[i][i] = i;
        dp[i][i] = freq[i];
        i++;
    }
    // divide the dp calculation among processes
    long long int k = 1;
    // for(long long int k = 1 ; k < n ; k++)
    while(k<n)
    {
        long long int ind_arr[n-k]; 
        long long int dp_distr[n-k];
        long long int new_dp[n-k];
        long long int root_send[n-k];
        long long int new_root[n-k];
        long long int ii=0;
        long long int i = 0;
        vector<vector<long long int> > temp_dp(n, vector<long long int> (n, 0)); 

        while (ii<(n-k))
        {
            ind_arr[ii] =LLONG_MAX/5;
            new_dp[ii] =LLONG_MAX/5;
            dp_distr[ii] =LLONG_MAX/5;
            ii++;
        }
        MPI_Bcast(&dp[0][0], n*n, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        vector<vector<long long int> > root_tmp(n, vector<long long int> (n, 0));
        MPI_Bcast(&root[0][0], n*n, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        long long int a =0;
        while(a<n){
            long long int b=0;
            while(b<n){
                long long int var = dp[a][b];
                temp_dp[a][b] = var;
                long long int var2 = root[a][b];
                root_tmp[a][b] = var2;
                b++;
            }
            a++;
        }
        long long int j = i+k;
        while(j<n){
            ind_arr[j-k] = j;
            j++;
        }
        long long int tmp = (n-k)/process_count;   
        MPI_Scatter(ind_arr, tmp, MPI_LONG_LONG_INT, dp_distr, (n-k)/process_count, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        dp_computation(k,temp_dp, root_tmp, dp_distr, root_send, prefix);
        MPI_Gather(dp_distr, tmp, MPI_LONG_LONG_INT, new_dp, (n-k)/process_count, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(root_send, (n-k)/process_count, MPI_LONG_LONG_INT, new_root, (n-k)/process_count, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        long long int fl;
        if(process_id == 0)
        {
            long long int bh1 = process_count*((n-k)/process_count);
            if(bh1 < (n-k)) {
                long long int i=process_count*((n-k)/process_count);
                while(i<(n-k)){
                    long long int val1 = i-process_count*((n-k)/process_count);
                    ind_arr[val1] = i+k;
                    i++;
                }   
                dp_computation(k,temp_dp, root_tmp, ind_arr, root_send, prefix);
                i = process_count*((n-k)/process_count);
                while(i<(n-k)){
                    long long int val2 = i-process_count*((n-k)/process_count);
                    long long int val3 = ind_arr[val2];
                    dp[i][i+k] = val3;
                    long long int val4 = root_send[val2] ;
                    root[i][i+k] = val4;
                    i++;
                }
            }
            long long int i=0;
            while(i<(process_count*((n-k)/process_count))){
                long long int val1 = new_dp[i];
                dp[i][i+k] = val1;
                long long int val2 = new_root[i];
                root[i][i+k] = val2;
                i++;
            }
        }
        k++;
    }

    // after finishing the matrix formation, just perform dfs traversal for parent array
    if(process_id == 0){
        vector<vector<long long int> > new_root(n, vector<long long int> (n, 0));
        cout << dp[0][n-1] << endl;
        long long int i=0;
        while(i<n){
            long long int j=0;   
            while(j<n){
                new_root[i][j] = root[i][j];
                j++;
            }
            i++;
        }
        vector<long long int> parents(n, 0);
        dfs(new_root, parents,new_root[0][n-1], 0, n-1);
        i=0;
        while(i<n){
            long long int val = key[parents[i]];
            parents[i] = val;
            i++;
        }
        long long int val2 = new_root[0][n-1];
        parents[val2] = 0;
        i=0;
        while(i<n){
            cout << parents[i] << " ";
            i++;
        }
        cout << endl;
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}