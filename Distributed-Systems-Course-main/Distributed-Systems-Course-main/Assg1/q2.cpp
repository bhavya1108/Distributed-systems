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
using namespace std;

map<char, pair<int,int> > dict_dir;
map< char, char> dict_rebound;
map< char, char> dict_collide;
int n, m, k, t;

pair<pair<int,int>, char> move_util(char dir, int x, int y) {
    int temp_x = x + dict_dir[dir].first;
    int temp_y = y + dict_dir[dir].second;
    if(temp_x < 0 or temp_x == n or temp_y < 0 or temp_y == m) {
        temp_x = x + dict_dir[dict_rebound[dir]].first;
        temp_y = y + dict_dir[dict_rebound[dir]].second;
        dir = dict_rebound[dir];
    }
    return make_pair(make_pair(temp_x, temp_y), dir);
}

// - let each process calc the next positions of each point and report points after each second
void move_points(vector<int> &x, vector<int> &y, vector<char> &dirn) {
    for(int i = 0 ; i < x.size() ; i++) {
        auto latest = move_util(dirn[i], x[i], y[i]);
        x[i] = latest.first.first;
        y[i] = latest.first.second;
        dirn[i] = latest.second;
    }
}

void print_points(vector<int> inp_x, vector<int> inp_y, vector<char> inp_d) {
    for(int i = 0 ; i < inp_x.size() ; i++)
        cout << inp_x[i] << " " << inp_y[i] << " " << inp_d[i] << endl;
}

// - If points are colliding, resolve this collision
vector<pair<int,int> > check_collision(vector<int> x, vector<int> y, vector<char> dirn) {
    map<pair<int,int>, vector<int> > vis;
    vector<pair<int,int> > indices;
    for(int i = 0 ; i < k ; i++) {
        vis[make_pair(x[i], y[i])].push_back(i);
    }
    for(auto i : vis) {
        if(i.second.size() == 2 and dirn[i.second[0]] == dict_rebound[dirn[i.second[1]]]) 
            indices.push_back(make_pair(i.second[0], i.second[1]));            
        
    }
    return indices;
}

int main(int argc, char *argv[])
{
    dict_dir['U'] = make_pair(0, 1);
    dict_dir['D'] = make_pair(0, -1);
    dict_dir['R'] = make_pair(1, 0);
    dict_dir['L'] = make_pair(-1, 0);
    dict_rebound['L'] = 'R';
    dict_rebound['R'] = 'L';
    dict_rebound['U'] = 'D';
    dict_rebound['D'] = 'U';
    dict_collide['L'] = 'D';
    dict_collide['R'] = 'U';
    dict_collide['U'] = 'L';
    dict_collide['D'] = 'R';
    MPI_Init(&argc, &argv);
    int size, rank;
    srand(time(NULL));
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // - do this till the specified amount of time
    if (rank == 0) {   
        cin >> n >> m >> k >> t;
        int step = k/size;
        int used_size = size;
        if(step == 0) {
            step = 1;
            used_size = k+1;
        }
        
        vector<int> info(5), inp_x(k, 0), inp_y(k, 0);
        info[0] = n;
        info[1] = m;
        info[2] = k;
        info[3] = t;
        info[4] = step;
        vector<char> inp_d(k, 'a');
        for(int i = 0 ; i < k ; i++)
            cin >> inp_x[i] >> inp_y[i] >> inp_d[i];

        for(int i = 1 ; i < used_size ; i += 1)
            MPI_Send(&info[0], 5, MPI_INT, i, 0, MPI_COMM_WORLD);
        info[4] = 0;
        for(int i = used_size ; i < size ; i++) 
            MPI_Send(&info[0], 5, MPI_INT, i, 0, MPI_COMM_WORLD);

        // - divide the k points to each process
        for(int j = 0 ; j < t ; j++) {
            for(int i = 1 ; i < used_size ; i += 1) {
                MPI_Send(&inp_x[(i-1)*step], step, MPI_INT, i, 1, MPI_COMM_WORLD);
                MPI_Send(&inp_y[(i-1)*step], step, MPI_INT, i, 2, MPI_COMM_WORLD);
                MPI_Send(&inp_d[(i-1)*step], step, MPI_CHAR, i, 3, MPI_COMM_WORLD);
            }
            int remaining_len = k - (size-1)*step;
            vector<int> inp_x_0, inp_y_0;
            vector<char> inp_d_0;
            if(k >= size){
                for(int i = (size-1)*step ; i < k ; i++) {
                    inp_x_0.push_back(inp_x[i]);
                    inp_y_0.push_back(inp_y[i]);
                    inp_d_0.push_back(inp_d[i]);
                }
                move_points(inp_x_0, inp_y_0, inp_d_0);
            }
            for (int i = 1; i < used_size; i++) {
                MPI_Recv(&inp_x[(i-1)*step], step, MPI_INT, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&inp_y[(i-1)*step], step, MPI_INT, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&inp_d[(i-1)*step], step, MPI_INT, i, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            int l = 0;
            if(k >= size) {
                for(int i = (size-1)*step ; i < k ; i++) {
                    inp_x[i] = inp_x_0[l];
                    inp_y[i] = inp_y_0[l];
                    inp_d[i] = inp_d_0[l];
                    l++;
                }
            }
            if(j == t-1)
                print_points(inp_x, inp_y, inp_d);
            auto indices_colliding = check_collision(inp_x, inp_y, inp_d);
            for(int i = 0 ; i < indices_colliding.size() ; i++) {
                inp_d[indices_colliding[i].first] = dict_collide[inp_d[indices_colliding[i].first]];
                inp_d[indices_colliding[i].second] = dict_collide[inp_d[indices_colliding[i].second]];
            }
            
        }
    }
    else {
        int step;
        vector<int> info(5);
        MPI_Recv(&info[0], 5, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        n = info[0];
        m = info[1];
        k = info[2];
        t = info[3];
        step = info[4];
        vector<int> inp_x(step, 0), inp_y(step, 0);
        vector<char> inp_d(step, 'a');
        if(step > 0) {
            for(int j = 0 ; j < t ; j++) {
                MPI_Recv(&inp_x[0], step, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&inp_y[0], step, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&inp_d[0], step, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                move_points(inp_x, inp_y, inp_d);
                MPI_Send(&inp_x[0], step, MPI_INT, 0, 4, MPI_COMM_WORLD);
                MPI_Send(&inp_y[0], step, MPI_INT, 0, 5, MPI_COMM_WORLD);
                MPI_Send(&inp_d[0], step, MPI_INT, 0, 6, MPI_COMM_WORLD);
            }
        }
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}