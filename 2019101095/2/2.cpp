#include <bits/stdc++.h>
#include <map>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

using namespace std;


int rows, columns, particle_count, time_limit;
map<pair<int, int>, vector<int> > particle_map;
map<char, pair<int, int> > movement_dict;

vector<pair<int, int> > collided_indices;


map<char, char> reflection_dict;
pair<pair<int, int>, char> move_particle(char direction, int x, int y) {
    int val1 = movement_dict[direction].first;
    int temp_x = x + val1;
    int val2 = movement_dict[direction].second;
    int temp_y = y + val2;
    if(temp_x < 0 || temp_x == rows || temp_y < 0 || temp_y == columns) {
        direction = reflection_dict[direction];
        temp_x = x + movement_dict[direction].first;
        temp_y = y + movement_dict[direction].second;
    }
    return {{temp_x,temp_y},direction};
}

void print_particles(vector<int> x_coords, vector<int> y_coords, vector<char> directions) {
     int sz = x_coords.size();
    for(int i = 0; i < sz; i++){
        cout << x_coords[i] << " ";
        cout << y_coords[i] << " ";
        cout << directions[i] << " ";
        cout<<endl;
    }
}
map<char, char> collision_dict;

vector<pair<int, int> > detect_collisions(int fl , vector<int> x_coords, vector<int> y_coords, vector<char> directions ) {
    fl=1;
    particle_map.clear();
    collided_indices.clear();
    for(int i = 0; i < particle_count; i++) {
        particle_map[{x_coords[i],y_coords[i]}].push_back(i);
    }
    for(auto i : particle_map) {
        int val = i.second.size();
        if(val == 2){
            if(directions[i.second[0]] == reflection_dict[directions[i.second[1]]]){
                 auto val2 = {i.second[0],i.second[1]};
                 collided_indices.push_back({i.second[0],i.second[1]});
            }
        }
    }
    return collided_indices;
}

void update_particles(vector<int> &x_coords, vector<int> &y_coords, vector<char> &directions) {
    int sz = x_coords.size();
    for(int i = 0; i < sz; i++) 
    {
        char dir = directions[i];
        int x =  x_coords[i];
        int y = y_coords[i];
        auto result = move_particle(dir,x,y);
        x_coords[i] = result.first.first;
        directions[i] = result.second;
        y_coords[i] = result.first.second;
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    srand(time(NULL));
    int num_procs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rank2 = rank;

    // movement_dict['U']= {0,1};
    movement_dict.insert({'U',{0,1}});
    movement_dict.insert({'D',{0,-1}});
    movement_dict.insert({'R',{1,0}});
    movement_dict.insert({'L',{-1,0}});

    reflection_dict.insert({'L','R'});
    reflection_dict.insert({'R','L'});
    reflection_dict.insert({'U','D'});
    reflection_dict.insert({'D','U'});

    collision_dict.insert({'L','D'});
    collision_dict.insert({'R','U'});
    collision_dict.insert({'U','L'});
    collision_dict.insert({'D','R'});

    if (rank == 0) 
    {   
        cin >> rows>> columns >> particle_count>>time_limit;
        int used_procs = num_procs;
        int step = particle_count/ num_procs;
        if (step == 0) {
            used_procs = particle_count+ 1;
            step = 1; 

        }
        
        vector<int> info = {rows, columns, particle_count, time_limit, step};
        vector<int> x(particle_count), y(particle_count);
        vector<char> dir(particle_count,'a');
        int i=0;
        while(i<particle_count){
            cin >> x[i]>>y[i]>>dir[i];
            i++;
        }
        int bh=0;
        i=1;
        while(i<used_procs){
            MPI_Send(&info[0], 5, MPI_INT, i, 0, MPI_COMM_WORLD);
            i++;
        }
        // cout<<bh<<endl;
        int bh2=0;
        info[4] = 0;
        i=used_procs;
        while(i<num_procs){
            MPI_Send(&info[0], 5, MPI_INT, i, 0, MPI_COMM_WORLD);
            i++;
        }
        // cout<<bh2<<endl;
        bh2=0;
        for (int j = 0; j < time_limit; j++) {
            int ans=0;
            int var1 = used_procs;
            int bh3=0;
            vector<int> x_zero, y_zero;
            int i=1;
            while(i<var1){
                int bh1=bh1+1;
                MPI_Send(&x[(i-1)*step], step, MPI_INT, i, 1, MPI_COMM_WORLD);
                bh2 = bh2+1;
                MPI_Send(&y[(i-1)*step], step, MPI_INT, i, 2, MPI_COMM_WORLD);
                bh3 = bh3+1;
                MPI_Send(&dir[(i-1)*step], step, MPI_CHAR, i, 3, MPI_COMM_WORLD);
                bh1++;
                i++;
            }
            // cout<<bh3<<endl;
            vector<char> dir_zero;

            int len_left = particle_count - (num_procs-1)*step;
            if(particle_count >= num_procs){
                for(int i = (num_procs-1)*step ; i < particle_count ; i++) {
                    x_zero.push_back(x[i]);
                }
                for(int i = (num_procs-1)*step ; i < particle_count ; i++) {
                    y_zero.push_back(y[i]);
                }for(int i = (num_procs-1)*step ; i < particle_count ; i++) {
                    dir_zero.push_back(dir[i]);
                }
                update_particles(x_zero, y_zero, dir_zero);
            }
            int bh4=0;
            i=1;
            while(i<used_procs){
                MPI_Recv(&x[(i-1)*step], step, MPI_INT, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&y[(i-1)*step], step, MPI_INT, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&dir[(i-1)*step], step, MPI_INT, i, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                i++;
            }
            // cout<<bh4<<endl;
            int l = 0;
            if(particle_count >= num_procs) {
                int var = (num_procs-1)*step;
                for(int i = var ; i < particle_count ; i++) {
                    x[i] = x_zero[l];
                    l++;
                }
                l=0;
                var =  (num_procs-1)*step;
                for(int i = var ; i < particle_count ; i++) {
                    y[i] = y_zero[l];
                    l++;
                } 
                l=0;
                var =  (num_procs-1)*step;
                for(int i = var ; i < particle_count ; i++) {
                    dir[i] = dir_zero[l];
                    l++;
                }
                l = 0;
                var =  (num_procs-1)*step;
            }

            if(j == time_limit-1){
               print_particles(x, y, dir);
               ans++;
            }
            auto indices_colliding =  detect_collisions(0,x, y, dir);
            int szz=indices_colliding.size();
            for(int i = 0 ; i < szz; i++) {
                auto vz = dir[indices_colliding[i].first];
                dir[indices_colliding[i].first] = collision_dict[vz];
                auto vz2 = dir[indices_colliding[i].second];
                dir[indices_colliding[i].second] = collision_dict[vz2];
            }
            
        }
    }
    else {
        vector<int> info(5);

        int step;
        MPI_Recv(&info[0], 5, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        step = info[4];
        vector<int> x(step, 0), y(step, 0);

        rows = info[0];
        columns = info[1];
        particle_count = info[2];
        vector<char> dir(step, 'a');

        time_limit = info[3];
        int flag =0;
        if(step > 0) {
            flag = 0;
            for(int j = 0 ; j < time_limit ; j++) {
                MPI_Recv(&x[0], step, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(flag==1) continue;
                MPI_Recv(&y[0], step, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&dir[0], step, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(flag==1) continue;
                update_particles(x, y, dir);
                MPI_Send(&x[0], step, MPI_INT, 0, 4, MPI_COMM_WORLD);
                MPI_Send(&y[0], step, MPI_INT, 0, 5, MPI_COMM_WORLD);
                if(flag==1) continue;
                MPI_Send(&dir[0], step, MPI_INT, 0, 6, MPI_COMM_WORLD);

            }
        }
        if(flag==1){
            flag = 0;
        }
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}

        

