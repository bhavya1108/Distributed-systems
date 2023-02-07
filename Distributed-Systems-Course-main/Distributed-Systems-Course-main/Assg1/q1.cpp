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
// #define double long double
int n, m, k, t;

int mandelbrot(double a, double b, int k, int rank) {
    double x = 0.0, y = 0.0;
    if(a*a + b*b > 4)
        return 0;
    for(int i = 0 ; i < k ; i++) {
        double temp = x;
        x = x*x - y*y + a;
        y = 2*temp*y + b;
        if(x*x + y*y > 4){
            return 0;
        }
    }
    if(x*x + y*y <= 4)
        return 1;
    return 0;
}

vector<vector<int> > iterate_points(vector<vector<double> > pts_x, vector<vector<double> > pts_y, vector<int> info, int rank) {
    vector<vector<int> > ans(pts_x.size(), vector<int> (pts_x[0].size(), 0));
    for(int i = 0 ; i < pts_x.size() ; i++)
        for(int j = 0 ; j < pts_x[0].size() ; j++)
            ans[i][j] = mandelbrot(pts_x[i][j], pts_y[i][j], info[2], rank);
    return ans;
}

void print_points(vector<vector<int> > ans) {
    for(int j = 0 ; j < ans[0].size() ; j++) {
        for(int i = 0 ; i < ans.size() ; i++) 
            cout << ans[i][j] << " ";
        cout << endl;
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int size, rank;
    srand(time(NULL));
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // - do this till the specified amount of time
    if (rank == 0) {   
        cin >> n >> m >> k;
        
        // make a grid of points
        vector<vector<double > > pts_x(n, vector<double>  (m, 0.0));
        vector<vector<double > > pts_y(n, vector<double>  (m, 0.0));
        double x_step = 2.5/(double)(n), y_step = 2.0/(double)(m);
        for(double i = 0 ; i < n ; i += 1) {
            double x = -1.5 + (i)*x_step;
            for(double j = 0 ; j < m ; j += 1) {
                pts_x[i][j] = x;
                pts_y[i][j] = -1.0 + (j)*y_step;
            }
        }
        
        int chunk = n/size, used_size = size;
        if(chunk == 0) {
            chunk = 1;
            used_size = n+1;
        }
        vector<int> info(3);
        info[0] = m;
        info[1] = chunk;
        info[2] = k;
        for(int i = 1 ; i < used_size ; i++)
            MPI_Send(&info[0], 3, MPI_INT, i, 0, MPI_COMM_WORLD);
        info[1] = 0;
        for(int i = used_size ; i < size ; i++)
            MPI_Send(&info[0], 3, MPI_INT, i, 0, MPI_COMM_WORLD);
        info[1] = chunk;
        // divide the elements of the grid  
        for(int i = 1 ; i < used_size ; i++) {
            for(int  j = 0 ; j < chunk ; j++) {
                MPI_Send(&pts_x[(i-1)*chunk + j][0], m, MPI_DOUBLE, i, j+1, MPI_COMM_WORLD);
                MPI_Send(&pts_y[(i-1)*chunk + j][0], m, MPI_DOUBLE, i, chunk + j+1, MPI_COMM_WORLD);
            }  
        }
        vector<vector<double> > pts1_x, pts1_y;
        if(n >= size) {
            pts1_x.resize(n - (used_size-1)*chunk);
            pts1_y.resize(n - (used_size-1)*chunk);
            for(int i = (used_size-1)*chunk ; i < n ; i++) {
                for(int j = 0 ; j < m ; j++) {
                    pts1_x[i-(used_size-1)*chunk].push_back(pts_x[i][j]);
                    pts1_y[i-(used_size-1)*chunk].push_back(pts_y[i][j]);
                }
            }
        }
        vector<vector<int> > anss;
        if(n >= size) {    
            anss = iterate_points(pts1_x, pts1_y, info, rank);
        }
        // recombine the results
        vector<vector<int> > ans(n, vector<int> (m, 0));
        

        for(int i = 1 ; i < used_size ; i++) {
            for(int j = 0 ; j < chunk ; j++)
                MPI_Recv(&ans[(i-1)*chunk + j][0], m, MPI_INT, i, j+10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if(n >= size) {    
            for(int i = (used_size-1)*chunk ; i < n ; i++) {
                for(int j = 0 ; j < m ; j++) {
                    ans[i][j] = anss[i-(used_size-1)*chunk][j];
                }
            }
        }
        print_points(ans);
    }
    else {
        vector<int> info(3);
        MPI_Recv(&info[0], 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(info[1] > 0) {
            int chunk = info[1];
            m = info[0];
            vector<vector<double> > pts_x(chunk, vector<double> (m, 0.0));
            vector<vector<double> > pts_y(chunk, vector<double> (m, 0.0));
            for(int i = 0 ; i < chunk ; i++) {
                MPI_Recv(&pts_x[i][0], m, MPI_DOUBLE, 0, i+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&pts_y[i][0], m, MPI_DOUBLE, 0, chunk + i+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            auto ans = iterate_points(pts_x, pts_y, info, rank);
            for(int i = 0 ; i < chunk ; i++) {
                MPI_Send(&ans[i][0], m, MPI_INT, 0, i+10, MPI_COMM_WORLD);
            }
        }
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}