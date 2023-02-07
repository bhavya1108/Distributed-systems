#include <bits/stdc++.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

using namespace std;
// #define double long double

int mandelbrot(complex<double> c, int k)
{
    complex<double> z(0.0, 0.0);

    for (int i = 0; i < k; i++)
    {
        z = z * z + c;
        if (abs(z) > 2)
        {
            return 0;
        }
    }
    return 1;
}

void fxn(complex<double> buf[], int n, int k, int rank)
{
    for (int i = 0; i < n; i++)
    {
        buf[i] = complex<double>((mandelbrot(buf[i], k)), 0);
    }
}

vector<vector<int>> iterate_points(vector<vector<complex<double>>> arr, int k)
{
    vector<vector<int>> ans(arr.size(), vector<int>(arr[0].size()));
    int ctr = 0;
    for (int i = 0; i < arr.size(); ++i)
    {
        for (int j = 0; j < arr[0].size(); ++j)
        {
            ans[i][j] = mandelbrot(arr[i][j], k);
            ctr++;
        }
    }
    return ans;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int size, rank;
    srand(time(NULL));
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n, m, k, chunk, used_size, padded_size;

    if (rank == 0)
    {
        cin >> n >> m >> k;
        // make a grid of points

        chunk = n / size;
        used_size = size;
        padded_size = ((n * m / size) + (n * m % size == 0 ? 0 : 1)) * size;
        if (chunk == 0)
        {
            chunk = 1;
            used_size = n + 1;
        }

        for (int i = 1; i < used_size; i++)
        {
            MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&m, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&chunk, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&k, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        for (int i = used_size; i < size; i++)
        {
            int tmp = 0;
            MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&m, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&tmp, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&k, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&n, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&m, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&chunk, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&k, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        padded_size = ((n * m / size) + (n * m % size == 0 ? 0 : 1)) * size;
    }

    complex<double> sendBuf[padded_size];
    complex<double> recvBuf[padded_size / size];
    vector<vector<complex<double>>> pts(n, vector<complex<double>>(m, (0.0, 0.0)));

    if (rank == 0)
    {
        double x_step = 2.5 / (double)(m - 1), y_step = 2.0 / (double)(n - 1);
        for (double i = 0; i < n; ++i)
        {
            for (double j = 0; j < m; ++j)
            {
                pts[i][j] = complex<double>(-1.5 + (j)*x_step, 1.0 - (i)*y_step);
            }
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                sendBuf[i * m + j] = pts[i][j];
            }
        }
    }

    MPI_Scatter(sendBuf, padded_size / size, MPI_C_DOUBLE_COMPLEX, recvBuf, padded_size / size, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    fxn(recvBuf, padded_size / size, k, rank);
    MPI_Gather(recvBuf, padded_size / size, MPI_C_DOUBLE_COMPLEX, sendBuf, padded_size / size, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        // cout << endl
            //  << size << endl;

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                cout << real(sendBuf[i * m + j]) << " ";
            }
            cout << endl;
        }
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}