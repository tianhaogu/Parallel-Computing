#include <mpi.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <iomanip>

using namespace std;

double func(double x) {
    double y = x;
    for (int i = 1; i <= 10; ++i) {
        y += (sin (x * i) / pow(2.0, i));
    }
    return y;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        cout << "Missing m or n!!!" << endl;
        exit(1);
    }

    MPI_Init(&argc, &argv);

    int p;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int m = atoi(argv[1]), n = atoi(argv[2]);
    vector<vector<double> > A(m, vector<double>(n, 0));

    int p_side = sqrt(p);
    int each_row_size = ceil((double)m / p_side), each_col_size = ceil((double)n / p_side);
    int curr_row_idx = floor(rank / p_side), curr_col_idx = floor(rank % p_side);
    int up_idx = each_row_size * curr_row_idx, left_idx = each_col_size * curr_col_idx;
    int down_idx = min(m, each_row_size * (curr_row_idx + 1)) - 1, right_idx = min(n, each_col_size * (curr_col_idx + 1)) - 1;

    int curr_row_size = down_idx - up_idx + 1, curr_col_size = right_idx - left_idx + 1;
    vector<vector<double> > localMatrix(curr_row_size + 2, vector<double>(curr_col_size + 2, 0.0));
    for (int i = 1; i <= curr_row_size; ++i) {
        for (int j = 1; j <= curr_col_size; ++j) {
            int org_i = curr_row_idx * each_row_size + (i - 1), org_j = curr_col_idx * each_col_size + (j - 1);
            localMatrix[i][j] = org_i * sin(org_i) + org_j * cos(org_j) + sqrt(org_i + org_j);
        }
    }
    vector<vector<double> > prevMatrix(localMatrix);
    
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time;
    if (rank == 0) {
        start_time = MPI_Wtime();
    }

    for (int iter = 0; iter < 10; ++iter) {
        MPI_Status status;
        // send upper to lower
        if (curr_row_idx < p_side - 1) {
            MPI_Send(&localMatrix[curr_row_size][1], curr_col_size, MPI_DOUBLE, rank + p_side, 0, MPI_COMM_WORLD);
        }

        // send lower to upper
        if (curr_row_idx > 0) {
            MPI_Send(&localMatrix[1][1], curr_col_size, MPI_DOUBLE, rank - p_side, 0, MPI_COMM_WORLD);
        }

        // receive from upper
        if (curr_row_idx > 0) {
            MPI_Recv(&localMatrix[0][1], curr_col_size, MPI_DOUBLE, rank - p_side, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        //receive from lower
        if (curr_row_idx < p_side - 1) {
            MPI_Recv(&localMatrix[curr_row_size + 1][1], curr_col_size, MPI_DOUBLE, rank + p_side, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        vector<double> border_col_send(curr_row_size);

        // send left to right
        if (curr_col_idx < p_side - 1) {
            for (int i = 1; i <= curr_row_size; ++i) {
                border_col_send[i-1] = localMatrix[i][curr_col_size];
            }
            MPI_Send(&border_col_send[0], curr_row_size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }

        // send right to left
        if (curr_col_idx > 0) {
            for (int i = 1; i <= curr_row_size; ++i) {
                border_col_send[i-1] = localMatrix[i][1];
            }
            MPI_Send(&border_col_send[0], curr_row_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }

        vector<double> border_col_recv(curr_row_size);
        
        // receive from left
        if (curr_col_idx > 0) {
            MPI_Recv(&border_col_recv[0], curr_row_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 1; i <= curr_row_size; ++i) {
                localMatrix[i][0] = border_col_recv[i-1];
            }
        }

        // receive from right
        if (curr_col_idx < p_side - 1) {
            MPI_Recv(&border_col_recv[0], curr_row_size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 1; i <= curr_row_size; ++i) {
                localMatrix[i][curr_col_size + 1] = border_col_recv[i-1];
            }
        }

        // calculation begins
        prevMatrix = localMatrix;
        for (int i = 1; i <= curr_row_size; ++i) {
            for (int j = 1; j <= curr_col_size; ++j) {
                if ((up_idx == 0 && i == 1) || (down_idx == m - 1 && i == curr_row_size)
                        || (left_idx == 0 && j == 1) || (right_idx == n - 1 && j == curr_col_size)) {
                    continue;
                }
                double z = (
                    func(prevMatrix[i-1][j]) + func(prevMatrix[i+1][j]) + func(prevMatrix[i][j-1]) + func(prevMatrix[i][j+1]) + func(prevMatrix[i][j])
                ) / 5;
                localMatrix[i][j] = max(-100.0, min(100.0, z));
            }
        }
    }

    double curr_sum_all = 0.0, curr_sum_square_all = 0.0;
    for (int i = 1; i <= curr_row_size; ++i) {
        for (int j = 1; j <= curr_col_size; ++j) {
            curr_sum_all += localMatrix[i][j];
            curr_sum_square_all += pow(localMatrix[i][j], 2);
        }
    }

    double end_time;
    if (rank == 0) {
        double sum_all = 0.0, sum_square_all = 0.0;
        for (int i = 1; i < p; ++i) {
            MPI_Status status;
            double sum_i = 0.0, sum_square_i = 0.0;
            MPI_Recv(&sum_i, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&sum_square_i, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum_all += sum_i;
            sum_square_all += sum_square_i;
        }
        sum_all += curr_sum_all;
        sum_square_all += curr_sum_square_all;
        cout << "Sum of entries for " << "m = " << m << ", n = " << n << ", and p = " << p << " is: " << sum_all << endl;
        cout << "Sum of square of entries for " << "m = " << m << ", n = " << n << ", and p = " << p << " is: " << sum_square_all << endl;
        end_time = MPI_Wtime();
        cout << "Elapsed time for " << "m = " << m << ", n = " << n << ", and p = " << p << " is: " << end_time << " seconds" << endl;
    }
    else {
        MPI_Send(&curr_sum_all, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&curr_sum_square_all, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
