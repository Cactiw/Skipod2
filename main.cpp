#include <iostream>
#include <vector>
#include <ctime>
#include <mpi.h>
#include <random>
#include <cmath>
#include <chrono>

struct MyPoint {
    double x;
    double y;
    double z;
};
std::ostream &operator<<(std::ostream &os, MyPoint const &m) {
    return os << "(" << m.x << ", " << m.y << ", " << m.z << ")";
}

double volume() {
    return 1/6.;
}
//double volume() {
//    return (1 - 0) * (1 - 0) * (1 - 0);
//}

double count_analytic_integral() {
    return log(2) / 2. - 5/16.;
}


bool check_point(double x, double y, double z) {
    return x + y + z < 1;
}

double count_function_value(double x, double y, double z) {
    if (!check_point(x, y, z)) {
        return 0;
    }
    return 1 / pow((1 + x + y + z), 3);
}

void generate_points(double * points, long long count, std::uniform_real_distribution<double> unif, std::default_random_engine re) {
    for (long long i = 0; i < count * 3; i += 3) {
        bool true_point = false;

        while(!true_point) {
            points[i] = unif(re);
            points[i + 1] = unif(re);
            points[i + 2] = unif(re);

            if (points[i + 2] + points[i + 1] + points[i] > 1) {
                true_point = false;
            } else {
                true_point = true;
            }
        }
    }
}

int main(int argc, char **argv) {
    int errCode;

    std::uniform_real_distribution<double> unif(0, 1);
    std::default_random_engine re;
    re.seed(std::chrono::system_clock::now().time_since_epoch().count());

    if ((errCode = MPI_Init(&argc, &argv)) != MPI_SUCCESS) {
        std::cerr << "MPI_Init error: code " << errCode << std::endl;
        MPI_Abort(MPI_COMM_WORLD, errCode);
    }
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    long long workers_count = size - 1;
    long long points_count = 10000000;
    long long points_single_worker_count = points_count / workers_count;
    points_count = points_single_worker_count * workers_count;

    if (rank == 0) {
        double* points;
        points = (double *) malloc(points_count * 3 * sizeof(double) + MPI_BSEND_OVERHEAD + 7);
        generate_points(points, points_count, unif, re);
        double startTime = MPI_Wtime();

        for (long long process_num = 0; process_num < workers_count; ++process_num) {
            MPI_Send(&points[process_num * points_single_worker_count * 3], points_single_worker_count * 3, MPI_DOUBLE, process_num + 1, 0, MPI_COMM_WORLD);
            std::cout << "Msg for process " << process_num << " sent." << std::endl;
        }
        double sum = 0;
        double received_sum;
        double time = 0, cur_time;
        for (long long process_num = 0; process_num < workers_count; ++process_num) {
            MPI_Recv(&received_sum, 1, MPI_DOUBLE, process_num + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::cout << "Got sum: " << received_sum << std::endl;
            sum += received_sum;
            MPI_Recv(&cur_time, 1, MPI_DOUBLE, process_num + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (cur_time > time) {
                time = cur_time;
            }
        }
        std::cout << "TOTAL SUM::: " << sum << std::endl;
        double result = volume() * sum / (double)points_count;
        std::cout << "==============\nTOTAL RESULT::: " << result << std::endl << "EPS:" <<
            count_analytic_integral() - result << "\nTIME: " << time << "\n==============" << std::endl;
        free(points);
    } else {
        double* points;
        points = (double *) malloc(points_single_worker_count * 3 * sizeof(double) + MPI_BSEND_OVERHEAD + 7);
        MPI_Recv(points, points_single_worker_count * 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        double startTime = MPI_Wtime();
        double sum = 0;
        for (long long i = 0; i < points_single_worker_count * 3; i += 3) {
            sum += count_function_value(points[i], points[i + 1], points[i + 2]);
        }
        double total_time = MPI_Wtime() - startTime;
        MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&total_time, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        free(points);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}