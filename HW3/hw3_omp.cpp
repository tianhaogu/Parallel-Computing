#include <iostream>
#include <omp.h>
#include <vector>
#include <stack>
#include <queue>
#include "f.h"

using namespace std;

int main() {
    double fa = f(a), fb = f(b);
    double interval_len = b - a;
    double max_val = max(fa, fb);
    queue<vector<double> > q;
    vector<double> vec{a, b, fa, fb};
    q.push(vec);
    double start_time = 0.0;
    while (interval_len > 1e-2) {
        vector<double> curr_interval = q.front();
        q.pop();
        double start = curr_interval[0];
        double end = curr_interval[1];
        double fstart = curr_interval[2];
        double fend = curr_interval[3];
        if (interval_len == b - a) {
            start_time = omp_get_wtime();
        }
        double mid = (start + end) / 2;
        double fmid = f(mid);
        double max_val_tmp = max(max(fstart, fend), fmid);
        max_val = max(max_val_tmp, max_val);
        if ((fstart + fmid + s * (mid - start)) / 2 >= max_val + epsilon) {
            vector<double> vec{start, mid, fstart, fmid};
            q.push(vec);
        }
        if ((fmid + fend + s * (end - mid)) / 2 >= max_val + epsilon) {
            vector<double> vec{mid, end, fmid, fend};
            q.push(vec);
        }
        interval_len /= 2;
    }

    #pragma omp parallel
    {
        #pragma omp single nowait
            cout << "Number of threads we use here is: " << omp_get_num_threads() << endl;
        while (!q.empty()) {
            stack<vector<double > > stk;
            vector<double> curr_interval;
            #pragma omp critical(get_value_from_queue)
            {
                if (!q.empty()) {
                    curr_interval = q.front();
                    q.pop();
                }
            }
            if (curr_interval.size() != 0) {
                stk.push(curr_interval);
            }
            while (!stk.empty()) {
                vector<double> original = stk.top();
                stk.pop();
                double start = original[0];
                double end = original[1];
                double fstart = original[2];
                double fend = original[3];
                double mid = (start + end) / 2;
                double fmid = f(mid);
                double max_val_tmp = max(max(fstart, fend), fmid);
                #pragma omp critical(update_value)
                {
                    max_val = max(max_val_tmp, max_val);
                }
                if ((fmid + fend + s * (end - mid)) / 2 >= max_val + epsilon) {
                    vector<double> vec{mid, end, fmid, fend};
                    stk.push(vec);
                }
                if ((fstart + fmid + s * (mid - start)) / 2 >= max_val + epsilon) {
                    vector<double> vec{start, mid, fstart, fmid};
                    stk.push(vec);
                }
            }
        }
    }
    double end_time = omp_get_wtime();
    cout << "Maximum of the function f between a and b is: " << max_val << endl;
    double duration = end_time - start_time;
    cout << "Total time getting the maximum is: " << duration << " seconds" << endl;
    return 0;
}
