#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

struct Point {
    int x;
    int y;
};

struct Trial {
    int nthRun;
    int trialNumber;
    double time;
};

struct Answer {
    double distance;
    int index1;
    int index2;
};


struct timespec bts, ets;

// const int MAX_RAND_VALUE = 32756;
const int MAX_RAND_VALUE = 50;
const int INT_MAX = 32756;

void PrintPoint(struct Point p) {
    printf("X: %d, Y: %d\n", p.x, p.y);
}

void PrintTrial(struct Trial t) {
    printf("nthRun: %d, trialNumber: %d, time: %f\n", t.nthRun, t.trialNumber, t.time);
}

double get_distance(struct Point a, struct Point b) {
    return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2));
}

bool containsPoint(struct Point P[], int currentIndex) {
    for (int i = 0; i < currentIndex; i++) {
        if (P[i].x == P[currentIndex].x && P[i].y == P[currentIndex].y) {
            return true;
        }
    }
    return false;
}

double diff_timespec(const struct timespec *time1, const struct timespec *time0) {
    return (time1->tv_sec - time0->tv_sec)
        + (time1->tv_nsec - time0->tv_nsec) / 1000000000.0;
}

void fillWithRandomPoints(struct Point P[], const int n) {
    for (int i = 0; i < n; i++) {
        do {
            // P[i].x = rand();
            // P[i].y = rand();
            P[i].x = rand() / (RAND_MAX / MAX_RAND_VALUE + 1);
            P[i].y = rand() / (RAND_MAX / MAX_RAND_VALUE + 1);
        } while (containsPoint(P, i));
    }
}

void int_swap(int* a, int* b) {
    const int temp = *a;
    *a = *b;
    *b = temp;
}

int x_partition(struct Point P[], int to_sort[], int low, int high) {
    const int pivot = P[to_sort[low]].x;
    int i = low;
    int j = high;

    while (i < j) {
        while (P[to_sort[i]].x <= pivot && i <= high - 1) {
            i++;
        }

        while (P[to_sort[j]].x > pivot && j >= low + 1) {
            j--;
        }

        if (i < j) {
            int_swap(&to_sort[i], &to_sort[j]);
        }
    }
    int_swap(&to_sort[low], &to_sort[j]);
    return j;
}

int y_partition(struct Point P[], int to_sort[], int low, int high) {
    const int pivot = P[to_sort[low]].y;
    int i = low;
    int j = high;

    while (i < j) {
        while (P[to_sort[i]].y <= pivot && i <= high - 1) {
            i++;
        }

        while (P[to_sort[j]].y > pivot && j >= low + 1) {
            j--;
        }

        if (i < j) {
            int_swap(&to_sort[i], &to_sort[j]);
        }
    }
    int_swap(&to_sort[low], &to_sort[j]);
    return j;
}

void quickSort(struct Point P[], int to_sort[], int low, int high, int (*f)(struct Point[], int[], int, int)) {
    if (low < high) {
        int partitionIndex = (*f)(P, to_sort, low, high);

        quickSort(P, to_sort, low, partitionIndex - 1, *f);
        quickSort(P, to_sort, partitionIndex + 1, high, *f);
    }
}

struct Answer DivideAndConquerClosestPoints_Rec(struct Point P[], int Pxi[], int Pyi[], int n) {
    struct Answer answer = { .distance = INT_MAX, .index1 = -1, .index2 = -2 };
    if (n <= 3) {
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                double d = get_distance(P[Pxi[i]], P[Pxi[j]]);
                if (d < answer.distance) {
                    answer.distance = d;
                    answer.index1 = Pxi[i];
                    answer.index2 = Pxi[j];
                }
            }
        }
        return answer;
    }
    
    int half = n / 2;
    if (half % 2 != 0) {
        half += 1;
    }
    int *Qx = malloc(sizeof(int)*(half));
    int *Qy = malloc(sizeof(int)*(half));
    int *Rx = malloc(sizeof(int)*(n-half));
    int *Ry = malloc(sizeof(int)*(n-half));
    
    for (int i = 0; i < half; i++) {
        Qx[i] = Pxi[i];
        Qy[i] = Pxi[i];
    }
    quickSort(P, Qy, 0, half, y_partition);
    
    for (int i = half; i < n; i++) {
        Rx[i] = Pxi[i];
        Ry[i] = Pxi[i];
    }
    quickSort(P, Ry, 0, n-half, y_partition);
    
    struct Answer q_answer = DivideAndConquerClosestPoints_Rec(P, Qx, Qy, half);
    struct Answer r_answer = DivideAndConquerClosestPoints_Rec(P, Rx, Ry, n-half);
    
    double delta = q_answer.distance;
    answer.distance = q_answer.distance;
    answer.index1 = q_answer.index1;
    answer.index2 = q_answer.index2;
    if (r_answer.distance < q_answer.distance) {
        delta = r_answer.distance;
        answer.distance = r_answer.distance;
        answer.index1 = r_answer.index1;
        answer.index2 = r_answer.index2;
    }
    
    int q_points = 1;
    for (int i = half-2; i >= 0; i--) {
        if ((P[Qx[half-1]].x - P[Qx[i]].x) < delta) {
            q_points += 1;
        }
        else {
            break;
        }
    }
    
    int r_points = 0;
    for (int i = 0; i < n-half; i++) {
        if ((P[Rx[i]].x - P[Qx[half-1]].x) < delta) {
            r_points += 1;
        }
        else {
            break;
        }
    }
    
    int *S = malloc(sizeof(int)*(q_points + r_points));
    for (int i = 0; i < q_points; i++) {
        S[i] = Qx[(half-1) - i];
    }
    int r_iteratror = 0;
    for (int i = q_points; i < q_points + r_points; i++) {
        S[i] = Rx[r_iteratror];
        r_iteratror += 1;
    }
    quickSort(P, S, 0, (q_points + r_points), y_partition);
    
    free(Qx);
    free(Qy);
    free(Rx);
    free(Ry);
    
    for (int i = 0; i < (q_points + r_points); i++) {
        for (int j = 1; j <= 15; j++) {
            double d = get_distance(P[S[i]], P[S[j]]);
            if (d < answer.distance) {
                answer.distance = d;
                answer.index1 = S[i];
                answer.index2 = S[j];
            }
        }
    }
    
    free(S);
    return answer;
}

struct Answer DivideAndConquerClosestPoints(struct Point P[], const int n) {
    int *Pxi = malloc(sizeof(int)*n);
    int *Pyi = malloc(sizeof(int)*n);
    for (int i = 0; i < n; i++) {
        Pxi[i] = i;
        Pyi[i] = i;
    }
    quickSort(P, Pxi, 0, n - 1, x_partition);
    quickSort(P, Pyi, 0, n - 1, y_partition);


    struct Answer answer = DivideAndConquerClosestPoints_Rec(P, Pxi, Pyi, n);
    free(Pxi);
    free(Pyi);
    return answer;
}

void CuteTest(struct Point P[], const int n) {
    int *Pxi = malloc(sizeof(int)*n);
    int *Pyi = malloc(sizeof(int)*n);
    for (int i = 0; i < n; i++) {
        Pxi[i] = i;
        Pyi[i] = i;
    }
    quickSort(P, Pxi, 0, n - 1, x_partition);
    quickSort(P, Pyi, 0, n - 1, y_partition);

    for (int i = 0; i < n; i++) {
        printf("X: %d\n", P[Pxi[i]].x);
    }
}

struct Answer BruteForceClosestPoints(struct Point P[], const int n) {
    struct Answer answer = { .distance = INT_MAX, .index1 = -1, .index2 = -1};
    for (int i = 0; i < n-1; i++) {
        for (int j = i + 1; j < n; j++) {
            double d = get_distance(P[i], P[j]);
            if (d < answer.distance) {
                answer.distance = d;
                answer.index1 = i;
                answer.index2 = j;

            }
        }
    }
    return answer;
}


int main(void) {
    srand(time(NULL));
    const int trials = 10;
    const int max_n = 1000;
    const int n_increment = 100;

    struct Trial *brute = malloc(sizeof(struct Trial) * trials * (max_n / n_increment));
    struct Trial *divide = malloc(sizeof(struct Trial) * trials * (max_n / n_increment));

    // main test loop
    int nth_run = 0;
    for (int n = n_increment; n <= max_n; n += n_increment) {
        for (int j = 0; j < trials; j++) {
            struct Point *P = malloc(sizeof(struct Point)*n);
            fillWithRandomPoints(P, n);
            clock_gettime(CLOCK_REALTIME, &bts);
            BruteForceClosestPoints(P, n);
            clock_gettime(CLOCK_REALTIME, &ets);

            brute[(nth_run * trials) + j].nthRun = nth_run;
            brute[(nth_run * trials) + j].trialNumber = j;
            brute[(nth_run * trials) + j].time = diff_timespec(&ets, &bts);


            clock_gettime(CLOCK_REALTIME, &bts);
            DivideAndConquerClosestPoints(P, n);
            clock_gettime(CLOCK_REALTIME, &ets);

            divide[(nth_run * trials) + j].nthRun = nth_run;
            divide[(nth_run * trials) + j].trialNumber = j;
            divide[(nth_run * trials) + j].time = diff_timespec(&ets, &bts);
            free(P);
        }
        nth_run += 1;
    }

    int brute_average[nth_run + 1];
    int divide_average[nth_run +1];

    // Calculate and Print Brute Averages
    printf("Brute Averages\n=====================\n");
    for (int i = 0; i < nth_run; i++) {
        double combined = 0;
        for (int j = 0; j < trials; j++) {
            combined += brute[(i * trials) + j].time;
        }
        printf("Number of items: %d, Average Time: %f\n", n_increment * (i + 1), combined / (double)trials);
    }

    printf("\n\n");

    // Calculate and Print Divide Averages
    printf("Divide Averages\n=====================\n");
    for (int i = 0; i < nth_run; i++) {
        double combined = 0;
        for (int j = 0; j < trials; j++) {
            combined += divide[(i * trials) + j].time;
        }
        printf("Number of items: %d, Average Time: %f\n", n_increment * (i + 1), combined / (double)trials);

    }
    return 0;

}
