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

void PrintPoint(struct Point p) {
    printf("X: %d, Y: %d\n", p.x, p.y);
}

void PrintTrial(struct Trial t) {
    printf("nthRun: %d, trialNumber: %d, time: %f\n", t.nthRun, t.trialNumber, t.time);
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
    struct Answer test = { .distance = 0, .index1 = 0, .index2 = 0 };
    return test;
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
            double d = sqrt(pow((P[i].x - P[j].x), 2) + pow((P[i].y - P[j].y), 2));
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
            CuteTest(P, n);
            return 0;
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


