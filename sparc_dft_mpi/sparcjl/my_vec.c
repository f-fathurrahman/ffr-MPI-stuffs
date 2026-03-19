// file: vec2.c
#include <stdlib.h>

typedef struct {
    int n;
    double *data;
} Vector;

void initialize(Vector *v, int n) {
    v->n = n;
    v->data = malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++)
        v->data[i] = i * 2.0;
}

double sum(Vector *v) {
    double s = 0.0;
    for (int i = 0; i < v->n; i++)
        s += v->data[i];
    return s;
}

void cleanup(Vector *v) {
    free(v->data);
}