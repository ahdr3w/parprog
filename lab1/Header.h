#pragma once
#include<iostream>
#include<math.h>
#include<mpi.h>
#include<time.h>

class problem{
public:
	int K, M;
	float X = 10, T = 10;
	float a = 3;
	float tau;
	float h;
	float alpha;
	float beta;
	problem(int K, int M) {
		this->K = K;
		this->M = M;
		this->tau = T / K;
		this->h = X / M;
		this->alpha = (1 / (2 * tau) - a / (2 * h)) / (1 / (2 * tau) - a / (2 * h));
		this->beta = 1 / (1 / (2 * tau) - a / (2 * h));
	}
	
	float func(float t, float x) {
		return (5*x+1)*cos(7*t);
	}
	float phi_func(float x) {
		return log(x+1);
	}
	float psi_func(float t) {
		return sin(14*t);
	}
};
float** init(problem* p);
float** par_alg(problem* p);
float** lin_alg(problem* p);
float   compare(problem* p, float** u1, float** u2);
void    print(problem* p, float** u);
void del(float** u,problem* p);