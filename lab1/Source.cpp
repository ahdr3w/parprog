#include"header.h"
using namespace std;

float** par_alg(problem* p, int commsize, int my_rank) {
	MPI_Status * status  = new MPI_Status;
	int str_num = 0;
	int str_point = 0;
	float** buf;
	float** answer = NULL;
	double t1, t2 = 0;
	t1 = MPI_Wtime();
	if (my_rank == 0) {
		int* str = new int[commsize];
		int* init_str = new int[commsize];
		for (int i = 0; i < commsize; i++) {
			str[i] = 0;
			init_str[i] = 0;
		}
		int idx = p->K;
		int count = commsize;
		while (idx > 0) {
			for (int i = 0; i < count; i++)
				str[i] += idx / count;
			idx %= count;
			count--;
		}
		for (int i = 1; i < commsize; i++)
			MPI_Send(&str[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		str_num = str[0];
		init_str[0] = 1;
		for (int i = 1; i < commsize; i++)
			init_str[i] = init_str[i - 1] + str[i - 1];
		str_point = init_str[0];
		for (int i = 1; i < commsize; i++) 
			MPI_Send(&init_str[i], 1, MPI_INT, i, 1, MPI_COMM_WORLD);
		delete[] str;
		delete[] init_str;
	}
	else {
		MPI_Recv(&str_num, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, status);
		MPI_Recv(&str_point, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, status);
	}

	buf = new float* [str_num + 1];
	for (int i = 0; i <= str_num; i++) {
		buf[i] = new float[p->M + 1];
		buf[i][0] = p->psi_func((str_point + i - 1) * p->tau);
	}
	if (my_rank == 0)
		for (int m = 0; m <= p->M; m++) 
			buf[0][m] = p->phi_func(m * p->h);
	int m = 1;
	for (int m = 1; m <= p->M; m++) {
		if (my_rank > 0)
			MPI_Recv(&buf[0][m], 1, MPI_FLOAT, my_rank - 1, 2, MPI_COMM_WORLD, status);
		for (int k = 1; k <= str_num; k++) 
			buf[k][m] = float(buf[k - 1][m - 1] + p->alpha * buf[k - 1][m] - p->alpha * buf[k][m - 1] + p->beta * p->func((str_point + k - 1 - (float)0.5) * p->tau, (m + (float)0.5) * p->h));
		if (my_rank < commsize - 1)
			MPI_Send(&buf[str_num][m], 1, MPI_FLOAT, my_rank + 1, 2, MPI_COMM_WORLD);
	}
	t2 = MPI_Wtime();
	if (my_rank == 0)
		cout << t2 - t1 << ", ";
	if (my_rank != 0) {
		MPI_Send(&str_num, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
		delete[] buf[0];
		for (int i = 1; i <= str_num; i++) {
			MPI_Send(&buf[i][0], p->M + 1, MPI_FLOAT, 0, 4, MPI_COMM_WORLD);
			delete[] buf[i];
		}
		delete[] buf;
	}
	else {
		answer = new float* [p->K + 1];
		for (int k = 0; k <= str_num; k++) {
			answer[k] = new float[p->M + 1];
			for (int m = 0; m <= p->M; m++)
				answer[k][m] = buf[k][m];
			delete[] buf[k];
		}
		delete[] buf;
		int idx = 0;
		int point = str_point + str_num;
		for (int i = 1; i < commsize; i++) {
			MPI_Recv(&idx, 1, MPI_INT, i, 3, MPI_COMM_WORLD, status);
			for (int k = 0; k < idx; k++) {
				answer[point + k] = new float[p->M + 1];
				MPI_Recv(&answer[point + k][0], p->M + 1, MPI_FLOAT, i, 4, MPI_COMM_WORLD, status);
			}
			point += idx;
		}
	}
	delete status;
	if (my_rank == 0)
		return answer;
	else
		return NULL;
}

float** lin_alg(problem* p) {
	float** u = init(p);
	for (int m = 1; m <= p->M; m++) {
		for (int k = 1; k <= p->K; k++)
			u[k][m] = u[k - 1][m - 1] + p->alpha * u[k - 1][m] - p->alpha * u[k][m - 1] + p->beta * p->func((k - (float)0.5) * p->tau, (m + (float)0.5) * p->h);
	}
	return u;
}

float** init(problem* p) {
	float** u = new float* [p->K + 1];
	for (int i = 0; i <= p->K; i++) {
		u[i] = new float[p->M+1];
	}
	for (int i = 0; i <= p->M; i++) {
		u[0][i] = p->phi_func(p->h * i);
	}
	for (int i = 0; i <= p->K; i++) {
		u[i][0] = p->psi_func(p->tau * i);
	}
	return u;
}

float compare(problem* p, float** u1, float** u2) {
	float max = 0;
	int idxi = -1;
	int idxj = -1;
	for (int k = 0; k <= p->K; k++) {
		for (int m = 0; m <= p->M; m++)
			if (abs(u1[k][m] - u2[k][m]) > max) {
				max = abs(u1[k][m] - u2[k][m]);
				idxi = k;
				idxj = m;
			}
	}
	cout << "[" << idxi << ", " << idxj << "]" << endl;
	return max;
}

void del(float** u, problem* p) {
	for (int i = 0; i <= p->K; i++)
		delete[] u[i];
	delete[] u;
}

void print(problem* p, float** u) {
	for (int k = 0; k <= p->K; k++)
		cout << k * p->tau << ", ";
	cout << endl;
	cout << endl;
	for (int m = 0; m <= p->M; m++)
		cout << m * p->h << ", ";
	cout << endl;
	cout << endl;

	for (int k = 0; k <= p->K; k++) {
		cout << "[";
		for (int m = 0; m <= p->M; m++)
			cout << u[k][m] << ", ";
		cout << "], ";
		cout << endl;
	}
}

int main(int argc, char* argv[]) {
	int commsize;
	int my_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);
	for (int k = 1; k <= 10; k++) {
		problem* p = new problem(k * 1000, 30000);
		float** u1 = par_alg(p, commsize, my_rank);
		if (my_rank == 0) {
		//	float** u2 = lin_alg(p);
		//	cout << "maxdif: " << compare(p, u1, u2) << " "; //ëèíåéíûé àëãîðèòì
		//	del(u2, p);
		}
		if (my_rank == 0) del(u1, p);
		delete p;
	}
	if (my_rank == 0) {
		cout << endl;
		for (int k = 1; k <= 10; k++)
			cout << k * 1000 << ", ";
	}
	MPI_Barrier(MPI_COMM_WORLD);
	problem* pr = new problem(1000, 3000);
	float** u3 = par_alg(pr, commsize, my_rank);
	if (my_rank == 0) {
		ofstream fout("datagraph.txt");
		for (int k = 0; k <= pr->K; k++) {
			fout << "[";
			for (int m = 0; m <= pr->M; m++)
				fout << u3[k][m] << ", ";
			fout << "]," << endl;
		}
		fout.close();
		del(u3, pr);
		delete pr;
	}
	MPI_Finalize();
	
}   
