#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#include <iostream>
#include <string>

#include <fstream>
#include <sstream>
#include <vector>
#include <array>
//#include <bits/stdc++.h>


#include <cuda.h>
using namespace std;

enum direction { d_down, d_right, none };

#define COORD std::pair<int, int>

//#define DEBUG

int iter = 0;

//////////////////////////////////////////////
//Auxiliary functions for preparing problem //
//////////////////////////////////////////////

void display_arr(int* arr, int n) {

	cout << "arr: ";

	for (int i = 0; i < n; i++) {
		cout << arr[i] << " ";
	}

	cout << endl;

}

void print_coords(COORD start, COORD end) {

	cout << "Start:" << start.first << "," << start.second << endl;
	cout << "End:" << end.first << "," << end.second << endl;

}

int find_length(COORD start, COORD end, direction dir) {

	if (dir == d_down)
		return end.first - start.first;
	if (dir == d_right)
		return end.second - start.second;

	return -1;
}

void convert_sol(int** mat, int**& sol_mat, int m, int n) {

	sol_mat = new int* [m]; //Rows
	for (int i = 0; i < m; i++) {
		sol_mat[i] = new int[n]; //Cols
	}

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			if (mat[i][j] == -2)
				sol_mat[i][j] = -2; //Empty value cell
			else
				sol_mat[i][j] = -1; //Hint or empty cell
		}
	}
}

void print_one_matrix(int** matrix, int m, int n) {
	std::cout << "Matrix: " << std::endl;
	for (int i = 0; i < m; i++) { //rows
		for (int j = 0; j < n; j++) { //cols
			std::cout << matrix[i][j] << "\t";
		}
		std::cout << "\n";
	}
}

///Auxiliary functions

struct sum {
	COORD start;
	COORD end;

	int hint;
	int dir;
	int length;
	int* arr;

	void print_sum() {
		cout << "############################" << endl;
		cout << "Creating sum with: " << endl;
		print_coords(start, end);
		cout << "Hint: " << hint << endl;
		cout << "Direction: " << dir << endl;
		cout << "Length: " << length << endl;
		cout << "############################" << endl;
	}

	sum(COORD _start, COORD _end, int _hint, direction _dir) :
		start(_start), end(_end), hint(_hint), dir(_dir)
	{
		length = find_length(_start, _end, _dir);
		arr = new int[length];
#ifdef DEBUG
		cout << "############################" << endl;
		cout << "Creating sum with: " << endl;
		print_coords(start, end);
		cout << "Hint: " << hint << endl;
		cout << "Direction: " << dir << endl;
		cout << "Length: " << length << endl;
		cout << "############################" << endl;
#endif
	}

	//~sum(){
	//delete arr;
	//}
};


COORD find_end(int** matrix, int m, int n, int i, int j, direction dir) { //0 down 1 right

	if (dir == d_right) {
		for (int jj = j + 1; jj < n; jj++) {
			if (matrix[i][jj] != -2 || jj == n - 1) {
				if (matrix[i][jj] == -2 && jj == n - 1)
					jj++;
				COORD END = COORD(i, jj);
				return END;
			}
		}
	}

	if (dir == d_down) {
		for (int ii = i + 1; ii < m; ii++) {
			if (matrix[ii][j] != -2 || ii == m - 1) {
				if (matrix[ii][j] == -2 && ii == m - 1)
					ii++;
				COORD END = COORD(ii, j);
				return END;
			}
		}
	}

}


vector<sum> get_sums(int** matrix, int m, int n) {

	vector<sum> sums;

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			int val = matrix[i][j];
			if (val != -1 && val != -2) {
				int hint = val;
				hint = hint / 10;

				if ((hint % 100) == 0) {
					hint = (int)(hint / 100);
					COORD START = COORD(i, j + 1);
					COORD END = find_end(matrix, m, n, i, j, d_right);
					sum _sum = sum(START, END, hint, d_right);
					sums.push_back(_sum);
				}

				else {
					int div = (int)(hint / 100);
					int rem = (int)(hint % 100);

					if (div == 0 && rem != 0) {
						COORD START = COORD(i + 1, j);
						COORD END = find_end(matrix, m, n, i, j, d_down);
						sum _sum = sum(START, END, rem, d_down);
						sums.push_back(_sum);
					}

					if (div != 0 && rem != 0) {
						COORD START1 = COORD(i + 1, j);
						COORD START2 = COORD(i, j + 1);
						COORD END1 = find_end(matrix, m, n, i, j, d_down);
						COORD END2 = find_end(matrix, m, n, i, j, d_right);
						sum _sum1 = sum(START1, END1, rem, d_down);
						sum _sum2 = sum(START2, END2, div, d_right);
						sums.push_back(_sum1);
						sums.push_back(_sum2);
					}
				}
			}


		}
	}
	return sums;
}


void read_matrix(int**& matrix, std::ifstream& afile, int m, int n) {

	matrix = new int* [m]; //rows

	for (int i = 0; i < m; i++) {
		matrix[i] = new int[n]; //cols
	}

	int val;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			afile >> val;
			matrix[i][j] = val;
		}
	}
}

void sol_to_file(int** mat, int** sol_mat, int m, int n) {

	string fname = "visualize.kakuro";
	ofstream to_write(fname);

	to_write << m << " " << n << "\n";

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (mat[i][j] != -2)
				to_write << mat[i][j] << " ";
			else
				to_write << sol_mat[i][j] << " ";
		}
		to_write << "\n";
	}

	to_write.close();
}

//////////////////////////////////////////////
//Auxiliary functions for preparing problem //
//////////////////////////////////////////////

///////////////////////////////////////////////////
//Auxiliary functions for preparing CUDA setting //
///////////////////////////////////////////////////

void flatten_sums(vector<sum> sums, int* h_sum_starts_x, int* h_sum_starts_y, int* h_sum_ends_x, int* h_sum_ends_y, int* h_sum_hints, int* h_sum_lengths, int* h_sum_dirs, int no_sums) {

	for (int i = 0; i < no_sums; i++) {

		h_sum_starts_x[i] = sums[i].start.first;
		h_sum_starts_y[i] = sums[i].start.second;

		h_sum_ends_x[i] = sums[i].end.first;
		h_sum_ends_y[i] = sums[i].end.second;

		h_sum_hints[i] = sums[i].hint;
		h_sum_lengths[i] = sums[i].length;

		h_sum_dirs[i] = sums[i].dir;
	}

}

void print_flattened(int* h_sum_starts_x, int* h_sum_starts_y, int* h_sum_ends_x, int* h_sum_ends_y, int* h_sum_hints, int* h_sum_lengths, int* h_sum_dirs, int no_sums) {

	cout << "###h_sum_starts_x: " << endl;
	for (int i = 0; i < no_sums; i++) {
		cout << h_sum_starts_x[i] << " ";
	}
	cout << endl;

	cout << "###h_sum_starts_y: " << endl;
	for (int i = 0; i < no_sums; i++) {
		cout << h_sum_starts_y[i] << " ";
	}
	cout << endl;

	cout << "###h_sum_ends_x: " << endl;
	for (int i = 0; i < no_sums; i++) {
		cout << h_sum_ends_x[i] << " ";
	}
	cout << endl;

	cout << "###h_sum_ends_y: " << endl;
	for (int i = 0; i < no_sums; i++) {
		cout << h_sum_ends_y[i] << " ";
	}
	cout << endl;

	cout << "###h_sum_hints: " << endl;
	for (int i = 0; i < no_sums; i++) {
		cout << h_sum_hints[i] << " ";
	}
	cout << endl;

	cout << "###h_sum_lengths: " << endl;
	for (int i = 0; i < no_sums; i++) {
		cout << h_sum_lengths[i] << " ";
	}
	cout << endl;

	cout << "###h_sum_dirs: " << endl;
	for (int i = 0; i < no_sums; i++) {
		cout << h_sum_dirs[i] << " ";
	}
	cout << endl;

}

void flatten_sol_mat(int** sol_mat, int* h_sol_mat, int m, int n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			h_sol_mat[i * n + j] = sol_mat[i][j];
		}
	}
}

void print_flattened_matrix(int* h_sol_mat, int m, int n) {

	cout << "###Flattened matrix: " << endl;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cout << h_sol_mat[i * n + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

///////////////////////////////////////////////////
//Auxiliary functions for preparing CUDA setting //
///////////////////////////////////////////////////


///////////////////
//CUDA FUNCTIONS //
///////////////////
__device__ void d_print_flattened_matrix(int* d_sol_mat, int m, int n) {
	printf("_device_ matrix: %d %d \n", m, n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			printf("%d ", d_sol_mat[i * n + j]);
		}
		printf("\n");
	}
	printf("\n");
}

__device__ int get_cord_to_change(int& current_ind, int last_index, int* mat_arr) {
	while (current_ind < last_index) {
		if (mat_arr[current_ind] == -2) {
			return current_ind;
		}
		current_ind++;
	}
	return -10; // why ? i dont know 
}

__device__ void copy_mat(int* from, int* to, int last_index) {
	for (int i = 0; i < last_index; i++) {
		to[i] = from[i];
	}
}

struct state_data {
	int x_cord;
	int y_cord;
	int val;
	__device__ state_data() : x_cord(-10), y_cord(-10), val(-10) {}
	__device__ state_data(int x, int y, int v) : x_cord(x), y_cord(y), val(v) {}
};

struct my_gpu_stack {
	int top;
	int maxSize;
	state_data* data;

	__device__ my_gpu_stack(int max) : top(-1), maxSize(max*10) {
		data = new state_data[maxSize];
	}

	__device__ ~my_gpu_stack() {
		delete[] data;
	}

	__device__ void push(state_data value) {
		if (top < maxSize - 1) {
			data[++top] = value;
		}
		else {
			printf("Stack Overflow\n");
		}
	}

	__device__ state_data pop() {
		if (top >= 0) {
			return data[top--];
		}
		else {
			printf("Stack Underflow\n");
			return state_data();
		}
	}

	__device__ bool is_empty() {
		return top == -1;
	}
};

__device__ bool check_a_spesific_sum_if_contains_duplicates(int* d_sol_mat, int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x, int* d_sum_ends_y, int* d_sum_hints, int* d_sum_lengths, int* d_sum_dirs, int no_sums, int spesific_sum_index, int m, int n) {
	int start_x = d_sum_starts_x[spesific_sum_index];
	int start_y = d_sum_starts_y[spesific_sum_index];
	int end_x = d_sum_ends_x[spesific_sum_index];
	int end_y = d_sum_ends_y[spesific_sum_index];
	int hint = d_sum_hints[spesific_sum_index];
	int length = d_sum_lengths[spesific_sum_index];
	int dir = d_sum_dirs[spesific_sum_index];
	bool has_duplicates = false;
	if (dir == 0) {
		for (int j = start_x; j <= end_x; j++) {
			for (int k = j + 1; k <= end_x; k++) {
				if (d_sol_mat[j * n + start_y] == d_sol_mat[k * n + start_y]) {
					has_duplicates = true;
					return has_duplicates;
				}
			}
		}
	}
	else {
		for (int j = start_y; j <= end_y; j++) {
			for (int k = j + 1; k <= end_y; k++) {
				if (d_sol_mat[start_x * n + j] == d_sol_mat[start_x * n + k]) {
					has_duplicates = true;
					return has_duplicates;
				}
			}
		}
	}
	return has_duplicates;
}

__device__ bool check_a_spesific_sum_if_contains_2(int* d_sol_mat, int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x, int* d_sum_ends_y, int* d_sum_hints, int* d_sum_lengths, int* d_sum_dirs, int no_sums, int spesific_sum_index, int m, int n) {
	int start_x = d_sum_starts_x[spesific_sum_index];
	int start_y = d_sum_starts_y[spesific_sum_index];
	int end_x = d_sum_ends_x[spesific_sum_index];
	int end_y = d_sum_ends_y[spesific_sum_index];
	int hint = d_sum_hints[spesific_sum_index];
	int length = d_sum_lengths[spesific_sum_index];
	int dir = d_sum_dirs[spesific_sum_index];
	bool does_contain_2 = false;
	if (dir == 0) {
		for (int j = start_x; j <= end_x; j++) {
			if (d_sol_mat[j * n + start_y] == -2) {
				does_contain_2 = true;
				return does_contain_2;
			}
		}
	}
	else {
		for (int j = start_y; j <= end_y; j++) {
			if (d_sol_mat[start_x * n + j] == -2) {
				does_contain_2 = true;
				return does_contain_2;
			}
		}
	}
	return does_contain_2;
}


__device__ int sum_a_spesific_sum(int* d_sol_mat, int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x, int* d_sum_ends_y, int* d_sum_hints, int* d_sum_lengths, int* d_sum_dirs, int no_sums, int spesific_sum_index, int m, int n) {
	int start_x = d_sum_starts_x[spesific_sum_index];
	int start_y = d_sum_starts_y[spesific_sum_index];
	int end_x = d_sum_ends_x[spesific_sum_index];
	int end_y = d_sum_ends_y[spesific_sum_index];
	int hint = d_sum_hints[spesific_sum_index];
	int length = d_sum_lengths[spesific_sum_index];
	int dir = d_sum_dirs[spesific_sum_index];
	int sum_of_sum = 0;
	if (dir == 0) {
		for (int j = start_x; j <= end_x; j++) {
			sum_of_sum += d_sol_mat[j * n + start_y];
		}
	}
	else {
		for (int j = start_y; j <= end_y; j++) {
			sum_of_sum += d_sol_mat[start_x * n + j];
		}
	}
	return sum_of_sum;
}

__device__ bool check_singe_sum_is_correct(int* d_sol_mat, int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x, int* d_sum_ends_y, int* d_sum_hints, int* d_sum_lengths, int* d_sum_dirs, int no_sums, int spesific_sum_index, int m, int n) {
	int sum_of_sum = sum_a_spesific_sum(d_sol_mat, d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints, d_sum_lengths, d_sum_dirs, no_sums, spesific_sum_index, m, n);
	int hint = d_sum_hints[spesific_sum_index];
	if (sum_of_sum != hint) {
		return false;
	}
	else {
		if (check_a_spesific_sum_if_contains_2(d_sol_mat, d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints, d_sum_lengths, d_sum_dirs, no_sums, spesific_sum_index, m, n)) {
			return false;
		}
		else {
			if (check_a_spesific_sum_if_contains_duplicates(d_sol_mat, d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints, d_sum_lengths, d_sum_dirs, no_sums, spesific_sum_index, m, n)) {
				return false;
			}
			else {
				return true;
			}
		}
	}
}

__device__ bool check_singe_sum_is_Viable(int* d_sol_mat, int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x, int* d_sum_ends_y, int* d_sum_hints, int* d_sum_lengths, int* d_sum_dirs, int no_sums, int spesific_sum_index, int m, int n) {
	int sum_of_sum = sum_a_spesific_sum(d_sol_mat, d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints, d_sum_lengths, d_sum_dirs, no_sums, spesific_sum_index, m, n);
	int hint = d_sum_hints[spesific_sum_index];
	if (sum_of_sum <= hint) {
		return true;
	}
	else {
		return false;
	}
}

__device__ bool check_all_the_sums_are_correct(int* d_sol_mat, int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x, int* d_sum_ends_y, int* d_sum_hints, int* d_sum_lengths, int* d_sum_dirs, int no_sums, int m, int n) {
	for (int i = 0; i < no_sums; i++) {
		if (!check_singe_sum_is_correct(d_sol_mat, d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints, d_sum_lengths, d_sum_dirs, no_sums, i, m, n)) {
			return false;
		}
	}
	return true;
}

// enum direction { d_down, d_right, none };
// for a given kakuro board make a lookup table which takes cordinate x,y and returns the two 
__device__ int* make_lookup_table(int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x, int* d_sum_ends_y, int* d_sum_hints, int* d_sum_lengths, int* d_sum_dirs, int no_sums, int m, int n) {
	int* lookup_table = new int[m * n * 2];
	for (int i = 0; i < m * n * 2; i++) {
		lookup_table[i] = -1;
	}
	for (int i = 0; i < no_sums; i++) {
		int start_x = d_sum_starts_x[i];
		int start_y = d_sum_starts_y[i];
		int end_x = d_sum_ends_x[i];
		int end_y = d_sum_ends_y[i];
		int dir = d_sum_dirs[i];
		if (dir == 1) { // if direction is down
			for (int j = start_y; j < end_y ; j++) {
				int index = (start_x * n + j) ;
				if (lookup_table[index] == -1) {
					lookup_table[index] = i;
				}
				else {
					lookup_table[index + m * n] = i; // tune this 
				}
			}
		}
		else { // if direction is right
			for (int j = start_x; j < end_x ; j++) {
				int index = (j * n + start_y);
				if (lookup_table[index] == -1) {
					lookup_table[index] = i;
				}
				else {
					lookup_table[index + m*n ] = i;
				}
			}
		}

	}
	return lookup_table;
}


__global__ void kakuro_kernel(
	int* d_sum_starts_x, int* d_sum_starts_y, int* d_sum_ends_x, int* d_sum_ends_y, int* d_sum_hints, int* d_sum_lengths, int* d_sum_dirs,
	int* d_sol_mat, int* d_perms,
	int* d_t_mats, int m, int n,
	int no_sums, volatile bool* solved) {
	//TO DO
	const int MAT_SIZE = m * n * sizeof(int);
	int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
	int* this_threads_mat = d_t_mats + thread_id * MAT_SIZE;
	my_gpu_stack stack(no_sums);
	copy_mat(d_sol_mat, this_threads_mat, MAT_SIZE / sizeof(int));
	__syncthreads();

	for (int i = 0; i < no_sums; i++) {
		printf(" hints : %d \n", d_sum_hints[i]);
	}

	d_print_flattened_matrix(this_threads_mat, m, n);
	int cur_ind = 0;


	// Initial state

	int tmp = get_cord_to_change(cur_ind, MAT_SIZE / sizeof(int), this_threads_mat);

	int i = 1;
	stack.push(state_data(tmp / m, tmp % n, i));
	
	printf("Hints :");
	for (int i = 0; i < no_sums; i++) {
		printf(" %d ", d_sum_hints[i]);
	}

	int* cord_to_sum_lookup = make_lookup_table(d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints, d_sum_lengths, d_sum_dirs, no_sums, m, n);

	printf("\n TABLe \n");	
	d_print_flattened_matrix(cord_to_sum_lookup, m, n);
	printf("\n TABLe2 \n");
	d_print_flattened_matrix(cord_to_sum_lookup + m * n, m, n);
	printf("\n TABLe \n");
	
	while (!stack.is_empty()) {
		d_print_flattened_matrix(this_threads_mat, m, n);


 		int tttmp = 0;
		int is_there_empthy_slot = get_cord_to_change(tttmp, MAT_SIZE / sizeof(int), this_threads_mat);
		
		// is at leaf node ?? 
		if (is_there_empthy_slot == -10) {						
			// TODO do stuff to check if solved and if not give the correct value to the stack 
			if (check_all_the_sums_are_correct(this_threads_mat, d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints, d_sum_lengths, d_sum_dirs, no_sums, m, n)) {
				*solved = true;
				return;
			}
			else {
				// TODO add viablitiy check here 
				state_data cur = stack.pop();
				i = cur.val;
				if (i < 10) {
					this_threads_mat[cur.x_cord * m + cur.y_cord] = i;
					i++; // next iterations value
					stack.push(state_data(cur.x_cord, cur.y_cord, i));
				}
				else {
					this_threads_mat[cur.x_cord * m + cur.y_cord] = -2;					
				}				
			} 
		}		
		else {
			state_data cur = stack.pop();
			i = cur.val;
			if (cur.val < 10) {
				this_threads_mat[cur.x_cord * m + cur.y_cord] = i;
				i++; // next iterations value

				// Pruning condition: Check if the related sums are still valid
				bool Viable = false;
				// get inpacted sums and check if they are valid
				int first_effected_sum_index = cord_to_sum_lookup[cur.x_cord * m + cur.y_cord];
				int second_effected_sum_index = cord_to_sum_lookup[cur.x_cord * m + cur.y_cord + m * n];
				printf(" %d  %d \n ", first_effected_sum_index, second_effected_sum_index);
				if (check_singe_sum_is_Viable(this_threads_mat, d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints, d_sum_lengths, d_sum_dirs, no_sums, first_effected_sum_index, m, n) 
					&& check_singe_sum_is_Viable(this_threads_mat, d_sum_starts_x,  d_sum_starts_y, d_sum_ends_x,  d_sum_ends_y,  d_sum_hints, d_sum_lengths,d_sum_dirs, no_sums, second_effected_sum_index,  m,  n) )
				{
					Viable = true;
				}

				if (Viable) {
					// Save the current state and push it back to the stack
					stack.push(state_data(cur.x_cord, cur.y_cord, i));

					int tttmp = 0;
					// Move on to the next state
					tmp = get_cord_to_change(tttmp, MAT_SIZE / sizeof(int), this_threads_mat);
					//next_cord = get<0>(tmp);
					//sums_index = get<1>(tmp);
					i = 1;
					stack.push(state_data(tmp / m, tmp % n, i));;
				}
				else {
					// Reset the current cell when backtracking
					this_threads_mat[cur.x_cord * m + cur.y_cord] = -2;
					//stack.push(state_data(next_cord, sums_index, i));
					stack.push(state_data(cur.x_cord, cur.y_cord, i));;

				}
			}
			else {
				// Reset the cell to its initial value when backtracking
				this_threads_mat[cur.x_cord * m + cur.y_cord] = -2;
			}
		}
		////////////















		//About volatile bool* solved:
		//You can get idea from https://stackoverflow.com/questions/12505750/how-can-a-global-function-return-a-value-or-break-out-like-c-c-does%5B/url%5D for how to break out of a CUDA kernel
		//You may or may not use it

	}
	

}
///////////////////
//CUDA FUNCTIONS //
///////////////////

int main(int argc, char** argv) {

	std::string filename(argv[1]);
	std::ifstream file;
	file.open(filename.c_str());

	int m, n;
	file >> m;
	file >> n;

	int** mat;
	read_matrix(mat, file, m, n);
	print_one_matrix(mat, m, n);

	int** sol_mat;
	convert_sol(mat, sol_mat, m, n);
	//print_one_matrix(sol_mat, m, n);

	vector<sum> sums = get_sums(mat, m, n);

	//CUDA
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	printf("==prop== Running on device: %d -- %s \n", 0, prop.name);
	printf("==prop== #of SM -- %d \n", prop.multiProcessorCount);
	printf("==prop== Max Threads Per Block: -- %d \n", prop.maxThreadsPerBlock);
	printf("==prop== Shared Memory Per Block: -- %zu bytes\n", prop.sharedMemPerBlock);
	printf("==prop== Total global memory: -- %zu bytes\n", prop.totalGlobalMem);
	printf("==prop== Clock rate: -- %d\n", prop.clockRate);
	printf("==prop== Compute capability: -- %d.%d\n", prop.major, prop.minor);



	int grid_dim = 1;//TO DO
	int block_dim = 1;//To DO

	int no_sums = sums.size();


	//Flattening sums and matrix
	int* h_sum_starts_x = new int[no_sums];
	int* h_sum_starts_y = new int[no_sums];
	int* h_sum_ends_x = new int[no_sums];
	int* h_sum_ends_y = new int[no_sums];
	int* h_sum_hints = new int[no_sums];
	int* h_sum_lengths = new int[no_sums];
	int* h_sum_dirs = new int[no_sums];


	int* h_perms = new int; // not sure what to do wiht this one here 

	flatten_sums(sums, h_sum_starts_x, h_sum_starts_y, h_sum_ends_x, h_sum_ends_y, h_sum_hints, h_sum_lengths, h_sum_dirs, no_sums);

	print_flattened(h_sum_starts_x, h_sum_starts_y, h_sum_ends_x, h_sum_ends_y, h_sum_hints, h_sum_lengths, h_sum_dirs, no_sums);

	int* h_sol_mat;
	h_sol_mat = new int[m * n];
	flatten_sol_mat(sol_mat, h_sol_mat, m, n);

	print_flattened_matrix(h_sol_mat, m, n);

	//Declare device pointers and copy data into device
	int* d_sum_starts_x, * d_sum_starts_y, * d_sum_ends_x, * d_sum_ends_y, * d_sum_hints, * d_sum_lengths, * d_sum_dirs, * d_sol_mat, * d_t_mats;

	int* d_perms;// not sure what to do wiht this one here 

	cudaMalloc(&d_sum_starts_x, no_sums * sizeof(int));
	cudaMalloc(&d_sum_starts_y, no_sums * sizeof(int));
	cudaMalloc(&d_sum_ends_x, no_sums * sizeof(int));
	cudaMalloc(&d_sum_ends_y, no_sums * sizeof(int));
	cudaMalloc(&d_sum_hints, no_sums * sizeof(int));
	cudaMalloc(&d_sum_lengths, no_sums * sizeof(int));
	cudaMalloc(&d_sum_dirs, no_sums * sizeof(int));
	cudaMalloc(&d_sol_mat, (m * n) * sizeof(int));
	cudaMalloc(&d_t_mats, (m * n * grid_dim * block_dim) * sizeof(int)); //Allocating invidual matrix for each GPU thread
	//You may use this array if you will implement a thread-wise solution

	cudaMalloc(&d_perms, sizeof(int));// not sure what to do wiht this one here 

	cudaMemcpy(d_sum_starts_x, h_sum_starts_x, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum_starts_y, h_sum_starts_y, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum_ends_x, h_sum_ends_x, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum_ends_y, h_sum_ends_y, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum_hints, h_sum_hints, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum_lengths, h_sum_lengths, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sum_dirs, h_sum_dirs, no_sums * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sol_mat, h_sol_mat, (m * n) * sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(d_perms, h_perms, sizeof(int), cudaMemcpyHostToDevice);// not sure what to do wiht this one here 


	bool* solved = new bool[1];
	*solved = false;
	bool* d_solved;

	cudaMalloc(&d_solved, sizeof(bool));
	cudaMemcpy(d_solved, solved, sizeof(bool), cudaMemcpyHostToDevice);


	kakuro_kernel << <grid_dim, block_dim >> > (d_sum_starts_x, d_sum_starts_y, d_sum_ends_x, d_sum_ends_y, d_sum_hints,
		d_sum_lengths, d_sum_dirs, d_sol_mat, d_perms, d_t_mats, m, n,
		no_sums, d_solved);
	cudaDeviceSynchronize();
	//CUDA

	cudaMemcpy(h_sol_mat, d_sol_mat, (m * n) * sizeof(int), cudaMemcpyDeviceToHost);


	print_flattened_matrix(h_sol_mat, m, n);
	// 
	//TO DO sol_mat_flattened_to_file(mat, d_sol_mat, m, n)
	//Similiar to sol_mat, use hints from mat and values from d_sol_mat

	for (int i = 0; i < n; i++) {
		delete mat[i];
		delete sol_mat[i];
	}

	delete mat;
	delete sol_mat;

	delete h_sum_starts_x;
	delete h_sum_starts_y;
	delete h_sum_ends_x;
	delete h_sum_ends_y;
	delete h_sum_hints;
	delete h_sum_lengths;
	delete h_sum_dirs;
	delete h_sol_mat;

	cudaFree(d_t_mats);
	cudaFree(d_sum_starts_x);
	cudaFree(d_sum_starts_y);
	cudaFree(d_sum_ends_x);
	cudaFree(d_sum_ends_y);
	cudaFree(d_sum_hints);
	cudaFree(d_sum_lengths);
	cudaFree(d_sum_dirs);
	cudaFree(d_sol_mat);


	return 0;
}
