/**
 * \file cuda_dijkstra.c
 * \brief Dokumentirana datoteka.
 *
 * Datoteka u kojoj su dokumentirane funkcije.
 */

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <iostream>
#include <time.h>

#define PRINT_THRESHOLD 10	// Number of vertices threshold to stop printing data
#define MAX_WEIGHT 100 // Max edge weight ([0, MAX_WEIGHT])


/**
 * \brief Checks if the number of vertices is smaller then the threshold allows and returns the result.
 *
 * Accepts a whole number that represents the number of vertices in a graph that's about to run and returns a bool value depending if it meets the requirement.
 *
 * \param V number of vertices
 * \return bool, true / false
 */
bool print_threshold_check(int V)
{
	if (V < PRINT_THRESHOLD)
	{
		return true;
	}
	return false;
}

/**
 * \brief Populates the graph (array representation) with weighted edges (macro MAX_LIMIT limits the edge weight).
 *
 * \param graph array representation of a graph
 * \param N num of vertices
 */
void createGraph(float *graph, int N)	
{
	int col;
	int row;

	srand((unsigned)time(0));

	for (col = 0; col < sqrt(N); col++)
	{
		for (row = 0; row < sqrt(N); row++)
		{
			if (col != row)
			{
				graph[(int)(row * sqrt(N)) + col] = rand() % MAX_WEIGHT; // assign random
				// copy the same value on the inverse position in the matrix
				graph[(int)(col * sqrt(N)) + row] = graph[(int)(row * sqrt(N)) + col];
			}
			else
			{
				graph[(int)(col * sqrt(N)) + row] = 0; // Points a vertex on itself, therefore we set weight to 0 beacuse no loops are allowed
			}
		}
	}
}

/**
 * \brief Print the graph (array representation) to console.
 *
 * \param graph array representation of a graph
 * \param size size of the array (numOfVertices^2)
 */
void printGraph(float *graph, int size)
{
	int index;
	printf("\Graf:\n");
	for (index = 0; index < size; index++)
	{
		if (((index + 1) % (int)sqrt(size)) == 0)
		{
			printf("%5.1f\n", graph[index]);
		}
		else
		{
			printf("%5.1f ", graph[index]);
		}
	}
	printf("\n");
}

/**
 * \brief Get's minimum distance from a "src" vertex to any other vertex in the graph.
 *
 * Pick the minimum distance vertex from the set of vertices not yet processed. "u" is always equal to "src" in the first iteration. 
 * 
 * \param dist resulting array of the graph that dijskta uses in it's iterations
 * \param sptSet graph array (numOfVertices^2) that contains bools if the vertex is included in SPT
 * \param V number of vertices
 * \return integer representing index of vertex
 */
int min_distance_cpu(float *dist, bool *sptSet, int V)
{
	// Initialize min value
	float min = INT_MAX;
	float min_index;

	// Find minimum distance
	for (int v = 0; v < V; v++)
	{
		if (!sptSet[v] && dist[v] <= min) 
		{
			min = dist[v];
			min_index = v;
		}
	}

	return min_index;
}

/**
 * \brief Function that run dijskta on a set of data and calculates SPT for one vertex.
 * 
 * \param graph (array) input graph with weighted edges and no loops
 * \param src (int) index of vertex
 * \param V number of vertices
 * \param result (array) resulting graph
 */
void dijkstra_cpu(float *graph, int src, int V, float *result)
{
	// sptSet[i] will be true if vertex i is included in shortest
	// path tree or shortest distance from src to i is finalized
	bool* sptSet = (bool*)malloc(V * sizeof(bool));

	// Initialize all distances as INFINITE and sptSet[] as false
	for (int i = 0; i < V; i++)
	{
		result[i] = INT_MAX;
		sptSet[i] = false;
	}

	// Distance of source vertex from itself is always 0
	result[src] = 0;

	// Find shortest path from src
	for (int count = 0; count < V-1; count++) 
	{
		// Pick the minimum distance vertex from the set of vertices not 
		// yet processed. "u" is always equal to "src" in the first iteration. 
		int u = min_distance_cpu(result, sptSet, V);

		// Mark the picked vertex as processed
		sptSet[u] = true;

		// Update dist value of the adjacent vertices of the picked vertex
		for (int v = 0; v < V; v++)
		{
			// Update result[v] only if is not in sptSet, there is an edge from 
			// u to v, and total weight of path from src to  v through u is 
			// smaller than current value of dist[v]
			if (
				!sptSet[v]
				&& graph[(u * V) + v] && result[u] != INT_MAX
				&& result[u] + graph[(u * V) + v] < result[v]
				)
			{
				result[v] = result[u] + graph[(u * V) + v];
			}
		}
	}
}

/**
 * \brief Untility function to print the solution of dijsktra run for a specific vertex
 * 
 * \param src (int) index of vertex
 * \param dist (array) resulting graph
 * \param V number of vertices
 */
void print_solution(int src, float *dist, int V)
{
	printf("\n Vrh	Udaljenost of izvora %d\n", src);
	// Loop and print the data from "src" to index
	for (int i = 0; i < V; i++) {
		printf("%d \t\t %.1f\n", i, dist[i]);
	}
}

/**
 * \brief Untility function to print the solution of dijsktra run for a specific vertex (whole array)
 * 
 * \param src (int) index of vertex
 * \param dist (array) resulting graph
 * \param V number of vertices
 * \param size suze of array (numOfVertices^2)
 */
void print_solution_interval(int src, float* dist, int V, int size)
{
	printf("\n Vrh	Udaljenost of izvora %d\n", src);
	// Loop and print the data from "src" to index
	for (int i = src; i < size; i+=V) {
		printf("%d \t\t %.1f\n", i, dist[i]);
	}
}

/**
 * \brief Populates an array with shortest paths for all vertices based on per-vertex dijsktra calucations.
 * 
 * \param src (int) index of vertex
 * \param dist (array) containing per vertex dijkstra solution
 * \param result array that (will) contain all shortest paths
 * \param V number of vertices
 * \param size suze of array (numOfVertices^2)
 */
void populate_result_array(int src, float *dist, float *result, int V, int size)
{
	for (int i = 0; i < V; i++)
	{
		//printf("\nPutting %5.1f at result array index %d", dist[i], src + (i * V));
		result[src+(i*V)] = dist[i];
	}
}

/**
 * \brief GPU Prepare the helper graphs for dijkstra algorithm.
 * 
 * Initializes GPU resulting graph with max distances that dijsktra algoithm requires.
 * Sets all indexes of "visited" graph to false - no vertecies are included yet in the SPT exept starting vertex 
 * 
 * \param result GPU resulting array graph
 * \param visited GPU array graph containing data on every vertex included-in-SPT state
 */
__global__ void gpu_setUpGraph(float* result, bool* visited) {

	// Initialize all distances as INFINITE and stpSet[] as false
	int index = threadIdx.x + blockIdx.x * blockDim.x;

	// Initially set all vertex not to have been visited
	visited[index] = false;

	// Set all distances to INFINITE (INT_MAX will do), but set distance of vertex to itself to 0
	if (index == ((blockDim.x * blockIdx.x) + blockIdx.x)) // "index" goes through every global threadId, and this matches it to the x*x place in the matrix representation
		result[index] = 0; // distance to itself is always 0

	else result[index] = INT_MAX;
}

/**
 * \brief GPU Performs dijkstra's algorithm for every vertice in the graph in separate cores.
 * 
 * \param graph GPU initial graph array with weighted edges and no loops
 * \param result GPU resulting array graph
 * \param visited GPU array graph containing data on every vertex included-in-SPT state
 * \param V number of vertices
 */
__global__ void gpu_dijkstra_threads(float* graph, float* result, bool* visited, int V) {

	// Find shortest path for all vertices
	for (int count = 0; count < V - 1; count++)
	{
		// Pick the minimum distance vertex from the set of vertices not
		// yet processed.
		int min = INT_MAX, u;
		for (int v = 0; v < V; v++)
			if (visited[(V * threadIdx.x) + v] == false && result[(V * threadIdx.x) + v] <= min)
				min = result[(V * threadIdx.x) + v], u = v;

		// Mark the picked vertex as processed
		visited[(V * threadIdx.x) + u] = true;

		// Update the wieght value 
		for (int v = 0; v < V; v++) {

			// Update only if is not in visited, there is an edge from 
			// u to v, and total weight of path from src to  v through u is 
			// smaller than current value
			if (!visited[(V * threadIdx.x) + v] && graph[(u * V) + v] && result[(V * threadIdx.x) + u] != INT_MAX
				&& result[(V * threadIdx.x) + u] + graph[(u * V) + v] < result[(V * threadIdx.x) + v])
				result[(V * threadIdx.x) + v] = result[(V * threadIdx.x) + u] + graph[(u * V) + v];
		}
	}
}

/**
 * \brief Takes two arrays and compares values on each index.
 * 
 * \param graph1 Graph 1
 * \param graph2 Graph 2
 * \param size size of arrays
 * 
 */
void compare_results(float *graph1, float *graph2, int size)
{
	for (int i = 0; i < size; i++)
	{
		if (graph1[i] != graph2[i]) 
		{
			printf("\n\n GRESKA:\n");
			printf("Vrijednost grafa na poziciji %d se ne podudaraju! : [%5.1f, %5.1f]", i, graph1[i], graph2[i]);
		}
	}
	printf("CPU i GPU matice se podudaraju.");
}

/**
 * \brief Returns the difference between two double numbers.
 * 
 * NOTE: Assumption is made that t2 is greater then t1
 * 
 * \param t1
 * \param t2
 * 
 */
double compare_times(float t1, float t2)
{
	return (t1-t2);
}

/**
 * \brief Main function.
 */
int main()
{
	// NOTE:
	// All printing threshold checks will be performed in the main function as to maintain printing functions reusability

	/**************************** TAKE USER INPUT *****************************/

	int* numOfVertices = (int*)std::malloc(sizeof(int));
	int* arrayLength = (int*) std::malloc(sizeof(int));

	// PROMPT USER FOR # OF VERTICES
	printf("\nUnesite broj vrhova: ");
	scanf("%d", numOfVertices);

	// WILL BE AN ARRAY REPRESENTATION OF A MATRIX
	*arrayLength = *numOfVertices * *numOfVertices;

	printf("Broj vrhova je %d, pa je velicina listne repreznetacije grafa %d\n", *numOfVertices, *arrayLength);

	// Store the print requirement bool
	bool* print_req = (bool*)std::malloc(sizeof(bool));
	*print_req = print_threshold_check(*numOfVertices);

	// Writing to console if the data will be printed or not
	if (*print_req)
	{
		printf("\nBroj vrhova je manji od limita printanja pa printam sve podatke matrica u konzolu.\n\n");
	}
	else
	{
		printf("\nBroj vrhova je veci od limita printanja pa ne printam podatke matrica u konzolu.\n\n");
	}

	// ============================== CPU ==============================

	// Variables to store time for dijkstra on CPU
	clock_t cpu_start, cpu_end;
	// Variable to store total time
	double cpu_total_time = 0;

	// Allocate CPU memory for initial and resulting graphs
	float* graph = (float*)malloc(*arrayLength * sizeof(float));
	float* result = (float*)malloc(*arrayLength * sizeof(float));
	float* result_graph = (float*)malloc(*arrayLength * sizeof(float));

	// Fill the graph with data
	createGraph(graph, *arrayLength);

	printf("\nGraf inicijaliziran i popunjen nasumičnim vrijednostima.");
	printf("\nTezine veza u grafu su u intervalu [%5.1f, %5.1f].\n", 0, (float)MAX_WEIGHT);

	// print the graph
	if (*print_req) 
	{
		printGraph(graph, *arrayLength);
	}

	if (*print_req)
	{
		printf("\nPrintam udaljenost od svakog vrha (CPU):\n");
	}

	// Run the dijsksta on every vertex and populate result array
	for (int vertex = 0; vertex < *numOfVertices; vertex++)
	{
		cpu_start = clock(); // start time of this iteration

		// dijkstra calculates shortest paths for vertex "vertex"
		dijkstra_cpu(graph, vertex, *numOfVertices, result);

		cpu_end = clock(); // end time of this iteration

		// Add time required for this dijsktra iteration to the total CPU time count
		cpu_total_time += (((double)cpu_end - (double)cpu_start) / CLOCKS_PER_SEC);

		if (*print_req)
		{
			// print the solution for that vertex
			print_solution(vertex, result, *numOfVertices);
		}

		// Fill the data for this dijsktra iteration in a full resulting graph
		populate_result_array(vertex, result, result_graph, *numOfVertices, *arrayLength);
	}

	// ============================== GPU ==============================

	// create GPU CUDA events to measure times
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// Initialize variable to store GPU resulting data
	float* gpu_result;

	// Allocate memory for device result array
	gpu_result = (float*)malloc(*arrayLength * sizeof(float));

	// Initialize GPU variables
	float* dev_graph;
	float* dev_result;
	bool* dev_visited;

	// Allocate memory for device variables
	cudaMalloc( (void**)&dev_graph, (*arrayLength * sizeof(float)) );
	cudaMalloc( (void**)&dev_result, (*arrayLength * sizeof(float)) );
	cudaMalloc( (void**)&dev_visited, (*arrayLength * sizeof(bool)) );

	// Copy CPU generated graph data to the device (GPU)
	cudaMemcpy(dev_graph, graph, (*arrayLength * sizeof(float)), cudaMemcpyHostToDevice);

	// Set up data on the GPU for dijkstra calculations
	gpu_setUpGraph << < *numOfVertices, *numOfVertices >> > (dev_result, dev_visited); // Every block executes one matrix

	cudaEventRecord(start); // Mark the event of start of GPU calculations

	// Perform dijstra on ALL vertices as src vertex using multiple threads
	//gpu_dijkstra_blocks << < *numOfVertices, 1 >> > (dev_graph, dev_result, dev_visited, *numOfVertices);
	gpu_dijkstra_threads << < 1, *numOfVertices >> > (dev_graph, dev_result, dev_visited, *numOfVertices);

	cudaEventRecord(stop); // Mark the event of end of GPU calcuations
	cudaEventSynchronize(stop);

	float gpu_total_time = 0; // stores total time required for GPU calucations
	cudaEventElapsedTime(&gpu_total_time, start, stop); // Calculates the time based on events

	// Copy result from GPU calculations to CPU (host)
	cudaMemcpy(gpu_result, dev_result, (*arrayLength * sizeof(float)), cudaMemcpyDeviceToHost);

	if (*print_req) 
	{
		// Printing by-vertex solutions
		printf("\nPrintam udaljenost od svakog vrha (GPU):\n");
		for (int v = 0; v < *numOfVertices; v++)
		{
			print_solution_interval(v, gpu_result, *numOfVertices, *arrayLength);
		}
	}

	if (*print_req)
	{
		// Printing resulting graph of CPU calculations
		printf("\nIspisujem rezultantnu matricu sa CPU (host):\n");
		printGraph(result_graph, *arrayLength);

		// Printing resulting graph of GPU calculations
		printf("\nIspisujem razultatnu marticu sa GPU (device):\n");
		printGraph(gpu_result, *arrayLength);
	}

	// Compare the two resulting arrays
	printf("\nUsporedujem...\n");
	compare_results(result_graph, gpu_result, *arrayLength);

	printf("\n\nIspisujem vremena:\n");

	printf("Potrebno vrijeme na CPU-u: %.10fs\n", cpu_total_time);
	printf("Potrebno vrijeme na GPU-u: %.10fs\n", gpu_total_time * .0001);

	printf("Vrijeme potrebno za izracun je %.10fs manje na GPU-u nego na CPU-u.", compare_times((double)cpu_total_time, gpu_total_time * 0.0001)); // Compare the times

	// Some more memory management
	//CPU
	free(numOfVertices);
	free(arrayLength);
	free(graph);
	free(result);
	free(result_graph);
	free(gpu_result);
	// GPU
	cudaFree(dev_graph);
	cudaFree(dev_result);
	cudaFree(dev_visited);

	return 0;
}

