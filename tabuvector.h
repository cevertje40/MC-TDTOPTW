#pragma once
#include <unordered_map>
#include <vector>
#include <iostream>
using namespace std;
class TabuVector
{
private:
	// For each tour, we maintain a vector of tabu tenures for each vertex
	unordered_map<int,vector<int>> tabu_vectors; // key = tour index, value = tabu vector for that tour
	int tabu_duration;  // Duration for which a vertex remains tabu
	int current_iteration;  // Keeps track of the current iteration of the search
	int num_vertices;  // Number of vertices in each tour

public:
	// Constructor: Initialize the tabu vectors for a given number of vertices
	TabuVector(int num_vertices, int tabu_duration)
		: tabu_duration(tabu_duration), current_iteration(0), num_vertices(num_vertices) {}

	// Move to the next iteration
	void nextIteration()
	{
		current_iteration++;
	}

	// Ensure that the tour has an initialized tabu vector
	void ensureTourInitialized(int tour_index)
	{
		if (tabu_vectors.find(tour_index) == tabu_vectors.end())
		{
			tabu_vectors[tour_index] = std::vector<int>(num_vertices, 0);  // Initialize the vector with zeros
		}
	}

	// Mark a vertex as tabu for a specific tour
	void addTabu(int vertex_index, int tour_index)
	{
		ensureTourInitialized(tour_index);
		tabu_vectors[tour_index][vertex_index] = current_iteration + tabu_duration;
	}

	// Check if a vertex is tabu for a specific tour
	bool isTabu(int vertex_index, int tour_index)
	{
		ensureTourInitialized(tour_index);
		return tabu_vectors[tour_index][vertex_index] > current_iteration;
	}

	// Print the current contents of the tabu vector for all tours
	void printTabuVector() const
	{
		std::cout << "Tabu Vectors:\n";
		for (const auto& tour : tabu_vectors)
		{
			int tour_index = tour.first;
			const auto& tabu_vector = tour.second;
			std::cout << "Tour " << tour_index << ":\n";
			for (size_t i = 0; i < tabu_vector.size(); ++i)
			{
				if (tabu_vector[i] > current_iteration)
				{
					std::cout << "  Vertex " << i << " is tabu until iteration " << tabu_vector[i] << "\n";
				}
			}
		}
	}
};
