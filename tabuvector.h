#pragma once
#include <unordered_map>
#include <vector>
#include <iostream>
using namespace std;
class TabuVector {
private:
    std::vector<std::vector<int>> tabu_vectors;  // A fixed-size vector of tabu vectors for each tour
    int tabu_duration;
    int current_iteration;
    int num_tours;
    int num_vertices;

public:
    // Constructor
    TabuVector(int num_tours, int num_vertices, int tabu_duration)
        : tabu_vectors(num_tours, std::vector<int>(num_vertices, 0)),
        tabu_duration(tabu_duration), current_iteration(0), num_tours(num_tours), num_vertices(num_vertices) {}

    // Move to the next iteration
    void nextIteration() {
        current_iteration++;
    }

    // Mark a vertex as tabu
    void addTabu(int vertex_index, int tour_index) {
        tabu_vectors[tour_index][vertex_index] = current_iteration + tabu_duration;
    }

    // Check if a vertex is tabu
    bool isTabu(int vertex_index, int tour_index) {
        return tabu_vectors[tour_index][vertex_index] > current_iteration;
    }

    void setDuration(int duration) { tabu_duration = duration; }

    //Print tabu vectors
    void printTabuVector() const {
        for (int tour_index = 0; tour_index < num_tours; ++tour_index) {
            for (int vertex_index = 0; vertex_index < num_vertices; ++vertex_index) {
                if (tabu_vectors[tour_index][vertex_index] > current_iteration) {
                    std::cout << "Vertex " << vertex_index << " in tour " << tour_index << " is tabu until " << tabu_vectors[tour_index][vertex_index] << "\n";
                }
            }
        }
    }
};
