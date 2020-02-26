#pragma once
#include <iostream>
#include <random>
#include <vector>
#include "Matrix.h"

struct Coordinate{
	int x, y,z;
};


class Brain {
protected:
	double Get_Square(double value) {
		return value * value;
	}
	double squared_3Point_distance(Coordinate first, Coordinate second) {
		return Get_Square(first.x - second.x) + Get_Square(first.y - second.y) + Get_Square(first.z - second.z);
	}



public:
	template <class Data_Type>
	std::vector<Data_Type> K_Means(const std::vector<Data_Type>& data, size_t k, size_t number_of_iterations,int number_of_points_struct = 3) {
		static std::random_device seed;
		static std::mt19937 random_number_generator(seed()); //merssene twisster
		std::uniform_int_distribution<size_t> indices(0, data.size() - 1);
		std::vector<Data_Type> means(k);
		
		for (auto& cluster : means) {
			cluster = data[indices(random_number_generator)];
		}

		std::vector<size_t> assignments(data.size());



		if (number_of_points_struct == 5) {

			for (size_t iteration = 0; iteration < number_of_iterations; ++iteration) {
				// Find assignments.
				for (size_t point = 0; point < data.size(); ++point) {
					double best_distance = std::numeric_limits<double>::max();
					size_t best_cluster = 0;
					for (size_t cluster = 0; cluster < k; ++cluster) {
						const double distance =
							squared_3Point_distance(data[point], means[cluster]);
						if (distance < best_distance) {
							best_distance = distance;
							best_cluster = cluster;
						}
					}
					assignments[point] = best_cluster;
				}

				// Sum up and count points for each cluster.
				std::vector<Data_Type> new_means(k);
				std::vector<size_t> counts(k, 0);


				for (size_t point = 0; point < data.size(); ++point) {
					const auto cluster = assignments[point];
					new_means[cluster].a += data[point].a;
					new_means[cluster].b += data[point].b;
					new_means[cluster].c += data[point].c;
					new_means[cluster].d += data[point].d;
					new_means[cluster].e += data[point].e;
					counts[cluster] += 1;
				}

				// Divide sums by counts to get new centroids.
				for (size_t cluster = 0; cluster < k; ++cluster) {
					// Turn 0/0 into 0/1 to avoid zero division.
					const auto count = std::max<size_t>(1, counts[cluster]);
					means[cluster].a = new_means[cluster].a / count;
					means[cluster].b = new_means[cluster].b / count;
					means[cluster].c = new_means[cluster].c / count;
					means[cluster].d = new_means[cluster].d / count;
					means[cluster].e = new_means[cluster].e / count;

				}
			}

		}


		if (number_of_points_struct == 4) {

			for (size_t iteration = 0; iteration < number_of_iterations; ++iteration) {
				// Find assignments.
				for (size_t point = 0; point < data.size(); ++point) {
					double best_distance = std::numeric_limits<double>::max();
					size_t best_cluster = 0;
					for (size_t cluster = 0; cluster < k; ++cluster) {
						const double distance =
							squared_3Point_distance(data[point], means[cluster]);
						if (distance < best_distance) {
							best_distance = distance;
							best_cluster = cluster;
						}
					}
					assignments[point] = best_cluster;
				}

				// Sum up and count points for each cluster.
				std::vector<Data_Type> new_means(k);
				std::vector<size_t> counts(k, 0);


				for (size_t point = 0; point < data.size(); ++point) {
					const auto cluster = assignments[point];
					new_means[cluster].a += data[point].a;
					new_means[cluster].b += data[point].b;
					new_means[cluster].c += data[point].c;
					new_means[cluster].d += data[point].d;
					counts[cluster] += 1;
				}

				// Divide sums by counts to get new centroids.
				for (size_t cluster = 0; cluster < k; ++cluster) {
					// Turn 0/0 into 0/1 to avoid zero division.
					const auto count = std::max<size_t>(1, counts[cluster]);
					means[cluster].a = new_means[cluster].a / count;
					means[cluster].b = new_means[cluster].b / count;
					means[cluster].c = new_means[cluster].c / count;
					means[cluster].d = new_means[cluster].d / count;

				}
			}

		}

		if (number_of_points_struct == 3) {

			for (size_t iteration = 0; iteration < number_of_iterations; ++iteration) {
				// Find assignments.
				for (size_t point = 0; point < data.size(); ++point) {
					double best_distance = std::numeric_limits<double>::max();
					size_t best_cluster = 0;
					for (size_t cluster = 0; cluster < k; ++cluster) {
						const double distance =
							squared_3Point_distance(data[point], means[cluster]);
						if (distance < best_distance) {
							best_distance = distance;
							best_cluster = cluster;
						}
					}
					assignments[point] = best_cluster;
				}

				// Sum up and count points for each cluster.
				std::vector<Data_Type> new_means(k);
				std::vector<size_t> counts(k, 0);


				for (size_t point = 0; point < data.size(); ++point) {
					const auto cluster = assignments[point];
					new_means[cluster].x += data[point].x;
					new_means[cluster].y += data[point].y;
					new_means[cluster].z += data[point].z;
					counts[cluster] += 1;
				}

				// Divide sums by counts to get new centroids.
				for (size_t cluster = 0; cluster < k; ++cluster) {
					// Turn 0/0 into 0/1 to avoid zero division.
					const auto count = std::max<size_t>(1, counts[cluster]);
					means[cluster].x = new_means[cluster].x / count;
					means[cluster].y = new_means[cluster].y / count;
					means[cluster].z = new_means[cluster].z / count;

				}
			}

		}
		if (number_of_points_struct == 2) {

			for (size_t iteration = 0; iteration < number_of_iterations; ++iteration) {
				// Find assignments.
				for (size_t point = 0; point < data.size(); ++point) {
					double best_distance = std::numeric_limits<double>::max();
					size_t best_cluster = 0;
					for (size_t cluster = 0; cluster < k; ++cluster) {
						const double distance =
							squared_3Point_distance(data[point], means[cluster]);
						if (distance < best_distance) {
							best_distance = distance;
							best_cluster = cluster;
						}
					}
					assignments[point] = best_cluster;
				}

				// Sum up and count points for each cluster.
				std::vector<Data_Type> new_means(k);
				std::vector<size_t> counts(k, 0);


				for (size_t point = 0; point < data.size(); ++point) {
					const auto cluster = assignments[point];
					new_means[cluster].x += data[point].x;
					new_means[cluster].y += data[point].y;
					counts[cluster] += 1;
				}

				// Divide sums by counts to get new centroids.
				for (size_t cluster = 0; cluster < k; ++cluster) {
					// Turn 0/0 into 0/1 to avoid zero division.
					const auto count = std::max<size_t>(1, counts[cluster]);
					means[cluster].x = new_means[cluster].x / count;
					means[cluster].y = new_means[cluster].y / count;

				}
			}

		}

		return means;

	}

};