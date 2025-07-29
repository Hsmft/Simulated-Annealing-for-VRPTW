# Simulated Annealing for the Vehicle Routing Problem with Time Windows (VRPTW)

![Language](https://img.shields.io/badge/language-C%2B%2B-blue.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains a C++ implementation of a **Simulated Annealing** metaheuristic designed to solve the **Vehicle Routing Problem with Time Windows (VRPTW)**. VRPTW is a well-known NP-hard optimization problem in logistics and transportation, and this project aims to find high-quality, near-optimal solutions by minimizing the number of vehicles used and the total travel distance.

---

## 📋 Project Overview

The Vehicle Routing Problem with Time Windows (VRPTW) involves finding the optimal set of routes for a fleet of vehicles to serve a set of customers with specific demands, service times, and time windows, starting from and returning to a central depot.

This solver implements the following workflow:

1.  **Instance Parsing:** Reads complex VRPTW problem instances from standard text files.
2.  **Initial Solution:** Generates a feasible starting solution using a **greedy nearest-neighbor heuristic**.
3.  **Optimization:** Employs a **Simulated Annealing** algorithm to iteratively improve the initial solution by exploring the solution space.
4.  **Neighborhood Operators:** The algorithm uses a variety of neighborhood "moves" to generate new candidate solutions at each iteration, including:
    * 2-Opt (intra-route)
    * Relocate Customer (inter-route)
    * Swap Customers (inter-route)
    * Insertion Move
    * Remove and Re-insert Route
5.  **Feasibility Checks:** Ensures that all generated routes adhere to vehicle capacity and customer time window constraints.

---

## 🛠️ Technologies Used

* **Language:** C++ (utilizing C++11 features like `<chrono>` and `<random>`)
* **Libraries:** C++ Standard Library only. No external optimization libraries were used.

---

## 🚀 How to Compile and Run

### Compilation
You can compile the source code using a standard C++ compiler like g++.

```bash
g++ -std=c++11 -o solver hw1.cpp
```

### Execution
The program is run from the command line with the following arguments:

```bash
./solver [instance-file-path] [max-execution-time] [max-evaluations]
```
* `instance-file-path`: The path to the problem instance file (e.g., `instances/toy.txt`).
* `max-execution-time`: The maximum run time in seconds. Use `0` for no time limit.
* `max-evaluations`: The maximum number of objective function evaluations. Use `0` for no limit.

**Example:**
```bash
./solver instances/C101.txt 60 0
```
This command runs the solver on the `C101.txt` instance for a maximum of 60 seconds. The solution, including the routes, number of vehicles, and total distance, will be printed to the console.

---

## 📄 License
This project is licensed under the MIT License.