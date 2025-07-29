#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <random>
#include <ctime>
#include <iomanip>
#include <chrono>

using namespace std;

// Structure to represent a customer with their details
struct Customer {
    int id;
    double x, y;
    int demand;
    double serviceTime;
    double earliest; 
    double latest; } ;

// Structure to represent a route with customer IDs, distances, and times
struct Route {
    vector<int> customerIds;
    double totalDistance;
    vector<double> arrivalTimes;  
    vector<double> departureTimes; 
    int load; // Total demand of customers in the route
    Route(): totalDistance(0.0), load(0) {}
    };

// Structure of problem data
struct ProblemData {
    vector<Customer> customers;
    int vehicleCapacity;
    int maxVehicles;
    double maxRouteDuration;
    };


// Euclidean distance between two customer
double distance(const Customer &a ,const Customer &b) 
{return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));}

// Reads problem instance
ProblemData readInstance(const string &filename) {
    ProblemData data;
    ifstream infile(filename);
    if (!infile) {cerr << "Cannot open file: " << filename << endl;
        exit(1);}
    string line;

    // Skip header lines and read vehicle info
    while (getline(infile, line)){
        if (line.find("CUST NO.") != string::npos) {
            getline(infile, line); // Skip the empty line
            break;}
        else if (line.find("NUMBER") != string::npos){
            getline(infile, line);
            istringstream iss(line);
            iss >> data.maxVehicles >> data.vehicleCapacity;}
    }
    // Read customer data
    int numCustomers = 0;
    while (getline(infile, line)) {
        istringstream issCust(line);
        Customer cust;
        if (issCust >> cust.id >> cust.x >> cust.y >> cust.demand >> cust.earliest >>  cust.latest >> cust.serviceTime) {
            data.customers.push_back(cust);
            numCustomers++;
        } else {break; // Stop if line is invalid
        }}
    data.maxRouteDuration = data.customers[0].latest; // Set from depot's due time
    infile.close();
    return data;
}

// Checks if a single route is feasible based on capacity, time windows, and duration
bool isRouteFeasible(const Route &route, const ProblemData &data) {
    int capacityUsed = 0;
    double currentTime = 0.0;
    int currentIndex = 0; // Start from depot (index 0)

    for (int custId : route.customerIds) {
        // Find customer index by ID
        auto it = find_if(data.customers.begin(), data.customers.end(), [custId](const Customer &c) {
            return c.id == custId; });
        if (it == data.customers.end())  {
             cerr << "Customer with id " << custId << " not found in isRouteFeasible." << endl; 
            return false; }
        int index = distance(data.customers.begin(), it);
        // Capacity check
        capacityUsed += data.customers[index].demand;
        if (capacityUsed > data.vehicleCapacity) {
            return false;}

        // Time window check
        double travelTime = distance(data.customers[currentIndex], data.customers[index]);
        double arrivalTime = currentTime + travelTime;
        double serviceStartTime = max(arrivalTime, data.customers[index].earliest);
        if (serviceStartTime > data.customers[index].latest) {
            return false;}
        currentTime = serviceStartTime + data.customers[index].serviceTime;
        currentIndex = index; // Update current index for next customer
        }
    // Return to depot
    double returnTravelTime = distance(data.customers[currentIndex], data.customers[0]);
    currentTime += returnTravelTime;
    if (currentTime > data.maxRouteDuration) {
        return false; }
    return true;
}

// Checks if the all routes is feasible
bool isSolutionFeasible(const vector<Route>& routes, const ProblemData& data) {
    // 1. Check feasibility of each route
    for (const auto& route : routes) {
        if (!isRouteFeasible(route, data)) {
        return false;  }}
    // 2. Check if all customers are visited exactly once
    vector<bool> visited(data.customers.size(), false);
    visited[0] = true; // Mark depot as visited
    for (const auto& route : routes) {
        for (int custId : route.customerIds) {
             auto it = find_if(data.customers.begin(), data.customers.end(), [custId](const Customer &c) { return c.id == custId;});
                int index = distance(data.customers.begin(), it);
            if (visited[index]) {
                return false; // Customer visited more than once
            }
            visited[index] = true;  } }
    // Check if all customers (except depot) are visited
    for (size_t i = 1; i < visited.size(); ++i) {
        if (!visited[i]) {
            return false; }  }
    return true;
}

// Finds the index of a customer in the customers list by their ID
int find_customer_index(const vector<Customer>& customers, int custId) {
    for (size_t i = 0; i < customers.size(); ++i) {
        if (customers[i].id == custId) {
            return i;
        }
    } return -1; // Not found 
}

// Constructs an initial solution using a greedy nearestneighbor 
vector<Route> constructInitialSolution(const ProblemData &data) {
    vector<Route> routes;
    vector<bool> visited(data.customers.size(), false);
    visited[0] = true; // Depot is always visited
    while (find(visited.begin() + 1, visited.end(), false) != visited.end()) {
        Route route;
        route.arrivalTimes.push_back(0.0); // Depot arrival time (always 0)
        double currentTime = 0.0;
        int currentId = 0; // Start from the depot
        int currentLoad = 0;
        while (true) {
            int nextCustomer = -1; // Initialize to -1
            double bestDist = numeric_limits<double>::max(); // Initialize to infinity
            for (size_t i = 1; i < data.customers.size(); i++) { // Start from 1 (skip depot)
                if (!visited[i]) {
                    double travelTime = distance(data.customers[currentId], data.customers[i]);
                    double arrivalTime = currentTime + travelTime;
                    // Correct time window check (wait if necessary)
                    if (arrivalTime < data.customers[i].earliest) {
                        arrivalTime = data.customers[i].earliest;}
                    if (arrivalTime <= data.customers[i].latest) {
                        if (travelTime < bestDist) {
                            bestDist = travelTime;
                            nextCustomer = i;} } }}
            if (nextCustomer == -1) break; // No more feasible customers
            // Calculate arrival and departure times *once*
            double arrivalTime = currentTime + distance(data.customers[currentId], data.customers[nextCustomer]);
            if (arrivalTime < data.customers[nextCustomer].earliest) {
                arrivalTime = data.customers[nextCustomer].earliest;}
            double departureTime = arrivalTime + data.customers[nextCustomer].serviceTime;
             int newLoad = currentLoad + data.customers[nextCustomer].demand; // Corrected line
            // check capacity 
            if (newLoad <= data.vehicleCapacity) { 
                Route tempRoute = route;
                tempRoute.customerIds.push_back(data.customers[nextCustomer].id);
                tempRoute.arrivalTimes.push_back(arrivalTime);
                tempRoute.departureTimes.push_back(departureTime);
                tempRoute.load = newLoad;
                //check full feasibility (time windows)
                if (isRouteFeasible(tempRoute, data)) {
                    // If feasible (both capacity and time windows), add to the actual route
                    route = tempRoute; // Copy the temporary route to the actual route
                    currentTime = departureTime; // Update current time
                    visited[nextCustomer] = true; // Mark as visited
                    currentId = data.customers[nextCustomer].id; 
                    currentLoad = newLoad; // Update load
                } else { break; // This customer is infeasible try a new route
                    }
            } else {   break; // Capacity exceeded try a new route
         }
        } // Inner while loop (building a single route)
        // Add return-to-depot information *only if the route is not empty*
        if (!route.customerIds.empty()) {
            int lastCustomerIndex = find_customer_index(data.customers, currentId); // Get the *index*
            double returnTime = currentTime + distance(data.customers[lastCustomerIndex], data.customers[0]);
            route.arrivalTimes.push_back(returnTime);
            route.departureTimes.push_back(returnTime);
            route.totalDistance += distance(data.customers[lastCustomerIndex], data.customers[0]);
            routes.push_back(route);}} // Outer while loop (building all routes)
    return routes;}



//objective function: number of vehicles and total distance
pair<double, double> objectiveFunction  (const vector<Route> &routes) {
    int vehicles = 0;
    double totalDistance = 0.0;
    for (const auto &route : routes) {
        if (!route.customerIds.empty()) { 
            vehicles++;
            totalDistance += route.totalDistance; }}
    return {vehicles, totalDistance}; }

// Compares two solutions to determine if the new one is better
// bool BetterSolution(const vector<Route> &newRoutes, const vector<Route> &currentRoutes) {
//     auto [vehiclesNew, distNew] = objectiveFunction(newRoutes);
//     auto [vehiclesCurrent, distCurrent] = objectiveFunction(currentRoutes);
//     double costNew = vehiclesNew * 9000 + distNew;
//     double costCurrent = vehiclesCurrent * 9000 + distCurrent;
//     return costNew < costCurrent;}
    
bool BetterSolution(const vector<Route> &newRoutes, const vector<Route> &currentRoutes) {
    auto [vehiclesNew, distNew] = objectiveFunction(newRoutes);
    auto [vehiclesCurrent, distCurrent] = objectiveFunction(currentRoutes);
    if (vehiclesNew < vehiclesCurrent) return true;
    if (vehiclesNew > vehiclesCurrent) return false;
    return distNew < distCurrent;
}

//prioritize distance
// bool BetterSolution(const vector<Route> &newRoutes, const vector<Route> &currentRoutes) {
//     auto [vehiclesNew, distNew] = objectiveFunction(newRoutes);
//     auto [vehiclesCurrent, distCurrent] = objectiveFunction(currentRoutes);
//     if (distNew < distCurrent) return true;
//     if (distNew > distCurrent) return false;
//     return vehiclesNew < vehiclesCurrent;
// }

// Function to update route times and load after a change
void updateRoute(Route& route, const ProblemData& data)  {
    route.arrivalTimes.clear();
    route.departureTimes.clear();
    route.load = 0;
    route.totalDistance = 0;
    if (route.customerIds.empty()) {
        return;}
    double currentTime = 0.0;
    route.arrivalTimes.push_back(currentTime); //start time with depot
    int currentIndex = 0; 
    for (int custId : route.customerIds) {
        if (custId < 0 || custId >= static_cast<int>(data.customers.size())) {
            cerr << "Invalid customer ID " << custId << " in updateRoute!" << endl;
              return;}
        double travelTime = distance(data.customers[currentIndex], data.customers[custId]);
        currentTime += travelTime;
        route.totalDistance += travelTime;
        if (currentTime < data.customers[custId].earliest) {
            currentTime = data.customers[custId].earliest;}
        route.arrivalTimes.push_back(currentTime);
        currentTime += data.customers[custId].serviceTime;
        route.departureTimes.push_back(currentTime);
        route.load += data.customers[custId].demand;
        currentIndex = custId ;}
    // return to depot
    double return_dist = distance(data.customers[currentIndex], data.customers[0]);
    route.totalDistance += return_dist;
    currentTime += return_dist;
    route.arrivalTimes.push_back(currentTime);
    route.departureTimes.push_back(currentTime);}

// Validates if a customer ID is within the valid range
bool isValidCustomerId(int custId, const ProblemData& data) {
    return custId >= 0 && custId < data.customers.size(); }

// Removes a route and redistributes its customers to other routes
void removeRouteMove(std::vector<Route>& routes, std::mt19937& rng, const ProblemData& data) {
    // At least one route must remain
    if (routes.size() <= 1) {
        return;}
    //select a route to remove randomly
    std::uniform_int_distribution<size_t> routeDist(0, routes.size() - 1);
    size_t worstRouteIdx = routeDist(rng);
    //Store customers of the selected route and remove the route
    std::vector<int> customersToRedistribute = routes[worstRouteIdx].customerIds;
    routes.erase(routes.begin() + worstRouteIdx);
    //  Redistribute customers randomly
    for (int customer : customersToRedistribute) {
        bool inserted = false;
        // Create and shuffle route indices
        std::vector<size_t> routeIndices(routes.size());
        std::iota(routeIndices.begin(), routeIndices.end(), 0);
        std::shuffle(routeIndices.begin(), routeIndices.end(), rng);
        for (size_t r : routeIndices) {
            if (inserted) break;
            // create and shuffle position indices
            std::vector<size_t> positions(routes[r].customerIds.size() + 1);
            std::iota(positions.begin(), positions.end(), 0);
            std::shuffle(positions.begin(), positions.end(), rng);
            for (size_t pos : positions) {
                routes[r].customerIds.insert(routes[r].customerIds.begin() + pos, customer);
                updateRoute(routes[r], data); // ypdate the route
                if (isRouteFeasible(routes[r], data)) {  //Check feasibility
                    inserted = true;
                    break;
                } else {
                    routes[r].customerIds.erase(routes[r].customerIds.begin() + pos);
                    updateRoute(routes[r], data);}}}
        // If the customer couldn't be inserted in any route, create a new route
        if (!inserted) {
            Route newRoute;
            newRoute.customerIds.push_back(customer);
            updateRoute(newRoute, data);
            if (isRouteFeasible(newRoute, data)) {
                routes.push_back(newRoute);}}}
    // Remove empty routes
    routes.erase(std::remove_if(routes.begin(), routes.end(),
                                [](const Route& r) { return r.customerIds.empty(); }),
                 routes.end());}

void twoOpt(std::vector<Route>& routes, std::mt19937& rng, const ProblemData& data) {
    //randomly select a route
    if (routes.empty()) return;
    std::uniform_int_distribution<size_t> routeDist(0, routes.size() - 1);
    size_t routeIdx = routeDist(rng);
    while (routes[routeIdx].customerIds.empty()) {
        routeIdx = routeDist(rng);}
    Route& route = routes[routeIdx];
    if (route.customerIds.size() < 4) return; //  At least 4 customers are needed for 2 Opt
    // randomlyy select two positions for swapping
    std::uniform_int_distribution<size_t> posDist(1, route.customerIds.size() - 3);
    size_t i = posDist(rng);
    size_t j = i + 1 + posDist(rng) % (route.customerIds.size() - i - 2);
    if (j >= route.customerIds.size()) j = route.customerIds.size() - 1;
    //Store the current route
    std::vector<int> originalCustomers = route.customerIds;
    // reverse the segment between i and j
    std::reverse(route.customerIds.begin() + i, route.customerIds.begin() + j + 1);
    updateRoute(route, data);
    //If the new route is infeasible, revert to the original
    if (!isRouteFeasible(route, data)) {
        route.customerIds = originalCustomers;
        updateRoute(route, data);}}

void relocateCustomer(std::vector<Route>& routes, std::mt19937& rng, const ProblemData& data) {
    //at least two routes are needed
    if (routes.size() < 2) {
        return;}
    //randomly select a source route
    std::uniform_int_distribution<size_t> routeDist(0, routes.size() - 1);
    size_t sourceRouteIdx = routeDist(rng);
    while (routes[sourceRouteIdx].customerIds.empty()) { // Ensure the route is not empty
        sourceRouteIdx = routeDist(rng);}
    std::uniform_int_distribution<size_t> custDist(0, routes[sourceRouteIdx].customerIds.size() - 1);
    size_t custPos = custDist(rng);
    int customer = routes[sourceRouteIdx].customerIds[custPos];
    // remove the customer from the source route
    routes[sourceRouteIdx].customerIds.erase(routes[sourceRouteIdx].customerIds.begin() + custPos);
    updateRoute(routes[sourceRouteIdx], data);
    //Randomly select a destination route (different from the source)
    size_t destRouteIdx = routeDist(rng);
    while (destRouteIdx == sourceRouteIdx) {
        destRouteIdx = routeDist(rng);}
    bool inserted = false;
    // Create and shuffle position indices
    std::vector<size_t> positions(routes[destRouteIdx].customerIds.size() + 1);
    std::iota(positions.begin(), positions.end(), 0);
    std::shuffle(positions.begin(), positions.end(), rng);
    // Try inserting the customer into the destination route
    for (size_t pos : positions) {
        routes[destRouteIdx].customerIds.insert(routes[destRouteIdx].customerIds.begin() + pos, customer);
        updateRoute(routes[destRouteIdx], data);
        if (isRouteFeasible(routes[destRouteIdx], data)) {
            inserted = true;
            break;
        } else {
            routes[destRouteIdx].customerIds.erase(routes[destRouteIdx].customerIds.begin() + pos);
            updateRoute(routes[destRouteIdx], data); }}
    // If insertion failed, try other routes
    if (!inserted) {
        std::vector<size_t> routeIndices(routes.size());
        std::iota(routeIndices.begin(), routeIndices.end(), 0);
        std::shuffle(routeIndices.begin(), routeIndices.end(), rng);
        for (size_t r : routeIndices) {
            if (r == sourceRouteIdx) continue; // Skip the source route
            if (inserted) break;
            std::vector<size_t> positions(routes[r].customerIds.size() + 1);
            std::iota(positions.begin(), positions.end(), 0);
            std::shuffle(positions.begin(), positions.end(), rng);
            for (size_t pos : positions) {
                routes[r].customerIds.insert(routes[r].customerIds.begin() + pos, customer);
                updateRoute(routes[r], data);
                if (isRouteFeasible(routes[r], data)) {
                    inserted = true;
                    break;
                } else {
                    routes[r].customerIds.erase(routes[r].customerIds.begin() + pos);
                    updateRoute(routes[r], data);
                }}}}
    // If the customer couldn't be inserted in any route, create a new route
    if (!inserted) {
        Route newRoute;
        newRoute.customerIds.push_back(customer);
        updateRoute(newRoute, data);
        if (isRouteFeasible(newRoute, data)) {
            routes.push_back(newRoute);}}
    routes.erase(std::remove_if(routes.begin(), routes.end(),
                                [](const Route& r) { return r.customerIds.empty(); }),
                    routes.end());}

void insertionMove(std::vector<Route>& routes, std::mt19937& rng, const ProblemData& data) {
    if (routes.empty()) {
        return;}
    std::uniform_int_distribution<size_t> routeDist(0, routes.size() - 1);
    size_t sourceRouteIdx = routeDist(rng);
    while (routes[sourceRouteIdx].customerIds.empty()) {
        sourceRouteIdx = routeDist(rng);}
    std::uniform_int_distribution<size_t> custDist(0, routes[sourceRouteIdx].customerIds.size() - 1);
    size_t custPos = custDist(rng);
    int customer = routes[sourceRouteIdx].customerIds[custPos];
    routes[sourceRouteIdx].customerIds.erase(routes[sourceRouteIdx].customerIds.begin() + custPos);
    updateRoute(routes[sourceRouteIdx], data);
    double minInsertionCost = std::numeric_limits<double>::max();
    size_t bestRouteIdx = -1;
    size_t bestPos = -1;
    for (size_t r = 0; r < routes.size(); ++r) {
        if (r == sourceRouteIdx) continue;
        for (size_t pos = 0; pos <= routes[r].customerIds.size(); ++pos) {
            double insertionCost = 0.0;
            if (routes[r].customerIds.empty()) {
                insertionCost = 2 * distance(data.customers[0], data.customers[customer]);
            } else {
                int prevCust = (pos == 0) ? 0 : routes[r].customerIds[pos - 1];
                int nextCust = (pos == routes[r].customerIds.size()) ? 0 : routes[r].customerIds[pos];
                insertionCost = distance(data.customers[prevCust], data.customers[customer]) +
                                distance(data.customers[customer], data.customers[nextCust]) -
                                distance(data.customers[prevCust], data.customers[nextCust]);}
            routes[r].customerIds.insert(routes[r].customerIds.begin() + pos, customer);
            updateRoute(routes[r], data);
            if (isRouteFeasible(routes[r], data) && insertionCost < minInsertionCost) {
                minInsertionCost = insertionCost;
                bestRouteIdx = r;
                bestPos = pos;}
            routes[r].customerIds.erase(routes[r].customerIds.begin() + pos);
            updateRoute(routes[r], data);}}
    bool inserted = false;
    if (bestRouteIdx != -1) {
        routes[bestRouteIdx].customerIds.insert(routes[bestRouteIdx].customerIds.begin() + bestPos, customer);
        updateRoute(routes[bestRouteIdx], data);
        inserted = true;}
    if (!inserted) {
        Route newRoute;
        newRoute.customerIds.push_back(customer);
        updateRoute(newRoute, data);
        if (isRouteFeasible(newRoute, data)) {
            routes.push_back(newRoute);}}
    routes.erase(std::remove_if(routes.begin(), routes.end(),
                                [](const Route& r) { return r.customerIds.empty(); }),
                 routes.end());}

void swapCustomers(std::vector<Route>& routes, std::mt19937& rng, const ProblemData& data) {
    // At least two routes are needed
    if (routes.size() < 2) {
        return;}
    //Randomly select two routes
    std::uniform_int_distribution<size_t> routeDist(0, routes.size() - 1);
    size_t route1Idx = routeDist(rng);
    while (routes[route1Idx].customerIds.empty()) { // Ensure the route is not empty
        route1Idx = routeDist(rng);}
    size_t route2Idx = routeDist(rng);
    while (route2Idx == route1Idx || routes[route2Idx].customerIds.empty()) {
        route2Idx = routeDist(rng);}
    //Randomly select one customer from each route
    std::uniform_int_distribution<size_t> custDist1(0, routes[route1Idx].customerIds.size() - 1);
    size_t pos1 = custDist1(rng);
    int customer1 = routes[route1Idx].customerIds[pos1];
    std::uniform_int_distribution<size_t> custDist2(0, routes[route2Idx].customerIds.size() - 1);
    size_t pos2 = custDist2(rng);
    int customer2 = routes[route2Idx].customerIds[pos2];
    // Swap the customers
    routes[route1Idx].customerIds[pos1] = customer2;
    routes[route2Idx].customerIds[pos2] = customer1;
    //Update the routes
    updateRoute(routes[route1Idx], data);
    updateRoute(routes[route2Idx], data);
    //  If the routes become infeasible, undo the swap
    if (!isRouteFeasible(routes[route1Idx], data) || !isRouteFeasible(routes[route2Idx], data)) {
        routes[route1Idx].customerIds[pos1] = customer1;
        routes[route2Idx].customerIds[pos2] = customer2;
        updateRoute(routes[route1Idx], data);
        updateRoute(routes[route2Idx], data);}
}

std::vector<Route> simulatedAnnealing(const ProblemData &data, std::vector<Route> initRoutes,
    double T0, double coolingRate, int iterPerTemp,
    int maxNoImprovement, int maxTime, int maxEvaluations,
    int ThresholdReheat, double reheatFactor) {
    std::vector<Route> currentRoutes = initRoutes;
    std::vector<Route> bestRoutes = initRoutes;
    if (!currentRoutes.empty()) {
        for (auto& route : currentRoutes) updateRoute(route, data);
        for (auto& route : bestRoutes) updateRoute(route, data);}
    double T = T0;
    int noImprovementCount = 0;
    auto startTime = std::chrono::steady_clock::now();
    int numEvaluations = 0;
    std::mt19937 rng(static_cast<unsigned>(std::time(nullptr)));
    // neighborhoods weights
    const std::vector<double> operator_weights = {1.0, 1.0, 1.0, 1.0, 1.0}; 
    std::discrete_distribution<int> opDist(operator_weights.begin(), operator_weights.end());
    std::cout << "Weights: {";
    for(size_t i = 0; i < operator_weights.size(); ++i) { std::cout << operator_weights[i] << (i == operator_weights.size() - 1 ? "" : ", "); }
    std::cout << "}" << std::endl;
    std::cout << "------------------------" << std::endl;
    while (T > 1e-4 && noImprovementCount < maxNoImprovement && (maxEvaluations == 0 || numEvaluations < maxEvaluations)) {
        bool improved = false;
        for (int iter = 0; iter < iterPerTemp; iter++) {
            auto currentTime = std::chrono::steady_clock::now();
            double elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
            if (maxTime > 0 && elapsedTime >= maxTime) { goto endLoop; }
            if (maxEvaluations > 0 && numEvaluations >= maxEvaluations) { goto endLoop; }

            std::vector<Route> neighborRoutes = currentRoutes; // Copy for passing to the operator
            int op = opDist(rng); // Weighted operator selection
            //neighborhoods operator 
            if (op == 0) { removeRouteMove(neighborRoutes, rng, data); }
            else if (op == 1) { relocateCustomer(neighborRoutes, rng, data); }
            else if (op == 2) { insertionMove(neighborRoutes, rng, data); }
            else if (op == 3) { swapCustomers(neighborRoutes, rng, data); }
            else { twoOpt(neighborRoutes, rng, data); }
            numEvaluations++;
            if (isSolutionFeasible(neighborRoutes, data)) {
                auto currentObj = objectiveFunction(currentRoutes);
                auto neighborObj = objectiveFunction(neighborRoutes);
                double deltaVehicles = neighborObj.first - currentObj.first;
                double deltaDistance = neighborObj.second - currentObj.second;
                bool isNeighborBetter = BetterSolution(neighborRoutes, currentRoutes);

                double deltaForAcceptance;
                if (deltaVehicles < 0) { deltaForAcceptance = -1.0; }
                else if (deltaVehicles > 0) { deltaForAcceptance = 10000.0; }
                else { deltaForAcceptance = deltaDistance; }
                if (isNeighborBetter || (T > 1e-9 && std::exp(-deltaForAcceptance / T) > std::uniform_real_distribution<double>(0.0, 1.0)(rng))) {
                    currentRoutes = neighborRoutes; 
                    if (BetterSolution(currentRoutes, bestRoutes)) {
                        bestRoutes = currentRoutes;
                        improved = true; }
}
            } // else: neighbor solution infeasible


            // Prioritize distance in acceptance criterion
            // double deltaForAcceptance;
            // if (deltaDistance < 0) { 
            //     deltaForAcceptance = deltaDistance; // Negative distance difference is good
            // } else if (deltaDistance > 0) { 
            //     deltaForAcceptance = deltaDistance + deltaVehicles * 0.01; // Add a small penalty for vehicles
            // } else { 
            //     deltaForAcceptance = deltaVehicles; // If distances are equal, compare vehicles
            // }
            // if (isNeighborBetter || (T > 1e-9 && std::exp(-deltaForAcceptance / T) > std::uniform_real_distribution<double>(0.0, 1.0)(rng))) {
            //     currentRoutes = neighborRoutes; 
            //     if (BetterSolution(currentRoutes, bestRoutes)) {
            //         bestRoutes = currentRoutes;
            //         improved = true; 
            //     }}}

        } // End of inner for loop
        T *= coolingRate;
        //update the no improvement counter
        if (improved) {
            noImprovementCount = 0;
        } else {
            noImprovementCount++;}
        if (noImprovementCount >= ThresholdReheat) {
            T = T0 * reheatFactor;
            noImprovementCount = 0;}
        if (numEvaluations % (iterPerTemp * 20) == 0 || improved) {
            auto [currentBestVeh, currentBestDist] = objectiveFunction(bestRoutes);
            printf("Temp: %9.3f | BestDist: %8.2f | Veh: %.0f | NoImpr: %3d | Eval: %8d\n",
                   T, currentBestDist, currentBestVeh, noImprovementCount, numEvaluations);
        }} //end of main while loop
endLoop:
    auto endTime = std::chrono::steady_clock::now();
    double finalTime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    std::cout << "\n-----------------------------------------" << std::endl;
    std::cout << "Simulated Annealing Finished!" << std::endl;
    auto [finalVehicles, finalDistance] = objectiveFunction(bestRoutes);
    std::cout << "\nBest Solution Found:" << std::endl;
    std::cout << "Number of Vehicles: " << finalVehicles << std::endl;
    std::cout << "Total Distance: " << std::fixed << std::setprecision(2) << finalDistance << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    return bestRoutes;}

int main(int argc, char* argv[]) {
    // Check arguments
    if (argc < 4) { cerr << endl; return 1; }
    string instanceFile = argv[1];
    int maxTime = 0;
    int maxEvaluations = 0;
    try {
        maxTime = stoi(argv[2]);
        maxEvaluations = stoi(argv[3]);
        if (maxTime < 0 || maxEvaluations < 0) { cerr << "error:max_time/max_evaluations cannot be negative." << endl; return 1; }
    } catch (const std::exception& e) { cerr << "Error parsing args: " << e.what() << endl; return 1; }
    cout << "reading instance file: " << instanceFile << ".." << endl;
    ProblemData data = readInstance(instanceFile);
    if (data.customers.empty()) { cerr << "Error reading instance." << endl; return 1; }
    cout << "Instance read successfully. Customers: " << data.customers.size() -1  << endl;

    //cout << "\nConstructing initial solution.." << endl;
    vector<Route> initRoutes = constructInitialSolution(data);
    // Check initial solution feasibility
    if (!initRoutes.empty() && !isSolutionFeasible(initRoutes, data)) { cerr << " Warning:Initial solution infeasible!" << endl; }
    else if (initRoutes.empty() && data.customers.size() > 1) { cerr << "Warning: Initial solution construction failed/empty." << endl; }

    // double T0 = 551.51;
    // double coolingRate = 0.95;
    // int iterPerTemp = 50;
    // int maxNoImprovement = 969;
    // int ThresholdReheat = 41;
    // double reheatFactor = 0.12;
    double T0 = 1089.14;
    double coolingRate = 0.94;
    int iterPerTemp = 50;
    int maxNoImprovement = 241;
    int ThresholdReheat = 30;
    double reheatFactor = 0.36;

    auto mainStartTime = std::chrono::steady_clock::now();
    // run the algorithm 
    vector<Route> bestRoutes = simulatedAnnealing(data, initRoutes, T0, coolingRate,
                iterPerTemp, maxNoImprovement, maxTime, maxEvaluations,
                ThresholdReheat, reheatFactor);
    auto mainEndTime = std::chrono::steady_clock::now();
    double executionTime = std::chrono::duration_cast<std::chrono::duration<double>>(mainEndTime - mainStartTime).count();
    // Check final solution feasibility 
    if (bestRoutes.empty() && data.customers.size() > 1) { cout << "Final solution is empty." << endl; }
    else if (!bestRoutes.empty() && isSolutionFeasible(bestRoutes, data)) { cout << "Final solution is FEASIBLE." << endl; }
    else if (!bestRoutes.empty()) { cout << " Final solution is INFEASIBLE!" << endl; }
    int routeNumber = 1;
    for (size_t i = 0; i < bestRoutes.size(); ++i) {
        //print non empty routes
        if (!bestRoutes[i].customerIds.empty()) {
            cout << "Route " << routeNumber++ << ":";
            for (int custId : bestRoutes[i].customerIds) {
                cout << " " << custId; }
            cout << endl; }}
    if (bestRoutes.empty() && data.customers.size() > 1) { cout << "(No routes in final solution)" << endl; }
    auto [finalVehicles_check, finalDistance_check] = objectiveFunction(bestRoutes);
    cout << "Vehicles: " << finalVehicles_check << endl;
    cout << "Distance: " << fixed << setprecision(2) << finalDistance_check << endl;
    cout << "Total Execution Time: " << fixed << setprecision(2) << executionTime << " seconds" << endl;
    return 0;  }