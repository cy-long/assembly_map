# assembly_map
This code randomly generates a ecological community represented by its interaaction matrix.

For the dynamical simulation version, it checks the feasibility and stability of each community, measuring the relative size of its feasibility domain. After that, it determines whether a given sub-community can be assembled from its smaller sub-communities, with the adding of a new species. Finally it illustrates all the possible transitons as described before.

For the structural analysis version, it calculates the Omega value of each sub-community, and the Omega_overlap value with their successions with possible transitions. The ratio of these two measurements are regarded a prediction of the transition possibility. This can have a similar illustration of the whole assembly process.

By this kind of backward induction, this code can generate a map for all the potential assembly paths wich leads to the final state observed. The dynamical part is currently working for only 3 species cases while the structural part can be extended to higher and arbitrary species numbers.
