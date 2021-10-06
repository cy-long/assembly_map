# assembly_map
This code randomly generates the interaction matrix of a community and all its sub-communities. It checks the feasibility and stability of each community, measuring the relative size of its feasibility domain. After that, it determines whether a given sub-community can be assembled from its smaller sub-communities, with the adding of a new species. Finally it illustrates all the possible transitons as described before.

By this kind of backward induction, this code can generate a map for all the potential assembly paths wich leads to the final state observed. It is currently working for only 3 species cases but can be extended to higher species numbers.
