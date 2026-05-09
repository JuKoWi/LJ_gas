- restructure code:
    - positions, velocity, forces, masses in flat list accesses only by index
    - construct list of atom types upon initialization
    - nested neighbor list
    - bonds as list of bond pairs, same for dihedrals
    - charges full pairs, later ewald
- parsing: 
    - read job, initialize system with job parameters
    - create atom type list (vector)
    - directly create system from geom.json and atom types

- calculate pressure
- solve two body problem with LJ-potential analytically
- write max force to file
- check bottleneck
- optimize force computation

What I should know:
- How to define struct and set the elements
- how to do for loop
- how to loop over elements of a vector


