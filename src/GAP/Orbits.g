PrintElementOrbits := function(G)
    local autGroup, elements, orbits, orbit, sortedOrbit;

    autGroup := AutomorphismGroup(G);
    elements := Elements(G);

    # Calculate the orbits of the elements under the action of the automorphism group
    orbits := Orbits(autGroup, elements, OnPoints);

    for orbit in orbits do
        sortedOrbit := ShallowCopy(orbit); # Original is immutable
        Sort(sortedOrbit);
        Print(List(sortedOrbit, x -> Position(elements, x)), "\n");
    od;
end;
