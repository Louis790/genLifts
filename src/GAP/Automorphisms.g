PrintAutomorphismMappings := function(G, maxSize)
    local autGroup, elements, i, j, mapping, element_index, aut, automorphisms, autCount;

    autGroup := AutomorphismGroup(G);
    elements := Elements(G);
    element_index := [];

    for i in [1..Length(elements)] do
        element_index[Position(elements, elements[i])] := i;
    od;

    autCount := Size(autGroup);
    # Print("Automorphisms: ", autCount, "\n");

    automorphisms := [];

    if autCount <= maxSize then
        automorphisms := AllAutomorphisms(G);
    else
        while Length(automorphisms) < maxSize do
            aut := Random(autGroup);
            if not aut in automorphisms then
                Add(automorphisms, aut);
            fi;
        od;
    fi;

    for i in [1..Length(automorphisms)] do
        aut := automorphisms[i];
        mapping := [];

        # Print("Automorphism ", i, ":\n");
        for j in [1..Length(elements)] do
            mapping[j] := Position(elements, Image(aut, elements[j]));
        od;
        Print(mapping, "\n");
    od;
end;

# G := SmallGroup(10, 1);
# PrintAutomorphismMappings(G);