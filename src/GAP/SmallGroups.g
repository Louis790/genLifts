Read("./src/GAP/Automorphisms.g");
Read("./src/GAP/Orbits.g");

lowerBound := 1;
upperBound := 50;
automorphismCount := 1000;

for order in [lowerBound..upperBound] do
  numGroups := NumberSmallGroups(order);

  for i in [1..numGroups] do
    G := SmallGroup(order, i);
    M := MultiplicationTable(G);

    Print("Order: ", order, ", Number: ", i, "\n");
    Print(M);
    Print("\nA\n");
    PrintAutomorphismMappings(G, automorphismCount);
    Print("\nO\n");
    PrintElementOrbits(G);
    Print("####");
    Print("\n");
  od;
od;

quit;