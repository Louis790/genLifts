Read("/home/veetoo/CLionProjects/masterproef/Scripts/GAP/Automorphisms.g");
Read("/home/veetoo/CLionProjects/masterproef/Scripts/GAP/Orbits.g");

lowerBound := 1;
upperBound := 100;
automorphismCount := 2000;

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

# Read("/home/veetoo/CLionProjects/masterproef/Scripts/GAP/SmallGroups.g");