


graph(
    Edge(state1 = "mu0",state2 = "mu0", penalty = 0, K = 3),
    Edge(state1 = "mu0",state2 = "Coll", penalty = 10, type = "std"),
    Edge(state1 = "Coll",state2 = "Coll", penalty = 0),
    Edge(state1 = "Coll",state2 = "mu0", penalty = 0, type = "std", K = 3),
    StartEnd(start = "mu0", end = c("mu0", "Coll")),
    Node(state = "mu0", min = 0, max = 0))


explore(graphReorder(mygraph))




######### ERROR

mygraph <- graphReorder(
  graph(
    Edge(1,1),
    Edge(0,1, decay = 0.5),
    Edge(2,2),
    Edge(2,1, decay = 0.3),
    StartEnd(start = 0, end = 1)
  ))

explore(mygraph)



mygraph <- graphReorder(
  graph(
    Edge(0,1, decay = 0.5),
    Edge(1,2),
    Edge(2,1, decay = 0.3),
    StartEnd(start = 0, end = 1)
  ))

explore(mygraph)



mygraph <- graphReorder(
  graph(
    Edge(0,0),
    Edge(0,1, decay = 0.5),
    Edge(1,2),
    StartEnd(start = 0, end = 1)
  ))

explore(mygraph)

######### OK

mygraph <- graphReorder(
  graph(
    Edge(0,1, type = "up", gap = 0.5),
    Edge(1,0, type ="down", gap = 0.5),
    Edge(0,0),
    Edge(1,1),
    StartEnd(start = 0, end = 0)
  ))

explore(mygraph)



