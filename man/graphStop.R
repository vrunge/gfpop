######### ERROR

mygraph <- graphAnalysis(
  graph(
    Edge(1,1),
    Edge(0,1, decay = 0.5),
    Edge(2,2),
    Edge(2,1, decay = 0.3),
    StartEnd(start = 0, end = 1)
  ))

explore(mygraph)



mygraph <- graphAnalysis(
  graph(
    Edge(0,1, decay = 0.5),
    Edge(1,2),
    Edge(2,1, decay = 0.3),
    StartEnd(start = 0, end = 1)
  ))

explore(mygraph)



mygraph <- graphAnalysis(
  graph(
    Edge(0,0),
    Edge(0,1, decay = 0.5),
    Edge(1,2),
    StartEnd(start = 0, end = 1)
  ))

explore(mygraph)

######### OK

mygraph <- graphAnalysis(
  graph(
    Edge(0,1, type = "up", gap = 0.5),
    Edge(1,0, type ="down", gap = 0.5),
    Edge(0,0),
    Edge(1,1),
    StartEnd(start = 0, end = 0)
  ))

explore(mygraph)



