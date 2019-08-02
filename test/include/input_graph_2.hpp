/**
 * A grig graph with three rows and five columns, 15 vertices and
 * 44 edges in total.
 * If a delegate partitiond graph is created with the parameter -d 4,
 * vertices 6, 7 and 8 are the delegates. The third value represents the
 * weight associated with the corresponding edge.
 * Same edges have same weight.
 */
std::vector<std::tuple<uint64_t, uint64_t, EdgeData>>
  grid_graph_sym_weighted_edges() {
  std::vector<std::tuple<uint64_t, uint64_t, EdgeData>> graph;
  graph.push_back(std::make_tuple(0, 1, 19));
  graph.push_back(std::make_tuple(0, 5, 3));
  graph.push_back(std::make_tuple(1, 0, 19));
  graph.push_back(std::make_tuple(1, 2, 10));
  graph.push_back(std::make_tuple(1, 6, 12));
  graph.push_back(std::make_tuple(2, 1, 10));
  graph.push_back(std::make_tuple(2, 3, 19));
  graph.push_back(std::make_tuple(2, 7, 20));
  graph.push_back(std::make_tuple(3, 2, 19));
  graph.push_back(std::make_tuple(3, 4, 1));
  graph.push_back(std::make_tuple(3, 8, 7));
  graph.push_back(std::make_tuple(4, 3, 1));
  graph.push_back(std::make_tuple(4, 9, 12));
  graph.push_back(std::make_tuple(5, 0, 3));
  graph.push_back(std::make_tuple(5, 6, 8));
  graph.push_back(std::make_tuple(5, 10, 5));
  graph.push_back(std::make_tuple(6, 1, 12));
  graph.push_back(std::make_tuple(6, 5, 8));
  graph.push_back(std::make_tuple(6, 7, 16));
  graph.push_back(std::make_tuple(6, 11, 3));
  graph.push_back(std::make_tuple(7, 2, 20));
  graph.push_back(std::make_tuple(7, 6, 16));
  graph.push_back(std::make_tuple(7, 8, 19));
  graph.push_back(std::make_tuple(7, 12, 15));
  graph.push_back(std::make_tuple(8, 3, 7));
  graph.push_back(std::make_tuple(8, 7, 19));
  graph.push_back(std::make_tuple(8, 9, 18));
  graph.push_back(std::make_tuple(8, 13, 1));
  graph.push_back(std::make_tuple(9, 4, 12));
  graph.push_back(std::make_tuple(9, 8, 18));
  graph.push_back(std::make_tuple(9, 14, 4));
  graph.push_back(std::make_tuple(10, 5, 5));
  graph.push_back(std::make_tuple(10, 11, 6));
  graph.push_back(std::make_tuple(11, 6, 3));
  graph.push_back(std::make_tuple(11, 10, 6));
  graph.push_back(std::make_tuple(11, 12, 4));
  graph.push_back(std::make_tuple(12, 7, 15));
  graph.push_back(std::make_tuple(12, 11, 4));
  graph.push_back(std::make_tuple(12, 13, 1));
  graph.push_back(std::make_tuple(13, 8, 1));
  graph.push_back(std::make_tuple(13, 12, 1));
  graph.push_back(std::make_tuple(13, 14, 8));
  graph.push_back(std::make_tuple(14, 9, 4));
  graph.push_back(std::make_tuple(14, 13, 8));
  return graph;
}

std::vector<uint64_t> grid_graph_weighted_edges_degree = 
{2, 3, 3, 3, 2,
 3, 4, 4, 4, 3,
 2, 3, 3, 3, 2
};

std::vector<uint64_t> grid_graph_weighted_edges_offset =
{0, 2, 5, 8, 11,
 13, 16, 20, 24, 28,
 31, 33, 36, 39, 42, 44
};

void print_adjacency_list(
  std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>& graph, 
  std::vector<uint64_t>& offset) {

  for (size_t i = 0; i < offset.size() - 1; ++i) { 
    std::cout << i << ": ";
    for(size_t j = offset[i]; j < offset[i+1] ; ++j) {
      std::cout << std::get<1>(graph[j]) << ", ";   
    }  
    std::cout << std::endl;      
  }
}
