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
  graph.push_back(std::make_tuple(0, 1, 1));
  graph.push_back(std::make_tuple(0, 5, 4));
  graph.push_back(std::make_tuple(1, 0, 1));
  graph.push_back(std::make_tuple(1, 2, 6));
  graph.push_back(std::make_tuple(1, 6, 6));
  graph.push_back(std::make_tuple(2, 1, 6));
  graph.push_back(std::make_tuple(2, 3, 8));
  graph.push_back(std::make_tuple(2, 7, 4));
  graph.push_back(std::make_tuple(3, 2, 8));
  graph.push_back(std::make_tuple(3, 4, 8));
  graph.push_back(std::make_tuple(3, 8, 8));
  graph.push_back(std::make_tuple(4, 3, 8));
  graph.push_back(std::make_tuple(4, 9, 6));
  graph.push_back(std::make_tuple(5, 0, 4));
  graph.push_back(std::make_tuple(5, 6, 8));
  graph.push_back(std::make_tuple(5, 10, 2));
  graph.push_back(std::make_tuple(6, 1, 6));
  graph.push_back(std::make_tuple(6, 5, 8));
  graph.push_back(std::make_tuple(6, 7, 4));
  graph.push_back(std::make_tuple(6, 11, 6));
  graph.push_back(std::make_tuple(7, 2, 4));
  graph.push_back(std::make_tuple(7, 6, 4));
  graph.push_back(std::make_tuple(7, 8, 7));
  graph.push_back(std::make_tuple(7, 12, 7));
  graph.push_back(std::make_tuple(8, 3, 8));
  graph.push_back(std::make_tuple(8, 7, 7));
  graph.push_back(std::make_tuple(8, 9, 9));
  graph.push_back(std::make_tuple(8, 13, 6));
  graph.push_back(std::make_tuple(9, 4, 6));
  graph.push_back(std::make_tuple(9, 8, 9));
  graph.push_back(std::make_tuple(9, 14, 6));
  graph.push_back(std::make_tuple(10, 5, 2));
  graph.push_back(std::make_tuple(10, 11, 3));
  graph.push_back(std::make_tuple(11, 6, 6));
  graph.push_back(std::make_tuple(11, 10, 3));
  graph.push_back(std::make_tuple(11, 12, 9));
  graph.push_back(std::make_tuple(12, 7, 7));
  graph.push_back(std::make_tuple(12, 11, 9));
  graph.push_back(std::make_tuple(12, 13, 8));
  graph.push_back(std::make_tuple(13, 8, 6));
  graph.push_back(std::make_tuple(13, 12, 8));
  graph.push_back(std::make_tuple(13, 14, 8));
  graph.push_back(std::make_tuple(14, 9, 6));
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
