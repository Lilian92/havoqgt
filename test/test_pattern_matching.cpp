#include <gtest/gtest.h>
#include <include/create_delegate_graph.hpp>
#include <include/havoqgt_setup.hpp>
#include <include/pattern_matching.hpp>
#include <include/input_graph_1.hpp>
#include <include/util.hpp>

namespace havoqgt { namespace test {

const std::string graph_unique_instance_name = "graph_obj";
const std::string edge_data_unique_instance_name = "graph_edge_data_obj";
const std::string input_graph_file_name = "/dev/shm/test_havoqgt_graph_6";
std::string pattern_input_filename = "../../../examples/prunejuice/rmat_log2_tree_pattern/12";

std::vector<std::tuple<uint64_t, uint64_t, EdgeData>> input_graph;
std::vector<std::tuple<uint64_t, uint64_t, EdgeData>> vec_global_edges; 
std::set<uint64_t> delegate_vertices;

void create_local_edge_list(graph_type& g, edge_data_t& edge_data_ptr, 
                            vloc_type vertex,  
                            std::vector<uint64_t>& edge_source, 
                            std::vector<uint64_t>& edge_target, 
                            std::vector<EdgeData>& edge_data) {
  for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); 
      ++eitr) {
    edge_source.push_back(g.locator_to_label(vertex));
    edge_target.push_back(g.locator_to_label(eitr.target()));
    //edge_data.push_back((uint64_t)eitr.edge_data());   
    edge_data.push_back((uint64_t)edge_data_ptr[eitr]); 
  }
}

bool test_pattern_matching() {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  const std::string output_filename = "./results/";
  const std::string vertex_metadata_filename = "/p/lustre1/an4/metadata/head";

  input_graph = grid_graph_sym_weighted_edges();

  //MPI_Barrier(MPI_COMM_WORLD);

  create_delegate_graph(input_graph, input_graph_file_name, 
                        graph_unique_instance_name, 
                        edge_data_unique_instance_name, mpi_rank);

  MPI_Barrier(MPI_COMM_WORLD);

  havoqgt::distributed_db ddb(havoqgt::db_open(), 
                              input_graph_file_name.c_str());
  graph_type *graph = ddb.get_segment_manager()->
    find<graph_type>(graph_unique_instance_name.c_str()).first;
  assert(graph != nullptr);

  edge_data_t* edge_data_ptr = ddb.get_segment_manager()->
    find<edge_data_t>(edge_data_unique_instance_name.c_str()).first;
  assert(edge_data_ptr != nullptr);

  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    std::cout << "MPI Rank: " << mpi_rank << " Graph Loaded Ready."
    << std::endl;
  }

  VertexMetadata vertex_metadata(*graph);
  generate_vertex_metadata(graph, vertex_metadata, vertex_metadata_filename);

  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    std::cout << "MPI Rank: " << mpi_rank << " Graph Metadata Generated."
    << std::endl;
  }

  if(mpi_rank == 0) {
    std::cout << "Pattern Matching Prunejuice ... " << std::endl;
  }

  size_t count_pj = pattern_matching_prunejuice(graph, vertex_metadata, edge_data_ptr,
          pattern_input_filename, output_filename, false, true);
  size_t count_seq = pattern_matching_seq(graph, vertex_metadata, edge_data_ptr,
          pattern_input_filename);

  return (count_pj == count_seq);
//  // build a local copy of the entire graph on each rank 
//  
//  std::vector<uint64_t> vec_local_edge_source;
//  std::vector<uint64_t> vec_global_edge_source;
//  std::vector<uint64_t> vec_local_edge_target;
//  std::vector<uint64_t> vec_global_edge_target;
//  std::vector<uint64_t> vec_local_edge_data;
//  std::vector<uint64_t> vec_global_edge_data;  
//
//  for (vitr_type vitr = g->vertices_begin(); vitr != g->vertices_end(); 
//       ++vitr) {
//    vloc_type vertex = *vitr;
//    if (vertex.is_delegate()) {
//      delegate_vertices.insert(g->locator_to_label(vertex));
//    }
//    create_local_edge_list(*g, *edge_data_ptr, vertex, vec_local_edge_source, 
//                           vec_local_edge_target, vec_local_edge_data); 
//  }
//
//  for (vitr_type vitr = g->delegate_vertices_begin(); 
//       vitr != g->delegate_vertices_end(); ++vitr) { 
//    vloc_type vertex = *vitr;
//    if (vertex.is_delegate()) {
//      delegate_vertices.insert(g->locator_to_label(vertex));
//    } 
//    create_local_edge_list(*g, *edge_data_ptr, vertex, vec_local_edge_source,
//                           vec_local_edge_target, vec_local_edge_data);
//  }
// 
//  MPI_Barrier(MPI_COMM_WORLD);
//
//
//  // gather global edges
//  mpi_all_gather(vec_local_edge_source, vec_global_edge_source, MPI_COMM_WORLD); 
//  mpi_all_gather(vec_local_edge_target, vec_global_edge_target, MPI_COMM_WORLD);
//  mpi_all_gather(vec_local_edge_data, vec_global_edge_data, MPI_COMM_WORLD);
//
//  MPI_Barrier(MPI_COMM_WORLD);
//
//  for (size_t i; i < vec_global_edge_source.size(); ++i) {
//    vec_global_edges.push_back(std::make_tuple(vec_global_edge_source[i], 
//                                               vec_global_edge_target[i], 
//                                               vec_global_edge_data[i]));
//  }
//
//  // sort edges
//  std::stable_sort(vec_global_edges.begin(),vec_global_edges.end(),
//    [](const std::tuple<uint64_t, uint64_t, uint64_t>& a,
//       const std::tuple<uint64_t, uint64_t, uint64_t>& b) -> bool {
//         return std::get<0>(a) < std::get<0>(b);
//       });
//  
//  for (size_t i = 0; i < grid_graph_weighted_edges_offset.size() - 1; ++i) {
//     size_t start = grid_graph_weighted_edges_offset[i];
//     size_t end = grid_graph_weighted_edges_offset[i+1];
//      
//     std::stable_sort(vec_global_edges.begin() + start, 
//                      vec_global_edges.begin() + end,
//       [](const std::tuple<uint64_t, uint64_t, uint64_t>& a,
//          const std::tuple<uint64_t, uint64_t, uint64_t>& b) -> bool {
//            return std::get<1>(a) < std::get<1>(b);
//          });    
//  } // sort edges  
//
//  MPI_Barrier(MPI_COMM_WORLD);
}

TEST(my_test, test_pattern_matching) {
  bool res = test_pattern_matching();
  MPI_Barrier(MPI_COMM_WORLD);
  EXPECT_TRUE(res);
}

}} //end namespace havoqgt::test

//mpi main for gteset
GTEST_API_ int main(int argc, char **argv) {
  // set up environment
  int mpi_rank(0), mpi_size(0), to_return;
  havoqgt::init(&argc, &argv);
  {
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  
  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    //havoqgt::get_environment().print();
    //print_system_info(false); 
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (argc > 1)
      havoqgt::test::pattern_input_filename = argv[1];
  // execute tests
  testing::InitGoogleTest(&argc, argv);
  to_return =  RUN_ALL_TESTS();
  }

  // delete the generated files

  ;

  return to_return;
}
