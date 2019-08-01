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
std::string pattern_input_filename = "../../../examples/prunejuice/rmat_log2_tree_pattern/13";

std::vector<std::tuple<uint64_t, uint64_t, EdgeData>> input_graph;
std::vector<std::tuple<uint64_t, uint64_t, EdgeData>> vec_global_edges; 
std::set<uint64_t> delegate_vertices;

bool test_pattern_matching() {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  const std::string output_filename = "./results/";
  const std::string vertex_metadata_filenames_head = "/p/lustre1/an4/metadata/head";
  const std::string vertex_metadata_filename = "/p/lustre1/an4/metadata/vertex_metadata_1";

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
  generate_vertex_metadata(graph, vertex_metadata, vertex_metadata_filenames_head);

  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    std::cout << "MPI Rank: " << mpi_rank << " Graph Metadata Generated."
    << std::endl;
  }

  if(mpi_rank == 0) {
    std::cout << "Pattern Matching Prunejuice ... " << std::endl;
  }

  OUTPUT_SUIT res_pj;
  size_t count_pj = pattern_matching_prunejuice(graph, vertex_metadata, edge_data_ptr,
          pattern_input_filename, output_filename, res_pj, false, true);
  res_pj.output("prunejuice");

  OUTPUT_SUIT res_seq;
  size_t count_seq = pattern_matching_seq(input_graph, vertex_metadata_filename, pattern_input_filename, res_seq);
  res_seq.output("seq");
  return (count_pj == count_seq);
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
