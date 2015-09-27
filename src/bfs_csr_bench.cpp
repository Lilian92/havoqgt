#include <stdint.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>      // std::stringstream
#include <iomanip>      // std::setfill, std::setw
#include <random>
#include <chrono>

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/segment_manager.hpp>

#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/parallel_edge_list_reader.hpp>
#include <havoqgt/graphstore/csr/csr_graph.hpp>
#include <havoqgt/graphstore/graph_traversal/bfs.hpp>
#include <havoqgt/graphstore/graphstore_utilities.hpp>

enum {
  kNumBFSLoop = 64
};


using mapped_file_type      = boost::interprocess::managed_mapped_file;
using segment_manager_type  = boost::interprocess::managed_mapped_file::segment_manager;

using index_type = uint64_t;
using vertex_type = uint64_t;
using graph_type = csr_graph_struct::csr_graph<havoqgt::parallel_edge_list_reader, index_type, vertex_type, segment_manager_type>;

template <typename gen_type, typename rnd_type>
void run_bfs(graph_type& graph, size_t max_vertex_id, size_t num_edges, gen_type& gen, rnd_type& dis)
{
  std::cout << "\n--- BFS ---" << std::endl;

  std::cout << "max_vertex_id:\t" << max_vertex_id << std::endl;
  std::cout << "num_edges:\t"     << num_edges     << std::endl;

#if USE_SRC_CANDIDATE_TBL
  bool* table = new bool[graph->num_vertices()];
  std::cout << "Allocated src_candidats_table: " << sizeof(bool) * graph->num_vertices() / (1ULL << 30) << " GB" << std::endl;
  find_startvertex_candidates(graph, table);
#endif

  for (int i = 0; i < kNumBFSLoop; ++i) {
    vertex_type src = dis(gen);
    std::cout << "BFS[" << i << "]: src=\t" << src << std::endl;

    graphstore::utility::print_time();
    bfs_sync<graph_type, vertex_type, false>(graph, src, max_vertex_id, num_edges);
    std::cout << "finish: ";
    graphstore::utility::print_time();

    std::cout << "\n" << std::endl;
  }
  std::cout << "BFS done." << std::endl;

}

std::string fname_graph_;
std::string fname_segmentfile_;
size_t segment_size_log2_ = 30;
std::vector<std::string> fname_edge_list_;
size_t max_vertex_id_ = 0;
size_t num_edges_ = 0;

void parse_options(int argc, char **argv)
{
  char c;
  while ((c = getopt (argc, argv, "g:s:f:e:v:m:")) != -1) {
    switch (c) {
      case 'g':
        fname_graph_ = optarg;
        break;
      case 'f':
        fname_segmentfile_ = optarg;
        break;
      case 's':
        segment_size_log2_ = boost::lexical_cast<size_t>(optarg);
        break;
      case 'v':
        max_vertex_id_ = boost::lexical_cast<size_t>(optarg);
        break;
      case 'm':
        num_edges_ = boost::lexical_cast<size_t>(optarg);
        break;
      case 'e':
        std::string fname(optarg);
        std::ifstream fin(fname);
        std::string line;
        if (!fin.is_open()) {
          std::cerr << fname << std::endl;
          HAVOQGT_ERROR_MSG("Unable to open a file");
        }
        while (std::getline(fin, line)) {
          fname_edge_list_.push_back(line);
        }
        break;
    }
  }
}


int main(int argc, char* argv[])
{
  graph_type* graph;

  mapped_file_type mapped_file = mapped_file_type(
                                   boost::interprocess::create_only,
                                   fname_segmentfile_.c_str(),
                                   1ULL << 30);
  segment_manager_type* segment_manager = mapped_file.get_segment_manager();

  if (!fname_graph_.empty()) {
    std::cout << "--- Constructing a graph from file ---" << std::endl;
    graph = new graph_type(fname_graph_, segment_manager);
    max_vertex_id_ = graph->num_vertices() - 1;
    num_edges_ = graph->num_edges();

  } else {
    std::cout << "\n--- Initializing edgelist ---" << std::endl;
    graphstore::utility::print_time();
    havoqgt::havoqgt_init(&argc, &argv);
    {
      int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
      int mpi_size = havoqgt::havoqgt_env()->world_comm().size();
      havoqgt::get_environment();
      if (mpi_rank == 0) {
        havoqgt::parallel_edge_list_reader edge_list(fname_edge_list_);
        std::cout << "\n--- Constructing csr graph ---" << std::endl;
        graphstore::utility::print_time();
        graph = new graph_type(edge_list, max_vertex_id_, num_edges_, segment_manager);
      }
    }

  }

  /// ---------- Graph Traversal --------------- ///

  // obtain a seed from the system clock:
  //  std::random_device rd; /// not implemented in gccc ?
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::cout << " use seed: " << seed << std::endl;
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<vertex_type> dis(0, max_vertex_id_);

  run_bfs(*graph, max_vertex_id_, num_edges_, gen, dis);

  delete graph;
  delete segment_manager;
}