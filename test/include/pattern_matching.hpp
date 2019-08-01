#include <include/havoqgt_setup.hpp>

#include <bitset>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <metadata/vertex_data_db.hpp>
#include <metadata/vertex_data_db_degree.hpp>

#include <prunejuice/template.hpp>
#include <prunejuice/non_local_constraint.hpp>
#include <prunejuice/temporal_constraint.hpp>
#include <prunejuice/algorithm_state.hpp>
#include <prunejuice/local_constraint_checking.hpp>
#include <prunejuice/non_local_constraint_checking_unique.hpp>
#include <prunejuice/non_local_constraint_checking_tds_batch.hpp> // TODO: optimize for batch_size = mpi_size

#include <math.h>

#define TP_ASYNC

template<typename T>
using SegmentAllocator = bip::allocator<T, segment_manager_t>;

template<typename T>
using DelegateGraphVertexDataSTDAllocator = graph_type::vertex_data
<T, std::allocator<T>>; 

template<typename T>
using DelegateGraphEdgeDataSTDAllocator = graph_type::edge_data
<T, std::allocator<T>>;  

typedef short EdgeData; 
typedef graph_type::edge_data<EdgeData, 
        bip::allocator<EdgeData, segment_manager_t>> edge_data_t;
typedef uint64_t Vertex;
typedef uint64_t Edge;
typedef uint64_t VertexData; // for string hash
typedef uint64_t VertexRankType; // TODO: delete
static constexpr size_t max_bit_vector_size = 16; // TODO:
static constexpr size_t max_template_vertex_count = 16;  
typedef std::bitset<max_bit_vector_size> BitSet; // TODO: rename to TemplateVertexSet
typedef BitSet TemplateVertexBitSet; 
typedef uint16_t TemplateVertexType; // TODO: rename to TemplateVertexBitSetToUint
typedef uint8_t Boolean; // TODO: replace all bool with Boolean?

typedef graph_type::vertex_data<VertexData, std::allocator<VertexData> > VertexMetadata;
typedef graph_type::vertex_data<Boolean, std::allocator<Boolean> > VertexActive; // TODO: solution_graph // TODO: you are mixing bool and uint!
typedef graph_type::vertex_data<TemplateVertexType, std::allocator<TemplateVertexType> > TemplateVertex; // TODO: solution_graph, rename to VertexTemplateVertexBitSetToUint
typedef std::unordered_map<Vertex, std::tuple<EdgeData, EdgeData>> VertexMinMaxMap;
typedef graph_type::vertex_data<VertexMinMaxMap, std::allocator<VertexMinMaxMap> > VertexMinMax;


typedef graph_type::vertex_data<uint64_t, std::allocator<uint64_t> > VertexIteration; // TODO: delete
typedef graph_type::vertex_data<VertexRankType, std::allocator<VertexRankType> > VertexRank; // TODO: delete

typedef prunejuice::vertex_state<Vertex, VertexData, BitSet, VertexMinMaxMap> VertexState;
typedef std::unordered_map<Vertex, VertexState> VertexStateMap; // TODO: solution_graph

typedef std::unordered_set<Vertex> VertexSet;  
typedef graph_type::vertex_data<VertexSet, std::allocator<VertexSet> > VertexSetCollection; 
typedef std::unordered_set<std::pair<Vertex, std::vector<EdgeData>>, boost::hash<std::pair<Vertex, std::vector<EdgeData>>> > VertexEdgeDataVectorSet;  
typedef graph_type::vertex_data<VertexEdgeDataVectorSet, std::allocator<VertexEdgeDataVectorSet> > VertexEdgeDataVectorSetCollection; 

typedef std::unordered_map<Vertex, uint8_t> VertexUint8Map; 
typedef graph_type::vertex_data<VertexUint8Map, std::allocator<VertexUint8Map> > VertexUint8MapCollection;    
typedef std::unordered_map<Vertex, std::pair<uint8_t, EdgeData>> VertexUint8EdgeDataMap; 
typedef graph_type::vertex_data<VertexUint8EdgeDataMap, std::allocator<VertexUint8EdgeDataMap> > VertexUint8EdgeDataMapCollection;    


typedef graph_type::edge_data<EdgeData, std::allocator<EdgeData> > EdgeMetadata;
typedef graph_type::edge_data<Boolean, std::allocator<Boolean> > EdgeActive; // TODO: solution_graph 

typedef std::vector<Boolean> VectorBoolean;
typedef prunejuice::pattern_graph_csr<Vertex, Edge, VertexData, 
        EdgeData> PatternGraph;
typedef pattern_nonlocal_constraint<Vertex, Edge, VertexData, EdgeData, PatternGraph>
PatternNonlocalConstraint;

typedef pattern_temporal_constraint<Vertex, Edge, VertexData, EdgeData, PatternGraph, PatternNonlocalConstraint, BitSet, VertexMinMaxMap>
PatternTemporalConstraint;
typedef struct output_suit{
    size_t patterns_found;
    size_t active_vertices_count;
    size_t active_edges_count;

    output_suit() : patterns_found(0),
    active_vertices_count(0),
    active_edges_count(0) {}

    void output(std::string program_name = std::string("")) {
        if (program_name != "")
            std::cout << program_name << " ";
        std::cout << "running results:" << std::endl;
        std::cout << "patterns_found: " << patterns_found << std::endl;
        std::cout << "active_vertices_count: " << active_vertices_count << std::endl;
        std::cout << "active_edges_count: " << active_edges_count << std::endl;
    }
}OUTPUT_SUIT;

void generate_vertex_metadata(graph_type * graph, VertexMetadata & vertex_metadata, std::string vertex_metadata_input = std::string()) {
  if (vertex_metadata_input.size() > 0) {
    vertex_data_db_nostdfs<graph_type, VertexMetadata, Vertex, VertexData>
      (graph, vertex_metadata, vertex_metadata_input, 10000);
  } else {
    vertex_data_db_degree<graph_type, VertexMetadata, Vertex, VertexData>
      (graph, vertex_metadata);
  }
}

size_t pattern_matching_prunejuice(graph_type * graph,
        VertexMetadata & vertex_metadata,
        edge_data_t * edge_data_ptr,
        const std::string pattern_input,
        const std::string result_output,
        OUTPUT_SUIT & res,
        bool enable_edge_matching,
        bool  enable_temporal_edge_matching,
        uint64_t tp_vertex_batch_size);

size_t pattern_matching_seq(const std::vector<std::tuple<uint64_t, uint64_t, EdgeData>> edge_list,
        const std::string vertex_metadata_filename,
        const std::string pattern_input_filename,
        OUTPUT_SUIT & res,
        bool enable_edge_matching,
        bool enable_temporal_edge_matching);

#include"pattern_matching_pj.hpp"
#include"pattern_matching_seq.hpp"
