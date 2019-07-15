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
typedef struct {
    size_t patterns_found;
    size_t active_vertices_count;
    size_t active_edges_count;

    void output(std::string program_name = std::string("")) {
        if (program_name != "")
            std::cout << program_name << " ";
        std::cout << "running results:" << std::endl;
        std::cout << "patterns_found: " << patterns_found << std::endl;
        std::cout << "active_vertices_count: " << active_vertices_count << std::endl;
        std::cout << "active_edges_count: " << active_edges_count << std::endl;
    }
}OUTPUT_SUIT;

size_t pattern_matching_prunejuice(graph_type * graph,
        VertexMetadata & vertex_metadata,
        edge_data_t * edge_data_ptr,
        const std::string pattern_input,
        const std::string result_output,
        OUTPUT_SUIT & res,
        bool enable_edge_matching = false,
        bool  enable_temporal_edge_matching = false,
        uint64_t tp_vertex_batch_size = 0) {

    int mpi_rank(0), mpi_size(0);

    {

        CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
        CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

        if (mpi_rank == 0) {
            std::cout << "MPI Initialized With " << mpi_size << " Ranks." << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);  

        /////////////////////////////////////////////////////////////////////////////

        // parse command line

        tp_vertex_batch_size = comm_world().size();

        std::string pattern_dir = pattern_input; 
        std::string result_dir = result_output;

        MPI_Barrier(MPI_COMM_WORLD);

        /////////////////////////////////////////////////////////////////////////////

        MPI_Barrier(MPI_COMM_WORLD);
        if (mpi_rank == 0) {
            std::cout << "Done Loading Graph." << std::endl;
        }

        /////////////////////////////////////////////////////////////////////////////

        // pattern matching
        {

            ////////////////////////////////////////////////////////////////////////////// 

            if(mpi_rank == 0) {
                std::cout << "Pattern Matching ... " << std::endl;
            }

            //////////////////////////////////////////////////////////////////////////////

            double time_start = MPI_Wtime();
            double time_end = MPI_Wtime();

            // per rank containers 
            VertexStateMap vertex_state_map; 

            // vertex containers 

            VertexActive vertex_active(*graph);
            TemplateVertex template_vertices(*graph);
            VertexMinMax vertexminmax_vertices(*graph);
            VertexUint8EdgeDataMapCollection vertex_active_edges_map(*graph);
            VertexSetCollection vertex_token_source_set(*graph); // per vertex set
            VertexEdgeDataVectorSetCollection token_source_edge_data_vertices(*graph);

            uint8_t vertex_rank; // TODO: dummy
            uint8_t vertex_iteration; // TODO: dummy  

            if(mpi_rank == 0) {
                std::cout << "Pattern Matching | Allocated Vertex and Edge Containers" 
                    << std::endl;
            }

            ///////////////////////////////////////////////////////////////////////////// 

            bool do_output_vertex_data = true; // TODO: ?

            size_t token_passing_algo = 0; // TODO: ? 

            MPI_Barrier(MPI_COMM_WORLD);

            ///////////////////////////////////////////////////////////////////////////// 

            time_start = MPI_Wtime();

            MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this?
            time_end = MPI_Wtime();
            if(mpi_rank == 0) {
                std::cout << "Pattern Matching Time | Vertex Data DB : " 
                    << time_end - time_start << std::endl;
            }

            if (do_output_vertex_data) {
                std::string vertex_data_filename = result_dir +
                    "/all_ranks_vertex_data/vertex_data_" + std::to_string(mpi_rank);
                std::ofstream vertex_data_file(vertex_data_filename, std::ofstream::out);

                for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
                        ++vitr) {
                    vloc_type vertex = *vitr;
                    vertex_data_file << mpi_rank << ", l, " << graph->locator_to_label(vertex) 
                        << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";  
                    std::cout << mpi_rank << ", l, " << graph->locator_to_label(vertex) 
                        << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";  
                } 	

                for(vitr_type vitr = graph->delegate_vertices_begin();
                        vitr != graph->delegate_vertices_end(); ++vitr) {
                    vloc_type vertex = *vitr;

                    if (vertex.is_delegate() && (graph->master(vertex) == mpi_rank)) {
                        vertex_data_file << mpi_rank << ", c, " << graph->locator_to_label(vertex) 
                            << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";
                        std::cout << mpi_rank << ", c, " << graph->locator_to_label(vertex) 
                            << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";
 
                    } else {	
                        vertex_data_file << mpi_rank << ", d, " << graph->locator_to_label(vertex) 
                            << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";  
                        std::cout << mpi_rank << ", d, " << graph->locator_to_label(vertex) 
                            << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";  
 
                    }
                }

                vertex_data_file.close();
            }

            ///////////////////////////////////////////////////////////////////////////// 

            // result
            std::string pattern_set_result_filename = result_dir + "/result_pattern_set";
            std::ofstream pattern_set_result_file;  
            if (mpi_rank == 0) {
                pattern_set_result_file = std::ofstream(pattern_set_result_filename, std::ofstream::out);
            }
            /////////////////////////////////////////////////////////////////////////////
            // TODO: setup pattern set 
            // a pattern set is a collection of directories containing pattern files 

            // TODO: code indentation

            // loop over pattern set
            for (size_t ps = 0; ps < 1; ps++) { // TODO: for now, only reading from pattern_dir/0 
                // beginning of the pattern set

                // setup pattern to search
                if(mpi_rank == 0) { 
                    std::cout << "Setting up Pattern [" << ps << "] ... " << std::endl;
                }

                // setup pattern - for local constraint checking 
                std::string pattern_input_filename = pattern_dir + "/pattern";

                PatternGraph pattern_graph(
                        pattern_input_filename + "_edge",
                        pattern_input_filename + "_vertex",
                        pattern_input_filename + "_vertex_data",
                        pattern_input_filename + "_edge_data",
                        pattern_input_filename + "_stat",
                        false, false); // TODO: improve

                PatternNonlocalConstraint ptrn_util_two(pattern_graph,
                        pattern_dir + "/pattern_nonlocal_constraint");

                PatternTemporalConstraint ptrn_temp_const(pattern_graph,
                        ptrn_util_two,
                        pattern_dir + "/pattern_temporal_constraint", enable_temporal_edge_matching);

                vertex_state_map.clear(); // important
                vertex_active.reset(true); // initially all vertices are active
                vertex_active_edges_map.clear(); // important
                vertexminmax_vertices.clear();
                vertex_token_source_set.clear(); // clear all the sets on all the vertices
                token_source_edge_data_vertices.clear();

                bool global_init_step = true; // TODO: Boolean 
                bool global_not_finished = false; // TODO: Boolean

                bool do_nonlocal_constraint_checking = true; // TODO: Boolean

                uint64_t global_itr_count = 0;
                uint64_t active_vertices_count = 0;
                uint64_t active_edges_count = 0;
                uint64_t message_count = 0; 

                // result
                std::string itr_result_filename = result_dir +
                    "/result_iteration"; 
                std::ofstream itr_result_file(itr_result_filename, std::ofstream::out);

                std::string step_result_filename = result_dir + 
                    "/result_step";
                std::ofstream step_result_file(step_result_filename, std::ofstream::out); 

                std::string superstep_result_filename = result_dir + 
                    "/result_superstep";
                std::ofstream superstep_result_file(superstep_result_filename, std::ofstream::out);

                std::string active_vertices_count_result_filename = result_dir +
                    "/all_ranks_active_vertices_count/active_vertices_" + std::to_string(mpi_rank); 
                std::ofstream active_vertices_count_result_file(active_vertices_count_result_filename, std::ofstream::out);

                std::string active_vertices_result_filename = result_dir + 
                    "/all_ranks_active_vertices/active_vertices_" + std::to_string(mpi_rank);
                std::ofstream active_vertices_result_file(active_vertices_result_filename, std::ofstream::out);

                std::string active_edges_count_result_filename = result_dir + 
                    "/all_ranks_active_edges_count/active_edges_" + std::to_string(mpi_rank);
                std::ofstream active_edges_count_result_file(active_edges_count_result_filename, std::ofstream::out);

                std::string active_edges_result_filename = result_dir + 
                    "/all_ranks_active_edges/active_edges_" + std::to_string(mpi_rank);
                std::ofstream active_edges_result_file(active_edges_result_filename, std::ofstream::out);

                std::string message_count_result_filename = result_dir + 
                    "/all_ranks_messages/messages_" + std::to_string(mpi_rank);
                std::ofstream message_count_result_file(message_count_result_filename, std::ofstream::out);

                MPI_Barrier(MPI_COMM_WORLD);

                double pattern_time_start = MPI_Wtime();

                /////////////////////////////////////////////////////////////////////////////

                if (mpi_rank == 0) {
                    std::cout << "Running Constraint Checking ..." << std::endl;
                }

                global_not_finished = false;

                double itr_time_start = MPI_Wtime();

                /////////////////////////////////////////////////////////////////////////////

                // mark inactive edges

                /////////////////////////////////////////////////////////////////////////////
                double label_propagation_time_start = MPI_Wtime();

#ifdef OUTPUT_EDGEDATA
                for(auto vertex = graph->vertices_begin(); vertex != graph->vertices_end(); vertex++)
                {
                    for(auto eitr = graph->edges_begin(*vertex); eitr != graph->edges_end(*vertex);
                            ++eitr) {
                        std::cout << (graph->locator_to_label(*vertex)) << " ";
                        std::cout << (graph->locator_to_label(eitr.target())) << " ";
                        std::cout << "edge data" << (uint64_t)(*edge_data_ptr)[eitr] << std::endl;
                    }
                }
#endif

                prunejuice::label_propagation_pattern_matching_bsp<Vertex, VertexData, EdgeData, edge_data_t,
                    graph_type, VertexMetadata, VertexStateMap, VertexActive, 
                    VertexUint8EdgeDataMapCollection, TemplateVertexBitSet, TemplateVertex, VertexMinMaxMap, VertexMinMax, PatternGraph, PatternTemporalConstraint>
                        (graph, *edge_data_ptr, enable_edge_matching, enable_temporal_edge_matching, vertex_metadata, vertex_state_map, vertex_active, 
                         vertex_active_edges_map, template_vertices, vertexminmax_vertices, pattern_graph, ptrn_temp_const, global_init_step, 
                         global_not_finished, global_itr_count, superstep_result_file, 
                         active_vertices_count_result_file, active_edges_count_result_file,
                         message_count_result_file);


                MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
                double label_propagation_time_end = MPI_Wtime();
                if(mpi_rank == 0) {
                    std::cout << "Pattern Matching Time | Local Constraint Checking : " 
                        << label_propagation_time_end - label_propagation_time_start << std::endl;
                }

                // result
                if(mpi_rank == 0) {
                    step_result_file << global_itr_count << ", LP, "
                        << (label_propagation_time_end - label_propagation_time_start) << "\n";
                }

                /////////////////////////////////////////////////////////////////////////////

                if (global_init_step) { // Important
                    global_init_step = false;
                } 

                global_not_finished = havoqgt::mpi_all_reduce(global_not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here 

                if(mpi_rank == 0) {
                    std::cout << "Pattern Matching | Global Finished Status : "; 
                    if (global_not_finished) { 
                        std::cout << "Continue" << std::endl;
                    } else {
                        std::cout << "Stop" << std::endl;
                    } 
                }

                bool global_active_vertex = true;
                if (vertex_state_map.size() < 1) {
                    global_active_vertex = false;
                }

                global_active_vertex = havoqgt::mpi_all_reduce(global_active_vertex, std::greater<uint8_t>(), MPI_COMM_WORLD); 
                MPI_Barrier(MPI_COMM_WORLD);

                global_not_finished = (global_not_finished & global_active_vertex);

                if(mpi_rank == 0) {
                    std::cout << "Pattern Matching | Global Active Vertex Status : ";
                    if (global_active_vertex) {
                        std::cout << "Active vertices left." << std::endl;
                    } else {
                        std::cout << "No active vertex left." << std::endl;
                    }
                }

                /////////////////////////////////////////////////////////////////////////////

                double token_passing_time_start = MPI_Wtime(); 

                if (ptrn_util_two.input_patterns.size() < 1) {
                    do_nonlocal_constraint_checking = false; 
                }   

                if (do_nonlocal_constraint_checking && global_not_finished) {

                    global_not_finished = false;  

                    VertexUint8Map token_source_map;

                    VectorBoolean pattern_found(ptrn_util_two.input_patterns.size(), 0);
                    VectorBoolean pattern_token_source_found(ptrn_util_two.input_patterns.size(), 0);

                    for (size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++) {

                        std::string paths_result_filename = result_dir +
                            "/all_ranks_subgraphs/subgraphs_" +
                            std::to_string(pl) + "_" + std::to_string(mpi_rank);
                        std::ofstream paths_result_file(paths_result_filename, std::ofstream::out);

                        bool token_source_deleted = false;

                        // setup pattern 
                        auto pattern_tp = std::get<0>(ptrn_util_two.input_patterns[pl]);
                        auto pattern_indices_tp = std::get<1>(ptrn_util_two.input_patterns[pl]);
                        auto pattern_cycle_length_tp = std::get<2>(ptrn_util_two.input_patterns[pl]);
                        auto pattern_valid_cycle_tp = std::get<3>(ptrn_util_two.input_patterns[pl]);
                        auto pattern_edge_data_tp = std::get<6>(ptrn_util_two.input_patterns[pl]);

                        auto pattern_selected_vertices_tp = 0; // TODO: remove

                        auto pattern_is_tds_tp = std::get<4>(ptrn_util_two.input_patterns[pl]); // boolean
                        auto pattern_interleave_label_propagation_tp = std::get<5>(ptrn_util_two.input_patterns[pl]); // boolean

                        auto pattern_enumeration_tp = ptrn_util_two.enumeration_patterns[pl]; 
                        auto pattern_aggregation_steps_tp = ptrn_util_two.aggregation_steps[pl]; 

                        auto pattern_selected_edges_tp = false; // boolean
                        auto pattern_mark_join_vertex_tp = false; // boolean
                        auto pattern_ignore_join_vertex_tp = false; // boolean  
                        size_t pattern_join_vertex_tp = 0; // TODO: 

                        message_count = 0;

                        if(mpi_rank == 0) {
                            std::cout << "Token Passing [" << pl << "] | Searching Subpattern : ";
                            PatternNonlocalConstraint::output_pattern(pattern_tp);
                            std::cout << "Token Passing [" << pl << "] | Vertices : ";
                            PatternNonlocalConstraint::output_pattern(pattern_indices_tp);
                            std::cout << "Token Passing [" << pl << "] | Arguments : " 
                                << pattern_cycle_length_tp << " " 
                                << pattern_valid_cycle_tp << " " 
                                << pattern_interleave_label_propagation_tp << " "
                                << pattern_selected_vertices_tp << std::endl; // Test

                            std::cout << "Token Passing [" << pl << "] | Enumeration Indices : ";
                            PatternNonlocalConstraint::output_pattern(pattern_enumeration_tp);
                            std::cout << "Token Passing [" << pl << "] | Agreegation Steps : TODO" << std::endl;
                        }

                        if (!pattern_selected_vertices_tp) {
                            token_source_map.clear(); // Important
                            vertex_token_source_set.clear(); // Important
                            token_source_edge_data_vertices.clear();
                        } else {
                            token_source_map.clear(); // Important       

                            for (vitr_type vitr = graph->vertices_begin(); 
                                    vitr != graph->vertices_end(); ++vitr) {
                                vloc_type vertex = *vitr;
                                if (vertex_active[vertex] && 
                                        (vertex_metadata[vertex] == pattern_tp[pattern_tp.size() - 1])) {
                                    continue; 
                                } else {
                                    vertex_token_source_set[vertex].clear();
                                    token_source_edge_data_vertices.clear();
                                }
                            } 	

                            for(vitr_type vitr = graph->delegate_vertices_begin();
                                    vitr != graph->delegate_vertices_end(); ++vitr) {
                                vloc_type vertex = *vitr;
                                if (vertex_active[vertex] &&
                                        (vertex_metadata[vertex] == pattern_tp[pattern_tp.size() - 1])) {
                                    continue;
                                } else {
                                    vertex_token_source_set[vertex].clear();
                                    token_source_edge_data_vertices.clear();
                                } 
                            }   

                        }

                        MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

                        time_start = MPI_Wtime();

#ifdef TP_ASYNC
                        if (pattern_is_tds_tp) {
                            prunejuice::token_passing_pattern_matching<graph_type, Vertex, Edge, VertexData, 
                                EdgeData, VertexMetadata, EdgeMetadata, VertexActive, 
                                VertexUint8EdgeDataMapCollection, TemplateVertex, VertexMinMax, VertexStateMap, PatternGraph, 
                                /*PatternUtilities*/PatternNonlocalConstraint, PatternTemporalConstraint, VertexUint8Map, VertexSetCollection, 
                                DelegateGraphVertexDataSTDAllocator, Boolean, BitSet>
                                    (graph, vertex_metadata, enable_edge_matching, enable_temporal_edge_matching, vertex_active, vertex_active_edges_map, 
                                     template_vertices, vertexminmax_vertices, vertex_state_map, pattern_graph, ptrn_util_two, ptrn_temp_const, pl,
                                     token_source_map, vertex_token_source_set, 
                                     pattern_found[pl], tp_vertex_batch_size, paths_result_file, message_count);
                        } else {     
                            prunejuice::token_passing_pattern_matching<graph_type, VertexMetadata, decltype(pattern_tp), decltype(pattern_indices_tp), decltype(pattern_edge_data_tp), uint8_t, PatternGraph,
                                PatternTemporalConstraint,
                                VertexStateMap, VertexUint8Map, edge_data_t, EdgeData,
                                VertexEdgeDataVectorSetCollection, VertexActive, TemplateVertex, VertexMinMax, VertexUint8EdgeDataMapCollection, BitSet>(graph, vertex_metadata, pattern_tp,
                                        pattern_indices_tp, pattern_edge_data_tp, vertex_rank, pattern_graph, ptrn_temp_const, pl, vertex_state_map,
                                        token_source_map, pattern_cycle_length_tp, pattern_valid_cycle_tp,
                                        pattern_found[pl], *edge_data_ptr, enable_edge_matching, enable_temporal_edge_matching, token_source_edge_data_vertices, vertex_active, 
                                        template_vertices, vertexminmax_vertices, vertex_active_edges_map, pattern_selected_vertices_tp, //);
                                pattern_selected_edges_tp, pattern_mark_join_vertex_tp,
                                pattern_ignore_join_vertex_tp, pattern_join_vertex_tp, message_count);
                        }
#endif

                        MPI_Barrier(MPI_COMM_WORLD);
                        time_end = MPI_Wtime();
                        if(mpi_rank == 0) {
                            std::cout << "Pattern Matching Time | Token Passing (Traversal) [" 
                                << pl << "] : " << time_end - time_start << std::endl;
                        }

#ifdef TP_ASYNC

                        uint64_t remove_count = 0;      
                        size_t token_source_pattern_indices_tp_index = 0; // TODO: this is confusing, update 

                        bool is_token_source_map_not_empty = havoqgt::mpi_all_reduce
                            (!token_source_map.empty(), std::greater<uint8_t>(), MPI_COMM_WORLD); // TODO: less did not work?
                        MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here

                        pattern_token_source_found[pl] = is_token_source_map_not_empty;

                        for (auto& s : token_source_map) {
                            if (!s.second) {
                                auto v_locator = graph->label_to_locator(s.first);
                                BitSet v_template_vertices(template_vertices[v_locator]);

                                if (v_template_vertices.none()) {
                                    continue;
                                } 

                                if (v_template_vertices.test(pattern_indices_tp[token_source_pattern_indices_tp_index])) {
                                    assert(pattern_indices_tp[token_source_pattern_indices_tp_index] < max_template_vertex_count); // Test
                                    v_template_vertices.reset(pattern_indices_tp[token_source_pattern_indices_tp_index]);
                                    template_vertices[v_locator] = v_template_vertices.to_ulong();
                                }

                                if (v_template_vertices.none()) {
                                    vertex_active[graph->label_to_locator(s.first)] = false;
                                }

                                if (!global_not_finished) {
                                    global_not_finished = true;
                                }

                                if (!token_source_deleted) {  
                                    token_source_deleted = true;
                                }

                            }  
                        }

                        MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here? // New

                        vertex_active.all_min_reduce(); 
                        MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

                        for(vitr_type vitr = graph->delegate_vertices_begin();
                                vitr != graph->delegate_vertices_end(); ++vitr) {
                            auto vertex = *vitr;
                            if (vertex.is_delegate() && (graph->master(vertex) == mpi_rank)) {
                                continue;  // skip the controller
                            }
                            else {
                                auto find_vertex = vertex_state_map.find(graph->locator_to_label(vertex));
                                if (find_vertex == vertex_state_map.end()) {
                                    template_vertices[vertex] = 0;
                                }
                            }
                        }
                        MPI_Barrier(MPI_COMM_WORLD);

                        template_vertices.all_max_reduce();
                        MPI_Barrier(MPI_COMM_WORLD);

                        for (auto& s : token_source_map) {
                            auto v_locator = graph->label_to_locator(s.first);
                            if (!vertex_active[v_locator]) {
                                auto find_vertex = 
                                    vertex_state_map.find(s.first);

                                if (find_vertex != vertex_state_map.end()) { 

                                    if (vertex_state_map.erase(s.first) < 1) { // s.first is the vertex
                                        std::cerr << "Error: failed to remove an element from the map."
                                            << std::endl;
                                    } else {
                                        remove_count++;
                                    }  

                                }
                            }
                        }

#endif // ifdef TP_ASYNC

                        time_end = MPI_Wtime();
                        if(mpi_rank == 0) {
                            std::cout << "Pattern Matching Time | Token Passing [" << pl 
                                << "] : " << time_end - time_start << std::endl;
                        }

                        // result
                        if(mpi_rank == 0) {
                            superstep_result_file << global_itr_count << ", TP, "
                                << pl << ", "
                                << time_end - time_start << "\n";
                        }

                        // Important : This may slow down things -only for presenting results
                        active_vertices_count = 0;
                        active_edges_count = 0;   

                        for (auto& v : vertex_state_map) {
                            auto v_locator = graph->label_to_locator(v.first);
                            if (v_locator.is_delegate() && (graph->master(v_locator) == mpi_rank)) {
                                active_vertices_count++;

                                // edges
                                active_edges_count+=vertex_active_edges_map[v_locator].size();   
                            } else if (!v_locator.is_delegate()) {
                                active_vertices_count++;

                                // edges
                                active_edges_count+=vertex_active_edges_map[v_locator].size();   
                            }
                        }

                        // vertices
                        active_vertices_count_result_file << global_itr_count << ", TP, "
                            << pl << ", "  
                            << active_vertices_count << "\n";  
                        res.active_vertices_count = active_vertices_count;

                        // edges   
                        active_edges_count_result_file << global_itr_count << ", TP, "
                            << pl << ", "
                            << active_edges_count << "\n";
                        res.active_edges_count = active_edges_count;

                        // messages
                        message_count_result_file << global_itr_count << ", TP, "
                            << pl << ", "
                            << message_count << "\n";

                        // result

                        MPI_Barrier(MPI_COMM_WORLD);

                        if (is_token_source_map_not_empty) {

                            havoqgt::mpi_all_reduce_inplace(pattern_found, std::greater<uint8_t>(), MPI_COMM_WORLD);
                            MPI_Barrier(MPI_COMM_WORLD);
                            if(mpi_rank == 0) {   
                                std::string s = pattern_found[pl] == 1 ? "True" : "False";
                                std::cout << "Token Passing [" << pl << "] | Found Subpattern : " << s << std::endl;
                            }

                            MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here? // New

                            token_source_deleted = havoqgt::mpi_all_reduce(token_source_deleted, std::greater<uint8_t>(), MPI_COMM_WORLD);
                            MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here 
                            if(mpi_rank == 0) {
                                std::cout << "Token Passing [" << pl << "] | Token Source Deleted Status : ";
                                if (token_source_deleted) {
                                    std::cout << "Deleted" << std::endl;
                                } else {
                                    std::cout << "Not Deleted" << std::endl;
                                }
                            } 

                        } else {
                            if(mpi_rank == 0) {
                                std::cout << "Token Passing [" << pl << "] | No Token Source Found" << std::endl;   
                            }
                        }

                        if (token_source_deleted && pattern_interleave_label_propagation_tp) {

                            bool global_not_finished_dummy = true; // TODO: do we need this?

                            // lable propagation   
                            label_propagation_time_start = MPI_Wtime();

                            prunejuice::label_propagation_pattern_matching_bsp<Vertex, VertexData, EdgeData, edge_data_t,
                                graph_type, VertexMetadata, VertexStateMap, VertexActive, 
                                VertexUint8EdgeDataMapCollection, TemplateVertexBitSet, TemplateVertex, VertexMinMaxMap, VertexMinMax, PatternGraph, PatternTemporalConstraint>
                                    (graph, *edge_data_ptr, enable_edge_matching, enable_temporal_edge_matching, vertex_metadata, vertex_state_map, vertex_active, 
                                     vertex_active_edges_map, template_vertices, vertexminmax_vertices, pattern_graph, ptrn_temp_const, global_init_step,
                                     global_not_finished, global_itr_count, superstep_result_file, 
                                     active_vertices_count_result_file, active_edges_count_result_file,
                                     message_count_result_file);  


                            MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
                            label_propagation_time_end = MPI_Wtime();
                            if(mpi_rank == 0) {
                                std::cout << "Fuzzy Pattern Matching Time | Local Constraint Checking (Interleaved) : " 
                                    << label_propagation_time_end - label_propagation_time_start << std::endl;
                            }

                            // result
                            if(mpi_rank == 0) {
                                step_result_file << global_itr_count << ", LP, "
                                    << (label_propagation_time_end - label_propagation_time_start) << "\n";
                            }

                            /////////////////////////////////////////////////////////////////////////////
                        } else {
                            if(mpi_rank == 0) {
                                std::cout << "Pattern Matching | Skipping Local Constraint Checking (Interleaved)." << std::endl;
                            }        
                        }

                        // result
                        paths_result_file.close();

                        MPI_Barrier(MPI_COMM_WORLD);

                    }

                    havoqgt::mpi_all_reduce_inplace(pattern_found, std::greater<uint8_t>(), MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?   
                    if(mpi_rank == 0) {   
                        for (size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++) {
                            if (pattern_token_source_found[pl]) {
                                std::string s = pattern_found[pl] == 1 ? "True" : "False";
                                std::cout << "Token Passing [" << pl << "] | Found Subpattern : " << s << std::endl;
                            } else {
                                std::cout << "Token Passing [" << pl << "] | No Token Source Found" << std::endl;
                            }
                        }
                    }

                    /////////////////////////////////////////////////////////////////////////////

                } else {
                    if(mpi_rank == 0) { 
                        std::cout << "Pattern Matching | Skipping Token Passing." << std::endl;
                    }
                    global_not_finished = false;
                }

                ///////////////////////////////////////////////////////////////////////////// 

                global_not_finished = havoqgt::mpi_all_reduce(global_not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here 

                if(mpi_rank == 0) {
                    std::cout << "Pattern Matching | Global Finished Status : ";
                    if (global_not_finished) {
                        std::cout << "Continue" << std::endl;
                    } else {
                        std::cout << "Stop" << std::endl;
                    }
                }

                /////////////////////////////////////////////////////////////////////////////

                double itr_time_end = MPI_Wtime(); 
                // result 
                if(mpi_rank == 0) {
                    itr_result_file << global_itr_count << ", "
                        << (itr_time_end - itr_time_start) << "\n";
                }

                global_itr_count++; //TODO: sum of LCC and NLCC iterations

                /////////////////////////////////////////////////////////////////////////////

                MPI_Barrier(MPI_COMM_WORLD);  
                double pattern_time_end = MPI_Wtime();
                if(mpi_rank == 0) {
                    std::cout << "Pattern Matching Time | Pattern [" << ps << "] : " 
                        << pattern_time_end - pattern_time_start << std::endl;
                }

                /////////////////////////////////////////////////////////////////////////////

                // result
                if(mpi_rank == 0) {  
                    pattern_set_result_file << ps << ", " 
                        << mpi_size << ", "
                        << global_itr_count << ", " 
                        << (pattern_time_end - pattern_time_start) << ", "
                        << pattern_graph.edge_count << ", " 
                        << pattern_graph.vertex_count << ", "
                        << ptrn_util_two.input_patterns.size() << "\n"; 
                }

                // Important : This may slow things down -only for presenting results

                for (auto& v : vertex_state_map) {
                    auto v_locator = graph->label_to_locator(v.first);
                    if (v_locator.is_delegate() && (graph->master(v_locator) == mpi_rank)) {
                        BitSet v_template_vertices(template_vertices[v_locator]); 
                        active_vertices_result_file << mpi_rank << ", "
                            << v.first << ", " 
                            << v.second.vertex_pattern_index << ", "
                            << vertex_metadata[v_locator] << ", "
                            << v_template_vertices << "\n";

                        // edges 
                        for (auto& n : vertex_active_edges_map[v_locator]) { 	 
                            active_edges_result_file << mpi_rank << ", "
                                << v.first << ", "
                                << n.first //<< ", "
                                //<< (size_t)n.second << ", "
                                //<< vertex_active_edges_map[v_locator].size() 
                                << "\n";
                        }

                    } else if (!v_locator.is_delegate()) {
                        BitSet v_template_vertices(template_vertices[v_locator]);
                        active_vertices_result_file << mpi_rank << ", "
                            << v.first << ", "
                            << v.second.vertex_pattern_index << ", "
                            << vertex_metadata[v_locator] << ", "
                            << v_template_vertices << "\n";

                        // edges
                        for (auto& n : vertex_active_edges_map[v_locator]) {  
                            active_edges_result_file << mpi_rank << ", "
                                << v.first << ", "
                                << n.first
                                << "\n"; 
                        } 

                    }
                }

                // close files
                itr_result_file.close();
                step_result_file.close();
                superstep_result_file.close();
                active_vertices_count_result_file.close(); 
                active_vertices_result_file.close();
                active_edges_count_result_file.close();
                active_edges_result_file.close();
                message_count_result_file.close();

            }

            if (mpi_rank == 0) {
                pattern_set_result_file.close(); // close file
            }

        } // pattern matching  

    } // havoqgt_init
    ;

    return 0;
} // pattern_matching_prunejuice

//find last position that is <= element
template<typename DataType>
bool binary_find(std::vector<DataType> & array, DataType element, size_t begin, size_t end, size_t & pos) {
    if (end - begin <= 1)
        return false;
    if (array[end - 1] <= element) {
        begin = end - 1;
        pos = end - 1;
        return true;
    }

    while (begin < end) {
        pos = (begin + end) / 2;
        if (array[pos] <= element && array[pos + 1] > element) {
            return true;
        }
        if (array[pos+1] == element) {
            pos += 1;
            return true;
        }
        if (array[pos] > element) {
            end = pos;
        }
        if (array[pos + 1] < element) {
            begin = pos + 1;
        }
    }
    return false;
}


typedef struct label_summary{
    std::set<VertexData> unique_vertex_data;
    std::set<std::tuple<VertexData, VertexData, EdgeData>> unique_edges;
    bool repeated_vertex_label;
    bool repeated_edge_label;

    label_summary() : unique_vertex_data(),
                      unique_edges(),
                      repeated_vertex_label(false),
                      repeated_edge_label(false) {}
} LABEL_SUMMARY;

typedef struct csr {
    size_t vertex_num;
    size_t edge_num;

    //vertex data
    std::vector<VertexData> vertex_data;

    //edge list and data
    std::vector<size_t> vertex_edgelist_begin_index;
    std::vector<size_t> edges;
    std::vector<EdgeData> edge_data;

    csr() : vertex_num(0),
            edge_num(0),
            vertex_edgelist_begin_index(),
            vertex_data(),
            edges(),
            edge_data() {}

    void clear() {
        vertex_num = 0;
        edge_num = 0;
        vertex_edgelist_begin_index.clear();
        vertex_data.clear();
        edges.clear();
        edge_data.clear();
    }
} CSR;

typedef struct csr_sort_by_vertex_label {
    size_t vertex_num;
    size_t edge_num;
    size_t vertex_label_type_num;

    //for vertex data originization
    std::vector<VertexData> to_vertex_data;
    std::map<VertexData, size_t> to_vertex_data_id;
    std::vector<size_t> each_label_begin_index;

    //edge list and data
    std::vector<size_t> vertex_edgelist_begin_index;
    std::vector<Vertex> edges;
    std::vector<EdgeData> edge_data;

    //begin_pos of each label for each vertex
    std::vector<size_t> begin_pos_for_each_label;

    csr_sort_by_vertex_label() : vertex_num(0),
            edge_num(0),
            to_vertex_data(),
            to_vertex_data_id(),
            each_label_begin_index(),
            vertex_edgelist_begin_index(),
            edges(),
            edge_data() {}

    void clear() {
        vertex_num = 0;
        edge_num = 0;
        to_vertex_data.clear();
        to_vertex_data_id.clear();
        vertex_edgelist_begin_index.clear();
        edges.clear();
        edge_data.clear();
    }

    VertexData vertex_data(Vertex id) {
        size_t begin = 0;
        size_t end = each_label_begin_index.size() - 1;
        assert(id>=0 && id<each_label_begin_index[end]);
        size_t pos;
        assert(binary_find(each_label_begin_index, id, begin, end, pos) && "wrong vertex id");
        return to_vertex_data[pos];
    }

    size_t vertex_data_to_vertex_data_id(VertexData data) {
        auto it = to_vertex_data_id.find(data);
        assert (it != to_vertex_data_id.end() && "not valid vertex data input for begin end vertex id checking");
        return it -> second;
    }

    void begin_end_vertex_id(VertexData data, size_t & begin, size_t & end) {
        size_t vertex_data_id = vertex_data_to_vertex_data_id(data);
        begin = each_label_begin_index[vertex_data_id];
        end = each_label_begin_index[vertex_data_id+1];
    }

    void init_vertex_label_org(LABEL_SUMMARY & ls) {
        assert(to_vertex_data.size() == 0);
        assert(to_vertex_data_id.size() == 0);

        size_t vertex_data_id = 0;
        for (auto item : ls.unique_vertex_data) {
            to_vertex_data.push_back(item);
            to_vertex_data_id.insert(std::make_pair(item, vertex_data_id));
            vertex_data_id++;
        }
    }

    bool find_neigbors_with_given_vertex_data(Vertex cur_vertex, VertexData vertex_data,
            size_t & begin_index, size_t & end_index) {
        //find vertex to be added to map next
        //get the vertex list in the first connected vertex
        size_t edge_list_begin = vertex_edgelist_begin_index[cur_vertex];
        size_t edge_list_end = vertex_edgelist_begin_index[cur_vertex + 1];

        if (edge_list_end - edge_list_begin < 1) {
            return false;
        }

        size_t seq_begin = (vertex_label_type_num - 1) * cur_vertex;
        size_t vertex_data_id = vertex_data_to_vertex_data_id(vertex_data);
        if (vertex_data_id == 0) {
            begin_index = edge_list_begin;
        } else
            begin_index = begin_pos_for_each_label[seq_begin + vertex_data_id - 1];

        if (vertex_data_id == vertex_label_type_num - 1) {
            end_index = edge_list_end;
        } else
            end_index = begin_pos_for_each_label[seq_begin + vertex_data_id];

        assert(begin_index <= end_index);

        if (begin_index == end_index)
            return false;

        return true;
    }

} CSR_SORT_BY_VERTEX_LABEL;

typedef struct id_map {
    std::map<Vertex, Vertex> to_new;
    std::vector<Vertex> to_orignal;

    id_map() : to_orignal(), to_new() {}

    bool to_new_id(Vertex original_id, Vertex & new_id) {
        auto it = to_new.find(original_id);
        if (it == to_new.end())
            return false;
        new_id = it->second;
        return true;
    }

    bool to_original_id(Vertex new_id, Vertex & original_id) {
        assert(new_id < to_orignal.size());
        original_id = to_orignal[new_id];
        return true;
    }

} ID_MAP;

typedef struct {
    size_t patterns_found;
} OUPUT_SUIT;

typedef struct seq_temporal_constraints{
    std::vector<std::vector<Edge>> smaller_than_constraints;
    std::vector<std::vector<Edge>> bigger_than_constraints;

    seq_temporal_constraints (
            PatternTemporalConstraint & temporal_constraints,
            bool enable_temporal_edge_matching) :
        smaller_than_constraints(), bigger_than_constraints() {
            if (!enable_temporal_edge_matching)
                return;

            for (size_t i=0; i<temporal_constraints.all_constraints.size(); i++) {
                smaller_than_constraints.push_back(std::vector<Edge>());
                bigger_than_constraints.push_back(std::vector<Edge>());
            }

            for (size_t i=0; i<temporal_constraints.all_constraints.size(); i++) {
                for (auto item : temporal_constraints.all_constraints[i]) {
                    smaller_than_constraints[i].push_back(item);
                    bigger_than_constraints[item].push_back(i);
                }
            }
    }

    //assume that edges are (a, b) and a < b
    //edge id are generated by sorting (a, b) and then using the pos of edge
    //as ID
    bool pass_temporal_checking(std::vector<EdgeData> & mapped_edge_data,
            EdgeData new_edge_data) {
        Edge current_edge_id_in_pattern = mapped_edge_data.size();
        //current edge should smaller than edge_data[item]
        for (auto item : smaller_than_constraints[current_edge_id_in_pattern]) {
            assert(item != current_edge_id_in_pattern);
            if (item > current_edge_id_in_pattern)
                break;

            if (mapped_edge_data[item] <= new_edge_data)
                return false;
        }

        //current edge should bigger than edge_data[item]
        for (auto item : bigger_than_constraints[current_edge_id_in_pattern]) {
            assert(item != current_edge_id_in_pattern);
            if (item > current_edge_id_in_pattern)
                break;

            if (mapped_edge_data[item] >= new_edge_data)
                return false;
        }
        return true;
    }

    bool pass_temporal_checking(std::vector<EdgeData> & mapped_edge_data,
            EdgeData new_edge_data, Edge current_edge_id_in_pattern) {
        //current edge should smaller than edge_data[item]
        size_t stored_edge_data_size = mapped_edge_data.size();
        for (auto item : smaller_than_constraints[current_edge_id_in_pattern]) {
            assert(item != current_edge_id_in_pattern);
            if (item > stored_edge_data_size)
                break;

            if (mapped_edge_data[item] <= new_edge_data)
                return false;
        }

        //current edge should bigger than edge_data[item]
        for (auto item : bigger_than_constraints[current_edge_id_in_pattern]) {
            assert(item != current_edge_id_in_pattern);
            if (item > stored_edge_data_size)
                break;

            if (mapped_edge_data[item] >= new_edge_data)
                return false;
        }
        return true;
    }

    bool pass_temporal_checking(std::vector<EdgeData> & mapped_edge_data,
            EdgeData new_edge_data, Edge current_edge_id_in_pattern,
            Edge extra_edge_id, EdgeData extra_edge_data) {
        //current edge should smaller than edge_data[item]
        size_t stored_edge_data_size = mapped_edge_data.size();
        size_t index = 0;
        for (; index < smaller_than_constraints[current_edge_id_in_pattern].size(); index++) {
            auto item = smaller_than_constraints[current_edge_id_in_pattern][index];
            assert(item != current_edge_id_in_pattern);
            if (item > stored_edge_data_size)
                break;

            if (mapped_edge_data[item] <= new_edge_data)
                return false;
        }
        for (; index < smaller_than_constraints[current_edge_id_in_pattern].size(); index++) {
            auto item = smaller_than_constraints[current_edge_id_in_pattern][index];
            assert(item != current_edge_id_in_pattern);
            if (item == extra_edge_id) {
                if (extra_edge_data <= new_edge_data)
                    return false;
            }
        }

        //current edge should bigger than edge_data[item]
        index = 0;
        for (; index < bigger_than_constraints[current_edge_id_in_pattern].size(); index++) {
            auto item = bigger_than_constraints[current_edge_id_in_pattern][index];
            assert(item != current_edge_id_in_pattern);
            if (item > stored_edge_data_size)
                break;

            if (mapped_edge_data[item] >= new_edge_data)
                return false;
        }
        for (; index < bigger_than_constraints[current_edge_id_in_pattern].size(); index++) {
            auto item = bigger_than_constraints[current_edge_id_in_pattern][index];
            assert(item != current_edge_id_in_pattern);
            if (item == extra_edge_id) {
                if (extra_edge_data >= new_edge_data)
                    return false;
            }
        }
        return true;
    }

    bool pass_temporal_checking(std::vector<EdgeData> & mapped_edge_data,
            EdgeData new_edge_data, Edge current_edge_id_in_pattern,
            size_t extra_edge_id_offset,
            std::vector<Edge> & extra_edge_ids,
            std::vector<EdgeData> & extra_edge_data,
            size_t extra_edges_num) {
        //current edge should smaller than edge_data[item]
        assert(extra_edge_ids.size() >= extra_edges_num);
        assert(extra_edge_data.size() >= extra_edges_num);
        size_t stored_edge_data_size = mapped_edge_data.size();
        size_t index = 0;
        for (; index < smaller_than_constraints[current_edge_id_in_pattern].size(); index++) {
            auto item = smaller_than_constraints[current_edge_id_in_pattern][index];
            assert(item != current_edge_id_in_pattern);
            if (item > stored_edge_data_size)
                break;

            if (mapped_edge_data[item] <= new_edge_data)
                return false;
        }
        size_t index_extra_edge_list = 0;
        for (; index < smaller_than_constraints[current_edge_id_in_pattern].size(); index++) {
            auto item = smaller_than_constraints[current_edge_id_in_pattern][index];
            assert(item != current_edge_id_in_pattern);
            while (item > extra_edge_ids[index_extra_edge_list] + extra_edge_id_offset
                    && index_extra_edge_list < extra_edges_num - 1) {
                index_extra_edge_list++;
            }
            if (item == extra_edge_ids[index_extra_edge_list] + extra_edge_id_offset) {
                if (extra_edge_data[index_extra_edge_list] <= new_edge_data)
                    return false;
            }
        }

        //current edge should bigger than edge_data[item]
        index = 0;
        for (; index < bigger_than_constraints[current_edge_id_in_pattern].size(); index++) {
            auto item = bigger_than_constraints[current_edge_id_in_pattern][index];
            assert(item != current_edge_id_in_pattern);
            if (item > stored_edge_data_size)
                break;

            if (mapped_edge_data[item] >= new_edge_data)
                return false;
        }
        index_extra_edge_list = 0;
        for (; index < bigger_than_constraints[current_edge_id_in_pattern].size(); index++) {
            auto item = bigger_than_constraints[current_edge_id_in_pattern][index];
            assert(item != current_edge_id_in_pattern);
            while (item > extra_edge_ids[index_extra_edge_list] + extra_edge_id_offset && index_extra_edge_list < extra_edges_num - 1) {
                index_extra_edge_list++;
            }
            if (item == extra_edge_ids[index_extra_edge_list] + extra_edge_id_offset) {
                if (extra_edge_data[index_extra_edge_list] >= new_edge_data)
                    return false;
            }
        }
        return true;
    }

} PatternTemporalConstraintSeq;

void analysis_pattern(PatternGraph & pattern, LABEL_SUMMARY & ls, bool enable_edge_matching = false) {
    size_t count = 0;
    ls.unique_vertex_data.clear();
    for (size_t v=0; v<pattern.vertex_count; v++) {
        VertexData cur = pattern.vertex_data[v];
        if (ls.unique_vertex_data.find(cur) == ls.unique_vertex_data.end()) {
            ls.unique_vertex_data.insert(cur);
            count++;
        } else {
            ls.repeated_vertex_label = true;
        }
    }
    for (size_t e=0; e<pattern.edge_count; e++) {
        Vertex src = std::get<0>(pattern.edge_list[e]);
        Vertex dst = std::get<1>(pattern.edge_list[e]);
        EdgeData edata = enable_edge_matching? pattern.edge_data[e] : 0;
        if (ls.unique_edges.find(std::make_tuple(pattern.vertex_data[src],
                        pattern.vertex_data[dst], edata))
                == ls.unique_edges.end())
            ls.unique_edges.insert(std::make_tuple(pattern.vertex_data[src],
                        pattern.vertex_data[dst], edata));
        else {
            ls.repeated_edge_label = true;
        }
    }
}

bool read_vertex_metadata(const std::string vertex_metadata_filename,
        CSR_SORT_BY_VERTEX_LABEL & graph,
        LABEL_SUMMARY & ls, ID_MAP & id_mapping) {
    assert(graph.vertex_data.size() == 0 && "vertex data is not cleared");
    std::ifstream vertex_metadata_stream(vertex_metadata_filename, std::ifstream::in);
    if (is_file_empty(vertex_metadata_stream)) {
        return false;
    }
    //construct new ids and map
    std::vector<std::vector<Vertex>> each_label_vertex(ls.unique_vertex_data.size(), std::vector<Vertex>());
    size_t vertex_number = 0;
    size_t num_valid_vertex = 0;
    size_t vertex;
    while (!vertex_metadata_stream.eof()) {
        VertexData data;
        vertex_metadata_stream >> vertex >> data;
        assert(vertex == vertex_number && "metadata file format error");
        vertex_number++;

        auto it = graph.to_vertex_data_id.find(data);
        if (it == graph.to_vertex_data_id.end())
            continue;

        each_label_vertex[it->second].push_back(vertex);
        num_valid_vertex++;
    }
    //construct id_mapping
    id_mapping.to_new.clear();
    id_mapping.to_orignal.resize(num_valid_vertex, vertex);
    size_t new_id_count = 0;
    assert(graph.each_label_begin_index.size() == 0);
    for (size_t i=0; i<ls.unique_vertex_data.size(); i++) {
        if (each_label_vertex[i].size() == 0)
            return false;
        graph.each_label_begin_index.push_back(new_id_count);
        for (auto item : each_label_vertex[i]) {
            id_mapping.to_new.insert(std::make_pair(item, new_id_count));
            id_mapping.to_orignal[new_id_count] = item;
            new_id_count++;
        }
    }
    graph.each_label_begin_index.push_back(new_id_count);

    vertex_metadata_stream.close();
    graph.vertex_num = num_valid_vertex;
    graph.vertex_label_type_num = ls.unique_vertex_data.size();

    return true;
}

bool valid_edge(Vertex org_src, Vertex org_dst, EdgeData edge_data,
        Vertex & new_src, Vertex & new_dst,
        CSR_SORT_BY_VERTEX_LABEL & graph, ID_MAP & id_mapping,
        LABEL_SUMMARY & ls,
        std::map<std::tuple<VertexData, VertexData, EdgeData>, bool> & edge_type_found,
        bool enable_edge_matching = false) {
    if (!id_mapping.to_new_id(org_src, new_src))
        return false;

    if (!id_mapping.to_new_id(org_dst, new_dst))
        return false;

    EdgeData edata = enable_edge_matching ? edge_data : 0;
    VertexData src_label = graph.vertex_data(new_src);
    VertexData dst_label = graph.vertex_data(new_dst);

    if (ls.unique_edges.find(std::make_tuple(src_label, dst_label, edata)) != ls.unique_edges.end()) {
        auto it = edge_type_found.find(std::make_tuple(src_label, dst_label, edata));
        assert(it != edge_type_found.end());
        it -> second = true;
        return true;
    } else {
        return false;
    }
}

bool sorted_edge_list_to_csr(const std::vector<std::tuple<uint64_t, uint64_t, EdgeData>> edge_list,
        CSR_SORT_BY_VERTEX_LABEL & graph,
        LABEL_SUMMARY & ls,
        ID_MAP & id_mapping,
        bool enable_edge_matching = false) {
    assert(graph.vertex_edgelist_begin_index.size() == 0 && "vertex_edgelist_begin_index are not cleared");
    assert(graph.edges.size() == 0 && "edges are not cleared");
    assert(graph.edge_data.size() == 0 && "edge data are not cleared");

    //assume edge list is sorted and vertex num is given in graph
    size_t count = 0;
    size_t valid_edge_count = 0;
    Vertex pre_source = 0;
    Vertex pre_dest = 0;
    std::vector<std::vector<std::pair<Vertex, EdgeData>>> new_edge_list(
            graph.vertex_num, std::vector<std::pair<Vertex, EdgeData>>());

    //record if all edge types needed are found, default as false and
    //update to true when found at valid_edge function.
    std::map<std::tuple<VertexData, VertexData, EdgeData>, bool> edge_type_found;
    for (auto const item : ls.unique_edges) {
        edge_type_found.insert(std::make_pair(item, false));
    }

    while (count < edge_list.size()) {
        Vertex cur_source = std::get<0>(edge_list[count]);
        Vertex cur_dest = std::get<1>(edge_list[count]);
        EdgeData cur_edge_data = std::get<2>(edge_list[count]);
        //check it's sorted
        assert(cur_source >= pre_source && "edge_list is not sorted");
        assert(((cur_source != pre_source) || cur_dest > pre_dest) && "edge_list is not sorted");
        count ++;
        pre_source = cur_source;
        pre_dest = cur_dest;

        Vertex new_source = 0;
        Vertex new_dest = 0;
        if (valid_edge(cur_source, cur_dest, cur_edge_data, new_source, new_dest,
                    graph, id_mapping, ls, edge_type_found, enable_edge_matching)) {
            //Add edge to new_edge_list
            new_edge_list[new_source].push_back(std::make_pair(new_dest, cur_edge_data));
            valid_edge_count++;
        }
    }
    graph.edge_num = valid_edge_count;

    //check if there is a edge type with sign false
    for (auto item : edge_type_found) {
        if (item.second == false)
            return false;
    }

    //constrcut graph
    assert(graph.vertex_edgelist_begin_index.size() == 0);
    assert(graph.edges.size() == 0);
    assert(graph.edge_data.size() == 0);
    size_t begin_index = 0;
    for (size_t i = 0; i<graph.vertex_num; i++) {
        graph.vertex_edgelist_begin_index.push_back(begin_index);
        sort(new_edge_list[i].begin(), new_edge_list[i].end());

        size_t vertex_data_id = 0;
        for (size_t j=0; j<new_edge_list[i].size(); j++, begin_index++) {
            Vertex dst = std::get<0>(new_edge_list[i][j]);
            graph.edges.push_back(dst);
            graph.edge_data.push_back(std::get<1>(new_edge_list[i][j]));

            while (dst >= graph.each_label_begin_index[vertex_data_id+1]) {
                graph.begin_pos_for_each_label.push_back(begin_index);
                vertex_data_id++;
                assert(vertex_data_id < graph.each_label_begin_index.size() - 1);
            }
        }
        while (vertex_data_id < graph.vertex_label_type_num - 1) {
            graph.begin_pos_for_each_label.push_back(begin_index);
            vertex_data_id++;
        }
    }
    graph.vertex_edgelist_begin_index.push_back(begin_index);

    return true;
}

bool generate_csr_sort_by_vertex_label(const std::vector<std::tuple<uint64_t, uint64_t, EdgeData>> edge_list,
        const std::string vertex_metadata_filename, CSR_SORT_BY_VERTEX_LABEL & graph,
        LABEL_SUMMARY & ls, ID_MAP & id_mapping, bool enable_edge_matching = false) {
    //assume edge list is sorted
    //every vertex's metadata is given in order
    graph.clear();
    //read vertex metadata, only keep lebals involved and reorder vertex
    //based on label
    graph.init_vertex_label_org(ls);
    if (!read_vertex_metadata(vertex_metadata_filename, graph, ls, id_mapping))
        //existing a needed vertex data not found in background graph
        return false;
    if (!sorted_edge_list_to_csr(edge_list, graph, ls, id_mapping, enable_edge_matching))
        //existing a kind of edge not found in background graph
        return false;
    return true;
}

void binary_search_set_intersection_with_constraints(
        size_t first,
        size_t second,
        std::vector<Vertex> & connected_vertices,
        std::vector<std::pair<size_t, size_t>> & begin_end_index,
        std::vector<Vertex> & possible_next_vertices,
        std::vector<std::vector<EdgeData>> & possible_next_vertices_edgedata,
        CSR_SORT_BY_VERTEX_LABEL & graph,
        PatternGraph & pattern,
        size_t pattern_neighbor_list_begin,
        std::vector<EdgeData> & mapped_edge_data,
        PatternTemporalConstraintSeq & temporal_constraints_seq,
        bool enable_edge_matching = false,
        bool enable_temporal_edge_matching = false) {
    assert(possible_next_vertices.size() == 0);
    assert(possible_next_vertices_edgedata.size() == 0);

    size_t first_begin = std::get<0>(begin_end_index[first]);
    size_t first_end = std::get<1>(begin_end_index[first]);
    size_t second_begin = std::get<0>(begin_end_index[second]);
    size_t second_end = std::get<1>(begin_end_index[second]);
    for (size_t i=first_begin; i<first_end; i++) {
        EdgeData edata;
        if (enable_edge_matching || enable_temporal_edge_matching) {
            edata = graph.edge_data[i];
            if (enable_edge_matching) {
                if (!(pattern.edge_data[pattern_neighbor_list_begin + first] == edata))
                    continue;
            }
            if (enable_temporal_edge_matching) {
                if (!temporal_constraints_seq.pass_temporal_checking(mapped_edge_data, edata, pattern_neighbor_list_begin + first))
                    continue;
            }
        }
        size_t pos;
        Vertex neig = graph.edges[i];
        if (binary_find(graph.edges, neig, second_begin, second_end, pos)) {
            if (graph.edges[pos] == neig) {
                possible_next_vertices.push_back(neig);
                if (enable_temporal_edge_matching)
                    possible_next_vertices_edgedata.push_back(std::vector<EdgeData>{edata, graph.edge_data[pos]});
            }
            second_begin = pos + 1;
        }
        if (second_end - second_begin <= 1)
            break;
    }
}

void recurrence_pattern_matching (std::vector<Vertex> mapping,
        std::vector<EdgeData> mapped_edge_data, size_t & number_of_patterns_found,
        CSR_SORT_BY_VERTEX_LABEL & graph, PatternGraph & pattern,
        PatternTemporalConstraintSeq & temporal_constraints_seq,
        std::vector<bool> & active_vertices, std::vector<bool> & active_edges,
        bool enable_edge_matching = false,
        bool enable_temporal_edge_matching = false) {
    //When we find a valid mapping
    if (mapping.size() == pattern.vertex_count) {
        for (Vertex src=0; src<pattern.vertex_count; src++) {
            Vertex src_background_graph = mapping[src];
            active_vertices[src_background_graph] = true;
            std::vector<Vertex> active_neighbors;
            for (size_t j=pattern.vertices[src]; j<pattern.vertices[src+1]; j++) {
                Vertex dst = pattern.edges[j];
                active_neighbors.push_back(mapping[dst]);
            }
            assert(set_edges_active(active_edges, graph,
                        src_background_graph, active_neighbors));
        }
        number_of_patterns_found++;
        return ;
    }

    size_t cur_vertex = mapping.size();
    assert(cur_vertex < pattern.vertex_count);
    //TODO: making sure of this in pattern analysis; for each vertex,
    //except the first vertex, it always connects to at least one vertex that has
    //id smaller that it.

    //Collect vertices that new vertex should be connected to
    size_t pattern_neighbor_list_begin = pattern.vertices[cur_vertex];
    size_t pattern_neighbor_list_end = pattern.vertices[cur_vertex + 1];
    size_t pattern_cur_index = pattern_neighbor_list_begin;
    std::vector<Vertex> connected_vertices;
    while ((pattern_cur_index < pattern_neighbor_list_end) && pattern.edges[pattern_cur_index] < cur_vertex) {
        connected_vertices.push_back(mapping[pattern.edges[pattern_cur_index]]);
        pattern_cur_index++;
    }
    size_t connected_vertices_num = connected_vertices.size();

    //Collect connected vertices' neighbor index with given vertex data
    std::vector<std::pair<size_t, size_t>> begin_end_index;
    for (size_t i=0; i<connected_vertices_num; i++) {
        size_t begin_index, end_index;
        if (!graph.find_neigbors_with_given_vertex_data(connected_vertices[i],
                    pattern.vertex_data[cur_vertex], begin_index, end_index))
            //connected vertex doesn't have
            //neighbors with required vertex data.
            return ;
        begin_end_index.push_back(std::make_pair(begin_index, end_index));
    }

    if (connected_vertices_num <= 1) {
        //just add that vertex's all neighbors that have the required
        //vertex data
        size_t begin_index, end_index;
        if (connected_vertices_num == 1) {
            begin_index = std::get<0>(begin_end_index[0]);
            end_index = std::get<1>(begin_end_index[0]);
        } else {
            graph.begin_end_vertex_id(pattern.vertex_data[cur_vertex],
                    begin_index, end_index);
        }
        for (size_t i=begin_index; i<end_index; i++) {
            //edge matching
            if (connected_vertices_num == 1 &&
                    (enable_edge_matching || enable_temporal_edge_matching)) {
                EdgeData edata = graph.edge_data[i];
                if (enable_edge_matching) {
                    if (!(pattern.edge_data[pattern_neighbor_list_begin] == edata))
                        continue;
                }
                if (enable_temporal_edge_matching) {
                    if (!temporal_constraints_seq.pass_temporal_checking(mapped_edge_data, edata))
                        continue;
                    mapped_edge_data.push_back(edata);
                }
            }
            //temporal edge
            mapping.push_back(graph.edges[i]);
            recurrence_pattern_matching(mapping, mapped_edge_data, number_of_patterns_found,
                    graph, pattern, temporal_constraints_seq, active_vertices, active_edges);
            mapping.pop_back();
            mapped_edge_data.pop_back();
        }
    } else {
        //reorder so to reduce binary search amount
        std::vector<size_t> to_new_order(connected_vertices_num, 0);
        std::vector<size_t> to_org_order(connected_vertices_num, 0);
        std::vector<std::pair<size_t, size_t>> neighbors_count;
        for (size_t i=0; i<connected_vertices_num; i++) {
            size_t count = std::get<1>(begin_end_index[i]) - std::get<0>(begin_end_index[i]);
            neighbors_count.push_back(std::make_pair(count, i));
        }
        sort(neighbors_count.begin(), neighbors_count.end());
        for (size_t i=0; i<connected_vertices_num; i++) {
            size_t org = std::get<1>(neighbors_count[i]);
            to_new_order[org] = i;
            to_org_order[i] = org;
        }
        //reorder so to reduce binary search amount

        //for the first two vertices to process
        std::vector<Vertex> possible_next_vertices;
        std::vector<std::vector<EdgeData>> possible_next_vertices_edgedata;
        binary_search_set_intersection_with_constraints(
                to_org_order[0],
                to_org_order[1],
                connected_vertices,
                begin_end_index,
                possible_next_vertices,
                possible_next_vertices_edgedata,
                graph,
                pattern,
                pattern_neighbor_list_begin,
                mapped_edge_data,
                temporal_constraints_seq,
                enable_edge_matching,
                enable_temporal_edge_matching);
        if (possible_next_vertices.size() == 0)
            return;
        std::vector<bool> pass_checking(possible_next_vertices.size(), true);
        size_t possible_next_vertices_num = possible_next_vertices.size();
        size_t num_mapped_edges = mapped_edge_data.size();
        //herizontal verification
        for (size_t cur_order = 2; cur_order < connected_vertices_num; cur_order++) {
            size_t cur_begin_index = std::get<0>(begin_end_index[to_org_order[cur_order]]);
            size_t cur_end_index = std::get<0>(begin_end_index[to_org_order[cur_order]]);
            for (size_t i=0; i<possible_next_vertices_num; i++) {
                if (pass_checking[i] == true) {
                    Vertex possible_vertex = possible_next_vertices[i];
                    size_t pos;
                    if (!binary_find(graph.edges, possible_vertex, cur_begin_index, cur_end_index, pos)) {
                        pass_checking[i] = false;
                        continue;
                    }
                    if (graph.edges[pos] == possible_vertex) {
                        if (enable_edge_matching || enable_temporal_edge_matching) {
                            EdgeData edata = graph.edge_data[i];
                            if (enable_edge_matching) {
                                if (!(pattern.edge_data[num_mapped_edges + to_org_order[cur_order]] == edata)) {
                                    pass_checking[i] = false;
                                    continue;
                                }
                            }
                            if (enable_temporal_edge_matching) {
                                if (!temporal_constraints_seq.pass_temporal_checking(mapped_edge_data,
                                            edata, num_mapped_edges + to_org_order[cur_order],
                                            num_mapped_edges, to_org_order,
                                            possible_next_vertices_edgedata[i], cur_order)) {
                                    pass_checking[i] = false;
                                    continue;
                                }
                                possible_next_vertices_edgedata[i].push_back(edata);
                            }
                        }
                    }
                    cur_begin_index = pos + 1;
                    if (cur_end_index - cur_begin_index <= 1)
                        break;
                }
            }
        }

        for (size_t i=0; i<possible_next_vertices_num; i++) {
            if (pass_checking[i] == true) {
                mapping.push_back(possible_next_vertices[i]);
                if (enable_temporal_edge_matching) {
                    mapped_edge_data.resize(num_mapped_edges + connected_vertices_num);
                    for (size_t j=0; j<connected_vertices_num; j++) {
                        mapped_edge_data[num_mapped_edges + to_org_order[j]] = possible_next_vertices_edgedata[i][j];
                    }
                }
                recurrence_pattern_matching(mapping, mapped_edge_data, number_of_patterns_found,
                        graph, pattern, temporal_constraints_seq, active_vertices, active_edges);
                mapping.pop_back();
                if (enable_temporal_edge_matching)
                    mapped_edge_data.resize(num_mapped_edges);
            }
        }
    }
}

size_t pattern_matching_seq(const std::vector<std::tuple<uint64_t, uint64_t, EdgeData>> edge_list,
        const std::string vertex_metadata_filename,
        const std::string pattern_input_filename,
        OUTPUT_SUIT & res,
        bool enable_edge_matching = false,
        bool enable_temporal_edge_matching = false) {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    if (mpi_rank != 0)
        return 0;

    //construct pattern graph and read temporal constraints
    PatternGraph pattern(
            pattern_input_filename + "_edge",
            pattern_input_filename + "_vertex",
            pattern_input_filename + "_vertex_data",
            pattern_input_filename + "_edge_data",
            pattern_input_filename + "_stat",
            false, false);
    PatternTemporalConstraint temporal_constraints(pattern,
                        pattern_input_filename + "/pattern_temporal_constraint", enable_temporal_edge_matching);
    PatternTemporalConstraintSeq temporal_constraints_seq(temporal_constraints, enable_temporal_edge_matching);

    LABEL_SUMMARY ls;
    analysis_pattern(pattern, ls, enable_edge_matching);

    //change edge_list to csr format with weight and metadata
    //while filtering out vertex with label not needed and edges not
    //needed
    CSR_SORT_BY_VERTEX_LABEL graph;
    ID_MAP id_mapping; //to recover graphs and compare correctness
    if (!generate_csr_sort_by_vertex_label(edge_list,  vertex_metadata_filename, graph, ls, id_mapping, enable_edge_matching)) {
        //in basic local checking, found it won't have match.
        res.patterns_found = 0;
        return 0;
    }

    //create active vertices/edges list
    std::vector<bool> active_vertices(graph.vertex_num, false);
    //TODO: optimize
    std::vector<bool> active_edges(graph.edge_num, false);
    size_t number_of_patterns_found = 0;
    std::vector<Vertex> mapping;
    std::vector<EdgeData> mapped_edge_data;
    recurrence_pattern_matching(mapping, mapped_edge_data, number_of_patterns_found,
            graph, pattern, temporal_constraints_seq,
            active_vertices, active_edges,
            enable_edge_matching, enable_temporal_edge_matching);
    res.patterns_found = number_of_patterns_found;
    res.active_vertices_count = 0;
    res.active_edges_count = 0;
    for (size_t i=0; i<graph.vertex_num; i++)
        if (active_vertices[i])
            res.active_vertices_count++;
    for (size_t i=0; i<graph.edge_num; i++)
        if (active_edges[i])
            res.active_edges_count++;
    //search matched pattern in a recurrence way and active vertices and
    //edges involved matches
    return 0;
}

bool set_edges_active(std::vector<bool> & active_edges,
        CSR_SORT_BY_VERTEX_LABEL & graph,
        Vertex src, std::vector<Vertex> & active_neighbors) {
    sort(active_neighbors.begin(), active_neighbors.end());
    size_t begin = graph.vertex_edgelist_begin_index[src];
    size_t end = graph.vertex_edgelist_begin_index[src+1];
    size_t cur = (begin + end) / 2;
    size_t count = 0;
    while (count < active_neighbors.size()) {
        size_t pos;
        binary_find(graph.edges, active_neighbors[count], begin, end, pos);
        assert (pos >= begin && pos<end && "position found is not valid");
        assert (graph.edges[pos] == active_neighbors[count] && "can't find the given neighbor");
        active_neighbors[pos] = true;
        begin = pos + 1;
        count++;
    }
    return true;
}

void generate_vertex_metadata(graph_type * graph, VertexMetadata & vertex_metadata, std::string vertex_metadata_input = std::string()) {
  if (vertex_metadata_input.size() > 0) {
    vertex_data_db_nostdfs<graph_type, VertexMetadata, Vertex, VertexData>
      (graph, vertex_metadata, vertex_metadata_input, 10000);
  } else {
    vertex_data_db_degree<graph_type, VertexMetadata, Vertex, VertexData>
      (graph, vertex_metadata);
  }
}


