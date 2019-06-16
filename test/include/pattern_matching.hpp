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

typedef int EdgeData; 
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

typedef graph_type::vertex_data<uint64_t, std::allocator<uint64_t> > VertexIteration; // TODO: delete
typedef graph_type::vertex_data<VertexRankType, std::allocator<VertexRankType> > VertexRank; // TODO: delete

typedef prunejuice::vertex_state<Vertex, VertexData, BitSet> VertexState;
typedef std::unordered_map<Vertex, VertexState> VertexStateMap; // TODO: solution_graph

typedef std::unordered_set<Vertex> VertexSet;  
typedef graph_type::vertex_data<VertexSet, std::allocator<VertexSet> > VertexSetCollection; 

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

size_t pattern_matching_prunejuice(graph_type * graph,
        VertexMetadata & vertex_metadata,
        edge_data_t * edge_data_ptr,
        const std::string pattern_input,
        const std::string result_output, bool enable_edge_matching = false,
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
            VertexUint8EdgeDataMapCollection vertex_active_edges_map(*graph);
            VertexSetCollection vertex_token_source_set(*graph); // per vertex set

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

                vertex_state_map.clear(); // important
                vertex_active.reset(true); // initially all vertices are active
                vertex_active_edges_map.clear(); // important
                vertex_token_source_set.clear(); // clear all the sets on all the vertices

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
                    VertexUint8EdgeDataMapCollection, TemplateVertexBitSet, TemplateVertex, PatternGraph>
                        (graph, *edge_data_ptr, enable_edge_matching, vertex_metadata, vertex_state_map, vertex_active, 
                         vertex_active_edges_map, template_vertices, pattern_graph, global_init_step, 
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
                                } 
                            }   

                        }

                        MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

                        time_start = MPI_Wtime();

#ifdef TP_ASYNC
                        if (pattern_is_tds_tp) {
                            prunejuice::token_passing_pattern_matching<graph_type, Vertex, Edge, VertexData, 
                                EdgeData, VertexMetadata, EdgeMetadata, VertexActive, 
                                VertexUint8EdgeDataMapCollection, TemplateVertex, VertexStateMap, PatternGraph, 
                                PatternNonlocalConstraint, VertexUint8Map, VertexSetCollection, 
                                DelegateGraphVertexDataSTDAllocator, Boolean, BitSet>
                                    (graph, vertex_metadata, enable_edge_matching, vertex_active, vertex_active_edges_map, 
                                     template_vertices, vertex_state_map, pattern_graph, ptrn_util_two, pl,
                                     token_source_map, vertex_token_source_set, 
                                     pattern_found[pl], tp_vertex_batch_size, paths_result_file, message_count);
                        } else {     
                            prunejuice::token_passing_pattern_matching<graph_type, VertexMetadata, decltype(pattern_tp),
                                decltype(pattern_indices_tp), decltype(pattern_edge_data_tp), uint8_t, PatternGraph,
                                VertexStateMap, VertexUint8Map, edge_data_t, EdgeData,
                                VertexSetCollection, VertexActive, TemplateVertex, VertexUint8EdgeDataMapCollection, BitSet>
                                    (graph, vertex_metadata, pattern_tp,
                                     pattern_indices_tp, pattern_edge_data_tp, vertex_rank, pattern_graph, vertex_state_map,
                                     token_source_map, pattern_cycle_length_tp, pattern_valid_cycle_tp,
                                     pattern_found[pl], *edge_data_ptr, enable_edge_matching, vertex_token_source_set, vertex_active, 
                                     template_vertices, vertex_active_edges_map, pattern_selected_vertices_tp,
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

                        // edges   
                        active_edges_count_result_file << global_itr_count << ", TP, "
                            << pl << ", "
                            << active_edges_count << "\n";

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
                                VertexUint8EdgeDataMapCollection, TemplateVertexBitSet, TemplateVertex, PatternGraph>
                                    (graph, *edge_data_ptr, enable_edge_matching, vertex_metadata, vertex_state_map, vertex_active, 
                                     vertex_active_edges_map, template_vertices, pattern_graph, global_init_step, 
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

size_t pattern_matching_seq(graph_type * graph,
        VertexMetadata & vertex_metadata,
        edge_data_t * edge_data_ptr,
        const std::string pattern_input) {
    //return expected value
    return 0;
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


