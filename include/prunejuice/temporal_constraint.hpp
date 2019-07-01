#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>  
#include <vector>
#include <tuple>
#include <limits>
#include <algorithm>

#include <boost/algorithm/string.hpp>

#include <prunejuice/file_utilities.hpp>
#include <prunejuice/utilities.hpp>

using namespace prunejuice::utilities;

template <typename Vertex, typename Edge, typename VertexData, typename EdgeData, typename PatternGraph, typename PatternNonlocalConstraint, typename BitSet, typename VertexMinMaxMap>
class pattern_temporal_constraint {
    public:
    pattern_temporal_constraint(PatternGraph& pattern_graph,
            PatternNonlocalConstraint& pattern_non_local_constraint,
            std::string pattern_temporal_constraint_filename, bool enable_edge_temporal_matching = true, bool _directed=false) :
        edge_count(0),
        vertex_count(0),
        all_constraints(),
        local_vertices(0),
        local_neighbors(),
        local_constraints(),
        non_local_constraints(),
        directed(_directed) {
        if (!enable_edge_temporal_matching) {
            for (size_t i; i<pattern_non_local_constraint.input_patterns.size(); i++)
                non_local_constraints.push_back(std::vector<NonLocalOPT>());
            return ;
        }

        if (directed) {
            //TODO
            std::cerr << "directed not supported yet" << std::endl;
            return ;
        }
        init(pattern_graph);
        //output_pattern_graph_info();
        read_all_constraints(pattern_graph, pattern_temporal_constraint_filename);
        //output_all_constraints();
        check_all_constraints(pattern_graph);
        generate_local_opt(pattern_graph); // recognize local and global constraints too
        generate_non_local_opt(pattern_graph, pattern_non_local_constraint);
        output_temporal_constraints();
    }

    ~pattern_temporal_constraint() {
        //TODO: clear memory
    }

    void output_temporal_constraints() {
        std::cout << "local temporal constraints: " << std::endl;
        for (auto cl : local_constraints) {
            std::cout << "In " << (cl.first).first << "'s neighborhood, " << (cl.first).second << ":" << std::endl;
            for (auto v : cl.second) {
                std::cout << ((v.second == true)? "smaller than " : "larger than ") <<  v.first << std::endl;
            }
        }
        std::cout << "BitSet:" << std::endl;
        std::cout << "local vertices: " << local_vertices.to_string() << std::endl;
        std::cout << "local neigbor:" << std::endl;
        for (auto neig : local_neighbors)
            std::cout << neig.to_string() << std::endl;

        std::cout << "global temporal constraints:" << std::endl;
        for (size_t i=0; i<non_local_constraints.size(); i++) {
            std::cout << "non local checking " << i << " :" << std::endl;
            auto gc = non_local_constraints[i];
            for (size_t s=0; s<gc.size(); s++) {
                std::cout << "step " << i << "(" << s << ")" << " :";
                auto opt = gc[s];
                opt.output_operation();
                std::cout << std::endl;
            }
        }
    }

    //call when cur vertex recevies a message 1
    bool do_local_checking(const BitSet cur) {
        return (cur & local_vertices).any();
    }

    //call before the process of each role of cur vertex
    bool do_local_checking(const Vertex cur, const BitSet neighbor) {
        return (neighbor & local_neighbors[cur]).any();
    }
    bool do_local_checking(const Vertex cur, const Vertex neighbor) {
        return local_neighbors[cur].test(neighbor);
    }

    //TODO: Jing another way of shrink checking and resetting is:
    //resetting: clear all the items
    //shrink: search last_itr_min_max in cur, and it couldn't find it then
    //put it in remove list and then remove items in the remove list after
    //updating loop.
    bool shrink_min_max_range(VertexMinMaxMap & last_itr_min_max, VertexMinMaxMap & cur_itr_min_max) {
        bool updated = false;
        for (auto item : cur_itr_min_max) {
            auto find_vertex = last_itr_min_max.find(item.first);
            if (find_vertex == last_itr_min_max.end()) {
                std::cerr << "Error: couldn't find cur min max in last iteration" << std::endl;
                return false;
            }
            if (std::get<0>(find_vertex -> second) < std::get<0>(item.second)) {
                std::get<0>(find_vertex -> second) = std::get<0>(item.second);
                updated = true;
            } else if (std::get<0>(find_vertex -> second) > std::get<0>(item.second)) {
                std::cerr << "Error: last iteration's min is bigger" << std::endl;
                return false;
            }
            if (std::get<1>(find_vertex -> second) > std::get<1>(item.second)) {
                std::get<1>(find_vertex -> second) = std::get<1>(item.second);
                updated = true;
            } else if (std::get<1>(find_vertex -> second) < std::get<1>(item.second)) {
                std::cerr << "Error: last iteration's max is smaller" << std::endl;
                return false;
            }

        }
        return updated;
    }

    void reset_min_max(VertexMinMaxMap & cur_itr_min_max) {
        for (auto & item : cur_itr_min_max) {
            std::get<0>(item.second) = max_edgedata_type();
            std::get<1>(item.second) = min_edgedata_type();
        }
    }
 
    void update_min_max(Vertex cur, Vertex neighbor, EdgeData edge_data, VertexMinMaxMap & min_max) {
        auto find = min_max.find(neighbor);
        if (find == min_max.end()) {
            auto insert_status = min_max.insert( {neighbor, std::tuple<EdgeData, EdgeData>(max_edgedata_type(), min_edgedata_type())} );
            if (!insert_status.second) {
                std::cerr << "Fail to insert neighbor-<min, max> pair" << std::endl;
                return ;
            }
            find = insert_status.first;
        }
        if (std::get<0>(find->second) > edge_data)
            std::get<0>(find->second) = edge_data;
        if (std::get<1>(find->second) < edge_data)
            std::get<1>(find->second) = edge_data;
    }

    bool compare_min_max(Vertex cur, Vertex neighbor, EdgeData edge_data, VertexMinMaxMap & min_max) {
        auto find_compare_list = local_constraints.find(std::make_pair(cur, neighbor));
        if (find_compare_list == local_constraints.end()) {
            std::cerr << "Fail to find compare list in the temporal pattern" << std::endl;
            return false;
        }

        for (auto compare_to : find_compare_list -> second) {
            Edge compare_to_edge = std::get<0>(compare_to);
            bool smaller_to = std::get<1>(compare_to);

            auto find_min_max = min_max.find(compare_to_edge);
            if (find_min_max == min_max.end()) {
                std::cerr << "Fail to find <min, max> from last round" << std::endl;
                return false;
            }

            if (smaller_to) {
                if (edge_data >= std::get<1>(find_min_max -> second))
                    return false;
            } else {
                if (edge_data <= std::get<0>(find_min_max -> second))
                    return false;
            }
        }

        return true;
    }

    struct NonLocalOPT {
        static constexpr size_t OPT_KIND = 2; // TODO:
        static constexpr size_t STORE_POS = 1;
        static constexpr size_t COMP_POS = 0;
        std::bitset<OPT_KIND> opt;
        std::vector<std::pair<Edge, bool>> compare_to_edges;
        size_t store_at;

        NonLocalOPT() : opt(0) {}
        NonLocalOPT(std::bitset<OPT_KIND> _opt,
                std::vector<std::pair<Vertex, bool>> _compare_to_edges) :
            opt(_opt), compare_to_edges() {
                if (store())
                    assert(compare_to_edges.size() == 0);
        }

        void set_store(size_t _store_at = 0) {
            opt.set(STORE_POS);
            store_at = _store_at;
        }

        void set_compare(std::vector<std::pair<Vertex, bool>> & _compare_to_edges) {
            opt.set(COMP_POS);
            compare_to_edges.copy(_compare_to_edges);
        }

        void set_compare() {
            opt.set(COMP_POS);
        }

        void insert_compare(Edge compare_to_edge, bool smaller_to) {
            compare_to_edges.push_back(std::make_pair(compare_to_edge, smaller_to));
        }

        bool store() {
            return opt.test(STORE_POS);
        }

        bool compare() {
            return opt.test(COMP_POS);
        }

        bool operate(std::vector<EdgeData> & stored, EdgeData cur_edge_data) {
            if (update_stored_edgedata(stored, cur_edge_data))
                return checking(stored, cur_edge_data);
            return false;
        }

        bool update_stored_edge_data(std::vector<EdgeData> & stored, EdgeData cur_edge_data) {
            if (store()) {
                size_t stored_num = stored.size();
                if (stored_num == store_at)
                    stored.push_back(cur_edge_data);
                else if (stored_num > store_at)
                    stored[store_at] = cur_edge_data;
                else {
                    std::cerr << "Error: storing from operator; there is missing storage" << std::endl;
                    return false;
                }

            }
            return true;
        }

        bool checking(const std::vector<EdgeData> & stored, EdgeData cur_edge_data) {
            if (compare()) {
                for (auto comp : compare_to_edges) {
#ifdef DEBUG_TEMPORAL
                    std::cout << "opt comp: " << comp.second << " " << stored[comp.first] << " " << cur_edge_data << std::endl;
                    if (comp.first >= stored.size()) {
                        std::cerr << "Error: invalid access of compared position" << std::endl;
                    }
#endif
                    EdgeData compared_to_edge_data = stored[comp.first];
                    if (comp.second && (cur_edge_data >= compared_to_edge_data))
                        return false;
                    if (!comp.second && (cur_edge_data <= compared_to_edge_data))
                        return false;
                }
            }
            return true;
        }

        void output_operation_kinds () {
            std::cout << "no operation" << std::endl;
            std::cout << "compare" << std::endl;
            std::cout << "store" << std::endl;
            std::cout << "store and compare" << std::endl;
        }

        void output_operation () {
            if (opt.none()) {
                std::cout << " no";
                return ;
            }

            if (store())
                std::cout << " store at: " << store_at << ";";
            if (compare()) {
                std::cout << " compare:";
                for (auto comp : compare_to_edges) {
                    std::cout << ((comp.second == true)? "smaller than " : "larger than ") \
                                            << comp.first << ";";
                }
            }
        }
    };

  private:

    bool directed;
    Edge edge_count;
    Vertex vertex_count;
    //all constraints
    std::vector<std::vector<Edge>> all_constraints;
    // For storing local temporal checking
    BitSet local_vertices;
    std::vector<BitSet> local_neighbors;
    std::unordered_map<std::pair<Vertex, Vertex>, std::vector<std::pair<Vertex, bool>>, boost::hash<std::pair<Vertex, Vertex>>> local_constraints;
    //TODO: global opt
    std::vector<std::vector<NonLocalOPT>> non_local_constraints;
    std::unordered_map<std::tuple<Vertex, Vertex>, Edge, boost::hash<std::tuple<Vertex, Vertex>>> edge_to_id_map;
    std::unordered_map<Edge, std::tuple<Vertex, Vertex>> id_to_edge_map;

    void init(PatternGraph & pattern_graph) {
        //get edge count from edge_ID
        std::vector<Edge> edge_id(pattern_graph.edge_ID);
        std::sort(edge_id.begin(), edge_id.end());
        auto ip = std::unique(edge_id.begin(), edge_id.end());
        edge_count = std::distance(edge_id.begin(), ip);
        edge_id.clear();

        vertex_count = pattern_graph.vertex_count;
        create_edge_to_pos_map(pattern_graph);
    }

    //if undirected, edge <a, b> and <b, a> have same edge id
    void read_all_constraints(PatternGraph& pattern_graph,
            std::string pattern_temporal_constraint_filename) {
        std::ifstream pattern_temporal_constraint_file
            (pattern_temporal_constraint_filename, std::ifstream::in);
        if (is_file_empty(pattern_temporal_constraint_file)) {
            return;
        }
        std::string line;
        all_constraints.clear();
        for (Edge i=0; i<edge_count; i++) {
            all_constraints.push_back(std::vector<Edge>());
        }

        Edge line_count = 0;
        while (std::getline(pattern_temporal_constraint_file, line)) {

            boost::trim(line);
            auto tokens = prunejuice::utilities::split<int>(line, ' ');
            assert(tokens.size() == edge_count);
            assert(line_count < edge_count);

            for (Edge i=0; i<edge_count; i++) {
                if (tokens[i] != 0) {
                    assert(i != line_count); //TODO: later we can put edge data here as some of the edges are required exact data
                    //edge line_count is smaller than edge i
                    all_constraints[line_count].push_back(i);
                }
            }
            line_count++;
        }
        pattern_temporal_constraint_file.close();
    }

    void output_all_constraints() {
        for (Edge i=0; i<edge_count; i++) {
            std::cout << "edge " << i << " smaller to:";
                for (auto e : all_constraints[i])
                    std::cout << e << " ";
            std::cout << std::endl;
        }
    }

    void check_all_constraints(PatternGraph& pattern_graph) {
        //TODO: checking cycle
        //TODO: remove redandent checking: e0 < e1, e1 < e6, e0 < e6, and
        //the last checking is redandent.
    }

    void generate_local_opt(PatternGraph& pattern_graph) {
        if (directed) {
            //TODO: for directed, direction info should be added too
            std::cerr << "generate local opt: directed is not supported yet" << std::endl;
            return ;
        }

        local_constraints.clear();
        reset_local_bitsets();
        for (Edge s = 0; s < edge_count; s++) {
            for(Edge g : all_constraints[s]) {
                auto [local_found, a, b, c] = find_local(s, g);
                if (local_found) {
                    insert_local_constraint(a, b, c, true);
                    insert_local_constraint(a, c, b, false);
                    update_local_bitset(a, b, c);
                }
            }
        }
    }

    void generate_non_local_opt(PatternGraph& pattern_graph,
            PatternNonlocalConstraint& pattern_non_local_constraint) {
        //TODO: Jing for compress, we should only check/store edge data
        //that will be compared with and remove those that won't be
        //accessed again.
        //Then it can help work aggregation too.
        auto constraints_embeded = alloc_embeded_sign();

        for (auto input_pattern : pattern_non_local_constraint.input_patterns) {
            std::vector<NonLocalOPT> opts;
            opts.push_back(NonLocalOPT()); //Operator at source
            std::unordered_map<Edge, Edge> visited_edges; //map between: edge, and index in opts

            std::vector<Vertex> path(std::get<1>(input_pattern));
            Vertex source = path[0];
            for (Edge i=1; i<path.size(); i++) {
                Vertex dest = path[i];
                //Checking if there is anything in the visited that cur <source, dest>, <dest, source> edge
                //needs to compare with
                opts.push_back(NonLocalOPT());
                update_opts(i, visited_edges, opts, source, dest, constraints_embeded, pattern_graph);
                source = dest;
            }
            compress_stored(opts);
            non_local_constraints.push_back(opts);
        }

        add_remain_constraints(constraints_embeded, pattern_graph, pattern_non_local_constraint);
    }

    void add_remain_constraints(std::vector<std::vector<bool>> & constraints_embeded,
            PatternGraph & pattern_graph,
            PatternNonlocalConstraint& pattern_non_local_constraint) {
        //pick out non local remained constrains
        std::vector<std::pair<Edge, Edge>> remains_constraints;
        for (Edge i=0; i<constraints_embeded.size(); i++) {
            for (Edge j=0; j<constraints_embeded[i].size(); j++) {
                if (constraints_embeded[i][j] == false) {
                    Edge smaller_edge = i;
                    Edge bigger_edge = all_constraints[i][j];
                    auto [local_found, a, b, c] = find_local(smaller_edge, bigger_edge);
                    if (!local_found)
                        remains_constraints.push_back(std::make_pair(smaller_edge, bigger_edge));
                }
            }
        }

        if (remains_constraints.size() == 0) {
            std::cout << "all temporal constraints are embeded." << std::endl;
            return ;
        } else {
            //TODO: Jing create new non-local constraints to cover remain
            //temporal constraints
            std::cerr << "there are temporal constraints not covered, and it's not supported yet" << std::endl;
            return ;
        }
    }

    void update_opts(Edge step_at, std::unordered_map<Edge, Edge> & visited_edges,
            std::vector<NonLocalOPT> & opts, Vertex source, Vertex dest,
            std::vector<std::vector<bool>> & constraints_embeded, PatternGraph & pattern_graph) {
        //mapping edge to its edge_id
        Edge cur_edge_id = get_edge_id(std::tuple<Vertex, Vertex>(source, dest));
        auto find_edge_in_visited = visited_edges.find(cur_edge_id);
        if (find_edge_in_visited != visited_edges.end()) {
            //If edge is already in visited, then pass
            return ;
        }

        //Check if current edge needs to be smaller than any of the
        //visited edges
        for (Edge i=0; i<all_constraints[cur_edge_id].size(); i++) {
            auto smaller_to = all_constraints[cur_edge_id][i];
            find_edge_in_visited = visited_edges.find(smaller_to);
            if (find_edge_in_visited != visited_edges.end()) {
                opts[find_edge_in_visited->second].set_store();
                opts[step_at].set_compare();
                opts[step_at].insert_compare(find_edge_in_visited->second, true);
                constraints_embeded[cur_edge_id][i] = true;
            }
        }
        //Check if current edge needs to be bigger than any of the visited
        //edges
        for (auto visited_edge : visited_edges) {
            for (Edge i = 0; i<all_constraints[visited_edge.first].size(); i++) {
                auto edge = all_constraints[visited_edge.first][i];
                if (edge == cur_edge_id) {
                    opts[visited_edge.second].set_store();
                    opts[step_at].set_compare();
                    opts[step_at].insert_compare(visited_edge.second, false);
                    //set temporal constrained embeded
                    constraints_embeded[visited_edge.first][i] = true;
                }
            }
        }

        //insert current edge as  visited
        visited_edges.insert( {cur_edge_id, step_at} );
    }

    void compress_stored(std::vector<NonLocalOPT> & opts) {
        std::vector<Edge> compressed_position(opts.size(), 0);
        Edge stored_at = 0;
        for (Edge i=1; i<opts.size(); i++) {
            if (opts[i].store()) {
                opts[i].set_store(stored_at);
                compressed_position[i] = stored_at;
                stored_at++;
            }
            if (opts[i].compare()) {
                for (auto & compare_to_edge : opts[i].compare_to_edges) {
                    Edge org = compare_to_edge.first;
                    compare_to_edge.first = compressed_position[org];
                }
            }
        }
        compressed_position.clear();
    }

    std::vector<std::vector<bool>> alloc_embeded_sign(bool initial=false) {
        std::vector<std::vector<bool>> embeded_sign;
        for (Edge i=0; i<edge_count; i++) {
            embeded_sign.push_back(std::vector<bool>(all_constraints[i].size(), initial));
        }
        return embeded_sign;
    }

    //Generate a map between edge id and vertex pair
    void create_edge_to_pos_map(PatternGraph & pattern_graph) {
        edge_to_id_map.clear();
        id_to_edge_map.clear();

        for (Edge eid=0; eid<(pattern_graph.edge_list).size(); eid++) {
            {
                auto insert_status = edge_to_id_map.insert( {(pattern_graph.edge_list)[eid], (pattern_graph.edge_ID)[eid]} );
                if (!insert_status.second) {
                    std::cerr << "Fail to insert edge-id pair" << std::endl;
                    return ;
                }
            }
            {
                auto insert_status = id_to_edge_map.insert( {(pattern_graph.edge_ID)[eid], (pattern_graph.edge_list)[eid]} );
                if (((insert_status.first) -> first) != (pattern_graph.edge_ID)[eid]) {
                    std::cerr << "Fail to insert id-edge pair" << std::endl;
                    return ;
                }
            }
        }
    }

    void output_pattern_graph_info() {
        std::cout << "Pattern graph info:" << std::endl;
        std::cout << "edge count: " << edge_count << std::endl;
        std::cout << "vertex count: " << vertex_count << std::endl;

        std::cout << "edge to id map:" << std::endl;
        for (auto item : edge_to_id_map) {
            std::cout << std::get<0>(item.first) << ", " << std::get<1>(item.first) << "->" << item.second << std::endl;
        }
        std::cout << "id to edge map:" << std::endl;
        for (auto item : id_to_edge_map) {
            std::cout << item.first<< "->"  << std::get<0>(item.second) << ", " << std::get<1>(item.second) << std::endl;
        }
    }

    Edge get_edge_id(std::tuple<Vertex, Vertex> src_dst) {
        auto find_id = edge_to_id_map.find(src_dst);
        if (find_id == edge_to_id_map.end()) {
            std::cerr << "couldn't find edge id in the edge_to_id_map for edge: " \
                << std::get<0>(src_dst) \
                << std::get<1>(src_dst) \
                << std::endl;
            return 0;
        }
        return find_id -> second;
    }

    std::tuple<Vertex, Vertex> get_edge_from_id(Edge id) {
        auto find_id = id_to_edge_map.find(id);
        if (find_id == id_to_edge_map.end()) {
            std::cerr << "couldn't find edge in the id_to_edge_map for id: " << id << std::endl;
            return std::tuple<Vertex, Vertex>();
        }
        return find_id -> second;
    }

    std::tuple<bool, Vertex, Vertex, Vertex> find_local(Edge s, Edge g) {
        auto [s_from, s_to] = get_edge_from_id(s);
        auto [g_from, g_to] = get_edge_from_id(g);

        Vertex a, b, c;
        bool local_found=false;
        if (s_from == g_from) {
            assert(s_to != g_to);
            a = s_from; b = s_to; c = g_to;
            local_found = true;
        }
        if (s_from == g_to) {
            assert(s_to != g_from);
            a = s_from; b = s_to; c = g_from;
            local_found = true;
        }
        if (s_to == g_from) {
            assert(s_from != g_to);
            a = s_to; b = s_from; c = g_to;
            local_found = true;
        }
        if (s_to == g_to) {
            assert(s_from != g_from);
            a = s_to; b = s_from; c = g_from;
            local_found = true;
        }

        return {local_found, a, b, c};
    }

    void insert_local_constraint(Vertex a, Vertex b, Vertex c, bool b_smaller_than_c) {
        auto found = local_constraints.find(std::make_pair(a, b));
        if (found ==  local_constraints.end()) {
            auto insert_status = local_constraints.insert( {std::make_pair(a, b), std::vector<std::pair<Vertex, bool>>()} );
            if (!insert_status.second) {
                std::cerr << "Fail to insert temporal constraint" << std::endl;
                return ;
            }
            found = insert_status.first;
        }
        (found -> second).push_back(std::make_pair(c, b_smaller_than_c));
    }

    void reset_local_bitsets() {
        local_vertices.reset();
        local_neighbors.clear();
        for (Vertex i=0; i<vertex_count; i++) {
            local_neighbors.push_back(BitSet());
        }
    }

    void update_local_bitset(Vertex a, Vertex b, Vertex c) {
        local_vertices.set(a);
        local_neighbors[a].set(b);
        local_neighbors[a].set(c);
    }

    EdgeData min_edgedata_type() {
        return std::numeric_limits<EdgeData>::min();
    }

    EdgeData max_edgedata_type() {
        return std::numeric_limits<EdgeData>::max();
    }
};
