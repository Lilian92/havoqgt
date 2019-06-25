#pragma once

#include <unordered_map>

namespace prunejuice {

template<typename Vertex, typename VertexData, typename BitSet, typename EdgeData>
class vertex_state {
  public :
    vertex_state() :
    vertex_pattern_index(0),
    is_active(false) {}
  
    BitSet template_vertices;
    BitSet template_neighbors;
    //std::unordered_map<VertexData, IntegralType>
    //  template_neighbor_metadata_count_map;
 
    size_t vertex_pattern_index; // TODO: dummy, to be removed
    bool is_active; 		  

    std::unordered_map<Vertex, std::tuple<EdgeData, EdgeData>> last_itr_min_max,
    std::unordered_map<Vertex, std::tuple<EdgeData, EdgeData>> cur_itr_min_max,
};
  
} // end namespace prunejuice 
