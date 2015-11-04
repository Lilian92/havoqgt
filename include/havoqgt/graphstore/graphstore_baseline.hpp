#ifndef GRAPHSTORE_BASELINE_HPP
#define GRAPHSTORE_BASELINE_HPP

#include <boost/unordered_map.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/graphstore/graphstore_utilities.hpp>


namespace graphstore {

template <typename vertex_type,
          typename vertex_property_type,
          typename edge_property_type,
          typename segment_manager_type>
class graphstore_baseline
{
 private:

  using size_type             = size_t;
  template <typename t1, typename t2>
  using pair_type             = graphstore::utility::packed_pair<t1, t2>;
  using allocator_type        = boost::interprocess::allocator<void, segment_manager_type>;


  /// adjacency list
  using edge_vec_element_type = pair_type<vertex_type, edge_property_type>;
  using vec_allocator_type    = boost::interprocess::allocator<edge_vec_element_type, segment_manager_type>;
  using edge_vec_type         = boost::interprocess::vector<edge_vec_element_type, vec_allocator_type>;

  /// vertex table
  using map_value_type        = pair_type<vertex_property_type, edge_vec_type>;
  using map_element_type      = std::pair<const vertex_type, map_value_type>;
  using map_allocator_type    = boost::interprocess::allocator<map_element_type, segment_manager_type>;
  using map_table_type        = boost::unordered_map<vertex_type,
                                                     map_value_type,
                                                     boost::hash<vertex_type>,
                                                     std::equal_to<vertex_type>,
                                                     map_allocator_type>;

  enum : size_t {
    kEmptyValue = std::numeric_limits<vertex_type>::max()
  };


 public:

  explicit graphstore_baseline(segment_manager_type* segment_manager) :
    m_allocator(segment_manager),
    m_map_table(m_allocator)
  {}


  bool insert_edge(const vertex_type& src, const vertex_type& trg, const edge_property_type& edge_property)
  {
    auto value = m_map_table.find(src);
    if (value == m_map_table.end()) { // new vertex
      const vertex_property_type dummy_vp = false;
      map_value_type map_value(dummy_vp,
                               edge_vec_type(1, edge_vec_element_type(trg, edge_property), m_allocator));
      m_map_table.insert(map_element_type(src, map_value));
    } else {
      edge_vec_type& edge_vec = value->second.second;
      edge_vec_element_type trg_edge(trg, edge_property);
      // --- find duplicated edge --- //
      for (const auto edge : edge_vec) {
        if (edge.first ==  trg) {
          return false; /// found duplicated edge
        }
      }
      /// --- find blank space and insert it into there if found --- ///
      for (auto edge : edge_vec) {
        if (edge.first == kEmptyValue) {
          edge.first = trg;
          edge.second = edge_property;
          return true;
        }
      }

      /// --- since there is no empty spase, add it at end --- ///
      edge_vec.push_back(edge_vec_element_type(trg, edge_property));
    }

    return true;
  }

  size_type erase_edge(vertex_type& src, vertex_type& trg)
  {
  }

  inline vertex_property_type& vertex_meta_data(const vertex_type& vertex)
  {
    return m_map_table.find(vertex)->second.first;
  }

  typename edge_vec_type::iterator adjacencylist(const vertex_type& src)
  {
    const auto itr = m_map_table.find(src);
    edge_vec_type& edge_vec = itr->second.second;
    return edge_vec.begin();
  }

  typename edge_vec_type::iterator adjacencylist_end(const vertex_type& src)
  {
    const auto itr = m_map_table.find(src);
    edge_vec_type& edge_vec = itr->second.second;
    return edge_vec.end();
  }


  void clear()
  {

  }

  void print_status()
  {

  }

  void fprint_all_elements(std::ofstream& of)
  {
  }

  allocator_type m_allocator;
  map_table_type m_map_table;
};

}
#endif // GRAPHSTORE_BASELINE_HPP
