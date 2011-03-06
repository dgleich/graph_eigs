/**
 * @file permutation.hpp
 * A simple code to sort a vector and return the permutation vector.
 */

/**
 * History
 * -------
 * :2011-03-05: Initial version based on overlapping/hypercluster.cc
 *   and the new sparfun/sf_perm.cc
 */ 

#include <vector>
#include <algorithm>

template <typename T, typename index_type = int>
class sort_order_comparison {
  const std::vector<T>& items;
public:
  sort_order_comparison(const std::vector<T>& i) : items(i) {}
  bool operator() (index_type i, index_type j) {
    return items[i] > items[j];
  }
};

/**
 * This function sorts in descending order
 * @param order this parameter is an output and will be initialized
 */
template <typename T>
void sort_permutation(const std::vector<T>& a, std::vector<int>& order) {
  order.resize(a.size());
  for (size_t i=0; i<a.size(); ++i) {
    order[i] = (int)i;
  }
  sort_order_comparison<T> comp(a);
  std::sort( order.begin(), order.end(), comp );
}
