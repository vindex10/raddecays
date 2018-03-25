#ifndef JSON_TYPES_HPP
#define JSON_TYPES_HPP

#include <fifo_map.hpp>
#include <nlohmann/json.hpp>

// Define ordered json read/writer
template<class K, class V, class dummy_compare, class A>
using ordered_map = nlohmann::fifo_map<K, V, nlohmann::fifo_map_compare<K>, A>;
using json = nlohmann::basic_json<ordered_map>;

#endif
