#include "motif_set_reducer.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <unordered_map>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/misc/interval_tree.h>

#include "triplex_pattern.hpp"

void reduce_motif_set(motif_set_t& output, motif_set_t& input)
{
    typedef typename seqan::Position<motif_t>::Type interval_value_t;
    typedef seqan::Graph<seqan::Undirected<void>> graph_t;
    typedef seqan::VertexDescriptor<graph_t>::Type vertex_descriptor_t;
    typedef seqan::IntervalAndCargo<interval_value_t, vertex_descriptor_t> interval_t;
    typedef typename seqan::Size<graph_t>::Type graph_size_t;
    typedef seqan::IntervalTree<interval_value_t, vertex_descriptor_t> interval_tree_t;
    typedef seqan::Iterator<graph_t, seqan::VertexIterator>::Type vertex_iterator_t;
    typedef seqan::String<vertex_descriptor_t> counts_t;
    typedef seqan::String<graph_size_t> component_t;
    typedef std::unordered_map<vertex_descriptor_t, motif_t> prop_map_t;
    typedef std::unordered_map<graph_size_t, motif_t> comp_map_t;

    if (input.empty()) {
        return;
    } else if (input.size() == 1) {
        output.push_back(std::move(input[0]));
        return;
    }

    graph_t parser;
    prop_map_t prop_map;
    std::vector<interval_t> intervals;

    prop_map.reserve(input.size());
    intervals.resize(input.size());

    unsigned int count = 0;
    for (auto& motif : input) {
        vertex_descriptor_t vtx = seqan::addVertex(parser);

        intervals[count].i1 = static_cast<interval_value_t>(seqan::beginPosition(motif));
        intervals[count].i2 = static_cast<interval_value_t>(seqan::endPosition(motif));
        intervals[count].cargo = vtx;

        prop_map.insert(std::make_pair<graph_size_t, motif_t>(vtx, std::move(motif)));

        count++;
    }

    counts_t tree_results;
    interval_tree_t tree(intervals);
    vertex_iterator_t vertex_it(parser);
    while (!seqan::atEnd(vertex_it)) {
        auto& motif = prop_map.find(*vertex_it)->second;

        seqan::findIntervals(tree_results,
                             tree,
                             seqan::beginPosition(motif),
                             seqan::endPosition(motif));
        for (auto& result : tree_results) {
            if (*vertex_it == result) {
                continue;
            }
            seqan::addEdge(parser, *vertex_it, result);
        }

        seqan::clear(tree_results);
        vertex_it++;
    }

    component_t components;
    graph_size_t num_components = seqan::connectedComponents(components, parser);

    comp_map_t comp_map;
    comp_map.reserve(num_components);

    seqan::goBegin(vertex_it);
    while (!seqan::atEnd(vertex_it)) {
        auto result_ptr = comp_map.find(seqan::getProperty(components, *vertex_it));

        if (result_ptr != comp_map.end()) {
            auto& motif = result_ptr->second;
            merge(motif, prop_map.find(*vertex_it)->second);
        } else {
            auto& motif = prop_map.find(*vertex_it)->second;
            graph_size_t motif_size = seqan::getProperty(components, *vertex_it);
            comp_map.insert(std::make_pair<graph_size_t, motif_t>(std::move(motif_size),
                                                                  std::move(motif)));
        }

        vertex_it++;
    }

    output.reserve(output.size() + num_components);
    std::transform(comp_map.begin(), comp_map.end(), std::back_inserter(output), second(comp_map));
}
