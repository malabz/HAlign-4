#include "Graph.hpp"
#include <limits>
#include <stdexcept>


utils::AdjacencyList::AdjacencyList(size_t nodes_number)
        : nodes(nodes_number)
        , _node_num(nodes_number)
// , _edge_num(0)
{}

void utils::AdjacencyList::add_edge(size_t from, size_t to, unsigned weight)
{
    nodes[from].emplace_back(to, weight);
    // ++_edge_num;
}

std::vector<size_t> utils::AdjacencyList::topological_sort() const
{
    size_t *indegree = new size_t[_node_num]();
    for (size_t i = 0; i != _node_num; ++i)
        for (auto edge : nodes[i])
            ++indegree[edge.to];

    std::vector<size_t> topological_sequence;
    topological_sequence.reserve(nodes.size());
    for (size_t pending = _node_num; pending; )
        for (size_t i = 0; i != _node_num; ++i)
            if (indegree[i] == 0)
            {
                indegree[i] = std::numeric_limits<size_t>::max();
                topological_sequence.emplace_back(i);
                --pending;

                for (auto edge : nodes[i]) --indegree[edge.to];
            }

    delete[] indegree;
    return topological_sequence;
}

std::vector<size_t> utils::AdjacencyList::get_longest_path() const
{
    std::vector<size_t> path;
    if (_node_num == 0) return path;

    std::vector<reverse_node_type> reverse_nodes(_node_num);
    for (size_t from = 0; from != _node_num; ++from)
        for (auto edge : nodes[from])
            reverse_nodes[edge.to].emplace_back(from, edge.weight);

    auto topological_sequence = topological_sort();

    auto distance    = new size_t[_node_num]();
    auto path_record = new size_t[_node_num];

    for (size_t to = 1; to != _node_num; ++to)
    {
        const size_t target = topological_sequence[to];
        for (auto edge : reverse_nodes[target])
        {
            unsigned challenger = distance[edge.from] + edge.weight;
            if (challenger > distance[target])
            {
                distance[target] = challenger;
                path_record[target] = edge.from;
            }
        }
    }

    size_t end_node = 0;
    for (size_t i = 1; i != _node_num; ++i)
        if (distance[end_node] < distance[i]) end_node = i;

    for (size_t i = end_node; i != 0; i = path_record[i])
        path.emplace_back(i);

    for (size_t i = 0, mid = path.size() >> 1; i != mid; ++i)
        std::swap(path[i], path[path.size() - i - 1]);

    // std::copy(distance, distance + node_num, std::ostream_iterator<size_t>(std::cout, ", ")); std::cout << '\n';
    // std::copy(path_record, path_record + node_num, std::ostream_iterator<size_t>(std::cout, ", ")); std::cout << '\n';
    // std::copy(path.cbegin(), path.cend(), std::ostream_iterator<size_t>(std::cout, ", ")); std::cout << '\n';

    delete[] distance;
    delete[] path_record;

    return path;
}

unsigned utils::AdjacencyList::get_weight(size_t from, size_t to) const noexcept
{
    for (auto edge : nodes[from])
        if (edge.to == to) return edge.weight;
}
