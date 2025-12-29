#pragma once

#include <vector>
#include <cstddef>
namespace utils
{

    class AdjacencyList
    {
    public:
        class edge_type
        {
        public:
            size_t to;
            unsigned weight;

            edge_type(size_t to, unsigned weight)
                    : to(to)
                    , weight(weight)
            {}
        };

        class reverse_edge_type
        {
        public:
            size_t from;
            unsigned weight;

            reverse_edge_type(size_t from, unsigned weight)
                    : from(from)
                    , weight(weight)
            {}
        };

        using node_type = std::vector<edge_type>;
        using reverse_node_type = std::vector<reverse_edge_type>;

        std::vector<node_type> nodes;

        explicit AdjacencyList(size_t nodes_number = 0);

        void add_edge(size_t from, size_t to, unsigned weight);
        unsigned get_weight(size_t from, size_t to) const noexcept;

        std::vector<size_t> get_longest_path() const;
        std::vector<size_t> topological_sort() const;

    private:
        size_t _node_num;
        // size_t _edge_num;
    };

    // assume that nodes with higher index cannot have an edge linking to a node with lower index
    template <typename MatrixType>
    std::vector<size_t> longest_path_of(const MatrixType& graph, size_t node_num)
    {
        auto distance    = new size_t[node_num]();
        auto path_record = new size_t[node_num];

        for (size_t to = 1; to != node_num; ++to)
            for (size_t station = 0; station != to; ++station)
                if (graph[station][to])
                {
                    unsigned curr = distance[station] + graph[station][to];
                    if (curr > distance[to])
                    {
                        distance[to] = curr;
                        path_record[to] = station;
                    }
                }

        size_t end_node = 0;
        for (size_t i = 1; i != node_num; ++i)
            if (distance[end_node] < distance[i]) end_node = i;

        std::vector<size_t> path;
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

}
