#include "parallel_push_relabel.h"
#include "push_relabel_fifo.h"
#include "command_line_parser.h"
#include "util.h"
#include "graph_loader.h"
#include "measure.h"
#include <fstream>
#include <ios>

std::istream &get_input_stream(std::ifstream &in_file,
                               std::optional<std::string> &file_path) {
    if (file_path.has_value()) {
        in_file.open(*file_path);
        return in_file;
    }
    return std::cin;
}

void print_result(const measurement_result &res, std::string_view solver,
                  std::string_view filename, std::size_t thr_cnt) {
    std::cout << "solver:\t\t" << solver << '\n'
              << "filename:\t" << filename << '\n'
              << "flow:\t\t" << res.max_flow << '\n'
              << "time read:\t" << res.time_read.count() << " ms\n"
              << "time init:\t" << res.time_init.count() << " ms\n"
              << "time solve:\t" << res.time_solve.count() << " ms\n"
              << "# of threads:\t" << thr_cnt << "\n";
}

auto load_graph_and_run(std::istream &is, solver solver,
                        std::size_t &thread_count) {
    auto get_graph_with_cached_edge = [](std::istream &is) {
        auto start = std::chrono::high_resolution_clock::now();
        auto [graph, source, sink] = load_graph(is);
        set_reverse_edge_cap(graph);
        auto end = std::chrono::high_resolution_clock::now();
        return std::make_tuple(
            graph, source, sink,
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start));
    };

    measurement_result result{};

    switch (solver) {
    case solver::prf: {
        auto [graph, source, sink, time_read] = get_graph_with_cached_edge(is);
        result = measure_sequential<push_relabel_fifo::max_flow_instance>(
            std::move(graph), source, sink);
        result.time_read = time_read;
        thread_count = 1;
        break;
    }
    case solver::ppr: {
        auto [graph, source, sink, time_read] = get_graph_with_cached_edge(is);
        result = measure_parallel<parallel_push_relabel::max_flow_instance>(
            std::move(graph), source, sink, thread_count);
        result.time_read = time_read;
        break;
    }
    default:
        throw std::logic_error("Unknown solver");
    }
    return result;
}

int main(int argc, char *argv[]) {
    std::ios_base::sync_with_stdio(false);
    command_line_parser parser;
    if (!parser.parse_arguments(argc, argv))
        return 1;
    auto solver = parser.get_solver();
    auto solver_str = parser.get_solver_str();
    auto file_path = parser.get_filename();
    auto thr_cnt = parser.getnthreads();

    std::ifstream in_file;
    auto &stream = get_input_stream(in_file, file_path);
    if (!stream) {
        std::cerr << "Unable to open file " << *file_path << '\n';
        return 1;
    }

    // get filename without full path, or stdin
    std::string filename = "stdin";
    if (file_path.has_value()) {
        filename = (*file_path).substr((*file_path).find_last_of("/\\") + 1);
    }

    // run
    auto result = load_graph_and_run(stream, solver, thr_cnt);
    print_result(result, solver_str, filename, thr_cnt);
    return 0;
}
