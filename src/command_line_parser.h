#ifndef MAXFLOW_COMMAND_LINE_PARSER_H
#define MAXFLOW_COMMAND_LINE_PARSER_H

#include "util.h"
#include <cstring>
#include <exception>
#include <iostream>
#include <omp.h>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>

void print_usage(std::string_view program_name) {
    std::cerr << "usage: " << program_name
              << " <solver> [-f <path>] [-p <number>]\n\n";
    std::cerr << "list of possible solvers: [prf, ppr]\n\n";
    std::cerr << "prf:\tpush-relabel algorithm with FIFO vertex selection\n"
              << "ppr:\tparallel push-relabel algorithm\n"
              << "[-f <path>]:\tspecify path to a maxflow problem instance in "
                 "DIMACS format. Reads from stdin if omitted.\n";
    std::cerr
        << "[-p <number>]:\tspecify max number of threads for parallel "
           "solvers. "
           "Ignored for sequential solvers. Default = number of hw threads.\n";
}

enum class solver { prf, ppr };

class command_line_parser {
    solver _solver;
    std::string _solver_str;
    std::string _input_filename;
    std::size_t nthreads = 0;

    std::unordered_map<std::string, solver> solver_map = {{"prf", solver::prf},
                                                          {"ppr", solver::ppr}};

  public:
    bool parse_arguments(int argc, char **argv) {
        if (argc < 2) {
            std::cerr << argv[0] << ": missing solver\n";
            print_usage(argv[0]);
            return false;
        }

        if (argc % 2 || argc > 6) {
            std::cerr << argv[0] << ": invalid argument count\n";
            print_usage(argv[0]);
            return false;
        }

        if (!strcmp(argv[1], "--help")) {
            print_usage(argv[0]);
            return false;
        }

        auto iter = solver_map.find(argv[1]);
        if (iter == std::end(solver_map)) {
            std::cerr << argv[0] << ": invalid solver\n";
            print_usage(argv[0]);
            return false;
        }
        _solver = iter->second;
        _solver_str = argv[1];

        for (int i = 2; i < argc; i += 2) {
            if (!strcmp(argv[i], "-f"))
                _input_filename = argv[i + 1];
            else if (!strcmp(argv[i], "-p")) {
                try {
                    nthreads = std::stoul(argv[i + 1]);
                    if (nthreads == 0)
                        throw std::invalid_argument(
                            "positive integer expected");
                } catch (...) {
                    std::cerr << argv[i + 1]
                              << ": invalid format, expected positive unsigned "
                                 "integer\n";
                    return false;
                }
            } else {
                std::cerr << "Unknown option: " << argv[i] << '\n';
                return false;
            }
        }

        return true;
    }

    std::optional<std::string> get_filename() const noexcept {
        return _input_filename.empty()
                   ? std::nullopt
                   : std::optional<std::string>{_input_filename};
    }

    std::size_t getnthreads() const noexcept {
        return nthreads == 0 ? omp_get_max_threads() : nthreads;
    }

    solver get_solver() const noexcept { return _solver; }

    std::string get_solver_str() const noexcept { return _solver_str; }
};

#endif // MAXFLOW_COMMAND_LINE_PARSER_H
