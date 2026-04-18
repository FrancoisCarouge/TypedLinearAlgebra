/* Typed Linear Algebra
Version 0.2.0
https://github.com/FrancoisCarouge/TypedLinearAlgebra

SPDX-License-Identifier: Unlicense

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <https://unlicense.org> */

//! @file
//! @brief Benchmark result visualization tool.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <matplot/matplot.h>

//! @brief Main entry point for benchmark result visualization.
//! @details Reads benchmark results from results.txt, parses the data,
//!          and generates a visualization plot showing performance metrics
//!          across different matrix sizes and methods. Saves the plot as
//!          plot.png.
//! @return EXIT_SUCCESS on successful execution, EXIT_FAILURE on error.
int main() {
#if defined(_WIN32)
  constexpr auto gnuplot_bin{R"(C:\Program Files\gnuplot\bin)"};
  const auto gnuplot_exe{std::filesystem::path{gnuplot_bin} / "gnuplot.exe"};
  if (std::filesystem::exists(gnuplot_exe)) {
    const auto current_path{std::getenv("PATH")};
    auto path{current_path != nullptr ? std::string{current_path}
                                      : std::string{}};
    const auto gnuplot_dir{std::string{gnuplot_bin}};
    if (path.find(gnuplot_dir) == std::string::npos) {
      if (!path.empty() && path.back() != ';') {
        path += ';';
      }
      path += gnuplot_dir;
      _putenv_s("PATH", path.c_str());
    }
  }
#endif

  std::ifstream file{"results.txt"};
  if (!file.is_open()) {
    return EXIT_FAILURE;
  }

  std::string line;
  std::map<std::string, std::vector<std::pair<double, double>>> series_data;

  auto trim{[](std::string &value) {
    value.erase(value.begin(),
                std::find_if(value.begin(), value.end(), [](unsigned char ch) {
                  return !std::isspace(ch);
                }));
    value.erase(std::find_if(value.rbegin(), value.rend(),
                             [](unsigned char ch) { return !std::isspace(ch); })
                    .base(),
                value.end());
  }};

  while (std::getline(file, line)) {
    if (line.empty()) {
      continue;
    }

    std::vector<std::string> tokens;
    std::size_t start{0};
    while (start < line.size()) {
      const auto pos{line.find('|', start)};
      std::string field;
      if (pos == std::string::npos) {
        field = line.substr(start);
        start = line.size();
      } else {
        field = line.substr(start, pos - start);
        start = pos + 1;
      }
      trim(field);
      if (!field.empty()) {
        tokens.push_back(std::move(field));
      }
    }

    if (tokens.size() >= 3) {
      const auto method{tokens[0]};
      auto size_str{tokens[1]};
      const auto time_str{tokens[2]};
      auto time{std::stod(time_str) * 1e9};

      const auto x_pos{size_str.find('x')};
      if (x_pos != std::string::npos) {
        size_str = size_str.substr(0, x_pos);
      }
      const auto size{std::stod(size_str)};

      series_data[method].emplace_back(size, time);
    }
  }

  auto figure{matplot::figure()};
  figure->size(1200, 800);

  std::vector<std::string> methods;
  methods.reserve(series_data.size());
  std::ranges::transform(series_data, std::back_inserter(methods),
                         [](const auto &point) { return point.first; });

  std::sort(methods.begin(), methods.end());

  for (const auto &method : methods) {
    auto points{series_data[method]};
    std::sort(points.begin(), points.end());

    std::vector<double> sorted_sizes;
    std::vector<double> sorted_times;
    sorted_sizes.reserve(points.size());
    sorted_times.reserve(points.size());
    std::ranges::transform(points, std::back_inserter(sorted_sizes),
                           [](const auto &point) { return point.first; });
    std::ranges::transform(points, std::back_inserter(sorted_times),
                           [](const auto &point) { return point.second; });

    auto line_plot{matplot::plot(sorted_sizes, sorted_times, "o-")};
    line_plot->display_name(method);
    line_plot->line_width(2.5);
    line_plot->marker_size(8);
    matplot::hold(matplot::on);
  }

  matplot::xlabel("Matrix Size NxN");
  matplot::ylabel("Time (ns)");
  matplot::title("Matrix-Matrix Product Benchmark");
  matplot::legend()->location(matplot::legend::general_alignment::bottomright);
  matplot::grid(matplot::on);
  matplot::gca()->x_axis().scale(matplot::axis_type::axis_scale::log);
  matplot::gca()->y_axis().scale(matplot::axis_type::axis_scale::log);
  matplot::xlim({1.0, 128.0});
  matplot::xticks({1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0});
  matplot::save("plot.png");

  return EXIT_SUCCESS;
}
