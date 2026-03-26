#pragma once

#include <cstdio>
#include <string>

namespace curvatureengine::app {

enum class CommandKind {
  Geodesics,
  LightGeodesics,
  Riemann,
  Metric,
  Shadow,
  Bssn,
};

struct CommandLine {
  CommandKind kind = CommandKind::Geodesics;
  float spin = 0.0f;
  std::string metric;
};

struct ParseResult {
  bool ok = false;
  bool print_usage = true;
  int exit_code = 1;
  CommandLine command{};
  std::string error;
};

ParseResult parse_command_line(int argc, char **argv);
void print_usage(std::FILE *stream);
int run_command(const CommandLine &command);

} // namespace curvatureengine::app
