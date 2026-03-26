#include "app/CommandLine.h"
#include "app/RuntimeState.h"

#include <cstdio>

int main(int argc, char **argv) {
  const auto parsed = curvatureengine::app::parse_command_line(argc, argv);
  if (!parsed.ok) {
    if (!parsed.error.empty()) {
      std::fprintf(stderr, "%s\n\n", parsed.error.c_str());
    }
    if (parsed.print_usage) {
      curvatureengine::app::print_usage(parsed.exit_code == 0 ? stdout : stderr);
    }
    return parsed.exit_code;
  }

  a = parsed.command.spin;
  return curvatureengine::app::run_command(parsed.command);
}
