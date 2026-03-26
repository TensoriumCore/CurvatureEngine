#include "app/CommandLine.h"

#include "app/Problems.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

namespace curvatureengine::app {
namespace {

bool parse_spin(const char *arg, float &spin) {
  if (arg == nullptr) {
    return false;
  }

  char *end = nullptr;
  const float parsed = std::strtof(arg, &end);
  if (end == arg || (end != nullptr && *end != '\0') || !std::isfinite(parsed)) {
    return false;
  }

  spin = parsed;
  return true;
}

bool is_supported_riemann_metric(const std::string &metric) {
  return metric == "schwarzschild" || metric == "kerr" ||
         metric == "kerr-newman" || metric == "kds";
}

} // namespace

void print_usage(std::FILE *stream) {
  std::fprintf(stream, "Usage:\n");
  std::fprintf(stream, "  CurvatureEngine -G <spin>\n");
  std::fprintf(stream, "  CurvatureEngine -L <spin>\n");
  std::fprintf(stream, "  CurvatureEngine -M <spin>\n");
  std::fprintf(stream, "  CurvatureEngine -S <spin>\n");
  std::fprintf(stream, "  CurvatureEngine -C <spin>\n");
  std::fprintf(stream, "  CurvatureEngine -R <spin> <metric>\n");
  std::fprintf(stream, "\nOptions:\n");
  std::fprintf(stream, "  -G  Timelike geodesic calculation\n");
  std::fprintf(stream, "  -L  Light geodesic calculation\n");
  std::fprintf(stream, "  -M  Metric tensor calculation\n");
  std::fprintf(stream, "  -S  Black hole shadow generation\n");
  std::fprintf(stream, "  -C  ADM/BSSN solver\n");
  std::fprintf(stream, "  -R  Riemann tensor calculation\n");
  std::fprintf(stream, "\nRiemann metrics:\n");
  std::fprintf(stream, "  schwarzschild (requires spin 0)\n");
  std::fprintf(stream, "  kerr\n");
  std::fprintf(stream, "  kerr-newman\n");
  std::fprintf(stream, "  kds\n");
}

ParseResult parse_command_line(int argc, char **argv) {
  ParseResult result{};

  if (argc <= 1) {
    result.exit_code = 1;
    result.error = "Missing command line arguments.";
    return result;
  }

  if (std::strcmp(argv[1], "-h") == 0 || std::strcmp(argv[1], "--help") == 0) {
    result.ok = false;
    result.print_usage = true;
    result.exit_code = 0;
    return result;
  }

  if (argc < 3) {
    result.exit_code = 1;
    result.error = "Missing spin argument.";
    return result;
  }

  if (!parse_spin(argv[2], result.command.spin)) {
    result.exit_code = 1;
    result.error = "Invalid spin value.";
    return result;
  }

  if (std::strcmp(argv[1], "-G") == 0) {
    result.command.kind = CommandKind::Geodesics;
  } else if (std::strcmp(argv[1], "-L") == 0) {
    result.command.kind = CommandKind::LightGeodesics;
  } else if (std::strcmp(argv[1], "-M") == 0) {
    result.command.kind = CommandKind::Metric;
  } else if (std::strcmp(argv[1], "-S") == 0) {
    result.command.kind = CommandKind::Shadow;
  } else if (std::strcmp(argv[1], "-C") == 0) {
    result.command.kind = CommandKind::Bssn;
  } else if (std::strcmp(argv[1], "-R") == 0) {
    result.command.kind = CommandKind::Riemann;
    if (argc < 4) {
      result.exit_code = 1;
      result.error = "Riemann mode requires a metric argument.";
      return result;
    }
    result.command.metric = argv[3];
    if (!is_supported_riemann_metric(result.command.metric)) {
      result.exit_code = 1;
      result.error = "Unsupported Riemann metric.";
      return result;
    }
    if (result.command.metric == "schwarzschild" &&
        result.command.spin != 0.0f) {
      result.exit_code = 1;
      result.error = "Schwarzschild requires spin 0.";
      return result;
    }
  } else {
    result.exit_code = 1;
    result.error = "Unknown option.";
    return result;
  }

  result.ok = true;
  result.print_usage = false;
  result.exit_code = 0;
  return result;
}

int run_command(const CommandLine &command) {
  switch (command.kind) {
  case CommandKind::Geodesics:
    return Geodesics_prob();
  case CommandKind::LightGeodesics:
    return light_geodesics_prob();
  case CommandKind::Riemann:
    return Riemann_tensor(command.metric.c_str());
  case CommandKind::Metric:
    return Metric_prob();
  case CommandKind::Shadow:
    return shadow_prob();
  case CommandKind::Bssn:
    return grid_setup();
  }

  return 1;
}

} // namespace curvatureengine::app
