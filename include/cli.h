#ifndef PAM_ST_CLI_H
#define PAM_ST_CLI_H

#include "structs.h"

#include <string>
#include <vector>

struct CliOptions {
    std::vector<std::string> inputs;  // one per sample
    std::string input;                // the inputs joined with ';', for reports
    std::string output_dir;
    std::string output_file;
    AnalysisConfig analysis;
};

void usage(const char* program);
CliOptions parse_args(int argc, char** argv);

#endif  // PAM_ST_CLI_H
