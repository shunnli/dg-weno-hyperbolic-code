#include "output_manager.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char **argv) {
    auto output_path = OutputManager::parse_output(argc, argv);
    const OutputManager om(output_path.value_or("../outputs"), "cpp-test");

    std::cout << "Output base: " << om.base_dir << "\n";
    std::cout << "Output: " << om.cur_dir << "\n";

    // --- Numerical data ---
    std::vector<double> x(200);
    std::vector<double> y(200);
    for (int i = 0; i < 200; i++) {
        x[i] = i * 2 * M_PI / 199.0;
        y[i] = std::sin(x[i]);
    }

    std::string data_file = om.get_file_path("data/sin_wave.txt");
    std::ofstream fout(data_file);
    for (int i = 0; i < 200; i++) fout << x[i] << " " << y[i] << "\n";
    fout.close();

    std::cout << "Save data to: " << data_file << "\n";

    // --- Log file ---
    std::string log_file = om.get_file_path("logs/log.txt");
    std::ofstream log(log_file);
    log << "CPP Demo\n";
    log << "Generated sine-wave data with 200 points.\n";
    log.close();

    std::cout << "Write log to: " << log_file << "\n";

    return 0;
}
