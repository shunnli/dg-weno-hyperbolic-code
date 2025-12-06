#pragma once
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

class OutputManager {
public:
    std::filesystem::path base_dir;
    std::filesystem::path cur_dir;

    OutputManager(const std::string &output_path,
                  const std::string &label_name) {
        namespace fs = std::filesystem;

        base_dir = fs::path(output_path);

        if (!fs::exists(base_dir)) {
            try {
                fs::create_directories(base_dir);
            }
            catch (const std::exception &e) {
                std::cerr << "[Warning] Failed to create base_dir: " << base_dir
                          << "\nReason: " << e.what() << "\n";
                base_dir = fs::current_path();
            }
        }

        // Timestamp
        auto ts = timestamp();
        std::string dirname = label_name.empty() ? ts : (label_name + "-" + ts);

        cur_dir = base_dir / dirname;

        try {
            fs::create_directories(cur_dir);
        }
        catch (const std::exception &e) {
            std::cerr << "[Warning] Failed to create cur_dir: " << cur_dir
                      << "\nReason: " << e.what()
                      << "\nFallback to ./timestamp\n";

            cur_dir = fs::current_path() / dirname;

            try {
                fs::create_directories(cur_dir);
            }
            catch (const std::exception &e2) {
                std::cerr << "[Warning] Failed fallback dir: " << cur_dir
                          << "\nReason: " << e2.what()
                          << "\nFallback to pwd()\n";
                cur_dir = fs::current_path();
            }
        }
    }

    explicit OutputManager(const std::string &output_path)
        : OutputManager(output_path, "") {}

    OutputManager() : OutputManager("outputs", "") {}

    // Get full file path and ensure parent directory exists
    std::filesystem::path get_file_path(const std::string &filename,
                                        bool make_sure = true) const {
        namespace fs = std::filesystem;

        if (filename.find("..") != std::string::npos) {
            std::cerr << "[Warning] filename contains '..'\n";
        }

        fs::path p = cur_dir / filename;

        if (make_sure) {
            fs::path parent = p.parent_path();
            if (!parent.empty() && !fs::exists(parent)) {
                try {
                    fs::create_directories(parent);
                }
                catch (const std::exception &e) {
                    std::cerr
                        << "[Warning] Cannot create parent dir: " << parent
                        << "\nReason: " << e.what() << "\n";
                    // fallback: place file inside cur_dir
                    p = cur_dir / fs::path(filename).filename();
                }
            }
        }

        return p;
    }

    // om / "abc.txt"  ->  get_file_path("abc.txt", true)
    std::filesystem::path operator/(const std::string &filename) const {
        return get_file_path(filename, true);
    }

    static std::optional<std::string> parse_output(int argc, char **argv) {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];

            // --output=name
            if (arg.rfind("--output=", 0) == 0) { return arg.substr(9); }

            // --output name
            if (arg == "--output") {
                if (i + 1 < argc) { return std::string(argv[i + 1]); }
                return std::nullopt;
            }
        }
        return std::nullopt;
    }

private:
    // Timestamp formatter: yyyyMMdd-HHmmss
    static std::string timestamp() {
        auto now = std::chrono::system_clock::now();
        std::time_t t = std::chrono::system_clock::to_time_t(now);

        std::tm tm{};
#ifdef _WIN32
        localtime_s(&tm, &t);
#else
        localtime_r(&t, &tm);
#endif
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y%m%d-%H%M%S");
        return oss.str();
    }
};
