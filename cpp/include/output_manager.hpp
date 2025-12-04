#pragma once

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <sstream>

namespace flux {
inline std::string parse_output_from_argv(int argc, char **argv) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        // --output=name
        if (arg.rfind("--output=", 0) == 0) { return arg.substr(9); }

        // --output name
        if (arg == "--output") {
            if (i + 1 < argc) { return std::string(argv[i + 1]); }
            return "";
        }
    }
    return "";
}

class OutputManager {
public:
    OutputManager(const std::string &user_path, const std::string &label_name) {
        m_base_dir = user_path.empty() ? "outputs" : user_path;
        std::filesystem::create_directories(m_base_dir);

        m_cur_dir = (std::filesystem::path(m_base_dir)
                     / generate_timestamp_dir(label_name))
                        .string();
        std::filesystem::create_directories(m_cur_dir);
    }

    std::filesystem::path get_file_path(const std::string &filename,
                                        bool make_sure_exist) const {
        // refuse dangerous paths
        if (filename.find("..") != std::string::npos) {
            throw std::runtime_error(
                "OutputManager::get_file_path: filename contains '..'");
        }

        std::filesystem::path p = std::filesystem::path(m_cur_dir) / filename;

        if (make_sure_exist) {
            // create parent directory
            try {
                std::filesystem::create_directories(p.parent_path());
            }
            catch (const std::exception &e) {
                throw std::runtime_error(
                    std::string("OutputManager::get_file_path: "
                                "cannot create parent directory: ")
                    + e.what());
            }
        }

        return p;
    }

    std::filesystem::path operator/(const std::string &filename) const {
        return get_file_path(filename, true);
    }

private:
    std::string m_base_dir;
    std::string m_cur_dir;

    static std::string generate_timestamp_dir(const std::string &label_name) {
        auto now = std::chrono::system_clock::now();
        auto t = std::chrono::system_clock::to_time_t(now);
        std::tm tm{};
#ifdef _WIN32
        localtime_s(&tm, &t);
#else
        localtime_r(&t, &tm);
#endif

        const std::string base_name =
            label_name.empty() ? "" : (label_name + '-');

        std::ostringstream oss;
        oss << base_name << std::put_time(&tm, "%Y%m%d-%H%M%S");
        return oss.str();
    }
};
}  // namespace flux
