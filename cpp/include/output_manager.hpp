#pragma once

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <ranges>
#include <sstream>
#include <string>

class OutputManager {
public:
    std::filesystem::path base_dir;
    std::filesystem::path cur_dir;

    OutputManager(const std::string &output_path,
                  const std::string &label_name) {
        namespace fs = std::filesystem;

        // base_dir: path-like, allow relative, allow disk prefix
        base_dir = fs::path(sanitize_path(output_path));

        if (!fs::exists(base_dir)) {
            try {
                fs::create_directories(base_dir);
            }
            catch (const std::exception &e) {
                std::cerr << "[Warning] Failed to create base_dir: " << base_dir
                          << "\nReason: " << e.what()
                          << "\nFallback to current_path()\n";
                base_dir = fs::current_path();
            }
        }

        base_dir = fs::absolute(base_dir);

        // dirname: must be a pure directory name
        std::string ts = timestamp();
        std::string dirname = sanitize_dirname(label_name);
        dirname = dirname.empty() ? ts : (dirname + "-" + ts);

        cur_dir = base_dir / dirname;

        try {
            fs::create_directories(cur_dir);
        }
        catch (const std::exception &e) {
            std::cerr << "[Warning] Failed to create cur_dir: " << cur_dir
                      << "\nReason: " << e.what()
                      << "\nFallback to current_path()/timestamp\n";

            cur_dir = fs::current_path() / ts;

            try {
                fs::create_directories(cur_dir);
            }
            catch (const std::exception &e2) {
                std::cerr << "[Warning] Failed fallback dir: " << cur_dir
                          << "\nReason: " << e2.what()
                          << "\nFallback to current_path()\n";
                cur_dir = fs::current_path();
            }
        }

        cur_dir = fs::absolute(cur_dir);
    }

    explicit OutputManager(const std::string &output_path)
        : OutputManager(output_path, "") {}

    OutputManager() : OutputManager("outputs", "") {}

    std::filesystem::path get_file_path(const std::string &filename,
                                        bool make_sure) const {
        namespace fs = std::filesystem;

        fs::path fname = fs::path(sanitize_filename(filename));
        fs::path p = (cur_dir / fname).lexically_normal();

        // safety: prevent escaping cur_dir
        if (!is_under_directory(p, cur_dir)) {
            std::cerr << "[Warning] filename escapes output directory: "
                      << filename << "\nFallback to basename\n";
            p = cur_dir / fname.filename();
        }

        if (make_sure) {
            fs::path parent = p.parent_path();
            if (!parent.empty() && !fs::exists(parent)) {
                try {
                    fs::create_directories(parent);
                }
                catch (const std::exception &e) {
                    std::cerr
                        << "[Warning] Cannot create parent dir: " << parent
                        << "\nReason: " << e.what()
                        << "\nFallback to cur_dir\n";
                    p = cur_dir / fname.filename();
                }
            }
        }

        return p;
    }

    std::filesystem::path get_file_path(const std::string &filename) const {
        return get_file_path(filename, true);
    }

    std::filesystem::path operator/(const std::string &filename) const {
        return get_file_path(filename, true);
    }

    static std::optional<std::string> parse_output(int argc, char **argv) {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];

            if (arg.starts_with("--output=")) { return arg.substr(9); }

            if (arg == "--output") {
                if (i + 1 < argc) { return std::string(argv[i + 1]); }
                return std::nullopt;
            }
        }
        return std::nullopt;
    }

private:
    static std::string trim(const std::string &s) {
        auto first = std::ranges::find_if(
            s, [](unsigned char c) { return !std::isspace(c); });
        auto last = std::ranges::find_if(
                        std::ranges::reverse_view(s),
                        [](unsigned char c) { return !std::isspace(c); })
                        .base();
        return (first < last) ? std::string(first, last) : std::string();
    }

    // For path-like input (output_path)
    static std::string sanitize_path(const std::string &s) {
        std::string out = trim(s);
#ifdef _WIN32
        std::ranges::replace(out, '/', '\\');
#endif
        return out;
    }

    // For filename (may include subdirs, but no illegal chars)
    static std::string sanitize_filename(const std::string &s) {
        std::string out = sanitize_path(s);
#ifdef _WIN32
        std::string invalid = "<>:\"|?*";
        for (char c : invalid) { std::ranges::replace(out, c, '_'); }
#endif
        return out;
    }

    // For directory name only (label_name)
    static std::string sanitize_dirname(const std::string &s) {
        std::string out = sanitize_filename(s);
        out.erase(std::ranges::remove(out, '\\').begin(), out.end());
        out.erase(std::ranges::remove(out, '/').begin(), out.end());
        out.erase(std::ranges::remove(out, '.').begin(), out.end());
        return out;
    }

    static bool is_under_directory(const std::filesystem::path &p,
                                   const std::filesystem::path &root) {
        auto p_it = p.begin();
        auto r_it = root.begin();

        for (; r_it != root.end(); ++r_it, ++p_it) {
            if (p_it == p.end() || *p_it != *r_it) { return false; }
        }
        return true;
    }

    // yyyyMMdd-HHmmss
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
