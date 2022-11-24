#include "algorithms/cfd/util/string_util.h"

// see ../algorithms/cfd/LICENSE

#include <algorithm>

// trim from start
[[maybe_unused]] std::string& LeftTrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
[[maybe_unused]] std::string& RightTrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
[[maybe_unused]] std::string& Trim(std::string &s) {
    return LeftTrim(RightTrim(s));
}

[[maybe_unused]] std::string Concat(int count, ...) {
    std::stringstream ss;
    va_list args;
    va_start(args, count);
    for (int i = 0; i < count; i++) {
        ss << va_arg(args, const char*);
    }
    va_end(args);
    return ss.str();
}

[[maybe_unused]] std::string ConcatCsv(int count, ...) {
    std::stringstream ss;
    va_list args;
    va_start(args, count);
    for (int i = 0; i < count; i++) {
        ss << va_arg(args, const char*);
        if (i < count - 1) {
            ss << ",";
        }
    }
    va_end(args);
    return ss.str();
}
