#include "ConfidenceInterval.h"

namespace util {

std::ostream &operator<<(std::ostream &ofs, ConfidenceInterval const& confidenceInterval) {
    return ofs << static_cast<std::string>(confidenceInterval);
}

} // namespace util

