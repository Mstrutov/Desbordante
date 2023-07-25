#include "util/config/thread_number/option.h"

#include <thread>

#include "util/config/names_and_descriptions.h"

namespace util::config {
using names::kThreads, descriptions::kDThreads;
extern const CommonOption<ThreadNumType> ThreadNumberOpt{
        kThreads, kDThreads, 0, [](auto &value) {
            if (value == 0) {
                value = std::thread::hardware_concurrency();
                if (value == 0) {
                    throw std::invalid_argument(
                            "Unable to detect number of concurrent threads supported by your "
                            "system. Please, specify it manually.");
                }
            }
        }};
}  // namespace util::config
