#include <filesystem>

#include <easylogging++.h>
#include <pybind11/pybind11.h>

#include "ac/bind_ac.h"
#include "ar/bind_ar.h"
#include "bind_main_classes.h"
#include "cfd/bind_cfd.h"
#include "data/bind_data.h"
#include "dd/bind_split.h"
#include "fd/bind_fd.h"
#include "fd/bind_fd_verification.h"
#include "gfd/bind_gfd_verification.h"
#include "ind/bind_ind.h"
#include "mfd/bind_mfd_verification.h"
#include "od/bind_od.h"
#include "statistics/bind_statistics.h"
#include "ucc/bind_ucc.h"
#include "ucc/bind_ucc_verification.h"

INITIALIZE_EASYLOGGINGPP

namespace python_bindings {

PYBIND11_MODULE(desbordante, module) {
    using namespace pybind11::literals;

    if (std::filesystem::exists("logging.conf")) {
        el::Loggers::configureFromGlobal("logging.conf");
    } else {
        using std::filesystem::is_empty, std::filesystem::remove;
        el::Configurations conf;
        conf.set(el::Level::Global, el::ConfigurationType::Enabled, "false");
        el::Loggers::reconfigureAllLoggers(conf);
        try {
            static constexpr auto kElppDefaultFile = "myeasylog.log";
            if (is_empty(kElppDefaultFile)) remove(kElppDefaultFile);
        } catch (std::filesystem::filesystem_error&) {
        }
    }

    for (auto bind_func :
         {BindMainClasses, BindDataModule, BindFd, BindCfd, BindAr, BindUcc, BindAc, BindOd,
          BindFdVerification, BindMfdVerification, BindUccVerification, BindStatistics, BindInd,
          BindGfdVerification, BindSplit}) {
        bind_func(module);
    }
}

}  // namespace python_bindings
