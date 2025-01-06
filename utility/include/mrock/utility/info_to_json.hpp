#pragma once
#include <nlohmann/json.hpp>
#include <string>

namespace mrock::utility {
    template <class info>
    nlohmann::json generate_json(std::string const& key_prefix) {
        return nlohmann::json {
            { key_prefix + "git_commit_version", std::string(info::GIT_COMMIT_VERSION) }, 
            { key_prefix + "git_commit_name",    std::string(info::GIT_COMMIT_NAME   ) }, 
            { key_prefix + "git_commit_date",    std::string(info::GIT_COMMIT_DATE   ) }, 
            { key_prefix + "make_date",          std::string(info::MAKE_DATE         ) }, 
            { key_prefix + "cxx_compiler",       std::string(info::CXX_COMPILER      ) }, 
            { key_prefix + "cxx_standard",       std::string(info::CXX_STANDARD      ) }, 
            { key_prefix + "project_version",    std::string(info::PROJECT_VERSION   ) }, 
            { key_prefix + "system_name",        std::string(info::SYSTEM_NAME       ) }, 
            { key_prefix + "system_version",     std::string(info::SYSTEM_VERSION    ) }, 
            { key_prefix + "hostname",           std::string(info::HOSTNAME          ) }
        };
    }
}