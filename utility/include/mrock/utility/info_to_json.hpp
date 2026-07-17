#ifndef MROCK_UTILITY_INCLUDE_MROCK_UTILITY_INFO_TO_JSON_HPP
#define MROCK_UTILITY_INCLUDE_MROCK_UTILITY_INFO_TO_JSON_HPP

#include <nlohmann/json.hpp>

#include <string>

#if defined(_MSVC_LANG)
#define MROCK_INFO_CXX_STANDARD_VALUE _MSVC_LANG
#else
#define MROCK_INFO_CXX_STANDARD_VALUE __cplusplus
#endif

namespace mrock::utility {
template <class info>
nlohmann::json generate_json(std::string const& key_prefix) {
    return nlohmann::json{{key_prefix + "git_commit_version", std::string(info::GIT_COMMIT_VERSION)},
                          {key_prefix + "git_commit_name", std::string(info::GIT_COMMIT_NAME)},
                          {key_prefix + "git_commit_date", std::string(info::GIT_COMMIT_DATE)},
                          {key_prefix + "make_date", std::string(info::MAKE_DATE)},
                          {key_prefix + "cxx_compiler", std::string(info::CXX_COMPILER)},
                          {key_prefix + "cxx_standard", std::to_string(MROCK_INFO_CXX_STANDARD_VALUE)},
                          {key_prefix + "project_version", std::string(info::PROJECT_VERSION)},
                          {key_prefix + "system_name", std::string(info::SYSTEM_NAME)},
                          {key_prefix + "system_version", std::string(info::SYSTEM_VERSION)},
                          {key_prefix + "hostname", std::string(info::HOSTNAME)}};
}
}  // namespace mrock::utility
#endif  // MROCK_UTILITY_INCLUDE_MROCK_UTILITY_INFO_TO_JSON_HPP
