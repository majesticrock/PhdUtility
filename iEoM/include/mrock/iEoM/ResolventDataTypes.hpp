#ifndef MROCK_IEOM_INCLUDE_MROCK_IEOM_RESOLVENTDATATYPES_HPP
#define MROCK_IEOM_INCLUDE_MROCK_IEOM_RESOLVENTDATATYPES_HPP
#include <array>
#include <iostream>
#include <string>
#include <vector>

#ifndef MROCK_IEOM_NO_NLOHMANN_JSON
#include <nlohmann/json.hpp>
#endif

namespace mrock::iEoM {
/**
 * @brief Stores Lanczos recurrence coefficients for one resolvent run.
 *
 * @tparam RealType Numeric scalar type for the recurrence coefficients.
 */
template <class RealType>
struct ResolventData {
    std::vector<RealType> a_i;
    std::vector<RealType> b_i;
};

/**
 * @brief Wraps multiple resolvent Lanczos results under a common name.
 *
 * @tparam RealType Numeric scalar type for the stored coefficients.
 */
template <class RealType>
struct ResolventDataWrapper {
    std::vector<ResolventData<RealType>> lanczos;
    std::string name;

    /**
     * @brief Default constructed empty wrapper.
     */
    ResolventDataWrapper() = default;
    /**
     * @brief Construct a wrapper with the given name.
     *
     * @param _name Name associated with this data wrapper.
     */
    ResolventDataWrapper(const std::string& _name) : name(_name){};

    /**
     * @brief Append a Lanczos result by moving it into the wrapper.
     *
     * @param data_point Lanczos coefficients to append.
     */
    void push_back(ResolventData<RealType>&& data_point) { lanczos.push_back(std::move(data_point)); };
    /**
     * @brief Append a Lanczos result by copying it into the wrapper.
     *
     * @param data_point Lanczos coefficients to append.
     */
    void push_back(const ResolventData<RealType>& data_point) { lanczos.push_back(data_point); };
};

/**
 * @brief Merge a named ResolventDataWrapper into an existing collection.
 *
 * If a wrapper with the same name already exists, its lanczos data are
 * appended. Otherwise, the new wrapper is added to the collection.
 *
 * @tparam RealType Numeric type of the resolvent data.
 * @param target Collection of resolvent data wrappers.
 * @param new_data New wrapper to merge.
 */
template <class RealType>
void join_data_wrapper(std::vector<ResolventDataWrapper<RealType>>& target,
                       ResolventDataWrapper<RealType> const& new_data) {
    for (auto& sub_target : target) {
        if (sub_target.name == new_data.name) {
            // Appends new_data.lanczos to sub_target.lanczos
            sub_target.lanczos.insert(sub_target.lanczos.end(), new_data.lanczos.begin(), new_data.lanczos.end());
            return;
        }
    }
    target.push_back(new_data);
}

/**
 * @brief Merge multiple resolvent wrappers into an existing collection.
 *
 * @tparam RealType Numeric type of the resolvent data.
 * @param target Collection of existing wrappers.
 * @param new_data New wrappers to merge.
 */
template <class RealType>
void join_data_wrapper(std::vector<ResolventDataWrapper<RealType>>& target,
                       std::vector<ResolventDataWrapper<RealType>> const& new_data) {
    if (target.empty()) {
        target = new_data;
        return;
    }
    for (const auto& _new : new_data) {
        join_data_wrapper(target, _new);
    }
}

/**
 * @brief Stream output operator for ResolventData.
 *
 * Writes the a_i and b_i coefficient sequences on separate lines.
 *
 * @tparam T Numeric type of the resolvent coefficients.
 * @param os Output stream.
 * @param data Resolvent data to stream.
 * @return Reference to the output stream.
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const ResolventData<T>& data) {
    for (const auto& elem : data.a_i) {
        os << elem << " ";
    }
    os << "\n";
    for (const auto& elem : data.b_i) {
        os << elem << " ";
    }
    os << "\n";
    return os;
}

/**
 * @brief Stores residual eigenvalue and Ritz vector information.
 *
 * @tparam RealType Numeric type for the residual data.
 * @tparam n Number of residual roots tracked.
 */
template <class RealType, int n>
struct ResidualInformation {
    std::array<RealType, n> eigenvalues{};
    std::array<std::vector<RealType>, n> eigenvectors{};
    std::array<RealType, n> weights{};
    std::array<RealType, n> residuals{};
    std::array<int, n> n_ghosts{};
    std::array<bool, n> converged{};

    ResidualInformation() {
        // Placeholder filling
        // Filling with 0 may produces errors if an actual eigenvalue is 0 and one asks whether
        // The eigenvalue is already contained within ResidualInformation
        eigenvalues.fill(std::numeric_limits<RealType>::max());
    }
};

/**
 * @brief Stores full diagonalization results for retained eigenvectors.
 *
 * @tparam RealType Numeric type for eigenvalues and eigenvectors.
 * @tparam n Number of retained eigenvectors.
 */
template <class RealType, int n>
struct FullDiagonalizationData {
    std::vector<RealType> eigenvalues;
    std::array<std::vector<RealType>, n> first_eigenvectors;
    std::vector<std::vector<RealType>> weights;
};

#ifndef MROCK_IEOM_NO_NLOHMANN_JSON
/**
 * @brief Serialize ResolventData to JSON.
 *
 * @tparam RealType Numeric type of the resolvent coefficients.
 * @param j JSON output object.
 * @param res_data Data to serialize.
 */
template <class RealType>
void to_json(nlohmann::json& j, const ResolventData<RealType>& res_data) {
    j = nlohmann::json{{"a_i", res_data.a_i}, {"b_i", res_data.b_i}};
}

/**
 * @brief Serialize a vector of ResolventDataWrapper objects to JSON.
 *
 * @tparam RealType Numeric type of the resolvent data.
 * @param j JSON output object.
 * @param vec_resolvent_data Collection of resolvent wrappers.
 */
template <class RealType>
void to_json(nlohmann::json& j, const std::vector<ResolventDataWrapper<RealType>>& vec_resolvent_data) {
    for (const auto& res : vec_resolvent_data) {
        j[res.name] = res.lanczos;
    }
}

/**
 * @brief Serialize ResidualInformation to JSON.
 *
 * @tparam RealType Numeric type of the residual data.
 * @tparam n Number of residual entries.
 * @param j JSON output object.
 * @param res_data Data to serialize.
 */
template <class RealType, int n>
void to_json(nlohmann::json& j, const ResidualInformation<RealType, n>& res_data) {
    j = nlohmann::json{
        {"eigenvalues", res_data.eigenvalues}, {"eigenvectors", res_data.eigenvectors}, {"weights", res_data.weights},
        {"residuals", res_data.residuals},     {"converged", res_data.converged},       {"n_ghosts", res_data.n_ghosts},
    };
}

/**
 * @brief Serialize FullDiagonalizationData to JSON.
 *
 * @tparam RealType Numeric type of the diagonalization data.
 * @tparam n Number of retained eigenvectors.
 * @param j JSON output object.
 * @param res_data Data to serialize.
 */
template <class RealType, int n>
void to_json(nlohmann::json& j, const FullDiagonalizationData<RealType, n>& res_data) {
    j = nlohmann::json{
        {"eigenvalues", res_data.eigenvalues},
        {"first_eigenvectors", res_data.first_eigenvectors},
        {"weights", res_data.weights},
    };
}
#endif
}  // namespace mrock::iEoM
#endif  // MROCK_IEOM_INCLUDE_MROCK_IEOM_RESOLVENTDATATYPES_HPP
