/**
 * @file AbstractTerm.hpp
 * @brief Defines the AbstractTerm structure, which serves as a parent to both \c Term and \c WickTerm.
 */

#pragma once

#include <vector>
#include <mrock/utility/Fractional.hpp>

#include "KroneckerDelta.hpp"
#include "KroneckerDeltaUtility.hpp"
#include "Coefficient.hpp"
#include "SumContainer.hpp"

namespace mrock::symbolic_operators {
    using IntFractional = mrock::utility::Fractional<int>;

    /**
	 * @class AbstractTerm
	 * @brief Serves as a parent to \c Term and \c WickTerm.
     * Defines and implements certain methods that are used by both classes as to avoid code duplication
     * 
     * @sa Term, WickTerm
     */
    template<class OperatorType>
    class AbstractTerm {
    public:
        IntFractional multiplicity; ///< Multiplicity of the term.
		SumContainer sums; ///< Sum container for the term. Contains e.g. \sum_{k,l} \sum_{sigma}
        std::vector<KroneckerDelta<Momentum>> delta_momenta; ///< Kronecker delta for momenta.
		std::vector<KroneckerDelta<Index>> delta_indizes; ///< Kronecker delta for indices.
        std::vector<OperatorType> operators; ///< Operators in the term, if empty the term is considered to contain the identiy operator
        std::vector<Coefficient> coefficients; ///< Coefficients of the term.

        /**
		 * @brief Default constructor.
		 */
        AbstractTerm() = default;

        /**
		 * @brief Constructs a Term with a summation over momenta and spins and multiple coefficients
		 * 
		 * @param _multiplicity The multiplicity of the term
		 * 
		 * @param _sums The sums
         * @param _coefficients The coefficients
		 */
        AbstractTerm(const IntFractional& _multiplicity, const SumContainer& _sums, const std::vector<OperatorType>& _operators, const std::vector<Coefficient>& coefficients)
            : multiplicity{_multiplicity}, sums{_sums}, operators{_operators}, coefficients{_coefficients}
        { };

        /**
		 * @brief Resolves the Kronecker deltas of the momenta in the term.
		 * @return True if successful, false otherwise.
		 */
		bool resolve_momentum_deltas();

		/**
		 * @brief Resolves the Kronecker deltas of the indizes in the term.
		 * @return True if successful, false otherwise.
		 */
		bool resolve_index_deltas();

		/**
		 * @brief Resolves the Kronecker deltas in the term ( calls \c resolve_momentum_deltas() and \c resolve_index_deltas() )
		 * @return True if successful, false otherwise.
		 */
		bool resolve_deltas();

        /**
		 * @brief Renames the sums in the term.
		 */
		void rename_sums();

        /**
		 * @brief Checks if the term is an identity term.
		 * 
		 * @return true if the term is an identity term.
		 * @return false otherwise.
		 */
		bool is_identity() const noexcept;

        /**
		 * @brief Inverts a momentum in the term.
		 * 
		 * @param what The momentum to invert.
		 */
		void invert_momentum(const MomentumSymbol::name_type what);

		/**
		 * @brief Inverts a momentum sum in the term.
		 * 
		 * @param what The momentum sum to invert.
		 */
		void invert_momentum_sum(const MomentumSymbol::name_type what);

		/**
		 * @brief Removes a momentum contribution from the term.
		 * 
		 * @param value The momentum value to remove.
		 */
		void remove_momentum_contribution(const MomentumSymbol::name_type value);
    };


    // Implementations
    template<class OperatorType>
    bool AbstractTerm<OperatorType>::resolve_momentum_deltas() 
	{
		for (auto delta_it = delta_momenta.begin(); delta_it != delta_momenta.end(); ) {
			delta_it->first -= delta_it->second;
			delta_it->second = Momentum();

			if (delta_it->first.momentum_list.empty() && delta_it->second.momentum_list.empty()) {
				// 0 = Q can never be achieved
				if (delta_it->first.add_Q != delta_it->second.add_Q) return false;
				// delta_(0,0) = 1
				delta_it = delta_momenta.erase(delta_it);
				continue;
			}
			
			MomentumSymbol resolve_to{ *(delta_it->first.begin()) };
			bool found_sum{};

			for (auto sum_it = sums.momenta.begin(); sum_it != sums.momenta.end(); ++sum_it) {
				const auto found_it = std::find_if(delta_it->first.begin(), delta_it->first.end(), [&sum_it](const MomentumSymbol& symbol) {
					return symbol.name == *sum_it;
				});
				if ( found_it != delta_it->first.end()) {
					resolve_to = *found_it;
					sums.momenta.erase(sum_it);
					found_sum = true;
					break;
				}
			}

			if (resolve_to.factor > 0) {
				delta_it->first.flip_momentum();
			}
			else {
				resolve_to.factor *= -1;
			}
			delta_it->second = Momentum(resolve_to);
			delta_it->first += delta_it->second;
			

			for (MomentumSymbol& symbol : delta_it->first) {
				assert(symbol.factor % delta_it->second.front().factor == 0);
				symbol.factor /= delta_it->second.front().factor;
			}

			for (auto& coeff : coefficients) {
				coeff.momenta.replace_occurances(delta_it->second.front().name, delta_it->first);
			}
			for (auto& op : operators) {
				op.momentum.replace_occurances(delta_it->second.front().name, delta_it->first);
			}
			for (auto delta_it2 = delta_momenta.begin(); delta_it2 != delta_momenta.end(); ++delta_it2) {
				if (delta_it2 == delta_it) continue;
				delta_it2->first.replace_occurances(delta_it->second.front().name, delta_it->first);
				delta_it2->second.replace_occurances(delta_it->second.front().name, delta_it->first);
			}

			if (found_sum) {
				delta_it = delta_momenta.erase(delta_it);
			}
			else {
				++delta_it;
			}
		}

		// Make sure that delta.first has always exactly one momentum symbol (or delta_{0,0})
		for (auto& delta : delta_momenta) {
			if (delta.first.empty() && !delta.second.empty()) {
				std::swap(delta.first, delta.second);
			} 
			if (delta.first.size() > 1U) {
				const Momentum shift = delta.first - Momentum(delta.first.front());
				delta.first -= shift;
				delta.second -= shift;
			}
			if (delta.first.front().factor < 0) {
				delta.first.flip_momentum();
				delta.second.flip_momentum();
			}
		}

		// Remove delta^2
		remove_delta_squared(this->delta_indizes);
		// Erase delta_k,k etc
		remove_delta_is_one(this->delta_indizes);

		return true;
	}


	template<class OperatorType>
    bool AbstractTerm<OperatorType>::resolve_index_deltas() 
	{
		if (is_always_zero(delta_indizes)) return false;

		for (auto delta_it = delta_indizes.begin(); delta_it != delta_indizes.end(); ) {
			Index to_resolve { Index::UndefinedIndex };
			Index change_to { Index::UndefinedIndex };
			bool found_sum{};
			auto sum_it = std::find_if( sums.spins.begin(), sums.spins.end(), [&delta_it](const Index& idx) {
				return idx == delta_it->first;
			});
			if (sum_it != sums.spins.end()) {
				to_resolve = delta_it->first;
				change_to = delta_it->second;
				found_sum = true;
			}
			else {
				sum_it = std::find_if( sums.spins.begin(), sums.spins.end(), [&delta_it](const Index& idx) {
					return idx == delta_it->second;
				});
				if (sum_it != sums.spins.end()) {
					to_resolve = delta_it->second;
					change_to = delta_it->first;
					found_sum = true;
				}
			}

			if (to_resolve == Index::UndefinedIndex) {
				if (is_mutable(delta_it->first)) {
					to_resolve = delta_it->first;
					change_to = delta_it->second;
				}
				else if (is_mutable(delta_it->second)) {
					to_resolve = delta_it->second;
					change_to = delta_it->first;
				}
				else if(delta_it->first != delta_it->second) {
					// Two differing, immutable indizes can never be equal
					return false;
				}
			}

			if(delta_it->first == delta_it->second) continue;

			for (auto& op : operators) {
				op.indizes.replace_index(to_resolve, change_to);
			}
			for (auto& coeff : coefficients) {
				coeff.indizes.replace_index(to_resolve, change_to);
			}

			for (auto delta_it2 = delta_indizes.begin(); delta_it2 != delta_indizes.end(); ++delta_it2) {
				if (delta_it2 == delta_it) continue;
				if (delta_it2->first == to_resolve) {
					delta_it2->first = change_to;
				}
				if (delta_it2->second == to_resolve) {
					delta_it2->second = change_to;
				}
			}

			if (found_sum) {
				sums.spins.erase(sum_it);
				delta_it = delta_indizes.erase(delta_it);
			}
			else {
				++delta_it;
			}
		}

		// Remove delta^2
		remove_delta_squared(this->delta_indizes);
		// Erase delta_k,k etc
		remove_delta_is_one(this->delta_indizes);

		return true;
	}

    template<class OperatorType>
    void AbstractTerm<OperatorType>::rename_sums()
	{
		constexpr int N_BUFFER = 11;
		constexpr MomentumSymbol::name_type name_list[N_BUFFER]   = { 'q', 'p', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z' };
		constexpr MomentumSymbol::name_type buffer_list[N_BUFFER] = { ':', ';', '|', '?', '!', '.', '-', '_', '+', '/', '=' };
		for (size_t i = 0U; i < sums.momenta.size(); ++i)
		{
			if (i >= N_BUFFER) {
				std::cerr << "More than " << N_BUFFER << "momenta, time to implement this..." << std::endl;
				break;
			}
			if (sums.momenta[i] == name_list[i]) continue;

			for (auto& op : operators) {
				op.momentum.replace_occurances(sums.momenta[i], Momentum(buffer_list[i]));
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replace_occurances(sums.momenta[i], Momentum(buffer_list[i]));
			}
			sums.momenta[i] = name_list[i];
		}

		for (size_t i = 0U; i < sums.momenta.size(); ++i)
		{
			for (auto& op : operators) {
				op.momentum.replace_occurances(buffer_list[i], Momentum(name_list[i]));
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replace_occurances(buffer_list[i], Momentum(name_list[i]));
			}
		}

		if (sums.spins.size() == 1U && sums.spins.front() == Index::SigmaPrime) {
			sums.spins.front() = Index::Sigma;
			for (auto& op : operators) {
				for (auto& index : op.indizes) {
					if (index == Index::SigmaPrime) index = Index::Sigma;
				}
			}
			for (auto& coeff : coefficients) {
				for (auto& index : coeff.indizes) {
					if (index == Index::SigmaPrime) index = Index::Sigma;
				}
			}
		}
	}

    template<class OperatorType>
    bool AbstractTerm<OperatorType>::is_identity() const noexcept 
    {
		return this->operators.empty();
	}

    template<class OperatorType>
    void AbstractTerm<OperatorType>::invert_momentum(const MomentumSymbol::name_type what) 
    {
		for (auto& coeff : coefficients) {
			coeff.invert_momentum(what);
		}
		for (auto& op : operators) {
			op.momentum.flip_single(what);
		}
	}

    template<class OperatorType>
    void AbstractTerm<OperatorType>::invert_momentum_sum(const MomentumSymbol::name_type what) 
    {
		if (std::find(sums.momenta.begin(), sums.momenta.end(), what) == sums.momenta.end()) {
			throw std::invalid_argument("You are trying to perform a sum transformation on a momentum that is not being summed over!");
		}
		invert_momentum(what);
	}

    template<class OperatorType>
    void AbstractTerm<OperatorType>::remove_momentum_contribution(const MomentumSymbol::name_type value) 
    {
	    for (auto& coeff : coefficients) {
	    	coeff.remove_momentum_contribution(value);
	    }
	    for (auto& op : operators) {
	    	op.remove_momentum_contribution(value);
	    }
	    for (auto& delta : delta_momenta) {
	    	delta.first.remove_contribution(value);
	    	delta.second.remove_contribution(value);
	    }
	    std::erase_if(sums.momenta.summations, [&](const MomentumSymbol::name_type sum_idx) { return sum_idx == value; });
	}


} // namespace mrock::symbolic_operators