#ifndef MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_ABSTRACTTERM_HPP
#define MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_ABSTRACTTERM_HPP
/**
 * @file AbstractTerm.hpp
 * @brief Defines the AbstractTerm structure, which serves as a parent to both \c Term and \c WickTerm.
 */

#include <vector>

#include "Fractional.hpp"
#include "KroneckerDelta.hpp"
#include "KroneckerDeltaUtility.hpp"
#include "Coefficient.hpp"
#include "SumContainer.hpp"

namespace mrock::symbolic_operators {
    /**
	 * @class AbstractTerm
	 * @brief Serves as a parent to \c Term and \c WickTerm.
     * Defines and implements certain methods that are used by both classes as to avoid code duplication
     * 
     * @sa Term, WickTerm
     */
    template<class tOperatorType>
    class AbstractTerm {
	protected:
		// Here, we define buffers to be used later on.
		// Later, we want to bring all terms to the same notation.
		// This will entail renaming sums, so that the first sum is always q, the second p, the third r, etc.
		// To avoid name clashes, we first rename all sums to a buffer name, e.g., the first sum is renamed to :, the second to ;, the third to |, etc.
		// As of now, this is limited to 11 sums (which I doubt will become a problem in the future).
		// If more is needed, the buffer_list and name_list need to be extended, and the N_BUFFER constant needs to be changed accordingly.
		constexpr static int N_BUFFER = 11;
		constexpr static MomentumSymbol::name_type name_list[N_BUFFER]   = { 'q', 'p', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z' };
		constexpr static MomentumSymbol::name_type buffer_list[N_BUFFER] = { ':', ';', '|', '?', '!', '.', '-', '_', '+', '/', '=' };

    public:
        IntFractional multiplicity; ///< Multiplicity of the term.
        std::vector<Coefficient> coefficients; ///< Coefficients of the term.
		SumContainer sums; ///< Sum container for the term. Contains e.g. \sum_{k,l} \sum_{sigma}
        std::vector<KroneckerDelta<Momentum>> delta_momenta; ///< Kronecker delta for momenta.
		std::vector<KroneckerDelta<Index>> delta_indizes; ///< Kronecker delta for indices.
        std::vector<tOperatorType> operators; ///< Operators in the term, if empty the term is considered to contain the identiy operator
       

        /**
		 * @brief Default constructor.
		 */
        AbstractTerm() = default;

        /**
		 * @brief Constructs a Term with a summation over momenta and spins and multiple coefficients and Kronecker deltas
		 * 
		 * @param _multiplicity The multiplicity of the term
		 * @param _sums The sums
         * @param _coefficients The coefficients
         * @param _operators The operators
         * @param _delta_momenta The Kronecker deltas for the momenta
         * @param _delta_indizes The Kronecker deltas for the indizes
		 */
        AbstractTerm(const IntFractional& _multiplicity, const std::vector<Coefficient>& _coefficients, const SumContainer& _sums, const std::vector<KroneckerDelta<Momentum>>& _delta_momenta, const std::vector<KroneckerDelta<Index>>& _delta_indizes, const std::vector<tOperatorType>& _operators)
            : multiplicity{_multiplicity}, coefficients{_coefficients}, sums{_sums}, delta_momenta{_delta_momenta}, delta_indizes{_delta_indizes}, operators{_operators}
        { };

        /**
		 * @brief Constructs a Term with a summation over momenta and spins and multiple coefficients
		 * 
		 * @param _multiplicity The multiplicity of the term
		 * @param _sums The sums
         * @param _coefficients The coefficients
         * @param _operators The operators
		 */
        AbstractTerm(const IntFractional& _multiplicity, const std::vector<Coefficient>& _coefficients, const SumContainer& _sums, const std::vector<tOperatorType>& _operators = std::vector<tOperatorType>())
            : multiplicity{_multiplicity}, coefficients{_coefficients}, sums{_sums}, operators{_operators}
        { };

        /**
		 * @brief Constructs a Term with a summation over momenta and spins and one coefficient
		 * 
		 * @param _multiplicity The multiplicity of the term
		 * @param _sums The sums
         * @param _coefficient The coefficient
         * @param _operators The operators
		 */
        AbstractTerm(const IntFractional& _multiplicity, const Coefficient& _coefficient, const SumContainer& _sums, const std::vector<tOperatorType>& _operators = std::vector<tOperatorType>())
            : multiplicity{_multiplicity}, coefficients(1, _coefficient), sums{_sums}, operators{_operators}
        { };

        /**
		 * @brief Constructs a Term with no summations
		 * 
		 * @param _multiplicity The multiplicity of the term
         * @param _coefficient The coefficient
         * @param _operators The operators
		 */
        AbstractTerm(const IntFractional& _multiplicity, const Coefficient& _coefficient, const std::vector<tOperatorType>& _operators = std::vector<tOperatorType>())
            : multiplicity{_multiplicity}, coefficients(1, _coefficient), operators{_operators}
        { };

        /**
		 * @brief Constructs a Term with no summations and no coefficient
		 * 
		 * @param _multiplicity The multiplicity of the term
         * @param _coefficient The coefficient
         * @param _operators The operators
		 */
        explicit AbstractTerm(const IntFractional& _multiplicity, const std::vector<tOperatorType>& _operators = std::vector<tOperatorType>())
            : multiplicity{_multiplicity}, operators{_operators}
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
		 * @brief Flips the sign of the term.
		 */
		void flip_sign();

		/**
		 * @brief Gets the operators in the term.
		 * @return The operators.
		 */
		const std::vector<tOperatorType>& get_operators() const;
        
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
    template<class tOperatorType>
    bool AbstractTerm<tOperatorType>::resolve_momentum_deltas() 
	{
		if (is_always_zero(delta_momenta)) return false;

		// Remove delta^2
		remove_delta_squared(this->delta_momenta);
		// Erase delta_k,k etc
		remove_delta_is_one(this->delta_momenta);

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

			// Try to find a momentum in the sums that is also in the delta
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

			delta_it->second = Momentum(resolve_to);
			delta_it->second.flip_momentum();
			delta_it->first += delta_it->second;
			if (delta_it->second.front().factor < 0) {
				delta_it->first.flip_momentum();
				delta_it->second.flip_momentum();
			}
			
			// Fractional momenta were never implemented, since they were never needed.
			// If they are ever needed, you unfortunately have to implement them yourself
			for (MomentumSymbol& symbol : delta_it->first) {
				assert(symbol.factor % delta_it->second.front().factor == 0);
				symbol.factor /= delta_it->second.front().factor;
			}

			// Replace set the delta everywhere, e.g., delta_{k,l+q} would replace each k with l+q
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
		remove_delta_squared(this->delta_momenta);
		// Erase delta_k,k etc
		remove_delta_is_one(this->delta_momenta);

		return true;
	}


	template<class tOperatorType>
    bool AbstractTerm<tOperatorType>::resolve_index_deltas() 
	{
		if (is_always_zero(delta_indizes)) return false;

		// Remove delta^2
		remove_delta_squared(this->delta_indizes);
		// Erase delta_k,k etc
		remove_delta_is_one(this->delta_indizes);

		for (auto delta_it = delta_indizes.begin(); delta_it != delta_indizes.end(); ) {
			Index to_resolve { Index::UndefinedIndex };
			Index change_to { Index::UndefinedIndex };
			bool found_sum{};

			// try to find a spin in the sums that is also in the delta
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

    template<class tOperatorType>
    void AbstractTerm<tOperatorType>::rename_sums()
	{
		for (std::size_t i = 0U; i < sums.momenta.size(); ++i)
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

		for (std::size_t i = 0U; i < sums.momenta.size(); ++i)
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
				op.indizes.replace_index(Index::SigmaPrime, Index::Sigma);
			}
			for (auto& coeff : coefficients) {
				coeff.indizes.replace_index(Index::SigmaPrime, Index::Sigma);
			}
		}
	}

    template<class tOperatorType>
    bool AbstractTerm<tOperatorType>::is_identity() const noexcept 
    {
		return this->operators.empty();
	}

    template<class tOperatorType>
    void AbstractTerm<tOperatorType>::flip_sign() {
		this->multiplicity *= -1;
	}

    template<class tOperatorType>
    const std::vector<tOperatorType>& AbstractTerm<tOperatorType>::get_operators() const {
		return this->operators;
	}

    template<class tOperatorType>
    void AbstractTerm<tOperatorType>::invert_momentum(const MomentumSymbol::name_type what) 
    {
		for (auto& coeff : coefficients) {
			coeff.invert_momentum(what);
		}
		for (auto& op : operators) {
			op.momentum.flip_single(what);
		}
	}

    template<class tOperatorType>
    void AbstractTerm<tOperatorType>::invert_momentum_sum(const MomentumSymbol::name_type what) 
    {
		if (std::find(sums.momenta.begin(), sums.momenta.end(), what) == sums.momenta.end()) {
			throw std::invalid_argument("You are trying to perform a sum transformation on a momentum that is not being summed over!");
		}
		invert_momentum(what);
	}

    template<class tOperatorType>
    void AbstractTerm<tOperatorType>::remove_momentum_contribution(const MomentumSymbol::name_type value) 
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
#endif  // MROCK_SYMBOLIC_OPERATORS_INCLUDE_MROCK_SYMBOLIC_OPERATORS_ABSTRACTTERM_HPP
