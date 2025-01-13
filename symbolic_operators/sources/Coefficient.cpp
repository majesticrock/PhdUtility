#include <mrock/symbolic_operators/Coefficient.hpp>
#include <mrock/utility/StringUtility.hpp>

namespace mrock::symbolic_operators {
	void Coefficient::invert_momentum(char what) {
		for (auto& mom : momenta) {
			if (this->inversion_symmetry && mom.size() == 1U) {
				// If the coefficient is translationally invariant, i.e., c(k) = c(-k),
				// we want to make sure that the momentum is always positive
				if (mom.front().first > 0) continue;
			}
			mom.flip_single(what);
		}
	}

	void Coefficient::use_symmetric_interaction_exchange() {
		if (this->is_symmetrized_interaction) {
			assert(this->momenta.size() == 3U);
			std::swap(this->momenta[0], this->momenta[1]);
			this->momenta.back().flip_momentum();
		}
	}

	void Coefficient::use_symmetric_interaction_inversion() {
		if (this->is_symmetrized_interaction) {
			assert(this->momenta.size() == 3U);
			for(auto& mom : momenta) {
				mom.flip_momentum();
			}
		}
	}

	void Coefficient::remove_momentum_contribution(char value)
	{
		for (auto& mom : momenta) {
			mom.remove_contribution(value);
		}
	}

    void Coefficient::apply_custom_symmetry()
    {
		if (custom_symmetry.has_value()) {
			custom_symmetry.value()(*this);
		}
    }

    Coefficient Coefficient::parse_string(const std::string &expression)
    {
		// Syntax:   name{Momentum_expression1,Momentum_expression1;index1,index2,...}

		Coefficient ret;
		ret.name = expression.substr(0U, expression.find('{'));
		std::vector<std::string> momentum_strings = mrock::utility::extract_elements(expression, '{', ';');
		std::vector<std::string> index_strings = mrock::utility::extract_elements(expression, ';', '}');

		ret.momenta.reserve(momentum_strings.size());
		for (const auto& arg : momentum_strings) {
			ret.momenta.push_back(Momentum(arg));
		}
		ret.indizes.reserve(index_strings.size());
		for (const auto& arg : index_strings) {
			ret.indizes.push_back(string_to_index.at(arg));
		}

		return ret;
	}

	std::ostream& operator<<(std::ostream& os, const Coefficient& coeff)
	{
		os << coeff.name;
		if (!coeff.indizes.empty()) {
			os << "_{ " << coeff.indizes << "}";
		}
		if (coeff.is_daggered) {
			os << "^*";
		}
		os << coeff.momenta << " ";
		return os;
	}
	std::ostream& operator<<(std::ostream& os, const std::vector<Coefficient>& coeffs) {
		for (std::vector<Coefficient>::const_iterator it = coeffs.begin(); it != coeffs.end(); ++it)
		{
			os << (*it) << " ";
		}
		return os;
	}

	Coefficient::Coefficient(const std::string& _name)
		: name(_name), momenta(), indizes(), Q_changes_sign(false), is_daggered(false) {}
	Coefficient::Coefficient(const std::string& _name, const Momentum& _momentum, const IndexWrapper& _indizes, bool _Q_changes_sign, bool _inversion_symmetry, bool _is_daggered)
		: name(_name), momenta(_momentum), indizes(_indizes), inversion_symmetry{_inversion_symmetry}, Q_changes_sign(_Q_changes_sign), is_daggered(_is_daggered) {}
	Coefficient::Coefficient(const std::string& _name, const Momentum& _momentum, bool _Q_changes_sign, bool _inversion_symmetry, bool _is_daggered)
		: name(_name), momenta(_momentum), indizes(), inversion_symmetry{_inversion_symmetry}, Q_changes_sign(_Q_changes_sign), is_daggered(_is_daggered) {}
	Coefficient::Coefficient(const std::string& _name, const MomentumList& _momenta, const IndexWrapper& _indizes, bool _Q_changes_sign, bool _inversion_symmetry, bool _is_daggered)
		: name(_name), momenta(_momenta), indizes(), inversion_symmetry{_inversion_symmetry}, Q_changes_sign(_Q_changes_sign), is_daggered(_is_daggered) { }

	Coefficient Coefficient::RealInversionSymmetric(const std::string& _name, const MomentumList& _momenta, 
		const std::optional<std::function<void(Coefficient&)>>& _custom_symmetry /* = std::nullopt */)
	{
		Coefficient ret(_name, _momenta, {}, true, false);
		ret.custom_symmetry = _custom_symmetry;
		return ret;
	}
	Coefficient Coefficient::RealInteraction(const std::string& _name, const MomentumList& _momenta, 
		const std::optional<std::function<void(Coefficient&)>>& _custom_symmetry /* = std::nullopt */)
	{
		assert(_momenta.size() == 3U);
		Coefficient ret(_name, _momenta, {}, false, false);
		ret.is_symmetrized_interaction = true;
		ret.custom_symmetry = _custom_symmetry;
		return ret;
	}
}