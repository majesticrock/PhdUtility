#include "HubbardDispersions.hpp"

namespace SymbolicOperators {
	std::vector<std::vector<Term>> HubbardDispersions::XP_basis() const
	{
		std::vector<std::vector<Term>> ret = Hubbard::XP_basis();
		for(auto& _v : ret){
			for(auto& v : _v){
				if(v.operators.front().isDaggered) {
					v.operators.front().momentum += Momentum('x');
				}
				else {
					v.operators.front().momentum += Momentum('x', -1);
				}
			}
		}
		return ret;
	}
	std::vector<std::vector<Term>> HubbardDispersions::STD_basis() const
	{
		std::vector<std::vector<Term>> ret = Hubbard::STD_basis();
		for(auto& _v : ret){
			for(auto& v : _v){
				if(v.operators.front().isDaggered) {
					v.operators.front().momentum += Momentum('x');
				}
				else {
					v.operators.front().momentum += Momentum('x', -1);
				}
			}
		}
		return ret;
	}

	std::string HubbardDispersions::get_subfolder() const
	{
		return "hubbard/dispersions/";
	}
}