#include <mrock/symbolic_operators/TermLoader.hpp>
#include <fstream>
#include <boost/archive/binary_iarchive.hpp>
#include <filesystem>

namespace mrock::symbolic_operators {
	void TermLoader::load(std::string const& folder, bool use_XP, int n_terms, int start_at/*=0*/) {
		M.resize(n_terms * n_terms);
		N.resize(n_terms * n_terms);
		const int name_offset = (start_at < 0) ? 0 : start_at;

		for (int i = 0; i < n_terms; i++)
		{
			for (int j = 0; j < n_terms; j++)
			{
				const std::string suffix = std::to_string(j + name_offset) + "_" + std::to_string(i + name_offset) + ".bin";
				{
					const auto filename = folder + (use_XP ? "XP_" : "") + "wick_M_" + suffix;
					if (!std::filesystem::exists(filename)) {
						throw std::runtime_error("Wick: FileNotFound: " + filename);
					}
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename, std::ios::binary);
					if(ifs.good()) {
						boost::archive::binary_iarchive ia(ifs);
						this->M[j * n_terms + i].clear();
						ia >> this->M[j * n_terms + i];
						ifs.close();
					}
					else {
						throw std::runtime_error("Inputstream for " + filename + " is bad!");
					}
				}
				{
					const auto filename = folder + (use_XP ? "XP_" : "") + "wick_N_" + suffix;
					if (!std::filesystem::exists(filename)) {
						throw std::runtime_error("Wick: FileNotFound: " + filename);
					}
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename, std::ios::binary);
					if(ifs.good()) {
						boost::archive::binary_iarchive ia(ifs);
						this->N[j * n_terms + i].clear();
						ia >> this->N[j * n_terms + i];
						ifs.close();
					}
					else {
						throw std::runtime_error("Inputstream for " + filename + " is bad!");
					}
				}
			}
		}
	}
}