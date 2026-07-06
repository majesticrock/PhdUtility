// Serialization support; these are the components of boost needed to use the program
// If you plan on using serialization using boost, you can simply include this header file
// instead of including the headers below individually.
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>