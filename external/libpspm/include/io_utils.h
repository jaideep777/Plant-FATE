#ifndef VECTOR_IO_UTILS_H_
#define VECTOR_IO_UTILS_H_

#include <ostream>
#include <istream>
#include <tuple>
#include <vector>
#include <array>
#include <iomanip>
#include <string_view>

// print a vector via ofstream
// prints: size | v1 v2 v3 ...
template <class T>
std::ostream& operator << (std::ostream &os, const std::vector<T> &v) {
	//os << std::setprecision(12);
	os << v.size() << " | ";
	for (const auto &x : v) {
		os << x << ' ';
	}
	// os << '\n'; // .FIXME: remove this newline and insert newline in every save() function
	return os;
}


// read a vector via ofstream
// reads size, "|", discards "|",
// then resizes the vector and reads values
template <class T>
std::istream& operator >> (std::istream &is, std::vector<T> &v) {
	int n; 
	std::string s;
	is >> n >> s; // s contains "|"
	// std::cout << "restoring vector: " << n << " " << s << std::endl;
	v.resize(n);
	for (int i=0; i<n; ++i) is >> v[i];
	return is;
}


// print an array via ofstream
// prints: size | v1 v2 v3 ...
template <class T, size_t n>
std::ostream& operator << (std::ostream &os, const std::array<T,n> &v) {
	//os << std::setprecision(12);
	for (const auto &x : v) {
		os << x << ' ';
	}
	// os << '\n'; // .FIXME: remove this newline and insert newline in every save() function
	return os;
}


// read an array via ofstream
// reads size, "|", discards "|",
// then resizes the vector and reads values
template <class T, size_t n>
std::istream& operator >> (std::istream &is, std::array<T,n> &v) {
	for (int i=0; i<n; ++i) is >> v[i];
	return is;
}


// print a tuple (requires C++17)
template<typename... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...> t) {
	apply([&](auto&&... args) {
		((os << args << " "), ...);
	}, t);
	return os;
}


// Convert an array to a vector
template <class T, size_t n>
inline std::vector<T> to_vector(const std::array<T,n>& a){
	return std::vector<T>(a.begin(), a.end());
}


// Get the type of an object as a string
template <typename T>
constexpr auto type_name() {
	std::string_view name, prefix, suffix;
#ifdef __clang__
	name = __PRETTY_FUNCTION__;
	prefix = "auto type_name() [T = ";
	suffix = "]";
#elif defined(__GNUC__)
	name = __PRETTY_FUNCTION__;
	prefix = "constexpr auto type_name() [with T = ";
	suffix = "]";
#elif defined(_MSC_VER)
	name = __FUNCSIG__;
	prefix = "auto __cdecl type_name<";
	suffix = ">(void)";
#endif
	name.remove_prefix(prefix.size());
	name.remove_suffix(suffix.size());
	return name;
}


#endif

