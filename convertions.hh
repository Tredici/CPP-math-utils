
#ifndef CONVERTIONS
#define CONVERTIONS

#include <string>
#include <stdexcept>

/**
 * This namespace contains a set of template
 * functions to cast std::string
 * to numbers.
 */
namespace math::convertions {

    template <typename T>
    T ston(const std::string& s) {
        return (T)s;
    }

    template <>
    long double ston(const std::string& s) {
        return std::stold(s);
    }

    template <>
    double ston(const std::string& s) {
        return std::stod(s);
    }

    template <>
    float ston(const std::string& s) {
        return std::stof(s);
    }

    template <>
    int ston(const std::string& s) {
        return std::stoi(s);
    }

    template <>
    unsigned int ston(const std::string& s) {
        auto cast = std::stol(s);
        if (cast != (unsigned int)cast) {
            using namespace std::literals;
            throw std::out_of_range("cannot cast '"s + s + "' to unsigned int"s);
        }
        return cast;
    }

    template <>
    long ston(const std::string& s) {
        return std::stol(s);
    }

    template <>
    unsigned long ston(const std::string& s) {
        return std::stoul(s);
    }

    template <>
    long long ston(const std::string& s) {
        return std::stoll(s);
    }

    template <>
    unsigned long long ston(const std::string& s) {
        return std::stoull(s);
    }

} // namespace math::convertions

#endif
