
#ifndef CONVERTIONS
#define CONVERTIONS

#include <string>
#include <stdexcept>

/**
 * This namespace contains a set of template
 * functions to cast std::string
 * to numbers.
 */
namespace math {
    namespace convertions {

        template <typename T>
        inline T ston(const std::string& s) {
            return (T)s;
        }

        template <>
        inline long double ston(const std::string& s) {
            return std::stold(s);
        }

        template <>
        inline double ston(const std::string& s) {
            return std::stod(s);
        }

        template <>
        inline float ston(const std::string& s) {
            return std::stof(s);
        }

        template <>
        inline int ston(const std::string& s) {
            return std::stoi(s);
        }

        template <>
        inline unsigned int ston(const std::string& s) {
            auto cast = std::stol(s);
            if (cast != (unsigned int)cast) {
                using namespace std::literals;
                throw std::out_of_range("cannot cast '"s + s + "' to unsigned int"s);
            }
            return cast;
        }

        template <>
        inline long ston(const std::string& s) {
            return std::stol(s);
        }

        template <>
        inline unsigned long ston(const std::string& s) {
            return std::stoul(s);
        }

        template <>
        inline long long ston(const std::string& s) {
            return std::stoll(s);
        }

        template <>
        inline unsigned long long ston(const std::string& s) {
            return std::stoull(s);
        }

    } // namespace convertions
} // namespace math

#endif
