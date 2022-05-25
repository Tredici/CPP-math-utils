
#ifndef CORRELATION
#define CORRELATION

#include <cmath>
#include <vector>
#include <stdexcept>
#include <string>

namespace math::statistics
{

    // This call willl be used to store partial
    // results of calculus of Pearson Correlation
    // Coefficient on large datasets:
    //  https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
    struct pcc_partial {
        pcc_partial() = default;
        pcc_partial(const pcc_partial&) = default;
        pcc_partial(pcc_partial&&) = default;

        auto& operator+=(const pcc_partial& p) {
            this->count += p.count;
            this->sum_1 += p.sum_1;
            this->sum_2 += p.sum_2;
            this->sum_1_squared += p.sum_1_squared;
            this->sum_2_squared += p.sum_2_squared;
            this->sum_prod += p.sum_prod;
            return *this;
        }

        auto operator+(const pcc_partial& p) const {
            auto ans = *this;
            ans += p;
            return ans;
        }

        // Calculate the Pearson Correlation Coefficient
        // with the accumulated data
        auto compute() const {
            return ( sum_prod - (sum_1*sum_2)/count )
                / sqrt( (sum_1_squared - (sum_1*sum_1 / count)) * (sum_2_squared - (sum_2*sum_2 / count)) );
        }

        long long count{};
        double sum_1{};
        double sum_2{};
        double sum_1_squared{};
        double sum_2_squared{};
        double sum_prod{};
    };

    template <typename T>
    pcc_partial pearson_correlation_coefficient(const std::vector<T>& v1, const std::vector<T>& v2) {
        using namespace std::literals;
        pcc_partial ans;
        if (v1.size() != v2.size()) {
            throw std::invalid_argument("Arguments must have the same length, found len(v1)="s + std::to_string(v1.size()) + ", len(v2)=" + std::to_string(v2.size()));
        }
        for (decltype(v1.size()) i{}; i!=v1.size(); ++i) {
            ans.sum_1 += v1[i];
            ans.sum_2 += v2[i];
            ans.sum_1_squared += v1[i]*v1[i];
            ans.sum_2_squared += v2[i]*v2[i];
            ans.sum_prod += v1[i]*v2[i];
        }
        ans.count = v1.size();
        return ans;
    }

} // namespace math::statistics


#endif
