
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
            if (count == 0) {
                return 0.0;
            }
            const auto num = ( sum_prod - (sum_1*sum_2)/count );
            const auto den = (sum_1_squared - (sum_1*sum_1 / count)) * (sum_2_squared - (sum_2*sum_2 / count));
            // check for div by 0
            return den ?  num / sqrt( den ) : 0;
        }

        // Fastest way to add two more elements
        auto& accumulate(double v_1, double v_2) {
            this->sum_1 += v_1;
            this->sum_2 += v_2;
            this->sum_1_squared += v_1*v_1;
            this->sum_2_squared += v_2*v_2;
            this->sum_prod += v_1*v_2;
            ++this->count;
            return *this;
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
            const auto v_1 = v1[i];
            const auto v_2 = v2[i];
            ans.sum_1 += v_1;
            ans.sum_2 += v_2;
            ans.sum_1_squared += v_1*v_1;
            ans.sum_2_squared += v_2*v_2;
            ans.sum_prod += v_1*v_2;
        }
        ans.count = v1.size();
        return ans;
    }

} // namespace math::statistics


#endif
