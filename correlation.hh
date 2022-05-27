
#ifndef CORRELATION
#define CORRELATION

#include <cmath>
#include <vector>
#include <stdexcept>
#include <string>
#include <map>
#include <utility>
#include <valarray>

namespace math::statistics
{

    // This call willl be used to store partial
    // results of calculus of Pearson Correlation
    // Coefficient on large datasets:
    //  https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
    template <typename T = double>
    struct pcc_partial {
        using value_type = T;

        pcc_partial() = default;
        pcc_partial(const pcc_partial<value_type>&) = default;
        pcc_partial(pcc_partial<value_type>&&) = default;

        auto& operator+=(const pcc_partial<value_type>& p) {
            this->count += p.count;
            this->sum_1 += p.sum_1;
            this->sum_2 += p.sum_2;
            this->sum_1_squared += p.sum_1_squared;
            this->sum_2_squared += p.sum_2_squared;
            this->sum_prod += p.sum_prod;
            return *this;
        }

        auto operator+(const pcc_partial<value_type>& p) const {
            auto ans = *this;
            ans += p;
            return ans;
        }

        // Calculate the Pearson Correlation Coefficient
        // with the accumulated data
        auto compute() const -> value_type {
            if (count == 0) {
                return 0;
            }
            const auto num = ( sum_prod - (sum_1*sum_2)/count );
            const auto den = (sum_1_squared - (sum_1*sum_1 / count)) * (sum_2_squared - (sum_2*sum_2 / count));
            // check for div by 0
            return den ?  num / std::sqrt( den ) : 0;
        }

        // Fastest way to add two more elements
        auto& accumulate(value_type v_1, value_type v_2) {
            this->sum_1 += v_1;
            this->sum_2 += v_2;
            this->sum_1_squared += v_1*v_1;
            this->sum_2_squared += v_2*v_2;
            this->sum_prod += v_1*v_2;
            ++this->count;
            return *this;
        }

        long long count{};
        value_type sum_1{};
        value_type sum_2{};
        value_type sum_1_squared{};
        value_type sum_2_squared{};
        value_type sum_prod{};
    };

    template <typename T, typename R = T>
    pcc_partial<R> pearson_correlation_coefficient(const std::vector<T>& v1, const std::vector<T>& v2) {
        using namespace std::literals;
        pcc_partial<R> ans;
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

    // This class has been conceived to easily
    // calculate PCC on all pairs of columns in
    // a large dataset
    template <typename T = double>
    class multicolumn_pcc_accumulator {
    public:
        using value_type = T;
    private:
        const int N;
        std::valarray<pcc_partial<value_type>> partials;
    public:
        multicolumn_pcc_accumulator(int N)
        : N{N}, partials((N-1)*N/2)
        {
            if (N < 2) {
                using namespace std::literals;
                std::invalid_argument("N must be at least 2, found "s + std::to_string(N));
            }
        }

        auto& accumulate(const std::vector<value_type>& row) {
            if (row.size() != N) {
                using namespace std::literals;
                throw std::runtime_error("Wrong number of columns received, expected "s + std::to_string(N) + " found "s + std::to_string(row.size()));
            }
            auto partial_iterator = std::begin(partials);
            for (int i{}; i!=N-1; ++i) {
                for (int j{i+1}; j!=N; ++j) {
                    partial_iterator->accumulate(row[i], row[j]);
                    ++partial_iterator;
                }
            }
            return *this;
        }

        auto& operator+=(const multicolumn_pcc_accumulator& o) {
            if (N != o.N) {
                using namespace std::literals;
                throw std::runtime_error("Size mismatch, this->N = "s + std::to_string(N) + ", other.N = "s + std::to_string(o.N));
            }
            partials += o.partials;
            return *this;
        }

        auto operator+(const multicolumn_pcc_accumulator& o) const {
            if (N != o.N) {
                using namespace std::literals;
                throw std::runtime_error("Size mismatch, this->N = "s + std::to_string(N) + ", other.N = "s + std::to_string(o.N));
            }
            auto ans = *this;
            ans += o;
            return ans;
        }

        // return a map containing all
        auto results() const {
            // use un'ordered map to sped up data access
            std::map<std::pair<int,int>,value_type> ans;
            auto partial_iterator = std::begin(partials);
            for (int position{}, i{}; i!=N-1; ++i) {
                for (int j{i+1}; j!=N; ++j, ++position) {
                    ans[std::make_pair(i,j)] = partial_iterator->compute();
                    ++partial_iterator;
                }
            }
            return ans;
        }
    };

} // namespace math::statistics


#endif
