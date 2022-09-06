
#ifndef CORRELATION
#define CORRELATION

#include <cmath>
#include <vector>
#include <stdexcept>
#include <string>
#include <map>
#include <utility>
#include <valarray>

namespace math
{
    namespace statistics
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

            pcc_partial& operator=(const pcc_partial<value_type>&) = default;
            pcc_partial& operator=(pcc_partial<value_type>&&) = default;

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
        inline pcc_partial<R> pearson_correlation_coefficient(const std::vector<T>& v1, const std::vector<T>& v2) {
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

        template <typename T, typename R = T>
        inline pcc_partial<R> pearson_correlation_coefficient(const T* v1, const T* v2, std::size_t size) {
            pcc_partial<R> ans;
            for (decltype(size) i{}; i!=size; ++i) {
                const auto v_1 = v1[i];
                const auto v_2 = v2[i];
                ans.sum_1 += v_1;
                ans.sum_2 += v_2;
                ans.sum_1_squared += v_1*v_1;
                ans.sum_2_squared += v_2*v_2;
                ans.sum_prod += v_1*v_2;
            }
            ans.count = size;
            return ans;
        }

        // to be used when data are organized by rows and column elements
        // are not stored consecutively. Scatter should be equal to size(row)
        template <typename T, typename R = T>
        inline pcc_partial<R> pearson_correlation_coefficient_scattered(const T* v1, const T* v2, std::size_t size, std::size_t scatter = 1) {
            pcc_partial<R> ans;
            if (scatter == 0) {
                throw std::invalid_argument("Scatter must be graeter than 0");
            }
            for (decltype(size) i{}; i!=size; ++i) {
                const auto v_1 = v1[i*scatter];
                const auto v_2 = v2[i*scatter];
                ans.sum_1 += v_1;
                ans.sum_2 += v_2;
                ans.sum_1_squared += v_1*v_1;
                ans.sum_2_squared += v_2*v_2;
                ans.sum_prod += v_1*v_2;
            }
            ans.count = size;
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
            // number of columns to analyze
            const int N;
            // it is unnecessary to calculate the mean and the
            // mean of the squared values of each column once per
            // pair of columns; it is sufficient to calculate it
            // just once per column.
            // for each column
            std::valarray<value_type> totals;           // sum of element in each column
            std::valarray<value_type> squared_totals;   // sum of squared element in each column
            // for each pair of columns - (0,1) (0,2) (0,3) (1,2) (1,3) (2,3)
            std::valarray<value_type> covariance_total;
            long long int count{};  // total number of rows
        public:
            multicolumn_pcc_accumulator(int N)
            : N{N}, totals(N), squared_totals(N), covariance_total((N-1)*N/2)
            {
                if (N < 2) {
                    using namespace std::literals;
                    std::invalid_argument("N must be at least 2, found "s + std::to_string(N));
                }
            }

            // count a single row
            auto& accumulate(const std::vector<value_type>& row) {
                if (row.size() != N) {
                    using namespace std::literals;
                    throw std::runtime_error("Wrong number of columns received, expected "s + std::to_string(N) + " found "s + std::to_string(row.size()));
                }
                auto covariance_iterator = std::begin(covariance_total);
                // first column
                for (int i{}; i!=N-1; ++i) {
                    const auto tmp = row[i];
                    totals[i] += tmp;
                    squared_totals[i] += tmp*tmp;
                    for (int j{i+1}; j!=N; ++j) {
                        *covariance_iterator += tmp*row[j];
                        ++covariance_iterator;
                    }
                }
                const auto tmp = row[N-1];
                totals[N-1] = tmp;
                squared_totals[N-1] = tmp*tmp;
                ++count;
                return *this;
            }

            // 
            auto& accumulate(
                const value_type* matrix,   // matrix containing the chunk to analyze
                std::size_t rows,           // rows in the matrix
                std::size_t cols,           // columns in the matrix
                std::size_t row_offset,     // offset between same index
                                            // elements of two adjacent rows
                std::size_t col_offset      // offset between same index
                                            // elements of two adjacent column
            ) {
                if (cols != (decltype(cols))N) {
                    using namespace std::literals;
                    throw std::runtime_error("Wrong number of columns received, expected "s + std::to_string(N) + " found "s + std::to_string(cols));
                }
                // offset between same index
                // elements of two adjacent column
                //const auto col_offset = /*sizeof(value_type) */ (by_columns ? rows : 1);
                // offset between same index
                // elements of two adjacent rows
                //const auto row_offset = /*sizeof(value_type) */ (by_columns ? 1 : cols);
                {
                    auto column = matrix; // pointer to first element of next column
                    // for each column
                    for (std::size_t c{}; c!=cols; ++c) {
                        // accumulators
                        value_type tot{};   // sum(row)
                        value_type tot2{};  // sum(row.^2)
                        auto row = column;    // point to first element of the column
                        // for each row
                        for (std::size_t r{}; r!=rows; ++r) {
                            const auto tmp = *row;
                            tot += tmp;
                            tot2 += tmp*tmp;
                            // point to next row
                            row += row_offset;
                        }
                        // store partials - once to limit memory deferencing
                        totals[c] += tot;
                        squared_totals[c] += tot2;
                        // point to next column
                        column += col_offset;
                    }
                }
                // for each couple of columns
                {
                    std::size_t couple_idx {};
                    // ptr to c1
                    auto column1 = matrix;
                    for (std::size_t c1{}; c1 != cols; ++c1) {
                        // ptr to c2
                        auto column2 = column1 + col_offset;
                        for (std::size_t c2{c1+1}; c2 != cols; ++c2) {
                            // accumulator
                            value_type cross{};
                            // for each couple (col1[i],col2[i])
                            auto row1 = column1;
                            auto row2 = column2;
                            for (std::size_t r{}; r!=rows; ++r) {
                                cross += (*row1)*(*row2);
                                // next pair
                                row1 += row_offset;
                                row2 += row_offset;
                            }
                            // store partial
                            covariance_total[couple_idx++] += cross;   
                            column2 += col_offset;
                        }
                        column1 += col_offset;
                    }
                }
                // +rows rows to be considered to calculate means
                count += rows;
                return *this;
            }

            auto& operator+=(const multicolumn_pcc_accumulator<T>& o) {
                if (N != o.N) {
                    using namespace std::literals;
                    throw std::runtime_error("Size mismatch, this->N = "s + std::to_string(N) + ", other.N = "s + std::to_string(o.N));
                }
                //partials += o.partials;
                totals += o.totals;
                squared_totals += o.squared_totals;
                covariance_total += o.covariance_total;
                count += o.count;
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

            // convert content to valarray of pcc_partial to adapt to
            // existing code.
            std::valarray<pcc_partial<value_type>> to_pcc_partial_valarray() const {
                std::valarray<pcc_partial<value_type>> ans(covariance_total.size());

                for (int couple_idx{}, c1{}; c1 != N; ++c1) {
                    const auto totals_c1 = totals[c1];
                    const auto squared_totals_c1 = squared_totals[c1];
                    for (int c2{c1+1}; c2 != N; ++c2) {
                        // store partial
                        pcc_partial<value_type> tmp;
                        tmp.count = count;
                        tmp.sum_1 = totals_c1;
                        tmp.sum_2 = totals[c2];
                        tmp.sum_1_squared = squared_totals_c1;
                        tmp.sum_2_squared = squared_totals[c2];
                        tmp.sum_prod = covariance_total[couple_idx];
                        ans[couple_idx++] = tmp;
                    }
                }

                return ans;
            }

            // return a map containing all
            auto results() const {
                auto partials = to_pcc_partial_valarray();
                // use unordered map to sped up data access
                std::map<std::pair<int,int>,value_type> ans;
                auto partial_iterator = std::begin(partials);
                for (int i{}; i!=N-1; ++i) {
                    for (int j{i+1}; j!=N; ++j) {
                        ans[std::make_pair(i,j)] = partial_iterator->compute();
                        ++partial_iterator;
                    }
                }
                return ans;
            }
        };

    } // namespace statistics
} // namespace math


#endif
