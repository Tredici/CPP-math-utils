
#ifndef COUPLE
#define COUPLE

#include <ostream>
#include <utility>
#include <optional>
#include <stdexcept>

namespace math::sets
{

    /**
     * Given a natural number n, this class is
     * used to generate all sets of two elements
     * both lower than n
     */
    class couple
    {
    private:
        const int n, limit = n*(n-1)/2;
        int i;
        int _1, _2;
        // if true no more inc() are possible
        bool _last{}, _finished{};

        static std::pair<int, int> pair(int n, int i) {
            // candidate supposing al pairs are ok
            std::pair<int, int> p{i/n, i%n};

            // first column and no overflow?
            if (p.first == 0 && p.second+1 < n) {
                return std::make_pair(p.first, p.second+1);
            } else if (p.first == 0 && p.second+1 == n) {
                return std::make_pair(1, 2);
            }

            // reduce problem with recursion
            // [0,1] for new base will be
            // translated to [p[0], p[0]+1]
            std::pair<int, int> base(p.first, p.first);
            // all points in the triangle
            // marked by [p[0], p[0]]
            // must be ignored, others must be
            // counted
            auto remaining = i - ((n-1)*p.first - (p.first-1)*p.first/2);
            auto p2 = pair(n-p.first, remaining);
            return std::make_pair(p2.first+base.first, p2.second+base.second);
        }
    public:
        couple(int n, int i = 0) : n{n}, i{i} {
            using namespace std::literals;
            if (n < 2) {
                throw std::invalid_argument("n must be at least 2, received "s + std::to_string(n));
            }
            if (i < 0 || limit <= i) {
                throw std::invalid_argument("i ("s + std::to_string(i) + ") must be in [0,"s + std::to_string(limit) + ")"s);
            }
            auto tmp = pair(n, i);
            _1 = tmp.first; _2 = tmp.second;
            if (n == 2) {
                _last = true;
            }
        }

        int first() const { return _1; }
        int second() const { return _2; }

        operator bool() const { return !_finished; }

        auto& inc() {
            if (_last) {
                _finished = true;
                return *this;
            }
            this->i++;
            if (_2+1 == n) {
                _2 = ++_1 + 1;
            } else {
                ++_2;
            }
            // reached last pair
            if (_2 == n-1 && _1 + 1 == _2) {
                _last = true;
            }
            return *this;
        }

        auto operator++() -> decltype(*this) {
            return inc();
        }

        auto operator++(int) {
            couple cp(*this);
            inc();
            return cp;
        }

        auto operator==(const couple& c) const {
            return this->n == c.n && this->i == c.i;
        }

        operator std::pair<int,int>() const {
            return as_pair();
        }

        std::pair<int,int> as_pair() const {
            return {_1,_2};
        }

        bool last() const { return _last; }
        bool finished() const { return _finished; }

        auto index() const {
            return this->i;
        }
    };

} // namespace math

inline std::ostream& operator<<(std::ostream& os, const math::sets::couple &p) {
    return os << "(" << p.first() << "," << p.second() << ")";
}

#endif

