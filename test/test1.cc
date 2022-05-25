#include "../modules/CPP-test-unit/tester.hh"
#include "../correlation.hh"

#include <vector>
#include <random>
#include <stdexcept>
#include <string>

using namespace std::literals;

tester t1([](){
    constexpr int N = 157;
    constexpr double expexted = 1.0;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(30,77);
    std::vector<double> v; v.reserve(N);
    for (int i{}; i!=N; ++i) {
        v.push_back(distribution(generator)); 
    }
    auto ans = math::statistics::pearson_correlation_coefficient(v, v).compute();
    if (ans != expexted) {
        throw std::runtime_error("Expected "s + std::to_string(expexted) + ", found "s + std::to_string(ans));
    }
});

tester t2([](){
    constexpr int N = 157;
    constexpr double expexted = -1.0;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(30,77);
    std::vector<double> v1, v2; v1.reserve(N); v2.reserve(N);
    for (int i{}; i!=N; ++i) {
        auto val = distribution(generator);
        v1.push_back(val); 
        v2.push_back(-val); 
    }
    auto ans = math::statistics::pearson_correlation_coefficient(v1, v2).compute();
    if (ans != expexted) {
        throw std::runtime_error("Expected "s + std::to_string(expexted) + ", found "s + std::to_string(ans));
    }
});


tester t3([](){
    constexpr int N = 157;
    constexpr double expexted = 0.0;
    std::vector<double> v; v.reserve(N);
    for (int i{}; i!=N; ++i) {
        v.push_back(0.0); 
    }
    auto ans = math::statistics::pearson_correlation_coefficient(v, v).compute();
    if (ans != expexted) {
        throw std::runtime_error("Expected "s + std::to_string(expexted) + ", found "s + std::to_string(ans));
    }
});

