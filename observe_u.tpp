#include <cmath>
#include <iostream>
#include <nlohmann/json.hpp>
#include <boost/numeric/odeint.hpp>


template <typename Eq, typename Obs>
void ObserveU<Eq, Obs>::operator() () {
    eq.initU(u, step);
    eq.initTu(stpr.prevdu, step);
    boost::numeric::odeint::integrate_adaptive(stpr, eq, u, rMin, rMax, step, obs);
}

template <typename Eq, typename Obs>
ObserveU<Eq, Obs>::ObserveU () {}

template <typename Eq, typename Observer>
void from_json(const nlohmann::json &j, ObserveU<Eq, Observer>& p) {
    try {
        p.eq = j.at("eq");
    } catch(nlohmann::json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(nlohmann::json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
    
    try {
        p.rMin = j.at("rMin").get<double>();
    } catch(nlohmann::json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(nlohmann::json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
    
    try {
        p.rMax = j.at("rMax").get<double>();
    } catch(nlohmann::json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(nlohmann::json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }

    try {
        p.step = j.at("step").get<double>();
    } catch(nlohmann::json::type_error& e) {
        std::cerr << e.what() << std::endl;
    } catch(nlohmann::json::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }
}
