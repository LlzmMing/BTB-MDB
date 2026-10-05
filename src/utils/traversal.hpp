#ifndef TRAVERSAL_HPP
#define TRAVERSAL_HPP

#include <iterator>

namespace traversal {

    template <typename T>
    class ReverseWrapper {
        T& external_iterable;
    public:
        explicit ReverseWrapper(T& iterable): external_iterable(iterable) {}
        auto begin() { return std::rbegin(external_iterable); }
        auto end() { return std::rend(external_iterable); }
    };

    template <typename T>
    static ReverseWrapper<T> reverse(T& container) {
        return ReverseWrapper<T>(container);
    }

}

#endif
