//=============================================================================
// Author: Joren Heit 
//=============================================================================

#ifndef MULTIARRAY_H
#define MULTIARRAY_H

#include <vector>
#include <algorithm> // std::fill

template <typename T, std::size_t D, template <typename...> class Container = std::vector>
class MultiArray
{
    Container<T> d_data;
    std::size_t d_dimensions[D];

public:
    // Constructor
    template <typename ... Args>
    MultiArray(Args ... args):
        d_data(product(args...)),
        d_dimensions{static_cast<std::size_t>(args)...}
    {
        static_assert(sizeof ... (args) == D, "Number of dimensions does not match number of constructor-args");
    }

    // Index operators
    template <typename ... Args>
    T &operator()(Args ... args)
    {
        return d_data[indexMap(static_cast<std::size_t>(args) ...)];
    }

    template <typename ... Args>
    T const &operator()(Args ... args) const
    {
        return d_data[indexMap(static_cast<std::size_t>(args) ...)];
    }

    // size information
    constexpr std::size_t dim() const
    {
        return D;
    }

    std::size_t size(std::size_t dim) const
    {
        return d_dimensions[dim];
    }

    
    // Range assignment
    MultiArray &assign(size_t const (&ranges)[D][2], T const &val)
    {
        Assign<0, decltype(*this)>::assign(*this, ranges, val);
        return *this;
    }

    MultiArray &assign(T const &val)
    {
        std::fill(d_data.begin(), d_data.end(), val);
        return *this;
    }

private:
    // Helper functions
    static std::size_t product();
    template <typename ... Tail>
    static std::size_t product(std::size_t head, Tail ... tail)
    {
        return (sizeof ... (Tail) == 0) 
            ? head 
            : head * product(tail ...);
    }

    template <std::size_t Count>
    std::size_t dimProduct() const
    {
        std::size_t ret = 1;
        for (std::size_t idx = D - Count; idx != D; ++idx)
            ret *= d_dimensions[idx];
        return ret;
    }


    std::size_t indexMap() const;
    template <typename ... Tail>
    std::size_t indexMap(std::size_t head, Tail ... tail) const
    {
        enum { count = sizeof ... (tail) };
        return head * dimProduct<count>() + ((count == 0) 
                                             ? 0
                                             : indexMap(tail ...));
        
    }

    // Range assignment
    template <size_t Level, typename Self>
    struct Assign
    {
        template <typename ... Indices>
        static void assign(Self &self, size_t const (&ranges)[D][2], T const &val, Indices ... indices)
        {
            for (size_t i = ranges[Level][0]; i != ranges[Level][1]; ++i)
                Assign<Level + 1, Self>::assign(self, ranges, val, indices..., i);
        }
    };

    template <typename Self>
    struct Assign<D, Self>
    {
        template <typename ... Indices>
        static void assign(Self &self, size_t const (&ranges)[D][2], T const &val, Indices ... indices)
        {
            self(indices ...) = val;
        }
    };
};

// convenient factory
template <typename T, typename ... Args>
inline MultiArray<T, sizeof ... (Args)> makeMultiArray(Args ... args)
{
    return MultiArray<T, sizeof ... (args)>(args ...);
}

#endif
