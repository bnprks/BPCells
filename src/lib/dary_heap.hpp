// Copyright Malte Skarupke 2020.
// Distributed under the Boost Software License, Version 1.0.
// (See http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <utility>
#include <cstdint>
#include <functional>

namespace dary_heap {
namespace dary_heap_helpers
{
template<int D>
uint64_t first_child_index(uint64_t index)
{
    return index * D + 1;
}
template<int D>
uint64_t last_child_index(uint64_t index)
{
    return index * D + D;
}
template<int D>
uint64_t last_grandchild_index(uint64_t index)
{
    return index * (D * D) + (D * D + D);
}
template<int D>
uint64_t parent_index(uint64_t index)
{
    return (index - 1) / D;
}
template<int D>
uint64_t grandparent_index(uint64_t index)
{
    return (index - (D + 1)) / (D * D);
}
template<int D>
uint64_t index_with_no_grandchild(uint64_t length)
{
    return grandparent_index<D>(length - 1) + 1;
}
template<int D, typename It, typename Compare>
inline It largest_child(It first_child_it, Compare && compare)
{
    if (D == 1) // constexpr
        return first_child_it;
    else if (D == 2) // constexpr
        return first_child_it + !!compare(first_child_it[0], first_child_it[1]);
    else
    {
        It first_half_largest = largest_child<D / 2>(first_child_it, compare);
        It second_half_largest = largest_child<D - D / 2>(first_child_it + D / 2, compare);
        return compare(*first_half_largest, *second_half_largest) ? second_half_largest : first_half_largest;
    }
}
template<int D, typename It, typename Compare>
It largest_child(It first_child_it, int num_children, Compare && compare)
{
    if (D == 2) // constexpr
        return first_child_it;
    else if (D == 3) // constexpr
    {
        if (num_children == 1)
            return first_child_it;
        else
            return first_child_it + !!compare(first_child_it[0], first_child_it[1]);
    }
    else if (D == 4) // constexpr
    {
        switch (num_children)
        {
        case 1: return first_child_it;
        case 2: return first_child_it + !!compare(first_child_it[0], first_child_it[1]);
        default:
            It largest = first_child_it + !!compare(first_child_it[0], first_child_it[1]);
            return compare(*largest, first_child_it[2]) ? first_child_it + 2 : largest;
        }
    }
    else
    {
        switch(num_children)
        {
        case 1: return first_child_it;
        case 2: return first_child_it + !!compare(first_child_it[0], first_child_it[1]);
        case 3:
        {
            It largest = first_child_it + !!compare(first_child_it[0], first_child_it[1]);
            return compare(*largest, first_child_it[2]) ? first_child_it + 2 : largest;
        }
        case 4:
        {
            It largest_first_half = first_child_it + !!compare(first_child_it[0], first_child_it[1]);
            It largest_second_half = first_child_it + 2 + !!compare(first_child_it[2], first_child_it[3]);
            return compare(*largest_first_half, *largest_second_half) ? largest_second_half : largest_first_half;
        }
        default:
            int half = num_children / 2;
            It first_half_largest = largest_child<D>(first_child_it, half, compare);
            It second_half_largest = largest_child<D>(first_child_it + half, num_children - half, compare);
            return compare(*first_half_largest, *second_half_largest) ? second_half_largest : first_half_largest;
        }
    }
}
}

template<int D, typename It, typename Compare>
void make_heap(It begin, It end, Compare && compare)
{
    using std::swap;
    uint64_t length = end - begin;
    if (length <= 1)
        return;
    uint64_t index = (length - 2) / D;
    // optimization: there can be only one item that has fewer than D children
    // handling that item up front simplifies the second loop a little, since
    // we know that all other items have two children
    int num_children_end = (length - 1) % D;
    if (num_children_end)
    {
        It largest_child = dary_heap_helpers::largest_child<D>(begin + dary_heap_helpers::first_child_index<D>(index), num_children_end, compare);
        if (compare(begin[index], *largest_child))
            swap(begin[index], *largest_child);
        if (index == 0)
            return;
        --index;
    }
    // optimization: half of all the items will have no grandchildren. this
    // simplifies the push_down function a lot, so we handle these items
    // first. we could then do another optimization where we know that
    // after the first half, the next quarter of items has grandchildren but
    // no great-grandchildren, but the code is already big enough
    if (index > 0)
    {
        uint64_t lowest_index_with_no_grandchildren = dary_heap_helpers::index_with_no_grandchild<D>(length);
        for (;;)
        {
            It largest_child = dary_heap_helpers::largest_child<D>(begin + dary_heap_helpers::first_child_index<D>(index), compare);
            if (compare(begin[index], *largest_child))
                swap(begin[index], *largest_child);
            if (index-- == lowest_index_with_no_grandchildren)
                break;
        }
    }
    for (;; --index)
    {
        typename std::iterator_traits<It>::value_type value = std::move(begin[index]);
        uint64_t move_down_index = index;
        for (;;)
        {
            uint64_t last_child_index = dary_heap_helpers::last_child_index<D>(move_down_index);
            uint64_t first_child_index = last_child_index - (D - 1);
            It largest_child = begin;
            if (last_child_index < length)
                largest_child = dary_heap_helpers::largest_child<D>(begin + first_child_index, compare);
            else if (first_child_index >= length)
                break;
            else
                largest_child = dary_heap_helpers::largest_child<D>(begin + first_child_index, length - first_child_index, compare);
            if (!compare(value, *largest_child))
                break;
            begin[move_down_index] = std::move(*largest_child);
            move_down_index = largest_child - begin;
        }
        begin[move_down_index] = std::move(value);
        if (index == 0)
            break;
    }
}
// template<int D, typename It>
// void make_heap(It begin, It end)
// {
//     make_heap<D>(begin, end, std::less<>{});
// }

template<int D, typename It, typename Compare>
bool is_heap(It begin, It end, Compare && compare)
{
    uint64_t length = end - begin;
    for (uint64_t i = 1; i < length; ++i)
    {
        uint64_t parent = dary_heap_helpers::parent_index<D>(i);
        if (compare(begin[parent], begin[i]))
            return false;
    }
    return true;
}
// template<int D, typename It>
// bool is_heap(It begin, It end)
// {
//     return is_heap<D>(begin, end, std::less<>{});
// }

template<int D, typename It, typename Compare>
void push_heap(It begin, It end, Compare && compare)
{
    typename std::iterator_traits<It>::value_type value = std::move(end[-1]);
    uint64_t index = (end - begin) - 1;
    while (index > 0)
    {
        uint64_t parent = dary_heap_helpers::parent_index<D>(index);
        if (!compare(begin[parent], value))
            break;
        begin[index] = std::move(begin[parent]);
        index = parent;
    }
    begin[index] = std::move(value);
}

// template<int D, typename It>
// void push_heap(It begin, It end)
// {
//     return push_heap<D>(begin, end, std::less<>{});
// }


template<int D, typename It, typename Compare>
void pop_heap(It begin, It end, Compare && compare)
{
    uint64_t length = (end - begin) - 1;
    typename std::iterator_traits<It>::value_type value = std::move(end[-1]);
    end[-1] = std::move(begin[0]);
    uint64_t index = 0;
    for (;;)
    {
        uint64_t last_child = dary_heap_helpers::last_child_index<D>(index);
        uint64_t first_child = last_child - (D - 1);
        if (last_child < length)
        {
            It largest_child = dary_heap_helpers::largest_child<D>(begin + first_child, compare);
            if (!compare(value, *largest_child))
                break;
            begin[index] = std::move(*largest_child);
            index = largest_child - begin;
        }
        else if (first_child < length)
        {
            It largest_child = dary_heap_helpers::largest_child<D>(begin + first_child, length - first_child, compare);
            if (compare(value, *largest_child))
            {
                begin[index] = std::move(*largest_child);
                index = largest_child - begin;
            }
            break;
        }
        else
            break;
    }
    begin[index] = std::move(value);
}
// template<int D, typename It>
// void pop_heap(It begin, It end)
// {
//     return pop_heap<D>(begin, end, std::less<>{});
// }

} // end namespace dary_heap