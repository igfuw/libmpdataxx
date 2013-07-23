#pragma once

#include <boost/exception/all.hpp>
#include <string>

// TODO: chapter on error handling in libmpdataxx-paper (let theuser know we're throwing boost exceptions)

namespace libmpdataxx
{
  struct error: virtual boost::exception, virtual std::exception 
  { 
    // fields
    std::string msg;

    // ctors
    error() {}
    error(const char* c_str) : msg(c_str) {}
    error(const std::string &str) : msg(str) {}

    // methods
    const char* what() const throw()
    {
      return msg.c_str();
    }
  };
}
