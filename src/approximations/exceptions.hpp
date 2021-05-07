/*!
    \file execptions.hpp
 */

#ifndef kappa_exceptions_hpp
#define kappa_exceptions_hpp

#include <exception>
#include <string>

namespace kappa {

class UnopenedFileException : public std::runtime_error {

 public:
  UnopenedFileException(std::string message):std::runtime_error(message) {}
};

class DataNotFoundException : public std::runtime_error {
 public:
  DataNotFoundException(std::string message):std::runtime_error(message) {}
};

class ModelParameterException : public std::runtime_error {
 public:
  ModelParameterException(std::string message):std::runtime_error(message) {}
};

class IncorrectValueException : public std::runtime_error {
 public:
  IncorrectValueException(std::string message):std::runtime_error(message) {}
};
}
#endif /* kappa_exceptions_hpp */
