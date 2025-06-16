// Copyright 2025 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
#include <cstdint>
#include <initializer_list>
#include <memory>
#include <string>
#include <vector>

namespace BPCells {

class Registry;
class AnalysisManager;

/// Data types for OpTree parameters
enum class ParamType : char {
    i32, /// int32_t
    i64, /// int64_t
    f32, /// float
    f64, /// double
};

/// Data types that can be produced/consumed in an OpTree
enum class OpType : char {
    Matrix,   /// Sparse matrix
    Fragment, /// Genomic fragments
};

class OpNode {
  class Param {
      protected:
        std::shared_ptr<const void> data_;
        std::string name_;
        ParamType type_;

      public:
        template <class T> bool is_type() const;

        template <class T> const T *get() const;

        std::string_view name() const;
    };

  public:
    std::string package;        /// Name of the package that defines this operation
    std::string op;             /// Operation name (i.e. log1p)
    std::vector<OpNode> inputs; /// Node inputs (zero or more)
    std::vector<Param> params;  /// Parameter values
};

class OpTree {
    OpNode root;
    Registry &registry;
    AnalysisManager &AM;
};

/// Specification of a parameter's length
class LengthSpec {
  public:
    /// Require a fixed length
    struct Fixed {uint32_t len;};
    /// Sources for data-dependent lengths
    enum LengthSource {
      Rows, /// number of matrix rows
      Cols, /// number of matrix cols
    };
    /// Match the length of some data-dependent attribute
    struct DataDependent {LengthSource type;};

};

/// Specification of a parameter's type, name, and error-checking metadata
template<typename T>
class ParamSpec {
  public:
    ParamSpec(const std::string &name);
    void set_length(const LengthSpec &length); 
};

/// Specification of an operator's type and parameters
class OpSpec {
  public:
    OpSpec(const std::string &name);

    void set_inputs(std::initializer_list<OpType>);
    void set_output(OpType);

    template <typename T>
    void add_param(const std::string &name, const ParamSpec<T> &param);
    void set_executor(void *(*fn)(const OpTree &tree));
};


class MatrixOpSpec : OpSpec {
  public:
    // Optionally add custom functions to infer types and dims
};

class FragmentOpSpec : OpSpec {
  // Optionally add custom functions to infer cell + chromosome counts
};


class Registry {
  public:
    void addOp(const OpSpec &op);

    template <typename T>
    T execute(const OpTree &tree);


};

} // namespace BPCells