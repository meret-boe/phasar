/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * Table.h
 *
 *  Created on: 07.11.2016
 *      Author: pdschbrt
 */

#ifndef PHASAR_UTILS_TABLE_H_
#define PHASAR_UTILS_TABLE_H_

#include <algorithm>
#include <ostream>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>

// we may wish to replace this by boost::multi_index at some point

namespace psr {

template <typename R, typename C, typename V> class Table {
private:
  std::unordered_map<R, std::unordered_map<C, V>> Tabl;

public:
  struct Cell {
    Cell() = default;
    Cell(R Row, C Col, V Val) : Row(Row), Col(Col), Val(Val) {}
    ~Cell() = default;
    Cell(const Cell &) = default;
    Cell &operator=(const Cell &) = default;
    Cell(Cell &&) noexcept = default;
    Cell &operator=(Cell &&) noexcept = default;

    R getRowKey() const { return Row; }
    C getColumnKey() const { return Col; }
    V getValue() const { return Val; }

    friend std::ostream &operator<<(std::ostream &Os, const Cell &Ce) {
      return Os << "Cell: " << Ce.Row << ", " << Ce.C << ", " << Ce.V;
    }
    friend bool operator<(const Cell &Lhs, const Cell &Rhs) {
      return std::tie(Lhs.Row, Lhs.Col, Lhs.Val) < std::tie(Rhs.Row, Rhs.Col, Rhs.Val);
    }
    friend bool operator==(const Cell &Lhs, const Cell &Rhs) {
      return std::tie(Lhs.Row, Lhs.Col, Lhs.Val) == std::tie(Rhs.Row, Rhs.Col, Rhs.Val);
    }

  private:
    R Row;
    C Col;
    V Val;
  };

  Table() = default;
  Table(const class Table &T) = default;
  class Table &operator=(const class Table &T) = default;
  Table(class Table &&T) noexcept = default;
  class Table &operator=(class Table &&T) noexcept = default;
  ~Table() = default;

  void insert(R Ro, C Co, V Va) {
    // Associates the specified value with the specified keys.
    Tabl[Ro][Co] = std::move(Va);
  }

  void insert(const Table &T) {
    Tabl.insert(T.Tabl.begin(), T.Tabl.end());
  }

  void clear() { Tabl.clear(); }

  [[nodiscard]] bool empty() const { return Tabl.empty(); }

  [[nodiscard]] size_t size() const { return Tabl.size(); }

  [[nodiscard]] std::set<Cell> cellSet() const {
    // Returns a set of all row key / column key / value triplets.
    std::set<Cell> S;
    for (const auto &M1 : Tabl) {
      for (const auto &M2 : M1.second) {
        S.emplace(M1.first, M2.first, M2.second);
      }
    }
    return S;
  }

  [[nodiscard]] std::vector<Cell> cellVec() const {
    // Returns a vector of all row key / column key / value triplets.
    std::vector<Cell> Va;
    for (const auto &M1 : Tabl) {
      for (const auto &M2 : M1.second) {
        Va.emplace_back(M1.first, M2.first, M2.second);
      }
    }
    return Va;
  }

  [[nodiscard]] std::unordered_map<R, V> column(C ColumnKey) const {
    // Returns a view of all mappings that have the given column key.
    std::unordered_map<R, V> Column;
    for (const auto &Row : Tabl) {
      if (Row.second.count(ColumnKey)) {
        Column[Row.first] = Row.second[ColumnKey];
      }
    }
    return Column;
  }

  [[nodiscard]] std::multiset<C> columnKeySet() const {
    // Returns a set of column keys that have one or more values in the table.
    std::multiset<C> Colkeys;
    for (const auto &M1 : Tabl) {
      for (const auto &M2 : M1.second) {
        Colkeys.insert(M2.first);
      }
    }
    return Colkeys;
  }

  [[nodiscard]] std::unordered_map<C, std::unordered_map<R, V>>
  columnMap() const {
    // Returns a view that associates each column key with the corresponding map
    // from row keys to values.
    std::unordered_map<C, std::unordered_map<R, V>> Columnmap;
    for (const auto &M1 : Tabl) {
      for (const auto &M2 : Tabl.second) {
        Columnmap[M2.first][M1.first] = M2.second;
      }
    }
    return Columnmap;
  }

  [[nodiscard]] bool contains(R RowKey, C ColumnKey) const {
    // Returns true if the table contains a mapping with the specified row and
    // column keys.
    if (auto RowIter = Tabl.find(RowKey); RowIter != Tabl.end()) {
      return RowIter->second.find(ColumnKey) != RowIter->second.end();
    }
    return false;
  }

  [[nodiscard]] bool containsColumn(C ColumnKey) const {
    // Returns true if the table contains a mapping with the specified column.
    for (const auto &M1 : Tabl) {
      if (M1.second.count(ColumnKey)) {
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] bool containsRow(R RowKey) const {
    // Returns true if the table contains a mapping with the specified row key.
    return Tabl.count(RowKey);
  }

  [[nodiscard]] bool containsValue(const V &Value) const {
    // Returns true if the table contains a mapping with the specified value.
    for (const auto &M1 : Tabl) {
      for (const auto &M2 : M1.second) {
        if (Value == M2.second) {
          return true;
        }
      }
    }
    return false;
  }

  [[nodiscard]] V &get(R RowKey, C ColumnKey) {
    // Returns the value corresponding to the given row and column keys, or null
    // if no such mapping exists.
    return Tabl[RowKey][ColumnKey];
  }

  V remove(R RowKey, C ColumnKey) {
    // Removes the mapping, if any, associated with the given keys.
    V Va = Tabl[RowKey][ColumnKey];
    Tabl[RowKey].erase(ColumnKey);
    return Va;
  }

  void remove(R RowKey) { Tabl.erase(RowKey); }

  [[nodiscard]] std::unordered_map<C, V> &row(R RowKey) {
    // Returns a view of all mappings that have the given row key.
    return Tabl[RowKey];
  }

  [[nodiscard]] std::multiset<R> rowKeySet() const {
    // Returns a set of row keys that have one or more values in the table.
    std::multiset<R> S;
    for (const auto &M1 : Tabl) {
      S.insert(M1.first);
    }
    return S;
  }

  [[nodiscard]] std::unordered_map<R, std::unordered_map<C, V>> rowMap() const {
    // Returns a view that associates each row key with the corresponding map
    // from column keys to values.
    return Tabl;
  }

  [[nodiscard]] std::multiset<V> values() const {
    // Returns a collection of all values, which may contain duplicates.
    std::multiset<V> S;
    for (const auto &M1 : Tabl) {
      for (const auto &M2 : M1.second) {
        S.insert(M2.second);
      }
    }
    return S;
  }

  friend bool operator==(const Table<R, C, V> &Lhs, const Table<R, C, V> &Rhs) {
    return Lhs.Tabl == Rhs.Tabl;
  }

  friend bool operator<(const Table<R, C, V> &Lhs, const Table<R, C, V> &Rhs) {
    return Lhs.Tabl < Rhs.Tabl;
  }

  friend std::ostream &operator<<(std::ostream &Os, const Table<R, C, V> &T) {
    for (const auto &M1 : T.Tabl) {
      for (const auto &M2 : M1.second) {
        Os << "< " << M1.first << " , " << M2.first << " , " << M2.second
           << " >\n";
      }
    }
    return Os;
  }
};

} // namespace psr

#endif
