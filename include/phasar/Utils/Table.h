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

template <typename RTy, typename CTy, typename VTy> class Table {
private:
  std::unordered_map<RTy, std::unordered_map<CTy, VTy>> Tabl;

public:
  struct Cell {
    Cell() = default;
    Cell(RTy Row, CTy Col, VTy Val) : Row(Row), Col(Col), Val(Val) {}
    ~Cell() = default;
    Cell(const Cell &) = default;
    Cell &operator=(const Cell &) = default;
    Cell(Cell &&) noexcept = default;
    Cell &operator=(Cell &&) noexcept = default;

    RTy getRowKey() const { return Row; }
    CTy getColumnKey() const { return Col; }
    VTy getValue() const { return Val; }

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
    RTy Row{};
    CTy Col{};
    VTy Val{};
  };

  Table() = default;
  Table(const class Table &T) = default;
  class Table &operator=(const class Table &T) = default;
  Table(class Table &&T) noexcept = default;
  class Table &operator=(class Table &&T) noexcept = default;
  ~Table() = default;

  void insert(RTy Ro, CTy Co, VTy Va) {
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

  [[nodiscard]] std::unordered_map<RTy, VTy> column(CTy ColumnKey) const {
    // Returns a view of all mappings that have the given column key.
    std::unordered_map<RTy, VTy> Column;
    for (const auto &Row : Tabl) {
      if (Row.second.count(ColumnKey)) {
        Column[Row.first] = Row.second[ColumnKey];
      }
    }
    return Column;
  }

  [[nodiscard]] std::multiset<CTy> columnKeySet() const {
    // Returns a set of column keys that have one or more values in the table.
    std::multiset<CTy> Colkeys;
    for (const auto &M1 : Tabl) {
      for (const auto &M2 : M1.second) {
        Colkeys.insert(M2.first);
      }
    }
    return Colkeys;
  }

  [[nodiscard]] std::unordered_map<CTy, std::unordered_map<RTy, VTy>>
  columnMap() const {
    // Returns a view that associates each column key with the corresponding map
    // from row keys to values.
    std::unordered_map<CTy, std::unordered_map<RTy, VTy>> Columnmap;
    for (const auto &M1 : Tabl) {
      for (const auto &M2 : Tabl.second) {
        Columnmap[M2.first][M1.first] = M2.second;
      }
    }
    return Columnmap;
  }

  [[nodiscard]] bool contains(RTy RowKey, CTy ColumnKey) const {
    // Returns true if the table contains a mapping with the specified row and
    // column keys.
    if (auto RowIter = Tabl.find(RowKey); RowIter != Tabl.end()) {
      return RowIter->second.find(ColumnKey) != RowIter->second.end();
    }
    return false;
  }

  [[nodiscard]] bool containsColumn(CTy ColumnKey) const {
    // Returns true if the table contains a mapping with the specified column.
    for (const auto &M1 : Tabl) {
      if (M1.second.count(ColumnKey)) {
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] bool containsRow(RTy RowKey) const {
    // Returns true if the table contains a mapping with the specified row key.
    return Tabl.count(RowKey);
  }

  [[nodiscard]] bool containsValue(const VTy &Value) const {
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

  [[nodiscard]] VTy &get(RTy RowKey, CTy ColumnKey) {
    // Returns the value corresponding to the given row and column keys, or null
    // if no such mapping exists.
    return Tabl[RowKey][ColumnKey];
  }

  VTy remove(RTy RowKey, CTy ColumnKey) {
    // Removes the mapping, if any, associated with the given keys.
    VTy Va = Tabl[RowKey][ColumnKey];
    Tabl[RowKey].erase(ColumnKey);
    return Va;
  }

  void remove(RTy RowKey) { Tabl.erase(RowKey); }

  [[nodiscard]] std::unordered_map<CTy, VTy> &row(RTy RowKey) {
    // Returns a view of all mappings that have the given row key.
    return Tabl[RowKey];
  }

  [[nodiscard]] std::multiset<RTy> rowKeySet() const {
    // Returns a set of row keys that have one or more values in the table.
    std::multiset<RTy> S;
    for (const auto &M1 : Tabl) {
      S.insert(M1.first);
    }
    return S;
  }

  [[nodiscard]] std::unordered_map<RTy, std::unordered_map<CTy, VTy>> rowMap() const {
    // Returns a view that associates each row key with the corresponding map
    // from column keys to values.
    return Tabl;
  }

  [[nodiscard]] std::multiset<VTy> values() const {
    // Returns a collection of all values, which may contain duplicates.
    std::multiset<VTy> S;
    for (const auto &M1 : Tabl) {
      for (const auto &M2 : M1.second) {
        S.insert(M2.second);
      }
    }
    return S;
  }

  friend bool operator==(const Table<RTy, CTy, VTy> &Lhs, const Table<RTy, CTy, VTy> &Rhs) {
    return Lhs.Tabl == Rhs.Tabl;
  }

  friend bool operator<(const Table<RTy, CTy, VTy> &Lhs, const Table<RTy, CTy, VTy> &Rhs) {
    return Lhs.Tabl < Rhs.Tabl;
  }

  friend std::ostream &operator<<(std::ostream &Os, const Table<RTy, CTy, VTy> &T) {
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
