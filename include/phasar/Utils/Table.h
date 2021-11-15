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
<<<<<<< HEAD
  std::unordered_map<R, std::unordered_map<C, V>> Tabl;
=======
  std::unordered_map<RTy, std::unordered_map<CTy, VTy>> Data;
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506

public:
  struct Cell {
    Cell() = default;
<<<<<<< HEAD
    Cell(R Row, C Col, V Val) : Row(Row), Col(Col), Val(Val) {}
=======
    Cell(RTy Row, CTy Col, VTy Val) : R(Row), C(Col), V(Val) {}
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
    ~Cell() = default;
    Cell(const Cell &) = default;
    Cell &operator=(const Cell &) = default;
    Cell(Cell &&) noexcept = default;
    Cell &operator=(Cell &&) noexcept = default;

<<<<<<< HEAD
    R getRowKey() const { return Row; }
    C getColumnKey() const { return Col; }
    V getValue() const { return Val; }

    friend std::ostream &operator<<(std::ostream &Os, const Cell &Ce) {
      return Os << "Cell: " << Ce.Row << ", " << Ce.C << ", " << Ce.V;
    }
    friend bool operator<(const Cell &Lhs, const Cell &Rhs) {
      return std::tie(Lhs.Row, Lhs.Col, Lhs.Val) < std::tie(Rhs.Row, Rhs.Col, Rhs.Val);
=======
    RTy getRowKey() const { return R; }
    CTy getColumnKey() const { return C; }
    VTy getValue() const { return V; }

    friend std::ostream &operator<<(std::ostream &Os, const Cell &Cel) {
      return Os << "Cell: " << Cel.R << ", " << Cel.C << ", " << Cel.V;
    }
    friend bool operator<(const Cell &Lhs, const Cell &Rhs) {
      return std::tie(Lhs.R, Lhs.Ce, Lhs.Va) < std::tie(Rhs.R, Rhs.C, Rhs.V);
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
    }
    friend bool operator==(const Cell &Lhs, const Cell &Rhs) {
      return std::tie(Lhs.Row, Lhs.Col, Lhs.Val) == std::tie(Rhs.Row, Rhs.Col, Rhs.Val);
    }

  private:
<<<<<<< HEAD
    R Row;
    C Col;
    V Val;
=======
    RTy R;
    CTy C;
    VTy V;
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
  };

  Table() = default;
  Table(const class Table &T) = default;
  class Table &operator=(const class Table &T) = default;
  Table(class Table &&T) noexcept = default;
  class Table &operator=(class Table &&T) noexcept = default;
  ~Table() = default;

<<<<<<< HEAD
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
=======
  void insert(RTy Row, CTy Cel, VTy Val) {
    // Associates the specified value with the specified keys.
    Data[Row][Cel] = std::move(Val);
  }

  void insert(const Table &T) { Data.insert(T.Data.begin(), T.Data.end()); }

  void clear() { Data.clear(); }

  [[nodiscard]] bool empty() const { return Data.empty(); }

  [[nodiscard]] size_t size() const { return Data.size(); }
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506

  [[nodiscard]] std::set<Cell> cellSet() const {
    // Returns a set of all row key / column key / value triplets.
    std::set<Cell> S;
<<<<<<< HEAD
    for (const auto &M1 : Tabl) {
=======
    for (const auto &M1 : Data) {
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
      for (const auto &M2 : M1.second) {
        S.emplace(M1.first, M2.first, M2.second);
      }
    }
    return S;
  }

  [[nodiscard]] std::vector<Cell> cellVec() const {
    // Returns a vector of all row key / column key / value triplets.
<<<<<<< HEAD
    std::vector<Cell> Va;
    for (const auto &M1 : Tabl) {
=======
    std::vector<Cell> V;
    for (const auto &M1 : Data) {
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
      for (const auto &M2 : M1.second) {
        Va.emplace_back(M1.first, M2.first, M2.second);
      }
    }
    return Va;
  }

  [[nodiscard]] std::unordered_map<RTy, VTy> column(CTy ColumnKey) const {
    // Returns a view of all mappings that have the given column key.
<<<<<<< HEAD
    std::unordered_map<R, V> Column;
    for (const auto &Row : Tabl) {
=======
    std::unordered_map<RTy, VTy> Column;
    for (const auto &Row : Data) {
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
      if (Row.second.count(ColumnKey)) {
        Column[Row.first] = Row.second[ColumnKey];
      }
    }
    return Column;
  }

  [[nodiscard]] std::multiset<CTy> columnKeySet() const {
    // Returns a set of column keys that have one or more values in the table.
<<<<<<< HEAD
    std::multiset<C> Colkeys;
    for (const auto &M1 : Tabl) {
=======
    std::multiset<CTy> Colkeys;
    for (const auto &M1 : Data) {
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
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
<<<<<<< HEAD
    std::unordered_map<C, std::unordered_map<R, V>> Columnmap;
    for (const auto &M1 : Tabl) {
      for (const auto &M2 : Tabl.second) {
=======
    std::unordered_map<CTy, std::unordered_map<RTy, VTy>> Columnmap;
    for (const auto &M1 : Data) {
      for (const auto &M2 : Data.second) {
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
        Columnmap[M2.first][M1.first] = M2.second;
      }
    }
    return Columnmap;
  }

  [[nodiscard]] bool contains(RTy RowKey, CTy ColumnKey) const {
    // Returns true if the table contains a mapping with the specified row and
    // column keys.
<<<<<<< HEAD
    if (auto RowIter = Tabl.find(RowKey); RowIter != Tabl.end()) {
=======
    if (auto RowIter = Data.find(RowKey); RowIter != Data.end()) {
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
      return RowIter->second.find(ColumnKey) != RowIter->second.end();
    }
    return false;
  }

  [[nodiscard]] bool containsColumn(CTy ColumnKey) const {
    // Returns true if the table contains a mapping with the specified column.
<<<<<<< HEAD
    for (const auto &M1 : Tabl) {
=======
    for (const auto &M1 : Data) {
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
      if (M1.second.count(ColumnKey)) {
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] bool containsRow(RTy RowKey) const {
    // Returns true if the table contains a mapping with the specified row key.
<<<<<<< HEAD
    return Tabl.count(RowKey);
=======
    return Data.count(RowKey);
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
  }

  [[nodiscard]] bool containsValue(const VTy &Value) const {
    // Returns true if the table contains a mapping with the specified value.
<<<<<<< HEAD
    for (const auto &M1 : Tabl) {
=======
    for (const auto &M1 : Data) {
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
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
<<<<<<< HEAD
    return Tabl[RowKey][ColumnKey];
=======
    return Data[RowKey][ColumnKey];
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
  }

  VTy remove(RTy RowKey, CTy ColumnKey) {
    // Removes the mapping, if any, associated with the given keys.
<<<<<<< HEAD
    V Va = Tabl[RowKey][ColumnKey];
    Tabl[RowKey].erase(ColumnKey);
    return Va;
  }

  void remove(R RowKey) { Tabl.erase(RowKey); }
=======
    VTy Val = Data[RowKey][ColumnKey];
    Data[RowKey].erase(ColumnKey);
    return Val;
  }

  void remove(RTy RowKey) { Data.erase(RowKey); }
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506

  [[nodiscard]] std::unordered_map<CTy, VTy> &row(RTy RowKey) {
    // Returns a view of all mappings that have the given row key.
<<<<<<< HEAD
    return Tabl[RowKey];
=======
    return Data[RowKey];
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
  }

  [[nodiscard]] std::multiset<RTy> rowKeySet() const {
    // Returns a set of row keys that have one or more values in the table.
<<<<<<< HEAD
    std::multiset<R> S;
    for (const auto &M1 : Tabl) {
=======
    std::multiset<RTy> S;
    for (const auto &M1 : Data) {
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
      S.insert(M1.first);
    }
    return S;
  }

  [[nodiscard]] std::unordered_map<RTy, std::unordered_map<CTy, VTy>> rowMap() const {
    // Returns a view that associates each row key with the corresponding map
    // from column keys to values.
<<<<<<< HEAD
    return Tabl;
=======
    return Data;
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
  }

  [[nodiscard]] std::multiset<VTy> values() const {
    // Returns a collection of all values, which may contain duplicates.
<<<<<<< HEAD
    std::multiset<V> S;
    for (const auto &M1 : Tabl) {
=======
    std::multiset<VTy> S;
    for (const auto &M1 : Data) {
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
      for (const auto &M2 : M1.second) {
        S.insert(M2.second);
      }
    }
    return S;
  }

<<<<<<< HEAD
  friend bool operator==(const Table<R, C, V> &Lhs, const Table<R, C, V> &Rhs) {
    return Lhs.Tabl == Rhs.Tabl;
  }

  friend bool operator<(const Table<R, C, V> &Lhs, const Table<R, C, V> &Rhs) {
    return Lhs.Tabl < Rhs.Tabl;
  }

  friend std::ostream &operator<<(std::ostream &Os, const Table<R, C, V> &T) {
    for (const auto &M1 : T.Tabl) {
=======
  friend bool operator==(const Table<RTy, CTy, VTy> &Lhs, const Table<RTy, CTy, VTy> &Rhs) {
    return Lhs.Data == Rhs.Data;
  }

  friend bool operator<(const Table<RTy, CTy, VTy> &Lhs, const Table<RTy, CTy, VTy> &Rhs) {
    return Lhs.Data < Rhs.Data;
  }

  friend std::ostream &operator<<(std::ostream &Os, const Table<RTy, CTy, VTy> &T) {
    for (const auto &M1 : T.Data) {
>>>>>>> b64c0176c1c39f7ad73feffb391fd6e22688d506
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
