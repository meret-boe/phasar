/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * MultiKeyTable.h
 *
 *  Created on: 07.02.2017
 *      Author: pdschbrt
 */

#ifndef PHASAR_UTILS_MULTIKEYTABLE_H_
#define PHASAR_UTILS_MULTIKEYTABLE_H_

#include <ostream>
#include <set>
#include <tuple>
#include <unordered_map>

namespace psr {

template <typename R, typename C, typename V> class MultiKeyTable {
private:
  std::unordered_multimap<R, std::unordered_multimap<C, V>> MultiKeyTableObj;

public:
  struct Cell {
    R Row;
    C Col;
    V Val;
    Cell(R Row, C Col, V Val) : Row(Row), Col(Col), Val(Val) {}
    R getRowKey() { return Row; }
    C getColumnKey() { return Col; }
    V getValue() { return Val; }
    friend std::ostream &operator<<(std::ostream &Os, const Cell &Cel) {
      return Os << "Cell: " << Cel.Row << ", " << Cel.Col << ", " << Cel.Val;
    }
    friend bool operator<(const Cell &Lhs, const Cell &Rhs) {
      return std::tie(Lhs.Row, Lhs.Col, Lhs.Val) < std::tie(Rhs.Row, Rhs.Col, Rhs.Val);
    }
    friend bool operator==(const Cell &Lhs, const Cell &Rhs) {
      return std::tie(Lhs.Row, Lhs.Col, Lhs.Val) == std::tie(Rhs.Row, Rhs.Col, Rhs.Val);
    }
  };

  MultiKeyTable() = default;

  ~MultiKeyTable() = default;

  MultiKeyTable(const MultiKeyTable &T) = default;

  MultiKeyTable(MultiKeyTable &&T)  noexcept = default;

  MultiKeyTable &operator=(const MultiKeyTable &T) = default;

  MultiKeyTable &operator=(MultiKeyTable &&T)  noexcept = default;

  V insert(R Ro, C Co, V Va) {
    // Associates the specified value with the specified keys.
    MultiKeyTableObj[Ro][Co] = Va;
    return Va;
  }

  void clear() { MultiKeyTableObj.clear(); }

  bool empty() { return MultiKeyTableObj.empty(); }

  size_t size() { return MultiKeyTableObj.size(); }

  std::multiset<Cell> cellSet() {
    // Returns a set of all row key / column key / value triplets.
    std::multiset<Cell> S;
    for (auto &M1 : MultiKeyTableObj) {
      for (auto &M2 : M1.second) {
        S.insert(Cell(M1.first, M2.first, M2.second));
      }
    }
    return S;
  }

  std::unordered_multimap<R, V> column(C ColumnKey) {
    // Returns a view of all mappings that have the given column key.
    std::unordered_multimap<R, V> Column;
    for (auto &Row : MultiKeyTableObj) {
      if (Row.second.count(ColumnKey)) {
        Column[Row.first] = Row.second[ColumnKey];
      }
    }
    return Column;
  }

  std::multiset<C> columnKeySet() {
    // Returns a set of column keys that have one or more values in the
    // table.
    std::multiset<C> Colkeys;
    for (auto &M1 : MultiKeyTableObj){
      for (auto &M2 : M1.second) {
        Colkeys.insert(M2.first);
      }
}
    return Colkeys;
  }

  std::unordered_multimap<C, std::unordered_multimap<R, V>> columnMap() {
    // Returns a view that associates each column key with the corresponding
    // map from row keys to values.
    std::unordered_multimap<C, std::unordered_multimap<R, V>> Columnmap;
    for (auto &M1 : MultiKeyTableObj) {
      for (auto &M2 : MultiKeyTableObj.second) {
        Columnmap[M2.first][M1.first] = M2.second;
      }
    }
    return Columnmap;
  }

  bool contains(R RowKey, C ColumnKey) {
    // Returns true if the table contains a mapping with the specified row
    // and column keys.
    if (MultiKeyTableObj.count(RowKey)) {
      return MultiKeyTableObj[RowKey].count(ColumnKey);
}
    return false;
  }

  bool containsColumn(C ColumnKey) {
    // Returns true if the table contains a mapping with the specified
    // column.
    for (auto &M1 : MultiKeyTableObj) {
      if (M1.second.count(ColumnKey)) {
        return true;
}
}
    return false;
  }

  bool containsRow(R RowKey) {
    // Returns true if the table contains a mapping with the specified row
    // key.
    return MultiKeyTableObj.count(RowKey);
  }

  bool containsValue(V Value) {
    // Returns true if the table contains a mapping with the specified
    // value.
    for (auto &M1 : MultiKeyTableObj) {
      for (auto &M2 : M1.second) {
        if (Value == M2.second) {
          return true;
}
}
}
    return false;
  }

  V get(R RowKey, C ColumnKey) {
    // Returns the value corresponding to the given row and column keys, or null
    // if no such mapping exists.
    // if (table.count(rowKey))
    //	if (table[rowKey].count(columnKey))
    return MultiKeyTableObj[RowKey][ColumnKey];
    // return V();
  }

  V remove(R RowKey, C ColumnKey) {
    // Removes the mapping, if any, associated with the given keys.
    V Va;
    if (contains(RowKey, ColumnKey)) {
      Va = MultiKeyTableObj[RowKey][ColumnKey];
      MultiKeyTableObj[RowKey].erase(ColumnKey);
    }
    return Va;
  }

  std::unordered_multimap<C, V> &row(R RowKey) {
    // Returns a view of all mappings that have the given row key.
    return MultiKeyTableObj[RowKey];
  }

  std::multiset<R> rowKeySet() {
    // Returns a set of row keys that have one or more values in the table.
    std::multiset<R> S;
    for (auto &M1 : MultiKeyTableObj) {
      S.insert(M1.first);
}
    return S;
  }

  std::unordered_multimap<R, std::unordered_multimap<C, V>> rowMap() {
    // Returns a view that associates each row key with the corresponding
    // map from column keys to values.
    return MultiKeyTableObj;
  }

  std::multiset<V> values() {
    // Returns a collection of all values, which may contain duplicates.
    std::multiset<V> S;
    for (auto &M1 : MultiKeyTableObj) {
      for (auto &M2 : M1.second) {
        S.insert(M2.second);
}
}
    return S;
  }

  friend bool operator==(const MultiKeyTable<R, C, V> &Lhs,
                         const MultiKeyTable<R, C, V> &Rhs) {
    return Lhs.MultiKeyTableObj == Rhs.MultiKeyTableObj;
  }

  friend bool operator<(const MultiKeyTable<R, C, V> &Lhs,
                        const MultiKeyTable<R, C, V> &Rhs) {
    return Lhs.MultiKeyTableObj < Rhs.MultiKeyTableObj;
  }

  friend std::ostream &operator<<(std::ostream &Os,
                                  const MultiKeyTable<R, C, V> &T) {
    for (auto &M1 : T.multi_key_table) {
      for (auto &M2 : M1.second) {
        Os << "< " << M1.first << " , " << M2.first << " , " << M2.second
           << " >\n";
}
}
    return Os;
  }
};

} // namespace psr

#endif
