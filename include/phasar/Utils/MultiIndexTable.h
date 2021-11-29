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
 *  Created on: 31.08.2016
 *      Author: pdschbrt
 */

#ifndef PHASAR_UTILS_MULTIINTEXTABLE_H_
#define PHASAR_UTILS_MULTIINTEXTABLE_H_

#include <ostream>
#include <string>

#include "boost/multi_index/composite_key.hpp"
#include "boost/multi_index/hashed_index.hpp"
#include "boost/multi_index/member.hpp"
#include "boost/multi_index/ordered_index.hpp"
#include "boost/multi_index_container.hpp"

namespace psr {

struct OrderedRowColKeyTag {};
struct HashedRowColKeyTag {};

template <typename R, typename C, typename V> struct MultiIndexTable {
  struct TableData {
    R Rowkey;
    C Columnkey;
    V Value;
    TableData(R Row, C Col, V Val)
        : Rowkey(Row), Columnkey(Col), Value(Val) {}
  };

  struct RowColKey
      : boost::multi_index::composite_key<
            TableData, BOOST_MULTI_INDEX_MEMBER(TableData, R, rowkey),
            BOOST_MULTI_INDEX_MEMBER(TableData, C, columnkey)> {};

  using InternTable = boost::multi_index_container<TableData, boost::multi_index::indexed_by<boost::multi_index::ordered_unique<boost::multi_index::tag<OrderedRowColKeyTag>, RowColKey>, boost::multi_index::hashed_unique<boost::multi_index::tag<HashedRowColKeyTag>, RowColKey>>>;

  // ordered indices
  using ordered_row_col_key_view_t = typename boost::multi_index::index<InternTable, OrderedRowColKeyTag>::type;
  using ordered_row_col_key_iterator_t = typename boost::multi_index::index<InternTable, OrderedRowColKeyTag>::type::const_iterator;

  // hashed indices
  using hashed_row_col_key_view_t = typename boost::multi_index::index<InternTable, HashedRowColKeyTag>::type;
  using hashed_row_col_key_iterator_t = typename boost::multi_index::index<InternTable, HashedRowColKeyTag>::type::const_iterator;

  // the indexed table containing instances of TableData
  InternTable IndexedTable;

  friend std::ostream &operator<<(std::ostream &Os, const InternTable & /*Itab*/) {
    return Os << "error: unsupported operation!";
  }
};

} // namespace psr

#endif
