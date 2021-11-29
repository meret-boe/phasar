/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * HexaStoreGraph.hh
 *
 *  Created on: 06.02.2017
 *      Author: pdschbrt
 */

#ifndef PHASAR_UTILS_HEXASTOREGRAPH_H_
#define PHASAR_UTILS_HEXASTOREGRAPH_H_

#include <iostream> // cerr, to suppress once it is not used anymore
#include <ostream>
#include <set>

//#include "phasar/DB/DBConn.h"
#include "phasar/Utils/Table.h"

namespace psr {

template <typename S, typename E, typename D> class HexaStoreGraph {
private:
  Table<S, E, D> Sed;
  Table<S, D, E> Sde;
  Table<E, S, D> Esd;
  Table<E, D, S> Eds;
  Table<D, S, E> Dse;
  Table<D, E, S> Des;

public:
  HexaStoreGraph() = default;

  ~HexaStoreGraph() = default;

  void insertEdge(S Src, E Edge, D Dest) {
    Sed.insert(Src, Edge, Dest);
    Sde.insert(Src, Dest, Edge);
    Esd.insert(Edge, Src, Dest);
    Eds.insert(Edge, Dest, Src);
    Dse.insert(Dest, Src, Edge);
    Des.insert(Dest, Edge, Src);
  }

  void removeEdge(S Src, E Edge, D Dest) {
    Sed.remove(Src, Edge);
    Sde.remove(Src, Dest);
    Esd.remove(Edge, Src);
    Eds.remove(Edge, Dest);
    Dse.remove(Dest, Src);
    Des.remove(Dest, Edge);
  }

  bool containsEdge(S Src, E Edge, D  /*Dest*/) { return Sed.contains(Src, Edge); }

  D getDbySE(S Src, E Edge) { return Sed.get(Src, Edge); }

  E getEbySD(S Src, D Dest) { return Sde.get(Src, Dest); }

  D getDbyES(E Edge, S Src) { return Esd.get(Edge, Src); }

  S getSbyED(E Edge, D Dest) { return Eds.get(Edge, Dest); }

  E getEbyDS(D Dest, S Src) { return Dse.get(Dest, Src); }

  S getSyDE(D Dest, E Edge) { return Des.get(Dest, Edge); }

  void clear() {
    Sed.clear();
    Sde.clear();
    Esd.clear();
    Eds.clear();
    Dse.clear();
    Des.clear();
  }

  bool empty() { return Sed.empty(); }

  std::set<typename Table<S, E, D>::Cell> tripleSet() { return Sed.cellSet(); }

  std::multiset<S> sourceSet() { return Sed.rowKeySet(); }

  std::multiset<E> edgeSet() { return Sed.columnKeySet(); }

  std::multiset<D> destinationSet() { return Sed.values(); }

  void loadHexaStoreGraphFromDB(const std::string & /*Tablename*/) {
    std::cerr << "Not implemented yet!" << std::endl;
  }

  void storeHexaStoreGraphToDB(const std::string & /*Tablename*/) {
    std::cerr << "Not implemented yet!" << std::endl;
  }

  friend bool operator==(const HexaStoreGraph<S, E, D> &Lhs,
                         const HexaStoreGraph<S, E, D> &Rhs) {
    return Lhs.Sed == Rhs.Sed;
  }

  friend bool operator<(const HexaStoreGraph<S, E, D> &Lhs,
                        const HexaStoreGraph<S, E, D> &Rhs) {
    return Lhs.Sed < Rhs.Sed;
  }

  friend std::ostream &operator<<(std::ostream &Os, const HexaStoreGraph &Hsg) {
    return Os << Hsg.Sed;
  }
};

} // namespace psr

#endif
