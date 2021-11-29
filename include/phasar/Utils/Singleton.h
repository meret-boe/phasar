/******************************************************************************
 * Copyright (c) 2017 Philipp Schubert.
 * All rights reserved. This program and the accompanying materials are made
 * available under the terms of LICENSE.txt.
 *
 * Contributors:
 *     Philipp Schubert and others
 *****************************************************************************/

/*
 * Singleton.h
 *
 *  Created on: 30.08.2016
 *      Author: pdschbrt
 */

#ifndef PHASAR_UTILS_SINGLETON_H_
#define PHASAR_UTILS_SINGLETON_H_

namespace psr {

template <typename T> class Singleton {
public:
  Singleton(const Singleton &S) = delete;
  Singleton(Singleton &&S) = delete;
  Singleton &operator=(const Singleton &S) = delete;
  Singleton &operator=(Singleton &&S) = delete;
  static T &instance() {
    static T Value;
    return Value;
  }

protected:
  Singleton() = default;
  ~Singleton() = default;
};

} // namespace psr

#endif
