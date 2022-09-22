<!--
SPDX-FileCopyrightText: Copyright (C) DUNE Project contributors, see file LICENSE.md in module root
SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
-->

DUNE-library
============

DUNE, the Distributed and Unified Numerics Environment is a modular toolbox
for solving partial differential equations with grid-based methods.

The main intention is to create slim interfaces allowing an efficient use of
legacy and/or new libraries. Using C++ techniques DUNE allows to use very
different implementation of the same concept (i.e. grid, solver, ...) under
a common interface with a very low overhead.

DUNE was designed with flexibility in mind. It supports easy discretization
using methods, like Finite Elements, Finite Volume and also Finite
Differences.

Dune-localfunctions
-------------------
This module provides interfaces and implementations for local
finite element ansatz spaces, i.e., shape functions, associated
interpolation functionals, and relation to grid entities.

More information
----------------

Check dune-common for more details concerning dependencies, known bugs,
license and installation.
