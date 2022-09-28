// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFINITEELEMENTTRAITS_HH
#define DUNE_LOCALFINITEELEMENTTRAITS_HH

namespace Dune {

  //! traits helper struct
  template<class LB, class LC, class LI>
  struct LocalFiniteElementTraits
  {
    /** \todo Please doc me !
     */
    typedef LB LocalBasisType;

    /** \todo Please doc me !
     */
    typedef LC LocalCoefficientsType;

    /** \todo Please doc me !
     */
    typedef LI LocalInterpolationType;
  };

}

#endif
