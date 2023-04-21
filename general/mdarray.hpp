// Copyright (c) 2010-2023, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_MDARRAY
#define MFEM_MDARRAY

#include "../config/config.hpp"

#include "array.hpp"
#include "mdspan.hpp"

namespace mfem
{

template<typename T, int N, typename Layout = MDLayoutLeft<N>>
struct MDArray : public MDSpan<mfem::Array<T>, T, N, Layout>
{
   using base_t = MDSpan<mfem::Array<T>, T, N, Layout>;

   using mfem::Array<T>::Read;
   using mfem::Array<T>::Write;
   using mfem::Array<T>::HostRead;
   using mfem::Array<T>::HostWrite;

   MDArray(): base_t() { }

   template <typename... Ts>
   MDArray(md_size_t n, Ts... args): MDArray(args...)
   {
      base_t::Setup(n, args...);
   }
};

} // namespace mfem

#endif // MFEM_MDARRAY
