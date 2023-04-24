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

#ifndef MFEM_MDVECTOR
#define MFEM_MDVECTOR

#include "../config/config.hpp"

#include "vector.hpp"
#include "general/mdspan.hpp"

namespace mfem
{

template<int N, typename Layout = MDLayoutLeft<N>>
struct MDVector : public MDSpan<mfem::Vector, N, Layout>
{
   using base_t = MDSpan<mfem::Vector, N, Layout>;

   using mfem::Vector::Read;
   using mfem::Vector::Write;
   using mfem::Vector::HostRead;
   using mfem::Vector::HostWrite;

   MDVector(): base_t() { }

   template <typename... Ts>
   MDVector(int n, Ts... args): MDVector(args...)
   {
      base_t::Setup(n, args...);
   }
};

} // namespace mfem

#endif // MFEM_MDVECTOR
