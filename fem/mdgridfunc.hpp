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

#ifndef MFEM_MDGRIDFUNC
#define MFEM_MDGRIDFUNC

#include "../config/config.hpp"

#include "fem/gridfunc.hpp"
#include "general/mdspan.hpp"

namespace mfem
{

template<md_size_t N, class Layout = MDLayoutLeft<N>>
class MDGridFunction : public MDSpan<GridFunction, double, N, Layout>
{
   using base_type = MDSpan<GridFunction, double, N, Layout>;
   using base_type::Nd;
   using base_type::Sd;

private:
   mutable md_size_t get_vdofs_offset = 0;

public:
   using GridFunction::Read;
   using GridFunction::Write;
   using GridFunction::HostRead;
   using GridFunction::HostWrite;

   MDGridFunction(): base_type() { }

   MDGridFunction(MDGridFunction&&) = delete;
   MDGridFunction(const MDGridFunction&) = delete;
   MDGridFunction& operator=(MDGridFunction&&) = delete;

   template <typename... Ts>
   MDGridFunction(FiniteElementSpace *ifes, Ts... args): MDGridFunction(args...)
   {
      this->GridFunction::fes = ifes;
      base_type::Setup(ifes->GetNDofs(), args...);
   }

   operator Vector&() { assert(false); }
   operator const Vector() const { assert(false); }

   template <typename... Ts>
   MDGridFunction(md_size_t n, Ts... args): MDGridFunction(args...)
   {
      base_type::Setup(n, args...);
   }

   template <md_size_t n = 1, typename... Ts>
   void GetVDofs(Array<int> &dofs, Ts... args) const
   {
      MFEM_VERIFY(dofs.Size() == Nd[n-1], "Error in difs size!");
      for (int s = 0; s < dofs.Size(); s++)
      {
         dofs[s] = get_vdofs_offset;
         dofs[s] += MDOffset<n,N,md_size_t,Ts...>::offset(Sd, s, args...);
      }
      get_vdofs_offset = 0;
   }

   template <md_size_t n = 1, typename... Ts>
   void GetVDofs(md_size_t dim, Ts&&... args) const
   {
      get_vdofs_offset += dim * Sd[n-1];
      MDGridFunction::GetVDofs<n+1>(std::forward<Ts>(args)...);
   }
};


} // namespace mfem

#endif // MFEM_MDGRIDFUNC
